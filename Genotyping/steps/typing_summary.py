"""Step: Summarise StrainPhlAn + MetaMLST status for every faecal sample."""

import pandas as pd
from pathlib import Path
from Bio import Phylo

from config_utils import get_output_dir
from steps.step import Step, run_cli
from steps.trim_stats import TrimStatsStep
from steps.strainphlan import StrainphlanStep, CollectStrainphlanStep
from steps.metamlst import MetaMLSTStep, MetaMLSTMergeStep


class TypingSummaryStep(Step):
    """Per-faecal-sample summary of StrainPhlAn and MetaMLST success/failure."""

    @classmethod
    def output_path(cls) -> Path:
        return get_output_dir("summary") / "typing_summary.csv"

    def _load_tree_leaves(self) -> set[str]:
        tree_path = CollectStrainphlanStep.output_path()
        if not tree_path.exists():
            self.print(f"Tree not found: {tree_path}")
            return set()
        tree = Phylo.read(str(tree_path), "newick")
        return {t.name for t in tree.get_terminals()}

    @staticmethod
    def _read_report(path) -> dict[str, int]:
        st_map: dict[str, int] = {}
        if not path.exists():
            return st_map
        with open(path) as f:
            next(f)
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    st_map[parts[2]] = int(parts[0])
        return st_map

    def _load_metamlst_reports(self) -> tuple[dict[str, int], dict[str, int]]:
        """Return (e1_map, e2_map) from merged MetaMLST reports."""
        merged = MetaMLSTMergeStep.output_path()
        e1 = self._read_report(merged / "escherichia1_report.txt")
        e2 = self._read_report(merged / "escherichia2_report.txt")
        return e1, e2

    def _strainphlan_failure_reason(self, row: pd.Series) -> str:
        run_dir = StrainphlanStep.run_folder(row)
        failed = run_dir / ".FAILED.txt"
        if failed.exists():
            return "per-sample failed: " + failed.read_text().strip()[:120]
        pkl = run_dir / "markers.pkl"
        if not pkl.exists():
            if not run_dir.exists():
                return "not run"
            return "markers.pkl missing"
        return "dropped by tree builder (insufficient markers)"

    def _metamlst_failure_reason(self, row: pd.Series) -> str:
        run_dir = MetaMLSTStep.run_folder(row)
        failed = run_dir / ".FAILED.txt"
        if failed.exists():
            return "per-sample failed: " + failed.read_text().strip()[:120]
        typing_dir = run_dir / "typing"
        if not typing_dir.exists():
            if not run_dir.exists():
                return "not run"
            return "typing dir missing"
        nfo_files = list(typing_dir.glob("*.nfo"))
        if not nfo_files:
            return "no .nfo produced (no E. coli detected)"
        species = set()
        for nfo in nfo_files:
            for line in nfo.read_text().splitlines():
                if line.strip():
                    species.add(line.split("\t")[0])
        return f"typed as other species: {', '.join(sorted(species))}"

    def _run(self) -> str:
        df = pd.read_csv(TrimStatsStep.output_path()[0])
        faecal = df[~df["is_isolate"]].copy()

        tree_leaves = self._load_tree_leaves()
        e1_map, e2_map = self._load_metamlst_reports()

        rows = []
        for _, row in faecal.iterrows():
            name = row["run_name"]
            fid = row["fid"]
            is_rep = row["is_rep"]
            bio_rep = row["bio_rep"]

            has_markers = (StrainphlanStep.run_folder(row) / "markers.pkl").exists()
            in_tree = name in tree_leaves
            has_e1 = name in e1_map
            has_e2 = name in e2_map
            mlst_scheme = "escherichia1" if has_e1 else ("escherichia2" if has_e2 else "")
            has_mlst_st = has_e1 or has_e2
            mlst_st = e1_map.get(name) or e2_map.get(name)

            sp_reason = ""
            if not in_tree:
                sp_reason = self._strainphlan_failure_reason(row)
            mm_reason = ""
            if not has_mlst_st:
                mm_reason = self._metamlst_failure_reason(row)

            in_pie_line = in_tree and not is_rep
            has_mlst_not_in_pie = has_mlst_st and not in_pie_line
            has_mlst_not_in_pie_and_not_rep = has_mlst_not_in_pie and not is_rep

            rows.append({
                "run_name": name,
                "fid": fid,
                "is_rep": is_rep,
                "bio_rep": bio_rep,
                "strainphlan_markers": has_markers,
                "strainphlan_in_tree": in_tree,
                "strainphlan_fail_reason": sp_reason,
                "metamlst_has_st": has_mlst_st,
                "metamlst_scheme": mlst_scheme,
                "metamlst_st": mlst_st,
                "metamlst_fail_reason": mm_reason,
                "in_pie_line_chart": in_pie_line,
                "has_mlst_not_in_pie": has_mlst_not_in_pie,
                "has_mlst_not_in_pie_and_not_rep": has_mlst_not_in_pie_and_not_rep,
            })

        out = pd.DataFrame(rows)
        out["metamlst_st"] = out["metamlst_st"].astype("Int64")
        out.to_csv(self.output_path(), index=False)

        n = len(out)
        n_nonrep = len(out[~out["is_rep"]])
        n_tree = out["strainphlan_in_tree"].sum()
        n_mlst = out["metamlst_has_st"].sum()
        n_gap = out["has_mlst_not_in_pie"].sum()
        summary = (
            f"{n} faecal samples ({n_nonrep} non-rep): "
            f"{n_tree} in tree, {n_mlst} with MLST ST, "
            f"{n_gap} have MLST but missing from pie-line chart"
        )
        self.print(summary)
        return summary


if __name__ == "__main__":
    run_cli(TypingSummaryStep)
