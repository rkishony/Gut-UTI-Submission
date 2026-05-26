#!/usr/bin/env python3
"""Step: Create UTI–faecal pairs table for read alignment."""

import numpy as np
import pandas as pd
from pathlib import Path

from config_utils import get_output_dir
from steps.trim_stats import TrimStatsStep
from steps.assembly_stats import AssemblyStatsStep
from steps.step import Step, run_cli

N_RANDOM_PATIENTS = 10
RANDOM_SEED = 42


MUST_CHOOSE_PAIRS: list[tuple[str, str]] = [
    ("U15A_hs1", "F109meta_ns"),
]


def _make_pair(uti_run: str, faecal_run: str, pair_type: str, match_num: int, bio_rep: int = 0) -> dict:
    return {
        "run_name": f"{uti_run}__{faecal_run}",
        "uti_run_name": uti_run,
        "faecal_run_name": faecal_run,
        "pair_type": pair_type,
        "match_num": match_num,
        "bio_rep": bio_rep,
    }


class CreatePairsStep(Step):
    def __init__(self, n_random_patients: int = N_RANDOM_PATIENTS, **kwargs):
        super().__init__(**kwargs)
        self.n_random_patients = n_random_patients

    @classmethod
    def output_path(cls) -> Path:
        return get_output_dir("summary") / "uti_faecal_pairs.csv"

    @classmethod
    def uti_status_path(cls) -> Path:
        return get_output_dir("summary") / "uti_status.csv"

    def _run(self) -> None:
        df = pd.read_csv(TrimStatsStep.output_path()[0])

        uti_df = df[df["is_isolate"] & ~df["is_faecal"]].copy()
        faecal = df[~df["is_isolate"] & df["is_faecal"]].copy()

        # Identify UTIs with no assembly
        assembly_stats = pd.read_csv(AssemblyStatsStep.output_path())
        utis_with_no_assembly = set(assembly_stats.loc[~assembly_stats["has_fasta"], "run_name"])

        # Primary faecal sample per fid: the non-replicate (is_rep=False)
        primary_faecal = faecal[~faecal["is_rep"]].copy()
        assert primary_faecal["fid"].is_unique, "Expected exactly one is_rep=False faecal per fid"
        assert primary_faecal["bio_rep"].isin([0, 1]).all(), "is_rep=False should have bio_rep 0 or 1"
        primary_faecal_by_fid = primary_faecal.set_index("fid")
        all_faecal_fids = sorted(primary_faecal["fid"])

        # bio_rep=2 faecal samples by fid
        biorep2_by_fid = (
            faecal[faecal["bio_rep"] == 2]
            .drop_duplicates("fid", keep="first")
            .set_index("fid")["run_name"]
            .to_dict()
        )

        # bio_rep lookup: faecal run_name -> bio_rep
        faecal_biorep = faecal.set_index("run_name")["bio_rep"].to_dict()

        # Build UTI status table
        uti_status = uti_df[["run_name", "fid"]].copy()
        uti_status["has_assembly"] = ~uti_status["run_name"].isin(utis_with_no_assembly)
        uti_status["has_faecal"] = uti_status["fid"].isin(primary_faecal_by_fid.index)
        uti_status["excluded"] = ~uti_status["has_assembly"] | ~uti_status["has_faecal"]
        uti_status.to_csv(self.uti_status_path(), index=False)

        rng = np.random.default_rng(RANDOM_SEED)
        pairs = []

        for _, uti_row in uti_status[~uti_status["excluded"]].iterrows():
            uti_fid = uti_row["fid"]
            uti_run = uti_row["run_name"]

            def _add_pair(faecal_run, fae_fid, pair_type, match_num):
                pairs.append(_make_pair(uti_run, faecal_run, pair_type, match_num,
                                        bio_rep=faecal_biorep.get(faecal_run, 0)))
                br2 = biorep2_by_fid.get(fae_fid)
                if br2 and br2 != faecal_run:
                    pairs.append(_make_pair(uti_run, br2, pair_type, match_num,
                                            bio_rep=faecal_biorep.get(br2, 2)))

            # 1. Matched
            _add_pair(primary_faecal_by_fid.loc[uti_fid, "run_name"], uti_fid, "matched", 0)

            # 2. Must-choose
            match_num = 1
            must_faecal_runs = {f for u, f in MUST_CHOOSE_PAIRS if u == uti_run}
            must_fids = set()
            for faecal_run_name in must_faecal_runs:
                match = faecal[faecal["run_name"] == faecal_run_name]
                fae_fid = match.iloc[0]["fid"] if not match.empty else None
                _add_pair(faecal_run_name, fae_fid, "must_choose", match_num)
                if fae_fid is not None:
                    must_fids.add(fae_fid)
                match_num += 1

            # 3. Random
            n_remaining = self.n_random_patients - len(must_faecal_runs)
            other_fids = [f for f in all_faecal_fids if f != uti_fid and f not in must_fids]
            if n_remaining > 0 and other_fids:
                chosen_fids = rng.choice(other_fids, size=min(n_remaining, len(other_fids)), replace=False)
                for chosen_fid in chosen_fids:
                    _add_pair(primary_faecal_by_fid.loc[chosen_fid, "run_name"], chosen_fid, "random", match_num)
                    match_num += 1

        pairs_df = pd.DataFrame(pairs)
        faecal_reads = faecal.set_index("run_name")["num_reads"]
        pairs_df["faecal_num_reads"] = pairs_df["faecal_run_name"].map(faecal_reads).fillna(0).astype(int)
        out = self.output_path()
        pairs_df.to_csv(out, index=False)

        n_no_asm = (~uti_status["has_assembly"]).sum()
        n_no_fae = (~uti_status["has_faecal"]).sum()
        type_counts = pairs_df["pair_type"].value_counts()
        counts_str = ", ".join(f"{t}: {c}" for t, c in type_counts.items())
        return (f"UTI status: {len(uti_status)} total, "
                f"{n_no_asm} no assembly, {n_no_fae} no faecal. "
                f"Pairs: {len(pairs_df)} ({counts_str})")


if __name__ == "__main__":
    run_cli(CreatePairsStep)
