"""Steps: StrainPhlAn per-sample marker extraction + collective tree building."""

import shutil

import pandas as pd
from pathlib import Path

from config import STRAINPHLAN_DEFAULT_CLADE
from config_utils import get_output_dir
from steps.extract_run_names import ExtractRunNamesStep
from steps.per_run_step import PerRunStep
from utils.file_utils import create_if_missing
from utils.metagenome_io import get_metagenome_reads
from steps.step import Step, run_cli
from tools.metaphlan import metaphlan_profile, sample2markers


class StrainphlanStep(PerRunStep):
    """Per-run: MetaPhlAn profiling + sample2markers → .pkl consensus markers."""

    output_dir_name = "strainphlan"
    output_paths = "markers.pkl"
    cores_per_run = 4
    clean_run_folder = False

    def run_input_paths(self, row: pd.Series) -> tuple:
        return get_metagenome_reads(row)

    def filter_df(self, df: pd.DataFrame) -> pd.DataFrame:
        return df[~df["is_isolate"]]

    def process_run(self, row: pd.Series) -> None:
        run_name = row["run_name"]
        r1, r2 = self.run_input_paths(row)
        run_dir = self.run_folder(row)

        sam = run_dir / f"{run_name}.sam.bz2"
        profile = run_dir / "profile.tsv"
        bowtie2 = run_dir / "bowtie2.bz2"

        nproc = self.threads_per_worker
        if not sam.exists():
            metaphlan_profile(r1, r2, profile, sam, bowtie2, nproc=nproc)

        sample2markers(sam, run_dir, nproc=nproc)

        # sample2markers names the PKL after the SAM stem; rename to standard name
        tool_pkl = run_dir / f"{run_name}.pkl"
        target_pkl = run_dir / "markers.pkl"
        if tool_pkl.exists() and tool_pkl != target_pkl:
            tool_pkl.rename(target_pkl)

        sam.unlink(missing_ok=True)
        bowtie2.unlink(missing_ok=True)
        for tmp_dir in run_dir.glob("tmp*"):
            if tmp_dir.is_dir():
                shutil.rmtree(tmp_dir, ignore_errors=True)


class CollectStrainphlanStep(Step):
    """Collective: run strainphlan on all consensus markers → phylogenetic tree."""

    allow_failure = True

    def __init__(self, csv_file: str | Path | None = None,
                 clade: str = STRAINPHLAN_DEFAULT_CLADE, **kwargs):
        super().__init__(**kwargs)
        self.csv_file = csv_file or ExtractRunNamesStep.output_path()
        self.clade = clade

    @classmethod
    def output_path(cls, clade: str = STRAINPHLAN_DEFAULT_CLADE) -> Path:
        return (get_output_dir("strainphlan") / "tree"
                / f"RAxML_bestTree.{clade}.StrainPhlAn4.tre")

    def _run(self) -> None:
        from tools.strainphlan import run_strainphlan

        df = pd.read_csv(self.csv_file)
        df = df[~df["is_isolate"]]
        df = df[~df["is_rep"]]

        output_dir = self.output_path(self.clade).parent

        # strainphlan uses os.path.basename() as sample ID, so each PKL
        # must have a unique filename.  Our per-run layout stores them all
        # as "markers.pkl".  Symlink into a flat staging dir with run_name.
        staging = create_if_missing(output_dir / "consensus_markers", clean=True)

        pkls = []
        missing = []
        for _, row in df.iterrows():
            src = StrainphlanStep.run_output_path(row)
            if not src.exists():
                missing.append(row['run_name'])
                continue
            dst = staging / f"{row['run_name']}.pkl"
            dst.symlink_to(src.resolve())
            pkls.append(dst)
        if missing:
            self.print(f"Skipping {len(missing)} samples with missing PKLs")

        self.print(f"Running strainphlan on {len(pkls)} samples, clade={self.clade}")
        try:
            tree = run_strainphlan(
                consensus_markers=pkls,
                output_dir=output_dir,
                clade=self.clade,
                nproc=self.total_cores,
            )
        finally:
            shutil.rmtree(staging, ignore_errors=True)
        return "Tree built"


if __name__ == "__main__":
    run_cli(StrainphlanStep, CollectStrainphlanStep)
