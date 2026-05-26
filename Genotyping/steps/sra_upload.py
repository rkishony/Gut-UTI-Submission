#!/usr/bin/env python3
"""Step: Upload trimmed reads to NCBI SRA."""

import pandas as pd
from pathlib import Path

import config
from config_utils import get_output_dir
from steps.per_run_step import PerRunStep
from steps.step import Step, run_cli
from utils.metagenome_io import get_metagenome_reads
from tools.sra_ftp import upload_to_sra_ftp


class SRAUploadStep(PerRunStep):
    """Upload trimmed FASTQs to NCBI SRA FTP preload area."""

    output_paths = "uploaded.txt"
    output_dir_name = "sra_upload"
    cores_per_run = 1

    def run_input_paths(self, row: pd.Series) -> tuple:
        return get_metagenome_reads(row)

    def setup(self) -> None:
        assert config.SRA_BIOPROJECT, "Set SRA_BIOPROJECT in config.py before uploading"
        assert config.SRA_FTP_USERNAME, "Set SRA_FTP_USERNAME in config.py"

    def process_run(self, row: pd.Series) -> None:
        run_name = row["run_name"]
        r1, r2 = self.run_input_paths(row)

        assert r1.exists(), f"Missing {r1}"
        assert r2.exists(), f"Missing {r2}"

        upload_to_sra_ftp([r1, r2], subfolder=run_name)

        receipt = self.run_output_path(row)
        receipt.write_text(f"Uploaded {r1.name} and {r2.name} to SRA FTP\n")
        self.print(f"uploaded {r1.name}, {r2.name}", row)


class SRAMetadataStep(Step):
    """Generate SRA metadata TSV for batch submission."""

    @classmethod
    def output_path(cls) -> Path:
        return get_output_dir("sra_upload") / "sra_metadata.tsv"

    def _run(self) -> None:
        df = pd.read_csv(
            next(get_output_dir("summary").glob("sequencing_runs.csv"))
        )

        rows = []
        for _, row in df.iterrows():
            run_name = row["run_name"]
            r1, r2 = get_metagenome_reads(row)
            rows.append({
                "bioproject_accession": config.SRA_BIOPROJECT,
                "biosample_accession": "",  # fill after BioSample registration
                "library_ID": run_name,
                "title": f"Paired-end sequencing of {run_name}",
                "library_strategy": "WGS",
                "library_source": "GENOMIC",
                "library_selection": "RANDOM",
                "library_layout": "paired",
                "platform": "ILLUMINA",
                "instrument_model": "Illumina NovaSeq 6000",
                "design_description": "Whole genome sequencing",
                "filetype": "fastq",
                "filename": r1.name,
                "filename2": r2.name,
            })

        meta = pd.DataFrame(rows)
        out = self.output_path()
        out.parent.mkdir(parents=True, exist_ok=True)
        meta.to_csv(out, sep="\t", index=False)
        return f"{len(meta)} entries"


if __name__ == "__main__":
    run_cli(SRAMetadataStep, SRAUploadStep)
