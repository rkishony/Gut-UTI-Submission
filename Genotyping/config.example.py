"""Example configuration for the Genotyping pipeline.

Copy this file to `config.py` and edit the paths below to point at your local
FASTQ downloads, tool binaries, and reference databases. `config.py` is not
tracked in git so each site keeps its own copy.
"""
from pathlib import Path

# ---------------------------------------------------------------------------
# Raw data folders
# ---------------------------------------------------------------------------
# One entry per sequencing batch. Keys are short batch tags that appear as a
# suffix in run names (e.g. `U32A_hs2`); values are the directories containing
# the raw paired-end FASTQ files for that batch.
SEQ_BATCHES_TO_PATHS = {
    "ns":  Path("/path/to/raw/novaseq_batch"),
    "hs1": Path("/path/to/raw/hiseq_batch_1"),
    "hs2": Path("/path/to/raw/hiseq_batch_2"),
}

# ---------------------------------------------------------------------------
# Sample sheets (relative to this file)
# ---------------------------------------------------------------------------
ROOT_DIR = Path(__file__).resolve().parent
UTI_ASSIGNMENT_CSV = ROOT_DIR / "inputs" / "UTIassembly_F_assignment.csv"
META_FOBTS_CSV     = ROOT_DIR / "inputs" / "meta_fobts.csv"

# ---------------------------------------------------------------------------
# Output directory
# ---------------------------------------------------------------------------
OUTPUT_DIR = Path("/path/to/output")

# ---------------------------------------------------------------------------
# External tool paths
# ---------------------------------------------------------------------------
# Absolute paths to each external binary/script. Adjust to match your install
# (e.g. conda envs, system packages).
FASTX_BIN = Path("/path/to/fastx/bin")

TOOLS = {
    "unicycler":         "/path/to/envs/unicycler/bin/unicycler",
    "samtools":          "/path/to/samtools",
    "bedtools":          "/path/to/bedtools",
    "bowtie2":           "/usr/bin/bowtie2",
    "bowtie2-build":     "/usr/bin/bowtie2-build",
    "cutadapt":          "/path/to/envs/marta_ms/bin/cutadapt",
    "seqtk":             "/path/to/envs/marta_ms/bin/seqtk",
    "metaphlan":         "/path/to/envs/metaphlan/bin/metaphlan",
    "sample2markers.py": "/path/to/envs/metaphlan/bin/sample2markers.py",
    "mash":              "/path/to/envs/marta_ms/bin/mash",
    "metamlst":          "/path/to/envs/metamlst/bin/metamlst.py",
    "metamlst-merge":    "/path/to/envs/metamlst/bin/metamlst-merge.py",
    "metamlst-python":   "/path/to/envs/metamlst/bin/python",
    "metamlst-bowtie2":  "/path/to/envs/metamlst/bin/bowtie2",
    "metamlst-samtools": "/path/to/envs/metamlst/bin/samtools",
    "resfinder":         "/path/to/resfinder/run_resfinder.py",
    "resfinder-python":  "/path/to/envs/resfinder/bin/python",
    "resfinder-kma":     "/path/to/envs/resfinder/bin/kma",
    "mlst":              "/path/to/mlst/bin/mlst",
    "strainphlan":       "/path/to/envs/metaphlan/bin/strainphlan",
}

# ---------------------------------------------------------------------------
# Reference databases
# ---------------------------------------------------------------------------
RESFINDER_DB   = Path("/path/to/resfinder/db_resfinder")
POINTFINDER_DB = Path("/path/to/resfinder/db_pointfinder")

STRAINPHLAN_DB_MARKERS      = Path("/path/to/strainphlan/db_markers/t__SGB10068.fna")
STRAINPHLAN_REF_GENOMES_DIR = Path("/path/to/strainphlan/reference_genomes")
STRAINPHLAN_DEFAULT_CLADE   = "t__SGB10068"

METAMLST_DB        = Path("/path/to/metamlst/metamlstDB_2022.db")
METAMLST_BT2_INDEX = Path("/path/to/metamlst/bowtie_index")

# ---------------------------------------------------------------------------
# Runtime settings
# ---------------------------------------------------------------------------
MAX_WORKERS = 10

# Metagenome downsampling based on coverage_factor:
#   if num_reads > THRESHOLD * coverage_factor  ->  downsample to TARGET * coverage_factor
# Set either to None to disable metagenome downsampling.
BASE_READS_PER_COVERAGE_THRESHOLD = 70_000
BASE_READS_PER_COVERAGE_TARGET    = 40_000

# Downsample reads before assembly if exceeding this threshold (read pairs).
# Set to None to disable downsampling.
MAX_READS_FOR_ASSEMBLY = 1_000_000

# ---------------------------------------------------------------------------
# SRA submission (fill in only when uploading)
# ---------------------------------------------------------------------------
SRA_BIOPROJECT    = None                       # e.g. "PRJNA123456"
SRA_FTP_HOST      = "ftp-private.ncbi.nlm.nih.gov"
SRA_FTP_SUBFOLDER = "uploads"                  # personal subfolder on SRA FTP
SRA_FTP_USERNAME  = None                       # NCBI account email
SRA_FTP_PASSWORD  = None                       # from NCBI submission portal

# ---------------------------------------------------------------------------
# Run selection
# ---------------------------------------------------------------------------
# None = process/list all runs; otherwise a list of run names, e.g.
# ['U32A_hs2', 'U32B_hs2'].
RUN_NAMES_TO_PROCESS = None
RUN_NAMES_TO_LIST    = None
