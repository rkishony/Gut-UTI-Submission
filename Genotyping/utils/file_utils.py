from pathlib import Path
from contextlib import contextmanager
import shutil
import subprocess
import tempfile

from tools.run_utils import run_cmd


@contextmanager
def decompress_if_needed(input_file: Path):
    """
    Context manager to decompress a fastq.gz file if needed.

    Args:
        input_file: Input file (may be .gz)

    Yields:
        Path to uncompressed file (temp file that will be cleaned up)
    """
    if input_file.suffix != '.gz':
        yield input_file
        return

    temp_file = Path(tempfile.mktemp(suffix='.fastq'))
    try:
        with open(temp_file, 'wb') as fh:
            subprocess.run(['gunzip', '-c', str(input_file)], stdout=fh, check=True)
        yield temp_file
    finally:
        temp_file.unlink(missing_ok=True)


def create_if_missing(path: Path | str, clean: bool = False) -> Path:
    """
    Create a directory if it doesn't exist.

    Args:
        path: Path to directory
    """
    if isinstance(path, str):
        path = Path(path)
    path = path.resolve()
    if clean and path.exists():
        shutil.rmtree(path, ignore_errors=True)
    path.mkdir(parents=True, exist_ok=True)
    return path


def count_reads_fastq(path: Path) -> int:
    """Count reads in a FASTQ file (4 lines per read, supports .gz)."""
    return count_lines(path) // 4


def count_lines(filepath: Path) -> int:
    """Count number of lines in a file using wc -l (supports .gz)."""
    with decompress_if_needed(filepath) as fp:
        result = run_cmd(['wc', '-l', str(fp)])
        return int(result.stdout.strip().split()[0])
