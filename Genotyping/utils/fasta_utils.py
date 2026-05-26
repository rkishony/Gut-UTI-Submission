import numpy as np
from pathlib import Path

from utils.file_utils import decompress_if_needed


def fastq_read_length_stats(fastq_path: Path, num_reads: int | None = None) -> dict:
    """Compute read length statistics from a fastq file (supports .gz)."""
    lengths = []
    with decompress_if_needed(fastq_path) as fp:
        with open(fp, 'r') as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    lengths.append(len(line.rstrip('\n')))
                if num_reads is not None and len(lengths) >= num_reads:
                    break
    arr = np.array(lengths, dtype=int)
    return {
        'n_reads': len(arr),
        'mean_len': round(float(arr.mean()), 1),
        'median_len': float(np.median(arr)),
        'min_len': int(arr.min()),
        'max_len': int(arr.max()),
    }
