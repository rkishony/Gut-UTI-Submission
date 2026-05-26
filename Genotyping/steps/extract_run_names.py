#!/usr/bin/env python3
"""Step 1: Extract and parse run names from raw sequencing folders."""

from pathlib import Path
import pandas as pd
import sys
import re

import config
from config import SEQ_BATCHES_TO_PATHS, UTI_ASSIGNMENT_CSV, META_FOBTS_CSV
from config_utils import get_output_dir
from steps.step import Step

RUN_NAME_PATTERN = re.compile(r'^([^_]+)_')
FILENAME = "sequencing_runs.csv"

FILTER_OUT_STARTS_WITH = ("Sample_Einat", "morb")


def _get_run_name(row: pd.Series) -> str:
    is_isolate = row['is_isolate']
    is_faecal = row['is_faecal']
    s = 'F' if is_faecal else 'U'
    if is_faecal:
        s += f"{row['fid']}"
    else:
        s += f"{row['uid']}"

    if is_isolate:
        assert row['isolate']
        s += f"{row['isolate']}"
    else:
        s += "meta"

    if not pd.isna(row['replicate']):
        s += f"_{row['replicate']}"

    s += f"_{row['batch']}"
    return s


def parse_run_name(run_name: str) -> dict:  # noqa: C901
    """Parse run name into structured components."""
    if run_name.startswith("Sample_"):
        run_name = run_name[7:]
    result = {
        'fid': None,
        'uid': None,
        'is_isolate': False,
        'is_faecal': False,
        'isolate': None,
        'replicate': None
    }

    # Faecal isolates: F29R12, F207B1, F12C1, etc.
    faecal_isolate_pattern = re.compile(r'^F(\d+)([A-Z]\d*)$')
    match = faecal_isolate_pattern.match(run_name)
    if match:
        result['fid'] = int(match.group(1))
        result['is_isolate'] = True
        result['is_faecal'] = True
        result['isolate'] = match.group(2)
        return result

    # Urinae isolates: U378A, U582A, etc.
    urinae_isolate_pattern = re.compile(r'^U(\d+)([A-Z])$')
    match = urinae_isolate_pattern.match(run_name)
    if match:
        result['uid'] = int(match.group(1))
        result['is_isolate'] = True
        result['is_faecal'] = False
        result['isolate'] = match.group(2)
        return result

    # Faecal scrapes: F120scrape, F277scrape, etc.

    faecal_scrape_pattern = re.compile(r'^F(\d+)scrape$')
    match = faecal_scrape_pattern.match(run_name)
    if match:
        result['fid'] = int(match.group(1))
        result['is_isolate'] = False
        result['is_faecal'] = True
        return result

    # Faecal scrapes with replicate: 1842spR, 183spR, 310dpR2, etc.
    scrape_replicate_pattern = re.compile(r'^(\d+)([a-zA-Z][a-zA-Z0-9]*)$')
    match = scrape_replicate_pattern.match(run_name)
    if match:
        result['fid'] = int(match.group(1))
        result['is_isolate'] = False
        result['is_faecal'] = True
        result['replicate'] = match.group(2)
        return result

    # Plain numbers: 97, 1524, etc. (assumed to be faecal scrapes)
    plain_number_pattern = re.compile(r'^\d+$')
    match = plain_number_pattern.match(run_name)
    if match:
        result['fid'] = int(run_name)
        result['is_isolate'] = False
        result['is_faecal'] = True
        return result

    return result


def extract_run_names() -> pd.DataFrame:
    """Extract run names from raw folders and return as DataFrame."""
    run_info: list[dict[str, str]] = []

    for batch, folder in SEQ_BATCHES_TO_PATHS.items():
        for item in folder.iterdir():
            if item.is_dir() and not item.name.startswith(FILTER_OUT_STARTS_WITH):
                r1_paths = list(item.glob("*_R1_*.fastq.gz"))
                assert len(r1_paths) == 1, f"Expected 1 R1 file in {item.name}, found {len(r1_paths)}"
                r1_path = r1_paths[0]
                r2_path = Path(str(r1_path).replace('_R1_', '_R2_'))
                assert r2_path.exists(), f"Missing R2 file for {r1_path.name}"
                raw_path = r1_path.name.replace('_R1_', '_R?_')
                run_info.append({
                    'raw_path': item.name + '/' + raw_path,
                    'raw_name': item.name, 'batch': batch, 'folder_name': batch,
                })

    df = pd.DataFrame(run_info)

    # Parse run names
    parsed_data = []
    for _, row in df.iterrows():
        parsed = parse_run_name(row['raw_name'])
        parsed['raw_name'] = row['raw_name']
        parsed['raw_path'] = row['raw_path']
        parsed['batch'] = row['batch']
        parsed_data.append(parsed)

    df_parsed = pd.DataFrame(parsed_data)
    df_parsed['fid'] = df_parsed['fid'].astype('Int64')
    df_parsed['uid'] = df_parsed['uid'].astype('Int64')

    # Fill in fid for UTI cases using the assignment CSV
    uti_map = pd.read_csv(UTI_ASSIGNMENT_CSV)
    uti_map = uti_map[uti_map['FOBT'] != '#N/A']
    uti_map = uti_map.dropna(subset=['Our UTI_ID', 'FOBT'])
    uid_to_fid = dict(zip(uti_map['Our UTI_ID'].astype(int), uti_map['FOBT'].astype(int)))
    assert all(df_parsed['uid'].isna() ^ df_parsed['fid'].isna())
    mask = df_parsed['uid'].notna()
    df_parsed.loc[mask, 'fid'] = df_parsed.loc[mask, 'uid'].map(uid_to_fid).astype('Int64')
    df_parsed['raw_path'] = df_parsed['raw_path'].astype(str)
    df_parsed['raw_name'] = df_parsed['raw_name'].astype(str)
    df_parsed['run_name'] = df_parsed.apply(lambda row: _get_run_name(row), axis=1)

    # Merge CFU and coverage_factor from meta_fobts for metagenome samples
    meta_fobts = pd.read_csv(META_FOBTS_CSV)
    meta_fobts = meta_fobts.rename(columns={'coverage factor': 'coverage_factor'})
    meta_fobts['FOBT'] = pd.to_numeric(meta_fobts['FOBT'], errors='coerce').astype('Int64')
    meta_fobts = meta_fobts.dropna(subset=['FOBT'])
    fid_to_cfu = dict(zip(meta_fobts['FOBT'], meta_fobts['CFU']))
    fid_to_cov = dict(zip(meta_fobts['FOBT'], meta_fobts['coverage_factor']))
    is_meta = (~df_parsed['is_isolate']) & df_parsed['is_faecal']
    df_parsed['CFU'] = pd.NA
    df_parsed['coverage_factor'] = pd.NA
    df_parsed.loc[is_meta, 'CFU'] = df_parsed.loc[is_meta, 'fid'].map(fid_to_cfu)
    df_parsed.loc[is_meta, 'coverage_factor'] = df_parsed.loc[is_meta, 'fid'].map(fid_to_cov)

    # assert that all run names are different:
    assert len(df_parsed['run_name'].unique()) == len(df_parsed), "Duplicate run names found"

    df_parsed.sort_values(
        by=['is_faecal', 'is_isolate', 'fid', 'uid', 'isolate', 'replicate', 'batch'],
        inplace=True,
    )

    return df_parsed


class ExtractRunNamesStep(Step):
    """Extract, parse, and compress run names from raw sequencing folders into a CSV."""

    @classmethod
    def output_path(cls) -> Path:
        return get_output_dir("summary") / FILENAME

    def _run(self) -> None:
        df = extract_run_names()
        if config.RUN_NAMES_TO_LIST is not None:
            df = df[df["run_name"].isin(config.RUN_NAMES_TO_LIST)]
        out = self.output_path()
        df.to_csv(out, index=False)
        return f"{len(df)} sequencing runs"


if __name__ == '__main__':
    force = "--force" in sys.argv
    ExtractRunNamesStep(force=force).run()
