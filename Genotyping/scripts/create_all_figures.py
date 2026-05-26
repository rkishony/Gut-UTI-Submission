import importlib
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPTS_DIR.parent
sys.path.insert(0, str(SCRIPTS_DIR))

from record_value import get_records, reset

FIGURE_SCRIPTS = [
    "fig_metagenome_depth",
    "fig_patient_trees",
    "fig_patient_trees_combined",
    "fig_uti_distance_to_closest_fecal",
    "fig_mlst_pie_lines",
    "fig_strainphlan_tree",
]


def main():
    reset()
    for script in FIGURE_SCRIPTS:
        print(f"\n{'='*60}\nRunning {script}\n{'='*60}")
        try:
            mod = importlib.import_module(script)
            mod.main()
        except Exception as e:
            print(f"WARNING: {script} failed: {e}")

    records = get_records()
    out = PROJECT_ROOT / "figures" / "record_values_python.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    records.to_csv(out, index=False)
    print(f"\nRecorded {len(records)} values -> {out}")


if __name__ == "__main__":
    main()
