"""Lightweight value recorder — Python equivalent of the MATLAB record_value."""

import pandas as pd

_records: list[dict] = []


def record_value(description: str, label: str, value, fmt: str | None = None):
    if isinstance(value, str):
        svalue = value
        num = float("nan")
    else:
        num = float(value)
        if fmt:
            svalue = fmt % value if "%" in fmt else f"{value:{fmt}}"
        elif num == int(num) and not (abs(num) == float("inf")):
            svalue = str(int(num))
        else:
            svalue = f"{num:.2g}"
        svalue = svalue.replace("e+0", "e+").replace("e-0", "e-")
    _records.append(
        dict(description=description, label=label, svalue=svalue, value=num)
    )
    print(f"*** {description:<35s} {label:>45s}: {svalue:>7s}")


def reset():
    _records.clear()


def get_records() -> pd.DataFrame:
    return pd.DataFrame(_records, columns=["description", "label", "svalue", "value"])
