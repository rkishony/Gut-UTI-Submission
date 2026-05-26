import re

# Color tags used in log messages: {cyan}, {magenta}, {bright_cyan}, {reset}, etc.
_ANSI_CODES = {
    "reset":         "0",
    "bold":          "1",
    "black":         "30",
    "red":           "31",
    "green":         "32",
    "yellow":        "33",
    "blue":          "34",
    "magenta":       "35",
    "cyan":          "36",
    "white":         "37",
    "bright_black":  "90",
    "bright_red":    "91",
    "bright_green":  "92",
    "bright_yellow": "93",
    "bright_blue":   "94",
    "bright_magenta":"95",
    "bright_cyan":   "96",
    "bright_white":  "97",
}

_TAG_RE = re.compile(r"\{(" + "|".join(_ANSI_CODES) + r")\}")


def _tags_to_ansi(msg: str) -> str:
    return _TAG_RE.sub(lambda m: f"\033[{_ANSI_CODES[m.group(1)]}m", msg)


def strip_tags(msg: str) -> str:
    return _TAG_RE.sub("", msg)


def print_color(*args, **kwargs):
    """Print with ANSI color tag support."""
    msg = " ".join(str(a) for a in args)
    print(_tags_to_ansi(msg), flush=True, **kwargs)
