"""Shared utilities for running external tools."""

import subprocess
import sys
import threading
from pathlib import Path
from typing import Union

from config import TOOLS

_thread_local = threading.local()


def set_run_stderr_path(path: Path) -> None:
    """Set per-thread stderr log path. File is only created if stderr is non-empty."""
    _thread_local.stderr_path = path
    _thread_local.stderr_buffer = []


def close_run_stderr_log() -> None:
    """Write collected stderr to file (if any) and reset."""
    buf = getattr(_thread_local, "stderr_buffer", [])
    path = getattr(_thread_local, "stderr_path", None)
    if buf and path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("".join(buf))
    _thread_local.stderr_buffer = []
    _thread_local.stderr_path = None


def _append_stderr(line: str) -> None:
    buf = getattr(_thread_local, "stderr_buffer", None)
    if buf is not None:
        buf.append(line)


def get_tool(name: str) -> Path:
    """Look up an external tool path by name."""
    if name not in TOOLS:
        raise KeyError(f"Unknown tool {name!r}. Add it to TOOLS in config.py.")
    return Path(TOOLS[name])


def run_cmd(cmd: Union[list, str], check: bool = True,
            stream: bool = False, prefix: str = "", **kwargs) -> subprocess.CompletedProcess:
    """Run a subprocess command.

    All list elements are stringified.  By default output is captured.

    Args:
        cmd: Command and arguments (list), or a shell string when shell=True.
        check: If True (default), raise RuntimeError on non-zero exit.
        stream: If True, stream stderr to console in real-time (for long-running commands).
        prefix: Optional prefix for streamed lines (e.g. "[strainphlan] ").
    """
    if isinstance(cmd, (list, tuple)):
        cmd = [str(c) for c in cmd]

    if stream:
        kwargs.pop("capture_output", None)
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, **kwargs)
        stderr_lines = []
        stdout_lines = []
        for line in proc.stderr:
            stderr_lines.append(line)
            print(f"{prefix}{line}", end="", file=sys.stderr, flush=True)
            _append_stderr(f"{prefix}{line}")
        stdout_data = proc.stdout.read()
        stdout_lines.append(stdout_data)
        proc.wait()
        result = subprocess.CompletedProcess(
            cmd, proc.returncode,
            stdout="".join(stdout_lines),
            stderr="".join(stderr_lines))
    else:
        result = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
        if result.stderr:
            for line in result.stderr.splitlines(keepends=True):
                _append_stderr(line)

    if check and result.returncode != 0:
        parts = filter(None, [result.stderr.strip(), result.stdout.strip()])
        msg = "\n".join(parts)
        raise RuntimeError(
            f"Command failed ({cmd[0] if isinstance(cmd, list) else cmd}): "
            f"returncode={result.returncode} {msg}"
        )
    return result
