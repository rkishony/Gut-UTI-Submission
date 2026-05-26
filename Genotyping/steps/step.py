"""Base classes for pipeline steps."""

import shutil
import sys
from abc import ABC, abstractmethod
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
import threading
import time
import traceback

import pandas as pd

import config
from config_utils import get_output_dir
from utils.file_utils import create_if_missing
from utils.tee_output import print_color, strip_tags

DEFAULT_TOTAL_CORES = 8

STEP_COLOR = "{cyan}"
RUN_COLOR = "{magenta}"


def _iter_paths(output: Path | Iterable[Path] | dict[str, Path] | None):
    """Yield paths from a single Path, iterable of Paths, dict of Paths, or None."""
    if output is None:
        return
    if isinstance(output, dict):
        yield from output.values()
    elif isinstance(output, Path):
        yield output
    else:
        yield from output


class Step(ABC):
    """
    Base class for all pipeline steps.

    Subclasses must implement:
        _run()  - do the work
    Subclasses may override:
        output_path()  - classmethod returning Path(s); used for skip logic
                         & by downstream steps to locate outputs
    """

    _name: str | None = None
    allow_failure: bool = False

    class _NameDescriptor:
        def __get__(self, obj, cls):
            return cls._name or cls.__name__.removesuffix("Step")

    name = _NameDescriptor()

    def __init__(self, force: bool = False,
                 total_cores: int = DEFAULT_TOTAL_CORES,
                 dry_run: bool = False,
                 **_kwargs):
        self.force = force
        self.total_cores = total_cores
        self.dry_run = dry_run
        self._log_file: Path | None = None

    @classmethod
    def output_path(cls) -> Path | list[Path] | None:
        """Return output path(s) for this step (used for skip logic).

        Other steps can call e.g. ``TrimStatsStep.output_path()``
        to locate upstream outputs without instantiating the step.
        """
        return None

    @classmethod
    def _iter_output_paths(cls):
        yield from _iter_paths(cls.output_path())

    def _output_exists(self) -> bool:
        paths = list(self._iter_output_paths())
        return bool(paths) and all(p.exists() for p in paths)

    def print(self, *args, time_stamp: bool = True, **kwargs) -> None:
        """Print with ``[StepName]`` prefix and write to step log."""
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if time_stamp:
            print_color(f"{STEP_COLOR}[{self.name}]{{reset}}", ts, *args, **kwargs)
        else:
            print_color(f"{STEP_COLOR}[{self.name}]{{reset}}", *args, **kwargs)
        plain = " ".join(str(a) for a in args)
        self._log(f"{ts}  {strip_tags(plain)}")

    def _log_path(self) -> Path:
        if self._log_file is None:
            stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            self._log_file = get_output_dir("logs") / f"{self.name}_{stamp}.log"
        return self._log_file

    def _log(self, msg: str) -> None:
        path = self._log_path()
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "a") as f:
            f.write(f"{msg}\n")

    @abstractmethod
    def _run(self) -> str | None:
        ...

    def run(self) -> None:
        if not self.force and self._output_exists():
            self.print("Skipped (output exists). Use force=True to rerun.")
            return
        for p in self._iter_output_paths():
            p.parent.mkdir(parents=True, exist_ok=True)
        self.print("STARTED")
        t0 = time.time()
        try:
            summary = self._run()
        except Exception:
            elapsed = time.time() - t0
            if not self.allow_failure:
                raise
            self.print(f"FAILED ({sys.exc_info()[1]}, {elapsed:.1f}s), continuing.\n{traceback.format_exc()}")
        else:
            elapsed = time.time() - t0
            msg = f"FINISHED ({elapsed:.1f}s)"
            if summary:
                msg += f" {summary}"
            self.print(msg)
            for p in self._iter_output_paths():
                if p.exists():
                    self.print(f"  → {p}")


class MultiStep(Step):
    """
    Step that iterates over rows of a CSV and processes each one.

    Subclasses must implement:
        process_run(row)   – do the work for one row
        run_folder(row)    – classmethod returning the absolute run directory
        output_paths       – class var with relative path(s) within run_folder
    Subclasses may override:
        setup()            – called once before the loop
        cores_per_run       – cores each worker needs (default 1);
                              max_workers = total_cores // cores_per_run
    """

    cores_per_run: int = 1
    clean_run_folder: bool = True
    sort_by: str | list[str] | None = None  # None do not sort
    ascending: bool | list[bool] | None = None  # None use default ascending for each column

    def __init__(self, csv_file: str | Path, **kwargs):
        super().__init__(**kwargs)
        self.csv_file = csv_file
        self.threads_per_worker = self.cores_per_run
        self.max_workers = max(1, self.total_cores // self.cores_per_run)
        self._lock = threading.Lock()

    @classmethod
    def run_folder(cls, row: pd.Series) -> Path:
        """Absolute path to the run's working directory."""
        raise NotImplementedError(f"{cls.__name__} must implement run_folder()")

    output_paths: str | list[str] | dict[str, str]

    @classmethod
    def run_output_path(cls, row: pd.Series) -> Path | Iterable[Path] | dict[str, Path]:
        """Absolute output path(s): run_folder(row) / output_paths."""
        folder = cls.run_folder(row)
        rel = cls.output_paths
        if isinstance(rel, dict):
            return {k: folder / v for k, v in rel.items()}
        if isinstance(rel, str):
            return folder / rel
        return [folder / p for p in rel]

    @classmethod
    def _iter_run_output_paths(cls, row: pd.Series):
        yield from _iter_paths(cls.run_output_path(row))

    def run_input_paths(self, row: pd.Series) -> Path | Iterable[Path] | dict[str, Path] | None:
        """Absolute input path(s) for a run. Override to enable staleness check."""
        return None

    def _iter_run_input_paths(self, row: pd.Series):
        yield from _iter_paths(self.run_input_paths(row))

    def _outputs_newer_than_inputs(self, row: pd.Series) -> bool:
        """True if all outputs are at least as new as all inputs."""
        input_paths = list(self._iter_run_input_paths(row))
        if not input_paths:
            return True
        output_paths = list(self._iter_run_output_paths(row))
        if not output_paths:
            return True
        newest_input = max(p.stat().st_mtime for p in input_paths if p.exists())
        oldest_output = min(p.stat().st_mtime for p in output_paths if p.exists())
        return oldest_output >= newest_input

    @abstractmethod
    def process_run(self, row: pd.Series) -> None:
        """Process a single run."""
        ...

    def setup(self) -> None:
        """Called once before iterating over runs."""
        pass

    def cleanup(self) -> None:
        """Called once after iterating over runs."""
        pass

    def _run_output_paths_satisfied(self, row: pd.Series) -> bool:
        """Return True if output paths indicate completion. Override for custom logic (e.g. any-of)."""
        paths = list(self._iter_run_output_paths(row))
        return bool(paths) and all(p.exists() for p in paths)

    def _run_output_exists(self, row: pd.Series) -> bool:
        if self._started_and_not_gracefully_completed(row):
            return False
        failed_path = self._status_path(row, "FAILED")
        if failed_path.exists() and not self._was_killed(row):
            return True  # real error — skip, don't retry
        if not self._run_output_paths_satisfied(row):
            return False
        if not self._outputs_newer_than_inputs(row):
            self.print("outputs stale (inputs newer), re-running", row)
            return False
        return True

    def _should_process(self, row: pd.Series) -> bool:
        """Return True if this row should be processed (respects RUN_NAMES_TO_PROCESS)."""
        if config.RUN_NAMES_TO_PROCESS is not None and row['run_name'] not in config.RUN_NAMES_TO_PROCESS:
            return False
        return True

    def _output_exists(self) -> bool:
        """Disable step-level skip; per-run skip handles it."""
        return False

    def filter_df(self, df: pd.DataFrame) -> pd.DataFrame:
        """Override to filter the runs DataFrame before processing."""
        return df

    def sort_df(self, df: pd.DataFrame) -> pd.DataFrame:
        """Override to sort the runs DataFrame before processing."""
        if self.sort_by is not None:
            kwargs = {} if self.ascending is None else {"ascending": self.ascending}
            df = df.sort_values(self.sort_by, **kwargs)
        return df

    def _status_path(self, row: pd.Series, status: str) -> Path:
        return self.run_folder(row) / f".{status}.txt"

    def _write_status(self, row: pd.Series, status: str, detail: str = "",
                      t0: float | None = None) -> None:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        msg = status
        if t0 is not None:
            msg += f" (elapsed: {time.time() - t0:.1f}s)"
        if detail:
            msg += f"\n{detail}"
        self._status_path(row, status).write_text(f"{ts} {msg}\n")
        self.print(msg, row)

    def print(self, msg: str, row: pd.Series | None = None):
        if row is not None:
            msg = f"{RUN_COLOR}{row['run_name']:<16}{{reset}} {msg}"
        super().print(msg)

    _STATUS_FILES = {".STARTED.txt", ".DONE.txt", ".FAILED.txt", ".STDERR.txt"}

    def _clean_run_folder(self, row: pd.Series) -> None:
        """Remove non-status files/dirs from the run folder."""
        folder = self.run_folder(row)
        if not folder.exists():
            return
        for item in folder.iterdir():
            if item.name in self._STATUS_FILES:
                continue
            if item.is_dir():
                shutil.rmtree(item)
            else:
                item.unlink()

    def _clear_status(self, row: pd.Series) -> None:
        for name in self._STATUS_FILES:
            path = self.run_folder(row) / name
            if path.exists():
                path.unlink()

    _KILL_PATTERNS = ("signal", "killed", "terminated", "returncode=-", "keyboardinterrupt")

    def _was_killed(self, row: pd.Series) -> bool:
        """True if .FAILED.txt exists and indicates a kill/signal rather than a real error."""
        failed_path = self._status_path(row, "FAILED")
        if not failed_path.exists():
            return False
        content = failed_path.read_text().lower()
        return any(p in content for p in self._KILL_PATTERNS)

    def _started_and_not_gracefully_completed(self, row: pd.Series) -> bool:
        started = self._status_path(row, "STARTED").exists()
        completed = self._status_path(row, "DONE").exists() or self._status_path(row, "FAILED").exists()
        killed = self._was_killed(row)
        return started and (not completed or killed)

    def _run(self) -> None:
        df = pd.read_csv(self.csv_file)
        df = self.sort_df(self.filter_df(df))
        total = len(df)

        self.setup()

        rows_to_process = []
        for _, row in df.iterrows():
            if not self._should_process(row):
                continue
            if not self.force and self._run_output_exists(row):
                continue
            rows_to_process.append(row)

        def _process(row: pd.Series) -> None:
            from tools.run_utils import set_run_stderr_path, close_run_stderr_log
            create_if_missing(self.run_folder(row), clean=self.clean_run_folder)
            self._clear_status(row)
            self._write_status(row, "STARTED")
            set_run_stderr_path(self._status_path(row, "STDERR"))
            t0 = time.time()
            try:
                self.process_run(row)
                self._write_status(row, "DONE", t0=t0)
            except Exception:
                self._write_status(row, "FAILED", t0=t0, detail=traceback.format_exc())
            finally:
                close_run_stderr_log()

        self.print(f"Processing {len(rows_to_process)}/{total} runs "
                   f"(max_workers={self.max_workers}, "
                   f"cores_per_run={self.cores_per_run})")

        if self.dry_run:
            for row in rows_to_process:
                self.print(f"[DRY RUN] would process", row)
            return f"[DRY RUN] {len(rows_to_process)}/{total} runs"

        if self.max_workers <= 1:
            for row in rows_to_process:
                _process(row)
        else:
            with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
                futures = {pool.submit(_process, row): row for row in rows_to_process}
                for future in as_completed(futures):
                    future.result()

        self.cleanup()
        return f"{len(rows_to_process)}/{total} runs"


def run_cli(*step_classes):
    """Common CLI entry point: parse --force, --total-cores, optional csv_file, then run each step."""
    import inspect
    import signal
    signal.signal(signal.SIGHUP, signal.SIG_IGN)  # survive SSH/terminal disconnect
    force = "--force" in sys.argv
    dry_run = "--dry-run" in sys.argv
    total_cores = DEFAULT_TOTAL_CORES
    positional = []
    for arg in sys.argv[1:]:
        if arg in ("--force", "--dry-run"):
            continue
        if arg.startswith("--total-cores="):
            total_cores = int(arg.split("=", 1)[1])
            continue
        positional.append(arg)
    csv_file = positional[0] if positional else None
    for cls in step_classes:
        params = inspect.signature(cls).parameters
        kwargs = {"force": force, "total_cores": total_cores, "dry_run": dry_run}
        if "csv_file" in params:
            kwargs["csv_file"] = csv_file
        cls(**kwargs).run()
