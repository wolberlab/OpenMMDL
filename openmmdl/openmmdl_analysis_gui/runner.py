"""Run ``openmmdl analysis`` as a subprocess and collect its console output."""

from __future__ import annotations

import json
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
import uuid
from pathlib import Path

from openmmdl.openmmdl_analysis_gui.options import build_argv, format_command

RUN_FILE = "analysis_gui_run.json"
MAX_LINES = 5000

_ANSI_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
_SPLIT_RE = re.compile(r"(\r\n|\n|\r)")
# A tqdm bar such as "Analyzing frames:  3%|3   | 3/99 [00:26<10:19,  6.45s/it]" plus what follows it.
_BAR_RE = re.compile(r"^(.*?\d+%\|[^|]*\|\s*\d+/\d+\s*\[[^\]]*\])(.*)$")


class AnalysisRun:
    """One analysis subprocess. Its log keeps finished lines plus the current tqdm bar."""

    def __init__(self, values, output_dir, copy_topology=True):
        self.id = uuid.uuid4().hex[:12]
        self.output_dir = Path(output_dir).resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.values = dict(values)
        if copy_topology and self.values.get("topology"):
            # The analysis rewrites its topology in place (water names, renumbering).
            source = Path(self.values["topology"])
            inputs = self.output_dir / "inputs"
            inputs.mkdir(exist_ok=True)
            target = inputs / source.name
            if source.resolve() != target.resolve():
                shutil.copy2(source, target)
            self.values["topology"] = str(target)
        self.argv = build_argv(self.values)
        self.command = format_command(self.argv)
        self.lines = []
        self.times = []
        self.dropped = 0
        self.progress = ""
        self.status = "running"
        self.returncode = None
        self.started = time.time()
        self.ended = None
        self._cancelled = False
        self._lock = threading.Lock()
        self._write_run_file()

        kwargs = {}
        if os.name == "nt":
            kwargs["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP
        else:
            kwargs["start_new_session"] = True
        # The analysis filters warnings only in its main process; its worker processes would
        # repeat library UserWarnings once per frame.
        env = dict(
            os.environ,
            PYTHONUNBUFFERED="1",
            PYTHONIOENCODING="utf-8",
            PYTHONWARNINGS="ignore",
            MPLBACKEND="Agg",
        )
        self._proc = subprocess.Popen(
            [sys.executable, "-m", "openmmdl.openmmdl_analysis.cli", *self.argv],
            cwd=self.output_dir,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env=env,
            **kwargs,
        )
        threading.Thread(target=self._pump, daemon=True).start()

    def _write_run_file(self):
        record = {
            "command": self.command,
            "argv": self.argv,
            "values": self.values,
            "started": self.started,
            "ended": self.ended,
            "status": self.status,
            "returncode": self.returncode,
        }
        (self.output_dir / RUN_FILE).write_text(json.dumps(record, indent=2), encoding="utf-8")

    def _add_line(self, text):
        text = _ANSI_RE.sub("", text).rstrip()
        if not text.strip():
            return
        self.lines.append(text)
        self.times.append(round(time.time() - self.started, 1))
        if len(self.lines) > MAX_LINES:
            del self.lines[0]
            del self.times[0]
            self.dropped += 1

    def _consume(self, buffer):
        parts = _SPLIT_RE.split(buffer)
        rest = parts.pop()
        with self._lock:
            for text, sep in zip(parts[::2], parts[1::2]):
                text = _ANSI_RE.sub("", text)
                # Worker processes can write right after a progress bar, on the same line.
                bar = _BAR_RE.match(text)
                if bar:
                    self.progress = bar.group(1).strip()
                    text = bar.group(2)
                if sep == "\r":
                    if text.strip():
                        self.progress = text.strip()
                elif text.strip():
                    self._add_line(text)
                    if not bar:
                        self.progress = ""
        return rest

    def _pump(self):
        buffer = ""
        stream = self._proc.stdout
        while True:
            chunk = stream.read1(4096)
            if not chunk:
                break
            buffer = self._consume(buffer + chunk.decode("utf-8", errors="replace"))
        returncode = self._proc.wait()
        with self._lock:
            self._add_line(buffer)
            self.progress = ""
            self.returncode = returncode
            self.ended = time.time()
            if self._cancelled:
                self.status = "cancelled"
            else:
                self.status = "finished" if returncode == 0 else "failed"
        self._write_run_file()

    def cancel(self):
        if self._proc.poll() is not None:
            return
        self._cancelled = True
        if os.name == "nt":
            subprocess.run(
                ["taskkill", "/F", "/T", "/PID", str(self._proc.pid)], capture_output=True
            )
        else:
            os.killpg(self._proc.pid, signal.SIGTERM)

    def snapshot(self, since=0):
        with self._lock:
            start = max(since - self.dropped, 0)
            return {
                "id": self.id,
                "status": self.status,
                "returncode": self.returncode,
                "command": self.command,
                "outputDir": str(self.output_dir),
                "elapsed": (self.ended or time.time()) - self.started,
                "values": self.values,
                "lines": self.lines[start:],
                "times": self.times[start:],
                "next": self.dropped + len(self.lines),
                "progress": self.progress,
            }
