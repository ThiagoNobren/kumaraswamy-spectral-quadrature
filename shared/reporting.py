from __future__ import annotations

import io
import sys
from datetime import datetime
from pathlib import Path
from typing import TextIO


class TeeStdout(TextIO):
    def __init__(self, original: TextIO, buffer: io.StringIO) -> None:
        self._original = original
        self._buffer = buffer

    def write(self, data: str) -> int:
        self._original.write(data)
        self._buffer.write(data)
        return len(data)

    def flush(self) -> None:
        self._original.flush()
        self._buffer.flush()

    def __getattr__(self, name: str):
        return getattr(self._original, name)


def write_terminal_md(body: str, report_dir: Path, prefix: str) -> Path:
    report_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = report_dir / f"{prefix}_{stamp}.md"
    path.write_text(
        f"# {prefix.replace('_', ' ')}\n\n```text\n{body}\n```\n",
        encoding="utf-8",
    )
    return path


def run_with_tee(main_fn, report_dir: Path, prefix: str) -> Path:
    buf = io.StringIO()
    old = sys.stdout
    sys.stdout = TeeStdout(old, buf)
    try:
        main_fn()
    finally:
        sys.stdout = old
    path = write_terminal_md(buf.getvalue(), report_dir, prefix)
    print(f"Relatorio Markdown gravado: {path}")
    return path
