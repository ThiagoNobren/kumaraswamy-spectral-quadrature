from __future__ import annotations

import sys
from pathlib import Path

_PKG = Path(__file__).resolve().parent
_ROOT = _PKG.parent
_SHARED = _ROOT / "shared"
for p in (_SHARED, _PKG):
    s = str(p)
    if s not in sys.path:
        sys.path.insert(0, s)
