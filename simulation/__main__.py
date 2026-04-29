"""Allow running as `python -m simulation.sweep_gelation` or `python simulation/sweep_gelation.py`."""

import sys
import os

# When running as a script (not via -m), ensure repo root is on sys.path.
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
