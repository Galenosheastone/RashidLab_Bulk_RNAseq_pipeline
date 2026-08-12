"""Compatibility entry point for Streamlit deployments.

The repository-level ``streamlit_app.py`` is the canonical app. This wrapper
keeps older Streamlit Cloud configurations that point at
``deseq2_enrich/app/streamlit_app.py`` on the same UI code path.
"""
from __future__ import annotations

from pathlib import Path
import runpy
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

runpy.run_path(str(ROOT / "streamlit_app.py"), run_name="__main__")
