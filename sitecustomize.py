"""Deployment hygiene for Streamlit Cloud.

Python imports ``sitecustomize`` automatically during interpreter startup when
this file is on ``sys.path``. Keep committed or persisted bytecode from taking
precedence over freshly pulled source files.
"""
from __future__ import annotations

import os
import shutil
import sys


sys.dont_write_bytecode = True

ROOT = os.path.dirname(os.path.abspath(__file__))
for dirpath, dirnames, _ in os.walk(ROOT):
    if "__pycache__" not in dirnames:
        continue
    shutil.rmtree(os.path.join(dirpath, "__pycache__"), ignore_errors=True)
    dirnames.remove("__pycache__")
