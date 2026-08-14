# Tests/test_imports.py

"""
Import smoke tests.

Every executable module in the repo must at minimum parse and import. This
guards against the class of breakage where source is pasted back from a
chat transcript carrying markdown fences, smart quotes, or flattened
indentation - all of which make a file look fine but fail to compile.
"""

import ast
import pathlib
import importlib

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent

PY_FILES = sorted(
    p for p in REPO_ROOT.rglob("*.py")
    if ".git" not in p.parts
)

# Modules importable without optional/heavy third-party extras.
IMPORTABLE_MODULES = [
    "CEED_universal_model",
    "simulation.convergence_model",
    "simulation.minimum_esm_code",
    "simulation.mhd_spatial_model",
    "simulation.unit_bridge",
    "simulation.hindcast",
    "experiments.run_mc",
    "Data.inputs",
]

# Curly double/single quotes, spelled as escapes so this file does not trip
# its own check below.
SMART_CHARS = "".join(chr(c) for c in (0x201C, 0x201D, 0x2018, 0x2019))
FENCE = "`" * 3


@pytest.mark.parametrize("path", PY_FILES, ids=lambda p: str(p.relative_to(REPO_ROOT)))
def test_file_parses(path):
    """Every .py file in the repo is syntactically valid Python."""
    source = path.read_text(encoding="utf-8")
    ast.parse(source, filename=str(path))


@pytest.mark.parametrize("path", PY_FILES, ids=lambda p: str(p.relative_to(REPO_ROOT)))
def test_no_markdown_artifacts(path):
    """No stray markdown fences or curly quotes in Python source."""
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        assert FENCE not in line, f"{path}:{lineno} contains a markdown code fence"
        for ch in SMART_CHARS:
            assert ch not in line, f"{path}:{lineno} contains smart quote {ch!r}"


@pytest.mark.parametrize("module", IMPORTABLE_MODULES)
def test_module_imports(module):
    """Core modules import cleanly."""
    importlib.import_module(module)
