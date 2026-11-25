import importlib.util

import pytest


@pytest.fixture
def no_pyarrow(monkeypatch):
    """
    Make pyarrow look uninstalled to delineator's parquet check.

    pyarrow is the optional extra delineator[parquet], so the check for it
    (settings._require_pyarrow) must be tested whether or not it is installed.
    """
    real_find_spec = importlib.util.find_spec

    def find_spec(name, *args, **kwargs):
        return None if name == "pyarrow" else real_find_spec(name, *args, **kwargs)

    monkeypatch.setattr(importlib.util, "find_spec", find_spec)
