import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


def pytest_configure(config):
    config.addinivalue_line("markers", "live: hits public web APIs (deselect with -m 'not live')")


@pytest.fixture(autouse=True)
def _isolated_cache(tmp_path, monkeypatch):
    """Each test gets its own cache dir unless it asks for the real one."""
    monkeypatch.setenv("ZEBRA_CACHE_DIR", str(tmp_path / "cache"))
    yield
