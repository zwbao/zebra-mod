import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


def pytest_configure(config):
    config.addinivalue_line("markers", "live: hits public web APIs; skipped unless ZEBRA_NETWORK_TESTS=1 "
                                       "or the run selects them with -m live")


def pytest_collection_modifyitems(config, items):
    """The default run is offline: a live test waits on a public API (Ensembl VEP, gnomAD …) and a slow
    or unreachable one hung the whole suite with no output. Opt in with ZEBRA_NETWORK_TESTS=1."""
    if os.environ.get("ZEBRA_NETWORK_TESTS") == "1" or "live" in (config.getoption("-m") or "").replace("not live", ""):
        return
    skip = pytest.mark.skip(reason="hits the network: set ZEBRA_NETWORK_TESTS=1 (or -m live) to run")
    for item in items:
        if item.get_closest_marker("live"):
            item.add_marker(skip)


@pytest.fixture(autouse=True)
def _isolated_cache(tmp_path, monkeypatch):
    """Each test gets its own cache dir unless it asks for the real one."""
    monkeypatch.setenv("ZEBRA_CACHE_DIR", str(tmp_path / "cache"))
    yield
