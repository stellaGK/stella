import pytest

# Allow us to pass the python version to the python tests
def pytest_addoption(parser):
    parser.addoption("--stella_version", action="store", default="master")
    parser.addoption("--include-slow-tests", action="store_true", default=False,
        help="Also run the tests marked as slow, e.g. the KBM and TAE benchmarks.")

# Register the 'slow' marker, used for long physics benchmarks
def pytest_configure(config):
    config.addinivalue_line("markers", "slow: slow test, only run with --include-slow-tests (make numerical-tests-slow)")

# Skip the slow tests, unless --include-slow-tests is given
def pytest_collection_modifyitems(config, items):
    if config.getoption("--include-slow-tests"): return
    skip_slow = pytest.mark.skip(reason="Slow test, run it with: make numerical-tests-slow")
    for item in items:
        if "slow" in item.keywords: item.add_marker(skip_slow)
