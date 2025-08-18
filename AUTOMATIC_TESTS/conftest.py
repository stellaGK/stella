
# Allow us to pass the python version to the python tests
def pytest_addoption(parser):
    parser.addoption("--stella_version", action="store", default="master")

