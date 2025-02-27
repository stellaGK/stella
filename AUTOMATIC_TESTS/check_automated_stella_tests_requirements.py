# This script just checks if we can import all the modules listed in
# requirements.txt

import importlib
import pathlib
import textwrap

requirements_file = pathlib.Path(__file__).parent / "requirements.txt"

with open(requirements_file, "r") as f:
    contents = f.readlines()

requirements = []
for line in contents:
    if line.startswith("#"):
        continue
    if line.replace(' ','').replace('\n','').replace('\t','')=='':
        continue
    module = line.split("=")[0]
    if 'pyrokinetics' in line: module = 'pyrokinetics'
    requirements.append(module.replace(">", "").replace("~", "").replace("\n", ""))

failed_modules = []
for requirement in requirements:
    try:
        importlib.import_module(requirement)
    except ImportError:
        failed_modules.append(requirement)

if failed_modules: 
    list_of_failed_modules = ", ".join(failed_modules)
    print(
        textwrap.dedent(
            f"""\
            Could not import: {list_of_failed_modules}
            Run:
                make create-test-virtualenv
                source AUTOMATIC_TESTS/venv/bin/activate
                pip install EXTERNALS/pyrokinetics/.
            to install and activate test dependencies
            """
        )
    )
    exit(1)
