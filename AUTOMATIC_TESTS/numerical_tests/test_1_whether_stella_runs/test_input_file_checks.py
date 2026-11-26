################################################################################
#                    Check the namelists of the input file                     #
################################################################################
# stella silently skips the namelists it does not read, so init_stella.f90 stops
# when the input file contains a namelist which is not in its list of namelists
# (check_for_unknown_namelists). Test that this list contains exactly the
# namelists which the code reads, and that stella stops on a misspelt namelist.
################################################################################

# Python modules
import pytest
import os, re
import pathlib
import subprocess

# Package to run stella
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# Global variables
input_filename = 'miller_nonlinear_CBC.in'
stella_code_directory = pathlib.Path(__file__).parent.parent.parent.parent / 'STELLA_CODE'

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#        Check whether the list of namelists matches the stella code           #
#-------------------------------------------------------------------------------
def test_whether_the_list_of_known_namelists_matches_the_code():
    '''The namelists in <known_namelists> of check_for_unknown_namelists() must be the namelists
    that the code looks up with input_unit_exist(), apart from the indexed &species_parameters_<i>.'''

    # Read <known_namelists> from init_stella.f90
    text = (stella_code_directory / 'init_stella.f90').read_text()
    block = re.search(r'known_namelists = \[character\(len=\d+\) ::(.*?)\]', text, re.S).group(1)
    known_namelists = set(re.findall(r"'([a-z0-9_]+)'", block))

    # Gather the namelists which the code looks up in the input file
    namelists_in_code = set()
    for path in list(stella_code_directory.rglob('*.f90')) + list(stella_code_directory.rglob('*.fpp')):
        namelists_in_code.update(re.findall(r"input_unit_exist\(\s*['\"]([a-z0-9_]+)['\"]", path.read_text(errors='ignore')))

    missing = sorted(namelists_in_code - known_namelists)
    extra = sorted(known_namelists - namelists_in_code)
    assert not missing, f'Add these namelists to <known_namelists> in init_stella.f90: {missing}'
    assert not extra, f'Remove these namelists from <known_namelists> in init_stella.f90: {extra}'
    print(f'  -->  The {len(known_namelists)} known namelists match the namelists read by the code.')
    return

#-------------------------------------------------------------------------------
#               Check whether stella stops on a misspelt namelist              #
#-------------------------------------------------------------------------------
@pytest.mark.parametrize('wrong_line', ['&initialize_distribution', '&Initialise_distribution', '&z_grid nzed = 24 /'])
def test_whether_stella_stops_on_an_unknown_namelist(tmp_path, stella_version, wrong_line):
    '''Replace the first line of a namelist by a misspelt name, a name with capitals, or a
    name followed by values, which stella would skip, and check that stella stops.'''
    if stella_version!='master': pytest.skip('Older stella versions do not check the namelists.')

    # Write an input file with the wrong namelist line instead of &initialise_distribution
    path_input_file = copy_input_file(input_filename, tmp_path)
    text = path_input_file.read_text()
    assert '\n&initialise_distribution\n' in text
    path_input_file.write_text(text.replace('\n&initialise_distribution\n', f'\n{wrong_line}\n'))

    # stella should stop before the simulation starts, and name the wrong line
    os.chdir(tmp_path)
    result = subprocess.run(['mpirun', '--oversubscribe', '-np', '2', get_stella_path(stella_version), input_filename],
        capture_output=True, text=True)
    output = result.stdout + result.stderr
    assert result.returncode != 0, f'stella did not stop on the namelist line "{wrong_line}".'
    assert f'The input file contains the line "{wrong_line}"' in output, f'stella did not name the wrong line:\n{output[-2000:]}'
    assert not (tmp_path / input_filename.replace('.in', '.out.nc')).exists(), 'stella ran the simulation.'
    print(f'  -->  stella stops on the namelist line "{wrong_line}".')
    return
