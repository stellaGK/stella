################################################################################
#              Reverse the magnetic field of a GVEC state and parameter file   #
################################################################################
# B = (Phi'/Jac) (...) in GVEC, so reversing the sign of the toroidal flux Phi
# (and with it the poloidal flux chi = iota Phi) reverses B, while the geometry
# of the flux surfaces, iota and the pressure stay the same.
#
# The python interface of GVEC (used by gvec_to_stella.py) takes the toroidal
# flux from the parameter file (phiedge), while its Fortran interface (used by
# stella's direct path) takes the profiles phi and chi stored in the state file,
# so both are reversed, which keeps the two files consistent.
################################################################################

import re


def reverse_magnetic_field_of_parameter_file(parameter_file, new_parameter_file):
    '''Write <parameter_file> with the opposite sign of phiedge.'''
    with open(parameter_file) as file:
        text = file.read()
    def reverse(match):
        value = match.group(2)
        return match.group(1) + (value[1:] if value.startswith('-') else '-' + value)
    text, n = re.subn(r'(?im)^(\s*phiedge\s*=\s*)([-+0-9.eEdD]+)', reverse, text)
    if n != 1:
        raise ValueError(f'Found {n} lines with phiedge in {parameter_file}, expected one.')
    with open(new_parameter_file, 'w') as file:
        file.write(text)


def reverse_magnetic_field_of_state_file(state_file, new_state_file):
    '''Write <state_file> with the opposite sign of the profiles phi and chi.'''
    with open(state_file) as file:
        lines = file.read().split('\n')
    start = next(i for i, line in enumerate(lines) if 'spos,phi,chi,iota,pressure' in line)
    end = next(i for i in range(start + 1, len(lines)) if lines[i].startswith('##'))
    for i in range(start + 1, end):
        values = [value.strip() for value in lines[i].split(',')]
        for column in (1, 2):
            values[column] = values[column][1:] if values[column].startswith('-') else '-' + values[column]
        lines[i] = '  ' + ',  '.join(values)
    with open(new_state_file, 'w') as file:
        file.write('\n'.join(lines))
