#!/usr/bin/env python3
################################################################################
#      Make the VMEC equilibria of test_2h_vmec_magnetic_field_sign            #
################################################################################
# Solves the four VMEC inputs of this folder with VMEC++ and writes
#   wout_tok_circ_Bplus.nc, wout_tok_circ_Bminus.nc,
#   wout_w7x_kjm_Bplus.nc,  wout_w7x_kjm_Bminus.nc
# in AUTOMATIC_TESTS/numerical_tests/test_2_geometry. The B+ and B- inputs only
# differ by the sign of PHIEDGE; iota is prescribed (tokamak) or the toroidal
# current is zero (W7-X), so the two equilibria have the same flux surfaces,
# iota and pressure, and the opposite magnetic field. The script checks that
# the B- file is the B+ file with the field quantities reversed.
#
# The W7-X inputs were written from wout_w7x_kjm.nc (boundary, axis, pressure
# and toroidal flux of its last iteration) by write_w7x_inputs() below, called
# with --write-w7x-inputs; the headers of the files were added by hand.
#
# Requirements (made with these versions on 25/11/2026):
#   pip install vmecpp==0.9.4 netCDF4
#
# Usage:
#   python make_vmec_files.py
################################################################################

import sys
import pathlib
import numpy as np
import netCDF4 as nc

here = pathlib.Path(__file__).resolve().parent
test_folder = here.parent
names = ['tok_circ', 'w7x_kjm']

# Quantities that keep or reverse their sign when the magnetic field is reversed
even_in_B = ['rmnc', 'zmns', 'lmns', 'bmnc', 'gmnc', 'iotaf', 'presf']
odd_in_B = ['phi', 'chi', 'bsupumnc', 'bsupvmnc', 'bsubumnc', 'bsubvmnc', 'buco', 'bvco']


def write_w7x_inputs():
    with nc.Dataset(test_folder / 'wout_w7x_kjm.nc') as wout:
        nfp = int(wout['nfp'][...])
        xm = np.asarray(wout['xm'][...]).astype(int)
        xn = (np.asarray(wout['xn'][...]) / nfp).astype(int)
        rbc, zbs = np.asarray(wout['rmnc'][-1]), np.asarray(wout['zmns'][-1])
        raxis, zaxis = np.asarray(wout['raxis_cc'][...]), np.asarray(wout['zaxis_cs'][...])
        am, phiedge = np.asarray(wout['am'][...]), float(wout['phi'][-1])
        mpol, ntor = int(wout['mpol'][...]), int(wout['ntor'][...])
    for sign, field in [(1, 'Bplus'), (-1, 'Bminus')]:
        lines = ['&INDATA', '  LFREEB = F', '  DELT = 0.9', '  NSTEP = 200', '  TCON0 = 2.0',
                 '  NS_ARRAY = 17 49', '  FTOL_ARRAY = 1.0e-10 1.0e-13', '  NITER_ARRAY = 2000 20000',
                 f'  NFP = {nfp}', '  LASYM = F', f'  MPOL = {mpol}', f'  NTOR = {ntor}',
                 f'  PHIEDGE = {sign*phiedge:.16e}', '  NCURR = 1', "  PCURR_TYPE = 'power_series'", '  AC = 0.0', '  CURTOR = 0.0',
                 "  PMASS_TYPE = 'power_series'", '  PRES_SCALE = 1.0', '  GAMMA = 0.0',
                 '  AM = ' + ' '.join(f'{x:.16e}' for x in am[:3]),
                 '  RAXIS_CC = ' + ' '.join(f'{x:.16e}' for x in raxis), '  ZAXIS_CS = ' + ' '.join(f'{x:.16e}' for x in zaxis)]
        for m, n, r, z in zip(xm, xn, rbc, zbs):
            if abs(r) > 1e-14: lines.append(f'  RBC({n},{m}) = {r:.16e}')
            if abs(z) > 1e-14: lines.append(f'  ZBS({n},{m}) = {z:.16e}')
        lines.append('/')
        with open(here / f'input.w7x_kjm_{field}', 'w') as file:
            file.write('\n'.join(lines) + '\n')


def solve_vmec_inputs():
    import vmecpp
    for name in names:
        for field in ['Bplus', 'Bminus']:
            output = vmecpp.run(vmecpp.VmecInput.from_file(here / f'input.{name}_{field}'), verbose=False)
            output.wout.save(str(test_folder / f'wout_{name}_{field}.nc'))


def check_vmec_files():
    for name in names:
        with nc.Dataset(test_folder / f'wout_{name}_Bplus.nc') as plus, nc.Dataset(test_folder / f'wout_{name}_Bminus.nc') as minus:
            for quantity in even_in_B + odd_in_B:
                sign = -1 if quantity in odd_in_B else 1
                a, b = np.asarray(plus[quantity][...]), sign * np.asarray(minus[quantity][...])
                difference = np.max(np.abs(a - b)) / max(np.max(np.abs(a)), 1e-300)
                if difference > 1e-12:
                    sys.exit(f'ERROR: {quantity} of wout_{name}_Bminus.nc is not {"-" if sign<0 else "+"}{quantity} of wout_{name}_Bplus.nc ({difference:.1e}).')
    print('The B- equilibria are the B+ equilibria with the magnetic field reversed.')


if __name__ == '__main__':
    if '--write-w7x-inputs' in sys.argv: write_w7x_inputs()
    solve_vmec_inputs()
    check_vmec_files()
