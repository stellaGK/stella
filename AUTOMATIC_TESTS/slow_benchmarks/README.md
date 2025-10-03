Automated stella benchmarks (slow)
==================================

Perform some simulations from which we know the expected (gamma,omega) or heat fluxes.

Various benchmarks of the stella gyrokinetic code have been performed. Their
published results are gathered here and are reproduced by the numerical tests.

## Benchmark of stella and GENE by A. Gonzalez 

The stella code has been benchmarked against GENE in:

[[2022 - Gonzalez - Electrostatic gyrokinetic simulations in Wendelstein 7-X geometry benchmark between the codes stella and GENE]]

The paper by A. Gonzalez contains 5 benchmark tests:

- **Test 1**: Linear scan in (kx,ky) of ITG turbulence, considering $a/L_{T_i} = 3$ and $a/L_{n_i} = 1$, and simulating kinetic ions and adiabatic electrons, performed in the flux tube centered around the bean-shaped cross-section of W7-X.
- **Test 2**: Linear scan in (kx) of ITG turbulence, considering $a/L_{T_i} = 3$ and $a/L_{n_i} = 1$, and simulating kinetic ions and adiabatic electrons, performed in the flux tube centered around the triangular-shaped cross-section of W7-X.
- **Test 3**: Linear scan in (ky) of TEM turbulence, considering $a/L_{T_s} = 0$ and $a/L_{n_s} = 3$, and simulating kinetic ions and electrons, performed in the flux tube centered around the bean-shaped cross-section of W7-X.
- **Test 4**: Rosenbluth-Hinton test which tests the zonal-flow relaxation in the bean flux tube, considering $a/L_{T_i} = 0$ and $a/L_{n_i} = 0$, performed in the flux tube centered around the bean-shaped cross-section of W7-X.
- **Test 5**: Nonlinear simulation of the ion heat flux, considering ITG turbulence with $a/L_{T_i} = 3$ and $a/L_{n_i} = 1$, and simulating kinetic ions and adiabatic electrons, performed in the flux tube centered around the bean-shaped cross-section of W7-X.

