# `field_equations` Directory

This directory contains the code for solving the main field equations used in stella.

## Equations Solved

The following equations are implemented:

---

### 1. Quasineutrality Equation  

$$
\sum_s Z_s n_s = 0
$$

where $Z_s$ is the charge and $n_s$ is the density of species $s$.
       
In terms of the distribution function this is:

$$
\sum_s Z_s \int \delta f_s d^3v = 0
$$

where $\delta f_s$ is the perturbed distribution function for species $s$, and the integral is performed at fixed particle position, $r$.
        
In the electromagnetic case, quasineutrality becomes

$$
\sum_s Z_s \int \langle g_s \rangle_{\mathbf{r}} + \frac{Z_s e}{T_s} F_0 \left( \langle \langle \chi \rangle_{R} \rangle_{r} - \phi \right) d^3v = 0
$$

where the angle brackets indicate the gyroaverage at fixed particle position, $r$, or gyrocenter position, $R$.  
Here $\chi = \phi - \mathbf{v} \cdot \mathbf{A}$ is the generalised potential, with $\phi$ the electrostatic potential, and $\mathbf{A}$ the magnetic vector potential. 

---

#### Normalisations

Quantities are normalised with respect to reference length scales, $L_r$, reference gyro-radius, $\rho_r$, reference temperature, $T_r$, and reference density, $n_r$, which are user-imposed.

$$
\tilde{\phi} = \frac{e}{T_r}\frac{L_r}{\rho_r} \phi
$$

$$
\tilde{A}_{\parallel} = \frac{L_r}{B_r\rho_r^2} A_{\parallel}
$$

$$
\delta \tilde{B}_{\parallel} = \frac{L_r}{B_r \rho_r} \delta B_{\parallel}
$$

Quasineturality is normalised by taking its product with the factor $\frac{L_r}{\rho_r e n_r}$.  
Any velocity integrals that are independent of $g_s$ can be evaluated, and the quasineutrality equation, in Fourier space, solved by **stella** is:

$$
\sum_s Z_s n_s \frac{2 B_0}{\sqrt{\pi}} \int d^3 v J_{0,s} \tilde{g}_s + \frac{Z_s^2 n_s}{T_s}(\Gamma_{0,s} - 1) \tilde{\phi} + \frac{Z_s n_s}{B_0} \Gamma_{1,s} \delta \tilde{B}_{\parallel} = 0
$$

---

### 2. Parallel Amperè's Law  

TO DO  

$$
-\nabla_\perp^2 A_\parallel = \mu_0 \sum_s j_{\parallel, s}
$$

where $A_\parallel$ is the parallel component of the vector potential and $j_{\parallel, s}$ is the parallel current.

---

### 3. Parallel Magnetic Field Equation  

TO DO  

$$
B_\parallel = \nabla \times \mathbf{A}
$$

where $B_\parallel$ is the parallel magnetic field fluctuation.
