# `field_equations` Directory

This directory contains the code for solving the main field equations used in stella.

## Equations Solved

The following equations are implemented:

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

---

### 1. Quasineutrality Equation  

Quasineutrality is defined as: 
$$
\sum_s Z_s \delta n_s = 0
$$

where $Z_s$ is the charge and $\delta n_s$ is the perturbed density of species $s$.
       
In terms of the distribution function this is:

$$
\sum_s Z_s \int \delta f_s d^3v = 0
$$

where $\delta f_s$ is the perturbed distribution function for species $s$, and the integral is performed at fixed particle position, $r$.
        
In the electromagnetic case, quasineutrality becomes

$$
\sum_s Z_s \int \langle g_s \rangle_{r} + \frac{Z_s e}{T_s} F_0 \left( \langle \langle \chi \rangle_{R} \rangle_{r} - \phi \right) d^3v = 0
$$

where the angle brackets indicate the gyroaverage at fixed particle position, $r$, or gyrocenter position, $R$.  
Here $\chi = \phi - \mathbf{v} \cdot \mathbf{A}$ is the generalised potential, with $\phi$ the electrostatic potential, and $\mathbf{A}$ the magnetic vector potential. 

Quasineturality is normalised by taking its product with the factor $\frac{L_r}{\rho_r e n_r}$.  
Any velocity integrals that are independent of $g_s$ can be evaluated, and the quasineutrality equation, in Fourier space, solved by **stella** is:

$$
\sum_s Z_s n_s \frac{2 B_0}{\sqrt{\pi}} \int d v_{\parallel} \int d \mu_s J_{0,s} \tilde{g}_s + \frac{Z_s^2 n_s}{T_s}(\Gamma_{0,s} - 1) \tilde{\phi} + \frac{Z_s n_s}{B_0} \Gamma_{1,s} \delta \tilde{B}_{\parallel} = 0
$$

---

### 2. Parallel Amperè's Law  

Amperè's law is $\nabla \times \delta \mathbf{B} = \frac{4 \pi}{c} \delta \mathbf{J}$. We can take parallel component to get 

$$
\nabla_{\perp}^2 A_{\parallel} = \frac{4 \pi}{c} \sum_{s} Z_{s}e \int d^3 v v_{\parallel} [\langle g_{s} \rangle_{r} + \frac{Z_{s} e}{T_{s} c} F_0 v_{\parallel} \langle \langle A_{\parallel} \rangle_{R} \rangle_{r} ] = 0.
$$ 

Parallel Amperè's Law is normalised by taking its product with the factor $\frac{L_r}{B_r}$. The equation solved in Fourier space is:

$$
\frac{\beta_r}{(k_{\perp}\rho_r)^2} \sum_s Z_s n_s v_{th} \frac{2B_0}{\sqrt{\pi}} \int d v_{\parallel} \int d \mu_s v_{\parallel} J_{0,s} g_{s} = \left[ 1+ \frac{\beta_r}{(k_{\perp}\rho_r)^2} \sum_s \frac{Z_s n_s}{m_s} \Gamma_{0,s} \right] A_{\parallel}
$$

---

### 3. Perpendicular Amperè's Law

Perpendicular Amperè's Law gives an equation for $\delta B_{\parallel}$: 

$$
8\pi \sum_s \frac{2\pi B_0}{m_s} \int d^3v \frac{J_1,s}{a_{s}} \mu_s g_{s} +\left[ 4\pi \sum_s \frac{Z_s e n_s}{B_0} \Gamma_{1,s} \phi \right] + \left[ 1+ 16\pi \sum_s \frac{n_s T_s}{B_0^2} \Gamma_{2,s}\delta B_{\parallel} \right] = 0.
$$ 

This is normalised by taking its product with the factor $\frac{L_r \rho_r}{B_r}$, to give:

$$
\nabla_{\perp}^2 \delta B_{\parallel} = \frac{2B_0 \beta_r}{\sqrt{\pi}} \sum_s Z_s n_s \int d v_{\parallel} \int d \mu_s \nabla \cdot \left\langle (g_{s} + \frac{Z_s}{T_{s}} \left\langle \chi_{s} \right\rangle_{r} e^{-v^2})(\hat{\mathbf{b}} \times \mathbf{v}_{\perp}) \right\rangle_{r} .
$$
