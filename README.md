This repository contains a modular Fortran 90 implementation for solving the 1D Convection-Diffusion equation. The project explores various spatial discretization schemes and high-order time integration.

## Physical Models
The code is designed to solve the general transport equation:
$$\frac{\partial \phi}{\partial t} = - c \frac{\partial \phi}{\partial x} + \nu \frac{\partial^2 \phi}{\partial x^2}$$

The repository is structured to handle three specific cases:
1. **Pure Convection** (Gaussian and Harmonic Initial COnditions)
2. **Pure Diffusion** (Gaussian Initial Condition)
3. **Convection-Diffusion** (Gaussian Initial Condition)

## Technical Features
- **Time Integration**: 3rd-order Runge-Kutta (RK3) scheme.
- **Spatial Schemes**: 
  - Central Schemes: 2nd-order (CS2) and 4th-order (CS4).
  - Upwind Schemes: 1st-order Forward (FW1) and Backward (BW1).
- **Boundary Conditions**: Periodic boundary conditions implemented for 1D domains.

## Project Structure
- `mod_param.f90`: Precision definitions (Double Precision).
- `mod_schemes.f90`: Implementation of spatial discretization (Convective and Diffusive terms).
- `mod_rk3.f90`: Third-order Runge-Kutta time-stepping logic.
- `mod_save.f90`: Utility for exporting results in `.txt` format for post-processing.
- `task1A.f90` / `task1B.f90`: Drivers for pure convection tests.
- `task2.f90`: Driver for pure diffusion tests.
- `task3.f90`: Driver for the complete convection-diffusion problem.

## Compilation
To compile the project, use a Fortran compiler (e.g., gfortran) following the dependency order:

```bash
gfortran -c mod_param.f90 mod_schemes.f90 mod_save.f90 mod_rk3.f90
gfortran task3.f90 mod_param.o mod_schemes.o mod_save.o mod_rk3.o -o task3
./task3
