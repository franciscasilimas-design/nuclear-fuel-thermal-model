# nuclear-fuel-thermal-model
Thermal modelling of a nuclear fuel rod (UO2) including nonlinear thermal conductivity.
# Thermal modelling of a nuclear fuel rod

This project presents a simplified thermal model of a nuclear fuel rod.

## Physical model

The model considers:

- Radial heat conduction
- Volumetric heat generation
- Three regions:
  - UO2 fuel pellet
  - Helium gap
  - Zirconium cladding
- Convective heat transfer with coolant

## Nonlinear model

The thermal conductivity of UO2 is temperature dependent:

k(T) = 1 / (a + bT)

The nonlinear heat equation is solved using:

- Kirchhoff transformation
- Numerical iterative solver in MATLAB

## Results

The comparison between linear and nonlinear models shows a difference of approximately **47 K** in the fuel centre temperature.

## Files

- `VF_Projet2026_Modélisation thermique simplifiée.pdf` : complete technical report
- `Projet_perso.m` : MATLAB implementation
- `Analyse_param.m` : parametric study
