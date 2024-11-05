# Code for 'Mathematical Modelling of Mechanotransduction via RhoA Signalling Pathways' 
Sofie Verhees<sup>1</sup>, Mariya Ptashnyk<sup>1</sup>, Chandrasekhar Venkataraman<sup>2</sup>

<sup>1</sup>Department of Mathematics, The Maxwell Institute for Mathematical Sciences, Heriot-Watt University, Edinburgh; <sup>2</sup>Department of Mathematics, University of Sussex, Brighton

## Abstract
We derive and analyse a mathematical model for mechanotransduction related to the Rho GTPase signalling pathway. The mathematical model addressed the bidirectional coupling between signalling processes and cell mechanics, regarded as linearly elastics. A numerical method based on bulk-surface finite elements is proposed for the approximation of the  coupled system of nonlinear reaction-diffusion equations, defined inside the cell and on the cell membrane, and equations of linear elasticity. Our simulation results illustrate novel emergent features such as the strong dependence of the dynamics on cell shape, a threshold like response to changes in substrate stiffness and  that a bidirectional coupling between cell mechanics and signalling processes can lead to robustness of cell deformation to larger changes in substrate stiffness, ensuring mechanical homeostasis as observed experimentally.

## Prerequisite Packages
Create an anaconda environment with
- conda 23.3.1
- conda-build 3.24.0
- Python 3.10.9.final.0

Then, you can use 'conda-forge' to install
- fenics-dolfin 2019.1.0
- fenics-ufl 2019.1.0

Additionally, install in the same environment with pip
- meshio 5.3.4
- numpy 1.24.3
- argparse
- matplotlib 3.7.1
- pyvista 0.44.1

## Code Structure
..

## Relate Code with Figures
...
