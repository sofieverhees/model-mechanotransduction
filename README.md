# Code for 'Mathematical Modelling of Mechanotransduction via RhoA Signalling Pathways' 
Sofie Verhees<sup>1</sup>, Mariya Ptashnyk<sup>1</sup>, Chandrasekhar Venkataraman<sup>2</sup>

<sup>1</sup>Department of Mathematics, The Maxwell Institute for Mathematical Sciences, Heriot-Watt University, Edinburgh; <sup>2</sup>Department of Mathematics, University of Sussex, Brighton


## Abstract
We derive and analyse a mathematical model for mechanotransduction related to the Rho GTPase signalling pathway. The mathematical model addressed the bidirectional coupling between signalling processes and cell mechanics, regarded as linearly elastics. A numerical method based on bulk-surface finite elements is proposed for the approximation of the  coupled system of nonlinear reaction-diffusion equations, defined inside the cell and on the cell membrane, and equations of linear elasticity. Our simulation results illustrate novel emergent features such as the strong dependence of the dynamics on cell shape, a threshold like response to changes in substrate stiffness and  that a bidirectional coupling between cell mechanics and signalling processes can lead to robustness of cell deformation to larger changes in substrate stiffness, ensuring mechanical homeostasis as observed experimentally.


## Prerequisite Packages
Create a [conda](https://docs.conda.io/en/latest/) environment with
- conda 23.3.1
- conda-build 3.24.0
- Python 3.10.9.final.0

Then, conda-forge can be used to [install fenics](https://fenics.readthedocs.io/en/latest/installation.html)
- fenics-dolfin 2019.1.0
- fenics-ufl 2019.1.0

Additionally, the next packages have to be installed (e.g. with pip)
- meshio 5.3.4
- numpy 1.24.3
- argparse 1.4.0
- matplotlib 3.7.1
- pyvista 0.44.1


## Code Structure
This code is organised in three different folders:
- 'meshes': includes mesh files (.geo and .msh), created with [gmsh](https://gmsh.info/) 4.11.1, for different cell shapes. These files should not be touched.
- 'signalling_model': includes all files needed to get results for the 'reduced model' without mechanics, as given in Appendix A.1.
- 'coupled_model': includes all files needed to get results for the model with bidirectional coupling between the chemistry and the mechanics.

Both '_model' folders include:
- 'loop.sh':
- 'main.py':
- 'summary_graphs.py':

For the coupled model, there is one extra file:
- 'tempstats.py':


## Relate Code to Figures
If the files are run as instructed, a folder 'results' will be created including the following subfolders
- 'simulations': includes simulation results (.vtu and .pvd), which are visualised with [Paraview](https://www.paraview.org/) 5.10.0-RC1 and correspond to Figures 1-4, 6-9, 11-14 & 16-19, and Figures 27 & 28 in Appendix A.5.
- 'temp': includes saved numpy arrays (.npy) with the discrete time derivative, L2-norm, minumum and maximum values of $\phi_a$, \phi_d, rho_a and u as well as other variables saved at each time step.
- 'figures': includes figures (.png) that correspond to Figures 5, 10, 15 & 20, Figures 22 & 24 in Appendix A.1, and Figure 25 in Appendix A.2.

