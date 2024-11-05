# Code for 'Mathematical Modelling of Mechanotransduction via RhoA Signalling Pathways' 
Sofie Verhees<sup>1</sup>, Mariya Ptashnyk<sup>1</sup>, Chandrasekhar Venkataraman<sup>2</sup>

<sup>1</sup>Department of Mathematics, The Maxwell Institute for Mathematical Sciences, Heriot-Watt University, Edinburgh; <sup>2</sup>Department of Mathematics, University of Sussex, Brighton


## Abstract
We derive and analyse a mathematical model for mechanotransduction related to the Rho GTPase signalling pathway. The mathematical model addressed the bidirectional coupling between signalling processes and cell mechanics, regarded as linearly elastics. A numerical method based on bulk-surface finite elements is proposed for the approximation of the  coupled system of nonlinear reaction-diffusion equations, defined inside the cell and on the cell membrane, and equations of linear elasticity. Our simulation results illustrate novel emergent features such as the strong dependence of the dynamics on cell shape, a threshold like response to changes in substrate stiffness and  that a bidirectional coupling between cell mechanics and signalling processes can lead to robustness of cell deformation to larger changes in substrate stiffness, ensuring mechanical homeostasis as observed experimentally.


## Prerequisite Packages
A [conda](https://docs.conda.io/en/latest/) environment is ceated with
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
- 'meshes': includes mesh files (.geo and .msh), created with [gmsh](https://gmsh.info/) 4.11.1, for different cell shapes. These files should not be touched to recreate any results.
- 'signalling_model': includes all files needed to recreate results for the 'reduced model' without mechanics, as given in Appendix A.1.
- 'coupled_model': includes all files needed to recreate any results for the model with bidirectional coupling between the chemistry and the mechanics.

Both '_model' folders include:
- 'loop.sh': bash file running 'main.py' for different parameters and versions of the model.
- 'main.py': creates simulations of the model, includes all main code.
- 'summary_graphs.py': given the saved data from 'main.py', creates graphs that plot $f(\phi_a)$, $div(u)$, $\phi_a$ and $\rho_a$ as a function of the susbstrate stiffness $E$, as can be found in Figures 5, 10, 15 & 20 (for the coupled model), and Figures 22 & 24 in Appendix A.1 (for the 'reduced model' without mechanics).

For the coupled model, there is one extra file:
- 'tempstats.py': given the saved data from 'main.py', creates graphs that plot $f(\phi_a)$, $div(u)$, $\phi_a$ and $\rho_a$ as a function of time, as can be found in Figure 25 in Appendix A.2.

> [!WARNING]
> Files 'summary_graphs.py' and 'tempstats.py' should only be run after running the corresponding 'loop.sh' as it requires the output from this programme.

### Different parameters/ versions
'loop.sh' contains different variables that determine which version of the model with which parameters is run.
- For 'signalling_model' without mechanics, the variables are:
  - `D1`: determines the diffusion coefficient of $\phi_d$ and $\phi_a$ (in x0.1 $\mu$ m<\sup>2<\sup>/s)
  - `twoDstim`: determines the kind of stimulus (2D, 2xD, 3D), as explained in the paper.
  - `E`: determines the substrate stiffness
- For the coupled model, the variables are:
  - `T`: gives the final time (in x10s).
  - `dt`: gives the time step (in x10s).
  - `mesh_name`: determines the shape of the cell, 'cell_substrate' being the radially symmetric cell, 'lamellipodium' the lamellipodium shape, and 'cell_substrate_empty-nucleus' the radially symmetric cell with a nucleus.
  - `coupling`: determines which coupling between the mechanics and the chemistry, where 1 -> $E_c=0.6$ & $C_1 = 0$; 2 -> $E_c=0.6$ & $C_1 = 1$; 3-> $E_c=f(\phi_a)$ & $C_1 = 0$; 4 -> $E_c=f(\phi_a)$ & $C_1 = 1$.
  - `C1`: determines the value of $C_1$
  - `E`: determines the substrate stiffness
  - `twoDstim`: determines where the substrate stiffness has effect with `yes` meaning only on the bottom boundary (2xD stimulus) and `no` meaning everywhere on the boundary (3D stimulus).
  - `partfixed`: determines which boundary condition is put on u with `yes` meaning no deformation in the vertical direction on the bottom of the cell (rigid substrate) and `no` meaning only the force boundary condition.

The resulting files will be named accordingly, e.g. for the 'reduced model' simululation results are named 'rhomodel_ _twoDstim_ _ _var_ _reduced_dt= _dt_ _T= _T_ _D= _D1_ _E_ E.pvd', where var is the variable $\phi_d$ for cd, $\phi_a$ for ca and $\rho_a$ for p; and for the coupled model simulations results are named 'coupled _coupling_ 2D_var_(partfixed/neumann)_m(0/1)_dt=.._T=.._k6=.._..E.pvd', where var is the variable \phi_d for cd, \phi_a for ca and \rho_a for p and coupling(1/2/3/4) determines the coupling between the mechanics and the chemistry as described before and if 2D is not there, it is the 3D stimulus and partfixed meaning u is fixed on the bottom boundary while neumann meaning only the force boundary condition and m1 means the radially symmtric cell shape while m0 means the lamellipodium cell shape


## Relate Results to Figures
If the files are run as instructed, a folder 'results' will be created including the following subfolders
- 'simulations': includes simulation results (.vtu and .pvd), which are visualised with [Paraview](https://www.paraview.org/) 5.10.0-RC1 and correspond to Figures 1-4, 6-9, 11-14 & 16-19, and Figures 27 & 28 in Appendix A.5.
- 'temp': includes saved numpy arrays (.npy) with the discrete time derivative, L2-norm, minumum and maximum values of $\phi_a$, $\phi_d$, $\rho_a$ and $u$ as well as other variables saved at each time step.
- 'figures': includes figures (.png) that correspond to Figures 5, 10, 15 & 20, Figures 22 & 24 in Appendix A.1, and Figure 25 in Appendix A.2.

