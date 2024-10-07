# mpsa_newmark
This repository contains everything needed to reproduce the results found in the
MPSA-Newmark paper.

That includes:
* Runscripts for the convergence and energy decay analyses.
* Runscripts for all simulation examples.
* Standardized model class setup for solving the elastic wave equation using PorePy.
  (https://github.com/pmgbergen/porepy).
* Utility material which is used in the various simulations.
* Additional material: We have included runscripts for separate space and time
convergence analyses of MPSA-Newmark with Dirichlet boundary conditions.

## How to use:
All scripts should be run from the mpsa_newmark directory. An example of how to run a
script after `cd` into the mpsa_newmark directory is as shown in the following:

`python convergence_analysis/runscript_energy_decay_vary_theta.py`

## Verification: Convergence and energy decay analyses
### Convergence analysis of MPSA-Newmark
These convergence analyses are performed with homogeneous Dirichlet conditions on a 3D
simplex grid. All convergence runscripts generate an output file which contains grid
size, number of cells, time step size, displacement error and traction error.

Note: These simulations are quite memory intensive and time consuming. Therefore we have
assigned some default parameters which lead to less demanding simulations. The parameters which reproduce the results in the article are detailed for each script:
* Convergence in space and time:
  * [runscript_space_time_convergence_dirichlet_boundaries](./convergence_analysis/runscript_space_time_convergence_dirichlet_boundaries.py)
    (this script also renders a figure of the convergence results).

    The default parameters run the space/time convergence analysis with two refinement
    levels. To reproduce the results which are found in the article:
    * Change the "levels" parameter in line 58 from 2 to 4.
    * Uncomment lines 108-128 to plot convergence slope triangles. 
* Convergence in space:
  * [runscript_space_convergence_dirichlet_boundaries](./convergence_analysis/runscript_space_convergence_dirichlet_boundaries.py)
    
    The default parameters run the space convergence analysis with two refinement
    levels. To produce convergence analysis results with four refinement levels instead
    of two:
      * Change the "levels" parameter in line 54 from 2 to 4.
* Convergence in time:
  * [runscript_time_convergence_dirichlet_boundaries](./convergence_analysis/runscript_time_convergence_dirichlet_boundaries.py) 

    The default parameters run the time convergence analysis with two refinement levels
    and a cell size coarser than ideal. The coarser grid is not sufficiently fine to
    observe the second order convergence in time. To reproduce results to observe second
    order convergence in time:
    * Change the "cell_size" from 0.1 to 0.03125 in line 46.
    * Change the "levels" parameter in line 53 from 2 to 4. 
    

All the runscripts utilize
[manufactured_solution_dynamic_3D](./convergence_analysis/convergence_analysis_models/manufactured_solution_dynamic_3D.py)
as the manufactured solution setup.

### Convergence analysis of MPSA-Newmark with absorbing boundaries
Convergence of the solution is performed in a quasi-1D setting. We have performed a
convergence analysis with successive refinemenet in both space and time. The script
generates a file with displacement and traction errors, and a convergence plot of the
results:
  * [runscript_space_time_convergence_absorbing_boundaries](./convergence_analysis/runscript_space_time_convergence_absorbing_boundaries.py)


The runscripts utilize
[model_convergence_ABC](./convergence_analysis/convergence_analysis_models/model_convergence_ABC.py)
as the model class setup. 

### Energy decay analysis of MPSA-Newmark with absorbing boundaries
The energy decay analysis is performed both for successive refinement of the grid, as
well as for varying wave incidence angles. In both cases the kinetic energy values for
each time step in each simulation are saved in files within the directory
[energy_values](./convergence_analysis/energy_values/).

Grid refinement:
* [runscript_energy_decay_space_refinement](./convergence_analysis/runscript_energy_decay_space_refinement.py).

  To reproduce the results in the article:
  * Change np.range(5, 7) with np.range(5, 10) in line 163.
  * Uncomment lines 205, 206 and 207 for plotting all of the refinements.


Varying the wave incidence angle, $\theta$:
* [runscript_energy_decay_vary_theta](./convergence_analysis/runscript_energy_decay_vary_theta.py).
  
  To reproduce the results in the article:
  * Change the cell size value from 0.1 to 0.015625 in line 118.


## Results: Simulation examples
Simulation example runscripts are found within a dedicated [example
runscripts](./example_runscripts/) directory.
* The simulation from Example 1.1, which considers a seismic source located inside an
  inner transversely isotropic domain, can be run by
  [runscript_example_1_1_source_in_inner_domain](./example_runscripts/runscript_example_1_1_source_in_inner_domain.py).

  * To reproduce the exact simulation in the article, change the cell_size from 0.1 to
    0.0125 in line 35.
* The simulation from Example 1.2, which considers a seismic source located outside an
  inner transversely isotropic domain, can be run by
  [runscript_example_1_2_source_in_outer_domain](./example_runscripts/runscript_example_1_2_source_in_outer_domain.py).

    * To reproduce the exact simulation in the article, change the cell_size from 0.1 to
      0.0125 in line 36.
* The simulation from Example 2, which considers a layered heterogeneous medium with an
  open fracture, is run by
  [runscript_example_2_heterogeneous_fractured_domain](./example_runscripts/runscript_example_2_heterogeneous_fractured_domain.py).

    * To reproduce the simulation in the article, change the cell_size from 0.25 to 0.0175
    in line 151. 
    
    Note: The simulation with a cell size of 0.0175 is time consuming and memory
    intensive. Be aware that it generates visualization files which require over 20GB of
    storage.
