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
A Docker image with the full environment (including PorePy and PETSc) for reproducing
the results in the article is available on Zenodo
[here](https://doi.org/10.5281/zenodo.13861514).

* All scripts should be run from the mpsa_newmark directory. An example of how to run a
script after `cd` into the mpsa_newmark directory is as shown in the following:

  `python convergence_analysis/runscript_energy_decay_vary_theta.py`

* There are certain boolean flags to adapt the simulations and outputs of the simulation
  runscripts:
  * `coarse`: Some of the simulations are quite time consuming and memory intensive. To
    ensure low default memory usage and low default wall clock time for running the
    runscripts, all simulation runscripts are assigned coarse default parameter values.
    Change the value of `coarse` from `True` to `False` to run simulations that match
    those detailed in the article.
  * `save_figure`: This flag is found in runscripts which allow for generating and
    saving figures. Change `save_figure` from `True` to `False` to not generate and save
    figures.
  * `limit_file_export`: Certain runscripts generate large amounts of vtu and pvd files.
    Default behavior is to generate files for each time step. This flag allows to only
    generate files for the time steps which are presented in the article. Change
    `limit_file_export` from `False` to `True` to limit the number of files exported.

Note that not all the runscripts have have all the flags. If a certain runscript does
not have access to a certain flag, it means that the flag is not needed for that
specific simulation.

## Verification: Convergence and energy decay analyses
### Convergence analysis of MPSA-Newmark
These convergence analyses are performed with homogeneous Dirichlet conditions on a 3D
simplex grid. All convergence runscripts generate an output file which contains grid
size, number of cells, time step size, displacement error and traction error.

* Convergence in space and time:
  * [runscript_space_time_convergence_dirichlet_boundaries](./convergence_analysis/runscript_space_time_convergence_dirichlet_boundaries.py)
    (this script can also save a figure of the convergence results).


* Convergence in space:
  * [runscript_space_convergence_dirichlet_boundaries](./convergence_analysis/runscript_space_convergence_dirichlet_boundaries.py)
    

* Convergence in time:
  * [runscript_time_convergence_dirichlet_boundaries](./convergence_analysis/runscript_time_convergence_dirichlet_boundaries.py) 

    The default parameters run the time convergence analysis with two refinement levels
    and a cell size coarser than ideal. The coarser grid is not sufficiently fine to
    observe the second order convergence in time. 

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


Varying the wave incidence angle, $\theta$:
* [runscript_energy_decay_vary_theta](./convergence_analysis/runscript_energy_decay_vary_theta.py).

## Results: Simulation examples
Simulation example runscripts are found within a dedicated [example
runscripts](./example_runscripts/) directory.


The simulation from Example 1.1, which considers a seismic source located inside an
inner transversely isotropic domain:
* [runscript_example_1_1_source_in_inner_domain](./example_runscripts/runscript_example_1_1_source_in_inner_domain.py).


The simulation from Example 1.2, which considers a seismic source located outside an
  inner transversely isotropic domain:
* [runscript_example_1_2_source_in_outer_domain](./example_runscripts/runscript_example_1_2_source_in_outer_domain.py).


The simulation from Example 2, which considers a layered heterogeneous medium with an
open fracture:
* [runscript_example_2_heterogeneous_fractured_domain](./example_runscripts/runscript_example_2_heterogeneous_fractured_domain.py).

    Note: The example 2 simulation with `coarse` = `False` is time consuming and memory
    intensive. Be aware that it generates visualization files which require over 20GB of
    storage if both `coarse = False` and `limit_file_export`= `False`.