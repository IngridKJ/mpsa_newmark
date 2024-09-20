# mpsa_newmark
This repository contains everything needed to run the simulation examples found in the
MPSA-Newmark paper.

That includes:
* Runscripts for the convergence and energy decay analyses.
* Runscripts for all simulation examples.
* Standardized model class setup for solving the elastic wave equation using PorePy
  (https://github.com/pmgbergen/porepy).
* Utility material which is used in the various simulations
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
* Convergence in space and time:
  * [runscript_space_time_convergence_dirichlet_boundaries](./convergence_analysis/runscript_space_time_convergence_dirichlet_boundaries.py)
    (this script also renders the convergence result figure).
* Convergence in space:
  * [runscript_space_convergence_dirichlet_boundaries](./convergence_analysis/runscript_space_convergence_dirichlet_boundaries.py)
* Convergence in time:
  * [runscript_time_convergence_dirichlet_boundaries](./convergence_analysis/runscript_time_convergence_dirichlet_boundaries.py) 

All the runscripts utilize
[manufactured_solution_dynamic_3D](./convergence_analysis/convergence_analysis_models/manufactured_solution_dynamic_3D.py)
as the manufactured solution setup.

### Convergence analysis of MPSA-Newmark with absorbing boundaries
Convergence of the solution is performed in a quasi-1D setting. We have performed a
convergence analysis with successive refinemenet in both space and time. The script also
renders a convergence plot of the results:
  * [runscript_space_time_convergence_absorbing_boundaries](./convergence_analysis/runscript_space_time_convergence_absorbing_boundaries.py)


The runscripts utilize
[model_convergence_ABC](./convergence_analysis/convergence_analysis_models/model_convergence_ABC.py)
as the model class setup. 

### Energy decay analysis of MPSA-Newmark with absorbing boundaries
The energy decay analysis is performed both for successive refinement 
of the grid, as well as for varying wave incidence angles. 

Rendering the figure for:
* Rendering the figure for successive grid refinement is done by running the script
[runscript_energy_decay_space_refinement](./convergence_analysis/runscript_energy_decay_space_refinement.py).


* Rendering the figure for changing the wave incidence angle, $\theta$, is done by running the script
[runscript_energy_decay_vary_theta](./convergence_analysis/runscript_energy_decay_vary_theta.py).


## Results: Simulation examples
Simulation example runscripts are found within [this](./example_runscripts/) directory.
* The simulation from Example 1.1, which considers a seismic source located inside an
  inner transversely isotropic domain, is run by
  [runscript_example_1_1_source_in_inner_domain](./example_runscripts/runscript_example_1_1_source_in_inner_domain.py).
* The simulation from Example 1.2, which considers a seismic source located outside an
  inner transversely isotropic domain, is run by
  [runscript_example_1_2_source_in_outer_domain](./example_runscripts/runscript_example_1_2_source_in_outer_domain.py).
* The simulation from Example 2, which considers a layered heterogeneous medium with an
  open fracture, is run by
  [runscript_example_2_heterogeneous_fractured_domain](./example_runscripts/runscript_example_2_heterogeneous_fractured_domain.py).


## Models
All the above scripts utilize a common model class for solving the elastic wave
equation. There are currently two model class setups:
* [elastic_wave_equation_abc](./models/elastic_wave_equation_abc.py) considers a general
  setup for solving the elastic wave equation with absorbing boundaries on all domain
  sides.
* [elastic_wave_equation_abc_linear](./models/elastic_wave_equation_abc_linear.py)
  inherits from [elastic_wave_equation_abc](./models/elastic_wave_equation_abc.py), but
  makes sure that the Jacobian is only assembled once. For linear problems it is not
  necessary to assemble the Jacobian every time step, as it is constant.