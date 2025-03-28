import sys
from typing import Union

import numpy as np
import porepy as pp
from porepy.models.momentum_balance import MomentumBalance

sys.path.append("../")
import time_derivatives
from solvers.solver_mixins import CustomSolverMixin
from utils import acceleration_velocity_displacement, body_force_function, u_v_a_wrap


class NamesAndConstants:
    def _is_nonlinear_problem(self) -> bool:
        """All problem setups in this repository are linear."""
        return False

    @property
    def beta(self) -> float:
        """Newmark time discretization parameter, beta.

        Discretization parameter which somehow represents how the acceleration is
        averaged over the coarse of one time step. The value of beta = 0.25 corresponds
        to the "Average acceleration method".

        """
        return 1 / 4

    @property
    def gamma(self) -> float:
        """Newmark time discretization parameter, gamma.

        Discretization parameter which somehow represents the amount of numerical
        damping (positive, ngeative or zero). The value of gamma = 0.5 corresponds to
        no numerical damping.

        """
        return 1 / 2

    @property
    def velocity_key(self) -> str:
        """Key/Name for the velocity variable/operator.

        Velocity is represented by a time dependent dense array with the name provided
        by this property, namely "velocity".

        """
        return "velocity"

    @property
    def acceleration_key(self) -> str:
        """Key/Name for acceleration variable/operator.

        Acceleration is represented by a time dependent dense array with the name
        provided by this property, namely "acceleration".

        """

        return "acceleration"

    @property
    def bc_values_mechanics_key(self) -> str:
        """Key for mechanical boundary conditions in the data dictionary."""
        return "bc_values_mechanics"

    def primary_wave_speed(self, is_scalar: bool = False) -> Union[float, np.ndarray]:
        """Primary wave speed.

        Speed of the compressive elastic waves.

        Parameters:
            is_scalar: Whether the primary wave speed should be scalar or not. Relevant
                for use with manufactured solutions which only supports a scalar wave speed for now.

        Returns:
            The value of the compressive elastic waves. Either scalar valued or the
            cell-center values, depending on the input parameter.

        """
        rho = self.solid.density
        if not is_scalar:
            return np.sqrt((self.lambda_vector + 2 * self.mu_vector) / rho)
        else:
            return np.sqrt(
                (self.solid.lame_lambda + 2 * self.solid.shear_modulus / rho)
            )

    def secondary_wave_speed(self, is_scalar: bool = False) -> Union[float, np.ndarray]:
        """Secondary wave speed.

        Speed of the shear elastic waves.

        Parameters:
            is_scalar: Whether the primary wave speed should be scalar or not. Relevant
                for use with manufactured solutions which only supports a scalar wave speed for now.

        Returns:
            The value of the shear elastic waves. Either scalar valued or the
            cell-center values, depending on the input parameter.

        """
        rho = self.solid.density
        if not is_scalar:
            return np.sqrt(self.mu_vector / rho)
        else:
            return np.sqrt(self.solid.shear_modulus / rho)


class BoundaryAndInitialConditions:
    @property
    def discrete_robin_weight_coefficient(self) -> float:
        """Additional coefficient for discrete Robin boundary conditions.

        After discretizing the time derivative in the absorbing boundary condition
        expressions, there might appear coefficients additional to those within the
        coefficient matrix. This model property assigns that coefficient.

        Returns:
            The Robin weight coefficient.

        """
        return 3 / 2

    def bc_type_mechanics(self, sd: pp.Grid) -> pp.BoundaryConditionVectorial:
        """Boundary condition type for the absorbing boundary condition model class.

        Assigns Robin boundaries to all subdomain boundaries by default. This also
        includes setting the Robin weight.

        Parameters:
            sd: The subdomain whose boundaries are to be set.

        Returns:
            The vectorial boundary condition operator for subdomain sd.

        """
        # Fetch boundary sides and assign type of boundary condition for the different
        # sides
        bounds = self.domain_boundary_sides(sd)
        bc = pp.BoundaryConditionVectorial(
            sd,
            bounds.north
            + bounds.south
            + bounds.east
            + bounds.west
            + bounds.bottom
            + bounds.top,
            "rob",
        )

        # Calling helper function for assigning the Robin weight
        self.assign_robin_weight(sd=sd, bc=bc)
        return bc

    def assign_robin_weight(
        self, sd: pp.Grid, bc: pp.BoundaryConditionVectorial
    ) -> None:
        """Assigns the Robin weight for Robin boundary conditions.

        Parameters:
            sd: The subdomain whose boundary conditions are to be defined.
            bc: The vectorial boundary condition object.

        """
        # Initiating the arrays for the Robin weight
        r_w = np.tile(np.eye(sd.dim), (1, sd.num_faces))
        value = np.reshape(r_w, (sd.dim, sd.dim, sd.num_faces), "F")

        # The Robin weight should only be assigned to subdomains of max dimension.
        if sd.dim != self.nd:
            bc.robin_weight = value
        else:
            # The coefficient matrix is constructed elsewhere, so we fetch it by making
            # this call.
            total_coefficient_matrix = self.total_coefficient_matrix(
                sd=sd,
            )

            # Discretizing the u_t term in the absorbing boundaries introduces an
            # additional weight to the u-term of the resulting Robin boundary condition
            # (see the method total_coefficient_matrix).
            # The next line includes that weight.
            total_coefficient_matrix *= self.discrete_robin_weight_coefficient

            # Fethcing all boundary faces for the domain and assigning the
            # total_coefficient_matrix to the boundary faces.
            boundary_faces = sd.get_boundary_faces()
            value[:, :, boundary_faces] *= total_coefficient_matrix.T

            # Finally setting the actual Robin weight.
            bc.robin_weight = value

    def bc_values_stress(self, boundary_grid: pp.BoundaryGrid) -> np.ndarray:
        """Method for assigning Robin boundary condition values (the right-hand side).

        Specifically, this method assigns the values corresponding to absorbing boundary
        conditions when we have used a second order approximation to u_t in:

            sigma * n + alpha * u_t = G

        See the method total_coefficient_matrix() for the expressions.

        Parameters:
            boundary_grid: The boundary grids on which to define boundary conditions.

        Returns:
            Array of boundary values.

        """
        if boundary_grid.dim != (self.nd - 1):
            return np.array([])

        data = self.mdg.boundary_grid_data(boundary_grid)
        sd = boundary_grid.parent
        boundary_faces = sd.get_boundary_faces()
        name = "boundary_displacement_values"

        if self.time_manager.time_index == 0:
            return self.initial_condition_bc(boundary_grid)

        displacement_values_0 = pp.get_solution_values(
            name=name, data=data, time_step_index=0
        )
        displacement_values_1 = pp.get_solution_values(
            name=name, data=data, time_step_index=1
        )

        # According to the expression for the absorbing boundaries, we have a
        # coefficient 2 in front of the values u_(n-1) and -0.5 in front of u_(n-2):
        displacement_values = 2 * displacement_values_0 - 0.5 * displacement_values_1

        # Transposing and reshaping displacement values to prepare for broadcasting
        displacement_values = displacement_values.reshape(
            boundary_grid.num_cells, self.nd, 1
        )

        # Assembling the vector representing the RHS of the Robin conditions
        total_coefficient_matrix = self.total_coefficient_matrix(sd=sd)
        robin_rhs = np.matmul(total_coefficient_matrix, displacement_values).squeeze(-1)
        robin_rhs *= sd.face_areas[boundary_faces][:, None]
        robin_rhs = robin_rhs.T
        return robin_rhs.ravel("F")

    def total_coefficient_matrix(self, sd: pp.Grid) -> np.ndarray:
        """Coefficient matrix for the absorbing boundary conditions.

        This method is used together with Robin boundary condition (absorbing boundary
        condition) assignment.

        The absorbing boundary conditions, in the continuous form, look like the
        following:

            sigma * n + D * u_t = 0,

        where u_t is the velocity and  D is a matrix containing material parameters
        and wave velocities.

        Approximate the time derivative by a second order backward difference:

            sigma * n + 3 / 2 * D_t * u_n = D_t * (2 * u_(n-1) - 0.5 * u_(n-2)),

        where D_t = 1/dt * D. Note the coefficient 3 / 2 (termed
        discrete_robin_weight_coefficient in the method assign_robin_weight)

        This method creates D_t, which is subsequently used both as a part of the weight
        in front of u_n and in the right hand side.

        Parameters:
            sd: The subdomain where the boundary conditions are assigned.

        Returns:
            A block array which is the sum of the normal and tangential component of the
            coefficient matrices.

        """
        # Fetching necessary grid related quantities
        boundary_cells = self.boundary_cells_of_grid
        boundary_faces = sd.get_boundary_faces()
        face_normals = sd.face_normals[:, boundary_faces][: self.nd]
        unitary_face_normals = face_normals / np.linalg.norm(
            face_normals, axis=0, keepdims=True
        )

        # This line of code does the same as the lines that are commented out
        # below:
        tensile_matrices = np.einsum(
            "ik,jk->kij", unitary_face_normals, unitary_face_normals
        )
        # tensile_matrices = np.array(
        #     [
        #         np.outer(normal_vector, normal_vector)
        #         for normal_vector in unitary_face_normals.T
        #     ]
        # )

        # Creating a block array of identity matrices. Subtracting the tensile matrix
        # block array from this provides us with the shear matrix block array.
        eye_block_array = np.tile(np.eye(self.nd), (len(tensile_matrices), 1, 1))
        shear_mat = eye_block_array - tensile_matrices

        # Scaling with the Robin weight value
        tensile_coeff = self.robin_weight_value(direction="tensile")[boundary_cells]
        shear_coeff = self.robin_weight_value(direction="shear")[boundary_cells]

        # Constructing the total coefficient matrix.
        tensile_matrix_with_coeff = tensile_matrices * tensile_coeff[:, None, None]
        shear_matrix_with_coeff = shear_mat * shear_coeff[:, None, None]
        total_coefficient_matrix = tensile_matrix_with_coeff + shear_matrix_with_coeff
        return total_coefficient_matrix

    def robin_weight_value(self, direction: str) -> float:
        """Robin weight value for the Robin (Absorbing) boundary conditions.

        The weight in the Robin (Absorbing) boundary conditions. The Robin weight is
        different for the tensile "direction" and the shear directions of a surface.
        Either shear or tensile Robin weight will be returned by this method. This
        depends on whether shear or tensile "direction" is chosen.

        Parameters:
            direction: Whether the boundary condition that uses the weight is the shear
                or tensile component of the displacement.

        Returns:
            The weight/coefficient for use in the Robin boundary conditions.

        """
        dt = self.time_manager.dt
        rho = self.solid.density
        mu = self.mu_vector

        if direction == "shear":
            value = np.sqrt(rho * mu) * 1 / dt
        elif direction == "tensile":
            lmbda = self.lambda_vector
            value = np.sqrt(rho * (lmbda + 2 * mu)) * 1 / dt
        return value

    def boundary_displacement(self, subdomains: list[pp.Grid]) -> pp.ad.Operator:
        """Method for reconstructing the boundary displacement.

        For now the only usage is within the "realm" of absorbing boundary conditions:
        we need the displacement values on the boundary at the previous time step.

        Note: This is for the pure mechanical problem. Modifications are needed when a
            coupling to fluid flow is introduced at some later point.

        Parameters:
            subdomains: List of subdomains. Should be of co-dimension 0.

        Returns:
            Ad operator representing the displacement on grid faces of the subdomains.

        """
        # Discretization
        discr = self.stress_discretization(subdomains)

        # Boundary conditions on external boundaries
        bc = self._combine_boundary_operators(  # type: ignore [call-arg]
            subdomains=subdomains,
            dirichlet_operator=self.displacement,
            neumann_operator=self.mechanical_stress,
            robin_operator=self.mechanical_stress,
            bc_type=self.bc_type_mechanics,
            dim=self.nd,
            name=self.bc_values_mechanics_key,
        )
        # Displacement
        displacement = self.displacement(subdomains)

        boundary_displacement = (
            discr.bound_displacement_cell() @ displacement
            + discr.bound_displacement_face() @ bc
        )

        boundary_displacement.set_name("boundary_displacement")
        return boundary_displacement

    def construct_and_save_boundary_displacement(
        self, boundary_grid: pp.BoundaryGrid
    ) -> None:
        """Construct and save boundary displacement for usage in boundary conditions.

        The absorbing boundary conditions require the previous displacement values on
        the domain boundary. This method makes sure that the displacement values two
        time steps back in time is stored on the appropriate time step solution index.
        Additionally it constructs the previous boundary displacement, and saves it to
        the current time step solution index (= 0).

        Paramters:
            boundary_grid: The grid that the displacement values are evaluated on.

        """
        data = self.mdg.boundary_grid_data(boundary_grid)
        name = "boundary_displacement_values"

        # The displacement value for the previous time step is constructed and the
        # one two time steps back in time is fetched from the dictionary.
        location = pp.TIME_STEP_SOLUTIONS
        pp.shift_solution_values(
            name=name,
            data=data,
            location=location,
            max_index=2,
        )

        sd = boundary_grid.parent
        displacement_boundary_operator = self.boundary_displacement([sd])
        displacement_values = self.equation_system.evaluate(
            displacement_boundary_operator
        )

        displacement_values_0 = boundary_grid.projection(self.nd) @ displacement_values
        pp.set_solution_values(
            name=name,
            values=displacement_values_0,
            data=data,
            time_step_index=0,
        )

    def initial_condition(self) -> None:
        """Assigning initial displacement, velocity and acceleration values."""
        super().initial_condition()
        for sd, data in self.mdg.subdomains(return_data=True, dim=self.nd):
            dofs = sd.num_cells

            initial_displacement = self.initial_displacement(dofs=dofs)
            initial_velocity = self.initial_velocity(dofs=dofs)
            initial_acceleration = self.initial_acceleration(dofs=dofs)

            pp.set_solution_values(
                name=self.displacement_variable,
                values=initial_displacement,
                data=data,
                time_step_index=0,
                iterate_index=0,
            )

            pp.set_solution_values(
                name=self.velocity_key,
                values=initial_velocity,
                data=data,
                time_step_index=0,
                iterate_index=0,
            )

            pp.set_solution_values(
                name=self.acceleration_key,
                values=initial_acceleration,
                data=data,
                time_step_index=0,
                iterate_index=0,
            )

    def initial_displacement(self, dofs: int) -> np.ndarray:
        """Cell-centered initial displacement values.

        Parameters:
            dofs: Number of degrees of freedom (typically cell number in the grid where
                the initial values are defined).

        Returns:
            An array with the initial displacement values.

        """
        sd = self.mdg.subdomains()[0]
        t = self.time_manager.time

        x = sd.cell_centers[0, :]
        y = sd.cell_centers[1, :]

        vals = np.zeros((self.nd, sd.num_cells))

        if self.nd == 2:
            displacement_function = u_v_a_wrap(self)

            vals[0] = displacement_function[0](x, y, t)
            vals[1] = displacement_function[1](x, y, t)

        elif self.nd == 3:
            z = sd.cell_centers[2, :]

            displacement_function = u_v_a_wrap(self, is_2D=False)

            vals[0] = displacement_function[0](x, y, z, t)
            vals[1] = displacement_function[1](x, y, z, t)
            vals[2] = displacement_function[2](x, y, z, t)

        return vals.ravel("F")

    def initial_velocity(self, dofs: int) -> np.ndarray:
        """Cell-centered initial velocity values.

        Parameters:
            dofs: Number of degrees of freedom (typically cell number in the grid where
                the initial values are defined).

        Returns:
            An array with the initial velocity values.

        """
        sd = self.mdg.subdomains()[0]
        t = self.time_manager.time

        x = sd.cell_centers[0, :]
        y = sd.cell_centers[1, :]

        vals = np.zeros((self.nd, sd.num_cells))

        if self.nd == 2:
            velocity_function = u_v_a_wrap(self, return_dt=True)

            vals[0] = velocity_function[0](x, y, t)
            vals[1] = velocity_function[1](x, y, t)

        elif self.nd == 3:
            z = sd.cell_centers[2, :]

            velocity_function = u_v_a_wrap(self, is_2D=False, return_dt=True)

            vals[0] = velocity_function[0](x, y, z, t)
            vals[1] = velocity_function[1](x, y, z, t)
            vals[2] = velocity_function[2](x, y, z, t)

        return vals.ravel("F")

    def initial_acceleration(self, dofs: int) -> np.ndarray:
        """Cell-centered initial acceleration values.

        Parameters:
            dofs: Number of degrees of freedom (typically cell number in the grid where
                the initial values are defined).

        Returns:
            An array with the initial acceleration values.

        """
        sd = self.mdg.subdomains()[0]
        t = self.time_manager.time

        x = sd.cell_centers[0, :]
        y = sd.cell_centers[1, :]

        vals = np.zeros((self.nd, sd.num_cells))

        if self.nd == 2:
            acceleration_function = u_v_a_wrap(self, return_ddt=True)

            vals[0] = acceleration_function[0](x, y, t)
            vals[1] = acceleration_function[1](x, y, t)

        elif self.nd == 3:
            z = sd.cell_centers[2, :]

            acceleration_function = u_v_a_wrap(self, is_2D=False, return_ddt=True)

            vals[0] = acceleration_function[0](x, y, z, t)
            vals[1] = acceleration_function[1](x, y, z, t)
            vals[2] = acceleration_function[2](x, y, z, t)

        return vals.ravel("F")

    def initial_condition_bc(self, bg: pp.BoundaryGrid) -> np.ndarray:
        """Sets the initial bc values for 0th and -1st time step in the data dictionary.

        We need initial bc values two time steps back in time. This method sets the
        values for 0th and -1st time step into the data dictionary. It also returns the
        0th time step values, which is needed when the bc_robin_values is called the
        first time (at initialization).

        Parameters:
            bg: Boundary grid whose boundary displacement value is to be set.

        Returns:
            An array with the initial displacement boundary values.

        """
        vals_0 = self.initial_condition_bc_0(bg=bg)
        vals_1 = self.initial_condition_bc_1(bg=bg)

        data = self.mdg.boundary_grid_data(bg)

        name = "boundary_displacement_values"
        # The values for the 0th and -1th time step are to be stored
        pp.set_solution_values(
            name=name,
            values=vals_1,
            data=data,
            time_step_index=1,
        )
        pp.set_solution_values(
            name=name,
            values=vals_0,
            data=data,
            time_step_index=0,
        )
        return vals_0

    def initial_condition_bc_0(self, bg: pp.BoundaryGrid) -> np.ndarray:
        """Initial boundary displacement values corresponding to time step 0."""
        return np.zeros((self.nd, bg.num_cells)).ravel("F")

    def initial_condition_bc_1(self, bg: pp.BoundaryGrid) -> np.ndarray:
        """Initial boundary displacement values corresponding to time step -1."""
        return np.zeros((self.nd, bg.num_cells)).ravel("F")


class DynamicMomentumBalanceEquations:
    def momentum_balance_equation(self, subdomains: list[pp.Grid]):
        """Dynamic momentum balance equation in the rock matrix.

        Parameters:
            subdomains: List of subdomains where the force balance is defined.

        Returns:
            Operator for the force balance equation in the matrix.

        """
        inertia_mass = self.inertia(subdomains)
        stress = pp.ad.Scalar(-1) * self.stress(subdomains)
        body_force = self.body_force(subdomains)

        equation = self.balance_equation(
            subdomains, inertia_mass, stress, body_force, dim=self.nd
        )
        equation.set_name("momentum_balance_equation")
        return equation

    def inertia(self, subdomains: list[pp.Grid]) -> pp.ad.Operator:
        """Inertia mass in the elastic wave equation/dynamic momentum balance equation.

        The elastic wave equation contains a term on the form M * u_tt (that is, M *
        acceleration term/inertia term). This method represents an operator for M.

        Parameters:
            subdomains: List of subdomains where the inertia mass is defined.

        Returns:
            Operator for the inertia mass.

        """
        mass_density = self.solid_density(subdomains)
        mass = self.volume_integral(mass_density, subdomains, dim=self.nd)
        mass.set_name("inertia_mass")
        return mass

    def balance_equation(
        self,
        subdomains: list[pp.Grid],
        inertia_mass: pp.ad.Operator,
        surface_term: pp.ad.Operator,
        source: pp.ad.Operator,
        dim: int,
    ) -> pp.ad.Operator:
        """Balance equation that combines an acceleration and surface term.

        The balance equation, namely the elastic wave equation, is given by
            d_tt(displacement) + div(surface_term) - source = 0.

        Parameters:
            subdomains: List of subdomains where the balance equation is defined.
            inertia_mass: Operator for the cell-wise mass of the acceleration term,
                integrated over the cells of the subdomains.
            surface_term: Operator for the surface term (e.g. flux, stress), integrated
                over the faces of the subdomains.
            source: Operator for the source term, integrated over the cells of the
                subdomains.
            dim: Spatial dimension of the balance equation.

        Returns:
            Operator for the balance equation.

        """
        div = pp.ad.Divergence(subdomains, dim=dim)

        # Fetch the necessary operators for creating the acceleration operator
        op = self.displacement(subdomains)
        dt_op = self.velocity_time_dep_array(subdomains)
        ddt_op = self.acceleration_time_dep_array(subdomains)

        # Create the acceleration operator:
        inertia_term = time_derivatives.inertia_term(
            model=self,
            op=op,
            dt_op=dt_op,
            ddt_op=ddt_op,
            time_step=pp.ad.Scalar(self.time_manager.dt),
        )

        return inertia_mass * inertia_term + div @ surface_term - source


class SolutionStrategyDynamicMomentumBalance:
    def prepare_simulation(self) -> None:
        """Run at the start of simulation. Used for initialization etc.

        This method overrides the original prepare_simulation. The reason for this is
        that we need to create the attributes boundary_cells_of_grid, lambda_vector and
        mu_vector. All of which are needed in the call to initial_condition as well as
        set_equations.

        """
        # Set the material and geometry of the problem. The geometry method must be
        # implemented in a ModelGeometry class.
        self.set_materials()
        self.set_geometry()

        self.set_vector_valued_mu_lambda()

        # Exporter initialization must be done after grid creation,
        # but prior to data initialization.
        self.initialize_data_saving()

        # Set variables, constitutive relations, discretizations and equations.
        # Order of operations is important here.
        self.set_equation_system_manager()
        self.create_variables()
        self.initial_condition()
        self.reset_state_from_file()
        self.set_equations()

        self.set_discretization_parameters()
        self.discretize()
        self._initialize_linear_solver()
        self.set_nonlinear_discretizations()

        # Export initial condition
        self.save_data_time_step()

    def set_vector_valued_mu_lambda(self) -> None:
        """Instantiate boundary cell information and assign vector valued mu and lambda.

        Allow mu and lambda (Lamé parameters) to be assigned as cell wise values instead
        of only scalars.

        Additionally, creating the coefficient matrix for the absorbing
        boundaries requires the attribute boundary_cells_of_grid to distinguish which
        cells are bordering the domain boundary. boundary_cells_of_grid is therefore defined as an attribute of the model class setup here.

        """
        sd = self.mdg.subdomains(dim=self.nd)[0]
        boundary_faces = sd.get_boundary_faces()
        self.boundary_cells_of_grid = sd.signs_and_cells_of_boundary_faces(
            faces=boundary_faces
        )[1]
        self.vector_valued_mu_lambda()

    def set_discretization_parameters(self) -> None:
        """Set discretization parameters for the simulation.

        Sets eta = 1/3 on all faces if it is a simplex grid. Default is to have 0 on the
        boundaries, but this causes divergence of the solution. 1/3 for all subfaces in
        the grid fixes this issue.

        """
        super().set_discretization_parameters()
        if self.params["grid_type"] == "simplex":
            num_subfaces = 0
            for sd, data in self.mdg.subdomains(return_data=True, dim=self.nd):
                subcell_topology = pp.fvutils.SubcellTopology(sd)
                num_subfaces += subcell_topology.num_subfno
                eta_values = np.ones(num_subfaces) * 1 / 3
                if sd.dim == self.nd:
                    pp.initialize_data(
                        sd,
                        data,
                        self.stress_keyword,
                        {
                            "mpsa_eta": eta_values,
                        },
                    )

    def velocity_time_dep_array(
        self, subdomains: list[pp.Grid]
    ) -> pp.ad.TimeDependentDenseArray:
        """Time dependent dense array for the velocity.

        Creates a time dependent dense array to represent the velocity, which is needed
        for the Newmark time discretization.

        Parameters:
            subdomains: List of subdomains on which to define the velocity.

        Returns:
            Operator representation of the acceleration.

        """
        return pp.ad.TimeDependentDenseArray(self.velocity_key, subdomains)

    def acceleration_time_dep_array(
        self, subdomains: list[pp.Grid]
    ) -> pp.ad.TimeDependentDenseArray:
        """Time dependent dense array for the acceleration.

        Creates a time dependent dense array to represent the acceleration, which is
        needed for the Newmark time discretization.

        Parameters:
            subdomains: List of subdomains on which to define the acceleration.

        Returns:
            Operator representation of the acceleration.

        """
        return pp.ad.TimeDependentDenseArray(self.acceleration_key, subdomains)

    def velocity_values(self, subdomain: pp.Grid) -> np.ndarray:
        """Update the velocity values at the end of each time step.

        The velocity values are updated once the system is solved for the cell-centered
        displacements (at the end of a time step). The values are updated by using the
        Newmark formula for velocity (see e.g. Dynamics of Structures by A. K. Chopra
        (pp. 175-176, 2014)).

        Parameters:
            subdomain: The subdomain the velocity is defined on.

        Returns:
            An array with the new velocity values.

        """
        data = self.mdg.subdomain_data(subdomain)
        dt = self.time_manager.dt

        beta = self.beta
        gamma = self.gamma

        (
            a_previous,
            v_previous,
            u_previous,
            u_current,
        ) = acceleration_velocity_displacement(model=self, data=data)

        v = (
            v_previous * (1 - gamma / beta)
            + a_previous * (1 - gamma - (gamma * (1 - 2 * beta)) / (2 * beta)) * dt
            + (u_current - u_previous) * gamma / (beta * dt)
        )
        return v

    def acceleration_values(self, subdomain: pp.Grid) -> np.ndarray:
        """Update the acceleration values at the end of each time step.

        See the method velocity_values for more extensive documentation.

        Parameters:
            subdomain: The subdomain the acceleration is defined on.

        Returns:
            An array with the new acceleration values.

        """
        data = self.mdg.subdomain_data(subdomain)
        dt = self.time_manager.dt

        beta = self.beta

        (
            a_previous,
            v_previous,
            u_previous,
            u_current,
        ) = acceleration_velocity_displacement(model=self, data=data)

        a = (
            (u_current - u_previous) / (dt**2 * beta)
            - v_previous / (dt * beta)
            - a_previous * (1 - 2 * beta) / (2 * beta)
        )
        return a

    def update_velocity_acceleration_time_dependent_ad_arrays(self) -> None:
        """Update the time dependent ad arrays for the velocity and acceleration.

        The new velocity and acceleration values (the value at the end of each time
        step) are set into the data dictionary.

        """
        for sd, data in self.mdg.subdomains(return_data=True, dim=self.nd):
            vals_acceleration = self.acceleration_values(sd)
            vals_velocity = self.velocity_values(sd)

            pp.set_solution_values(
                name=self.velocity_key,
                values=vals_velocity,
                data=data,
                time_step_index=0,
            )
            pp.set_solution_values(
                name=self.acceleration_key,
                values=vals_acceleration,
                data=data,
                time_step_index=0,
            )

            pp.set_solution_values(
                name=self.velocity_key,
                values=vals_velocity,
                data=data,
                iterate_index=0,
            )
            pp.set_solution_values(
                name=self.acceleration_key,
                values=vals_acceleration,
                data=data,
                iterate_index=0,
            )

    def after_nonlinear_convergence(self, iteration_counter: int = 1) -> None:
        """Method to be called after every non-linear iteration.

        The method update_velocity_acceleration_time_dependent_ad_arrays needs to be
        called at the end of each time step. This is not within PorePy itself, and
        therefore this method is overriding the default after_nonlinear_convergence
        method.

        Parameters:
            solution: The new solution, as computed by the non-linear solver.
            errors: The error in the solution, as computed by the non-linear solver.
            iteration_counter: The number of iterations performed by the non-linear
                solver.

        """
        solution = self.equation_system.get_variable_values(iterate_index=0)

        self.update_velocity_acceleration_time_dependent_ad_arrays()

        self.equation_system.shift_time_step_values()
        self.equation_system.set_variable_values(
            values=solution, time_step_index=0, additive=False
        )
        if self.time_manager.time_index >= 1:
            bg = self.mdg.boundaries(dim=self.nd - 1)[0]
            self.construct_and_save_boundary_displacement(boundary_grid=bg)

        self.convergence_status = True
        self.save_data_time_step()


class ConstitutiveLawsDynamicMomentumBalance:
    def vector_valued_mu_lambda(self) -> None:
        """Vector representation of mu and lambda.

        Cell-wise representation of the mu and lambda quantities in the rock matrix.
        Both are saved as an attribute of the model class setup.

        """
        subdomain = self.mdg.subdomains(dim=self.nd)[0]

        self.lambda_vector = self.solid.lame_lambda * np.ones(subdomain.num_cells)
        self.mu_vector = self.solid.shear_modulus * np.ones(subdomain.num_cells)

    def stiffness_tensor(self, subdomain: pp.Grid) -> pp.FourthOrderTensor:
        """Stiffness tensor [Pa].

        Overriding the stiffness_tensor method to accomodate vector valued mu and
        lambda. Vector valued mu and lambda is treated as default in the model class for
        MPSA-Newmark with absorbing boundaries.

        Parameters:
            subdomain: Subdomain where the stiffness tensor is defined.

        Returns:
            Cell-wise stiffness tensor in SI units.

        """
        return pp.FourthOrderTensor(self.mu_vector, self.lambda_vector)


class TimeDependentSourceTerm:
    def body_force(self, subdomains: list[pp.Grid]) -> pp.ad.Operator:
        """Time dependent dense array for the body force term.

        Creates a time dependent dense array to represent the body force/source.

        Parameters:
            subdomains: List of subdomains on which to define the body force.

        Returns:
            Operator representation of the body force.

        """
        external_sources = pp.ad.TimeDependentDenseArray(
            name="source_mechanics",
            domains=subdomains,
        )
        return external_sources

    def before_nonlinear_loop(self) -> None:
        """Update the time dependent mechanics source before every nonlinear loop.

        The default behaviour is to save source values which correspond to a source
        which is constructed from a known analytical reference solution. In the case of
        no analytical reference solution being chosen (this is done in the params
        dictionary), the source is just zero.

        """
        super().before_nonlinear_loop()

        sd = self.mdg.subdomains(dim=self.nd)[0]
        data = self.mdg.subdomain_data(sd)
        t = self.time_manager.time

        # Mechanics source
        if self.nd == 2:
            mechanics_source_function = body_force_function(self)
        elif self.nd == 3:
            mechanics_source_function = body_force_function(self, is_2D=False)

        mechanics_source_values = self.evaluate_mechanics_source(
            f=mechanics_source_function, sd=sd, t=t
        )
        pp.set_solution_values(
            name="source_mechanics",
            values=mechanics_source_values,
            data=data,
            iterate_index=0,
        )

    def evaluate_mechanics_source(self, f: list, sd: pp.Grid, t: float) -> np.ndarray:
        """Computes the values for the body force.

        The method computes the source values returned by the source value function (f)
        integrated over the cell. The function is evaluated at the cell centers.

        Parameters:
            f: Function expression for the source term. It depends on time and space for
                the source term. It is represented as a list, where the first list
                component corresponds to the first vector component of the source,
                second list component corresponds to the second vector component, and so
                on.
            sd: Subdomain where the source term is defined.
            t: Current time in the time-stepping.

        Returns:
            An array of source values.

        """
        cell_volume = sd.cell_volumes
        vals = np.zeros((self.nd, sd.num_cells))

        x = sd.cell_centers[0, :]
        y = sd.cell_centers[1, :]

        if self.nd == 2:
            x_val = f[0](x, y, t)
            y_val = f[1](x, y, t)

        elif self.nd == 3:
            z = sd.cell_centers[2, :]

            x_val = f[0](x, y, z, t)
            y_val = f[1](x, y, z, t)
            z_val = f[2](x, y, z, t)

            vals[2] = z_val * cell_volume

        vals[0] = x_val * cell_volume
        vals[1] = y_val * cell_volume
        return vals.ravel("F")


class DynamicMomentumBalanceABC(
    NamesAndConstants,
    CustomSolverMixin,
    BoundaryAndInitialConditions,
    DynamicMomentumBalanceEquations,
    ConstitutiveLawsDynamicMomentumBalance,
    TimeDependentSourceTerm,
    SolutionStrategyDynamicMomentumBalance,
    MomentumBalance,
):
    """Full model class for the dynamic momentum balance with absorbing boundaries."""
