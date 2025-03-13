"""Model class setup for the convergence analysis of MPSA-Newmark with absorbing
boundaries.

All off-diagonal/shear components of the stiffness tensor are discarded.

"""

import sys

import numpy as np
import porepy as pp

sys.path.append("../../")

from .model_convergence_ABC import ABCModel
from porepy.applications.convergence_analysis import ConvergenceAnalysis
import sympy as sym


class BoundaryConditions:
    def bc_values_displacement(self, bg: pp.BoundaryGrid) -> np.ndarray:
        """Method for setting Dirichlet boundary values.

        Sets a time dependent sine condition in the x-direction of the western boundary.
        Zero elsewhere.

        Parameters:
            bg: Boundary grid whose boundary displacement value is to be set.

        Returns:
            An array with the displacement boundary values at time t.

        """
        values = np.zeros((self.nd, bg.num_cells))
        bounds = self.domain_boundary_sides(bg)
        t = self.time_manager.time

        xmin = self.domain.bounding_box["xmin"]

        # Time dependent sine Dirichlet condition
        bc_left, _ = self.heterogeneous_analytical_solution()
        values[0][bounds.west] += np.ones(len(values[0][bounds.west])) * bc_left[0](
            xmin, t
        )

        return values.ravel("F")

    def initial_condition_value_function(
        self, bg: pp.BoundaryGrid, t: float
    ) -> np.ndarray:
        """Initial values for the absorbing boundary."""
        sd = bg.parent

        x = sd.face_centers[0, :]

        boundary_sides = self.domain_boundary_sides(sd)
        inds_east = np.where(boundary_sides.east)[0]

        bc_vals = np.zeros((sd.dim, sd.num_faces))

        _, displacement_function_left = self.heterogeneous_analytical_solution()

        # East
        bc_vals[0, :][inds_east] = displacement_function_left[0](x[inds_east], t)

        bc_vals = bg.projection(self.nd) @ bc_vals.ravel("F")
        return bc_vals


class InitialConditions:
    def heterogeneous_analytical_solution(
        self, return_dt=False, return_ddt=False, lambdify=True
    ):
        """Compute the analytical solution and its time derivatives."""
        x, t = sym.symbols("x t")

        L = self.heterogeneity_location
        cp = self.primary_wave_speed(is_scalar=False)

        heterogeneity_factor = self.heterogeneity_factor
        if heterogeneity_factor >= 1.0:
            left_speed = min(cp)
            right_speed = max(cp)
        else:
            left_speed = max(cp)
            right_speed = min(cp)

        u_left = sym.sin(t - (x - L) / left_speed) - (right_speed - left_speed) / (
            right_speed + left_speed
        ) * sym.sin(t + (x - L) / left_speed)
        u_right = (
            (2 * left_speed)
            / (right_speed + left_speed)
            * sym.sin(t - (x - L) / right_speed)
        )

        # Compute derivatives based on function arguments
        if return_dt:
            u_left, u_right = sym.diff(u_left, t), sym.diff(u_right, t)
        elif return_ddt:
            u_left, u_right = sym.diff(u_left, t, 2), sym.diff(u_right, t, 2)

        if lambdify:
            return [sym.lambdify((x, t), u_left, "numpy"), 0], [
                sym.lambdify((x, t), u_right, "numpy"),
                0,
            ]
        else:
            return x, t, u_left, u_right

    def _compute_initial_condition(self, return_dt=False, return_ddt=False):
        """Helper function to compute displacement, velocity or acceleration."""
        sd = self.mdg.subdomains()[0]
        x = sd.cell_centers[0, :]
        t = self.time_manager.time

        L = self.heterogeneity_location
        left_layer = x <= L
        right_layer = x > L

        vals = np.zeros((self.nd, sd.num_cells))

        left_solution, right_solution = self.heterogeneous_analytical_solution(
            return_dt=return_dt, return_ddt=return_ddt
        )

        vals[0, left_layer] = left_solution[0](x[left_layer], t)
        vals[0, right_layer] = right_solution[0](x[right_layer], t)
        return vals.ravel("F")

    def initial_displacement(self, dofs):
        """Compute the initial displacement."""
        return self._compute_initial_condition()

    def initial_velocity(self, dofs):
        """Compute the initial velocity."""
        return self._compute_initial_condition(return_dt=True)

    def initial_acceleration(self, dofs):
        """Compute the initial acceleration."""
        return self._compute_initial_condition(return_ddt=True)


class ConstitutiveLawsAndSource:
    def evaluate_mechanics_source(self, f: list, sd: pp.Grid, t: float) -> np.ndarray:
        vals = np.zeros((self.nd, sd.num_cells))
        return vals.ravel("F")

    def vector_valued_mu_lambda(self) -> None:
        """Setting a vertically layered medium."""
        subdomain = self.mdg.subdomains(dim=self.nd)[0]
        x = subdomain.cell_centers[0, :]

        lmbda1 = self.solid.lame_lambda
        mu1 = self.solid.shear_modulus

        lmbda2 = self.solid.lame_lambda * self.heterogeneity_factor
        mu2 = self.solid.shear_modulus * self.heterogeneity_factor

        lmbda_vec = np.ones(subdomain.num_cells)
        mu_vec = np.ones(subdomain.num_cells)

        left_layer = x <= self.heterogeneity_location
        right_layer = x > self.heterogeneity_location

        lmbda_vec[left_layer] *= lmbda1
        mu_vec[left_layer] *= mu1

        lmbda_vec[right_layer] *= lmbda2
        mu_vec[right_layer] *= mu2

        self.mu_vector = mu_vec
        self.lambda_vector = lmbda_vec


class ExactHeterogeneousSigmaAndForce:
    def exact_heterogeneous_sigma(self, u, lam, mu, x) -> list:
        """Representation of the exact stress tensor given a displacement field.

        Parameters:
            u: The sympy representation of the exact displacement.
            lam: The first Lamé parameter.
            mu: The second Lamé parameter, also called shear modulus.
            x: The sympy representation of the x-coordinate.

        Returns:
            A list which represents the sympy expression of the exact stress tensor.

        """
        y = sym.symbols("y")

        u = [u, 0]
        # Exact gradient of u and transpose of gradient of u
        grad_u = [
            [sym.diff(u[0], x), sym.diff(u[0], y)],
            [sym.diff(u[1], x), sym.diff(u[1], y)],
        ]

        grad_u_T = [[grad_u[0][0], grad_u[1][0]], [grad_u[0][1], grad_u[1][1]]]

        # Trace of gradient of u, in the linear algebra sense
        trace_grad_u = grad_u[0][0] + grad_u[1][1]

        # Exact strain (\epsilon(u))
        strain = 0.5 * np.array(
            [
                [grad_u[0][0] + grad_u_T[0][0], grad_u[0][1] + grad_u_T[0][1]],
                [grad_u[1][0] + grad_u_T[1][0], grad_u[1][1] + grad_u_T[1][1]],
            ]
        )

        # Exact stress tensor (\sigma(\epsilon(u)))
        sigma = [
            [2 * mu * strain[0][0] + lam * trace_grad_u, 2 * mu * strain[0][1]],
            [2 * mu * strain[1][0], 2 * mu * strain[1][1] + lam * trace_grad_u],
        ]
        return sigma

    def evaluate_exact_force(self, sd, time, sigma, inds, force_array) -> np.ndarray:
        """Evaluate exact elastic force at the face centers for certain face indices.

        This method is, as it is implemented now, called more than once. This is because
        the exact elastic force values may be different in two or more parts of the
        subdomain. In the case of two parts of the subdomain, the method fills half the
        `force_array` in the first call, and then the other half in the second call. The
        filling of the array is done index wise, determined by `inds`.

        Parameters:
            sd: Subdomain grid.
            time: Time in seconds.
            sigma: Exact stress tensor.
            inds: The face indices which we are computing the force at.
            force_array: Either empty or semi empty array of shape (self.nd, sd.
                num_faces) which we are filling with the force values. This is done by
                indices in `inds`.

        Returns:
            Array of `shape=(self.nd, sd.num_faces)` containing the exact elastic
            force at the face centers of `inds` for the given `time`.

        Notes:
            - The returned elastic force is _not_ given in PorePy's flattened vector
              format. Thus, it may be necessary to flatten it at a later point.
            - Recall that force = (stress dot_prod unit_normal) * face_area.

        """
        # Symbolic variables
        x, y, t = sym.symbols("x y t")

        fc = sd.face_centers[:, inds].squeeze()
        fn = sd.face_normals[:, inds].squeeze()

        # Lambdify expression
        sigma_total_fun = [
            [
                sym.lambdify((x, y, t), sigma[0][0], "numpy"),
                sym.lambdify((x, y, t), sigma[0][1], "numpy"),
            ],
            [
                sym.lambdify((x, y, t), sigma[1][0], "numpy"),
                sym.lambdify((x, y, t), sigma[1][1], "numpy"),
            ],
        ]

        # Face-centered elastic force
        force_total_fc: list[np.ndarray] = [
            # (sigma_xx * n_x + sigma_xy * n_y) * face_area
            sigma_total_fun[0][0](fc[0], fc[1], time) * fn[0]
            + sigma_total_fun[0][1](fc[0], fc[1], time) * fn[1],
            # (sigma_yx * n_x + sigma_yy * n_y) * face_area
            sigma_total_fun[1][0](fc[0], fc[1], time) * fn[0]
            + sigma_total_fun[1][1](fc[0], fc[1], time) * fn[1],
        ]

        # Insert values into force_array at the indices given by inds
        for i, force_component in enumerate(force_total_fc):
            # Update the appropriate row of force_array
            force_array[i, inds] = force_component

        return force_array

    def evaluate_exact_heterogeneous_force(self, sd) -> np.ndarray:
        """Evaluate the exact heterogeneous force in the entire domain.

        The domain is split in 2: One left and one right part, where the exact force may
        be different in each part of the domain. This method handles computing the exact
        force values for the entire domain.

        Parameters:
            sd: The subdomain grid where the forces are to be evaluated.

        Returns:
            A flattened array of the exact force values in the entire domain.

        """
        x, _, u_left, u_right = self.heterogeneous_analytical_solution(lambdify=False)

        mu_lambda_values = {
            "left": (self.solid.lame_lambda, self.solid.shear_modulus),
            "right": (
                self.solid.lame_lambda * self.heterogeneity_factor,
                self.solid.shear_modulus * self.heterogeneity_factor,
            ),
        }

        sigma = {}
        for side, (lam, mu) in mu_lambda_values.items():
            u = u_left if side == "left" else u_right
            sigma[side] = self.exact_heterogeneous_sigma(u, lam, mu, x)

        fc_x = sd.face_centers[0, :]

        empty_force_array = np.zeros((self.nd, sd.num_faces))

        inds_left = np.where(fc_x < self.heterogeneity_location)
        inds_right = np.where(fc_x >= self.heterogeneity_location)

        semi_full_force_array = self.evaluate_exact_force(
            sd, self.time_manager.time, sigma["left"], inds_left, empty_force_array
        )
        full_force_array = self.evaluate_exact_force(
            sd,
            self.time_manager.time,
            sigma["right"],
            inds_right,
            semi_full_force_array,
        )
        full_force_array: np.ndarray = np.asarray(full_force_array).ravel("F")
        return full_force_array


class ComputeErrorsHeterogeneity:
    def compute_and_save_errors(self, filename: str) -> None:
        sd = self.mdg.subdomains(dim=self.nd)[0]

        x = sd.cell_centers[0, :]
        L = self.heterogeneity_location

        left_solution, right_solution = self.heterogeneous_analytical_solution()
        left_layer = x <= L
        right_layer = x > L

        vals = np.zeros((self.nd, sd.num_cells))
        vals[0, left_layer] = left_solution[0](x[left_layer], self.time_manager.time)
        vals[0, right_layer] = right_solution[0](x[right_layer], self.time_manager.time)

        displacement_ad = self.displacement([sd])
        u_approximate = self.equation_system.evaluate(displacement_ad)
        exact_displacement = vals.ravel("F")

        exact_force = self.evaluate_exact_heterogeneous_force(sd=sd)
        force_ad = self.stress([sd])
        approx_force = self.equation_system.evaluate(force_ad)

        error_displacement = ConvergenceAnalysis.lp_error(
            grid=sd,
            true_array=exact_displacement,
            approx_array=u_approximate,
            is_scalar=False,
            is_cc=True,
            relative=True,
        )
        error_traction = ConvergenceAnalysis.lp_error(
            grid=sd,
            true_array=exact_force,
            approx_array=approx_force,
            is_scalar=False,
            is_cc=False,
            relative=True,
        )
        with open(filename, "a") as file:
            file.write(
                f"{sd.num_cells}, {self.time_manager.time_index}, {error_displacement}, {error_traction}\n"
            )


class ABCModelHeterogeneous(
    BoundaryConditions,
    ConstitutiveLawsAndSource,
    InitialConditions,
    ExactHeterogeneousSigmaAndForce,
    ComputeErrorsHeterogeneity,
    ABCModel,
): ...
