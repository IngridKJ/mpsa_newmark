"""File for anisotropic stiffness tensor mixin."""

import numpy as np
import porepy as pp
from utils.utility_functions import (
    use_constraints_for_inner_domain_cells,
    create_stiffness_tensor_basis,
)


class TransverselyIsotropicTensorMixin:
    """Mixin which defines a tensor with (possibly) heterogeneous properties.

    The tensor may be assigned to be anisotropic in the regions identified by
    use_constraints_for_inner_domain_cells(). The default parameter values for the
    anisotropic region is however isotropic, so this tensor may also be used in fully
    isotropic cases.

    NOTE: The utility function use_constraints_for_inner_domain_cells() used hereing
    assumes existence of "constraints" within the domain. For simplex grids that means
    constraints which influence the meshing, while for Cartesian grids it means polygons
    which are just used to determine what is part of the "inner" domain and what is part
    of the "outer" domain. The constriants are set by set_fractures() and the polygons
    are set by set_polygons(). See the function use_constraints_for_inner_domain_cells()
    and its documentation for further details.

    """

    def stiffness_tensor(self, subdomain: pp.Grid) -> pp.FourthOrderTensor:
        """Compute the stiffness tensor for a given subdomain.

        Parameters:
            subdomain: The subdomain we compute the stiffness tensor for.

        Returns:
            The stiffness tensor.

        """
        # Fetch inner domain indices
        inner_cell_indices = use_constraints_for_inner_domain_cells(
            model=self,
            sd=subdomain,
        )

        # Preparing basis arrays for inner and outer domains
        inner = np.zeros(subdomain.num_cells)
        inner[inner_cell_indices] = 1

        outer = np.ones(subdomain.num_cells)
        outer = outer - inner

        # Compute the stiffness tensors for each term using the extracted constants
        stiffness_matrices = create_stiffness_tensor_basis(
            lambda_val=1,
            lambda_parallel=1,
            lambda_perpendicular=1,
            mu_parallel=1,
            mu_perpendicular=1,
            n=self.params.get("symmetry_axis", [0, 0, 1]),
        )

        # Extract stiffness matrices for each component
        lmbda_mat = stiffness_matrices["lambda"]
        lambda_parallel_mat = stiffness_matrices["lambda_parallel"]
        lambda_orthogonal_mat = stiffness_matrices["lambda_perpendicular"]
        mu_parallel_mat = stiffness_matrices["mu_parallel"]
        mu_orthogonal_mat = stiffness_matrices["mu_perpendicular"]

        # Extract individual constants from the anisotropy constants dictionary. The
        # default here is set according to an isotropic tensor.
        anisotropy_constants = self.params.get(
            "anisotropy_constants",
            {
                "mu_parallel": self.solid.shear_modulus,
                "mu_orthogonal": self.solid.shear_modulus,
                "lambda_parallel": 0.0,
                "lambda_orthogonal": 0.0,
                "volumetric_compr_lambda": self.solid.lame_lambda,
            },
        )

        volumetric_compr_lambda = anisotropy_constants["volumetric_compr_lambda"]
        mu_parallel = anisotropy_constants["mu_parallel"]
        mu_orthogonal = anisotropy_constants["mu_orthogonal"]
        lambda_parallel = anisotropy_constants["lambda_parallel"]
        lambda_orthogonal = anisotropy_constants["lambda_orthogonal"]

        # Standard material values: assigned to the outer domain
        lmbda = self.solid.lame_lambda * outer
        mu = self.solid.shear_modulus * outer

        # Values for inner domain with anisotropic constants
        mu_parallel_inner = mu_parallel * inner
        mu_orthogonal_inner = mu_orthogonal * inner
        volumetric_compr_lambda_inner = volumetric_compr_lambda * inner
        lambda_parallel_inner = lambda_parallel * inner
        lambda_orthogonal_inner = lambda_orthogonal * inner

        # Create the final stiffness tensor
        stiffness_tensor = pp.FourthOrderTensor(
            mu=mu,
            lmbda=lmbda,
            other_fields={
                "mu_parallel": (mu_parallel_mat, mu_parallel_inner),
                "mu_orthogonal": (mu_orthogonal_mat, mu_orthogonal_inner),
                "volumetric_compr_lambda": (lmbda_mat, volumetric_compr_lambda_inner),
                "lambda_parallel": (lambda_parallel_mat, lambda_parallel_inner),
                "lambda_orthogonal": (lambda_orthogonal_mat, lambda_orthogonal_inner),
            },
        )
        return stiffness_tensor
