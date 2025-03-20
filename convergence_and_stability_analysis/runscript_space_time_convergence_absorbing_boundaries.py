import os
import sys

sys.path.append("../")

from itertools import product

import matplotlib.pyplot as plt
import numpy as np
import porepy as pp
import run_models.run_linear_model as rlm
from plotting.plot_utils import draw_multiple_loglog_slopes

from convergence_and_stability_analysis.analysis_models.model_convergence_ABC import (
    ABCModel,
)

# Prepare path for generated output files
folder_name = "convergence_analysis_results"
header = "num_cells, num_time_steps, displacement_error, traction_error\n"
script_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(script_dir, folder_name)
os.makedirs(output_dir, exist_ok=True)

# Run three refinement levels (coarse = True) or five refinement levels (coarse =
# False), and choose whether the figure should be saved (True) or not (False).
coarse = True
save_figure = True


class ConvergenceAnalysisHeterogeneity(ABCModel):
    def meshing_arguments(self) -> dict:
        cell_size = self.units.convert_units(0.125 / 2 ** (self.refinement), "m")
        mesh_args: dict[str, float] = {"cell_size": cell_size}
        return mesh_args

    def set_geometry(self):
        """Perturb all internal boundary nodes randomly to ensure an unstructured
        grid."""

        # Choose a seed for reproducibility.
        np.random.seed(42)
        super().set_geometry()

        sd = self.mdg.subdomains()[0]
        h = self.meshing_arguments()["cell_size"]

        inds = sd.get_internal_nodes()
        inds_not_constraint = np.where(sd.nodes[0, inds] != self.heterogeneity_location)
        inds = inds[inds_not_constraint]

        perturbation = 0.1 * h
        signs = np.random.choice([-1, 1], size=len(inds))
        sd.nodes[:2, inds] += perturbation * signs

        sd.compute_geometry()

    def write_pvd_and_vtu(self) -> None:
        """Override method such that pvd and vtu files are not created."""
        self.data_to_export()


heterogeneity_coefficients = [
    # Homogeneous case
    1,
    # # Heterogeneous case
    1 / 2**4,
    # # Heterogeneous case
    1 / 2**8,
]
anisotropy_coefficients = [
    # Isotropic case
    0,
    # Anisotropic case
    1e1,
]

for heterogeneity_factor_index in range(0, len(heterogeneity_coefficients)):
    r_h = heterogeneity_coefficients[heterogeneity_factor_index]
    for r_a in anisotropy_coefficients:
        h_lambda_ind = anisotropy_coefficients.index(r_a)
        filename = f"errors_heterogeneity_{str(heterogeneity_factor_index)}_mu_lam_{str(h_lambda_ind)}.txt"

        filename = os.path.join(output_dir, filename)

        refinements = np.arange(2, 5) if coarse else np.arange(2, 7)
        for refinement_coefficient in refinements:
            if refinement_coefficient == 2:
                with open(filename, "w") as file:
                    file.write(header)
            tf = 15.0
            time_steps = 15 * (2**refinement_coefficient)
            dt = tf / time_steps

            time_manager = pp.TimeManager(
                schedule=[0.0, tf],
                dt_init=dt,
                constant_dt=True,
            )

            lame_lambda, shear_modulus = 0.01, 0.01

            solid_constants = pp.SolidConstants(
                lame_lambda=lame_lambda, shear_modulus=shear_modulus
            )
            material_constants = {"solid": solid_constants}

            anisotropy_constants = {
                "mu_parallel": shear_modulus,
                "mu_orthogonal": shear_modulus,
                "lambda_parallel": 0.0,
                "lambda_orthogonal": lame_lambda * r_a,
                "volumetric_compr_lambda": lame_lambda,
            }

            params = {
                "time_manager": time_manager,
                "grid_type": "simplex",
                "progressbars": True,
                "heterogeneity_factor": r_h,
                "heterogeneity_location": 0.5,
                "material_constants": material_constants,
                "anisotropy_constants": anisotropy_constants,
                "symmetry_axis": [0, 1, 0],
                "meshing_kwargs": {"constraints": [0, 1, 2, 3]},
                "run_type": "vertical_anisotropy",
            }

            model = ConvergenceAnalysisHeterogeneity(params)
            model.refinement = refinement_coefficient
            model.filename_path = filename
            rlm.run_linear_model(model, params)


# Plotting from here and down:
if save_figure:
    # Set paths for modules and data directories
    folder_name = "convergence_analysis_results"
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, folder_name)
    os.makedirs(output_dir, exist_ok=True)

    # Define legend labels using a dictionary
    label_dict = {
        (0, 0, 0): r"$r_h = 1, r_{a} = 0$",
        (1, 0, 0): r"$r_h = 2^{-4}, r_{a} = 0$",
        (2, 0, 0): r"$r_h = 2^{-8}, r_{a} = 0$",
        (0, 1, 1): r"$r_h = 1, r_{a} = 10$",
        (1, 1, 1): r"$r_h = 2^{-4}, r_{a} = 10$",
        (2, 1, 1): r"$r_h = 2^{-8}, r_{a} = 10$",
    }

    # Get list of heterogeneity error files
    heterogeneity_files = [
        f
        for f in os.listdir(output_dir)
        if f.startswith("errors_heterogeneity") and f.endswith(".txt")
    ]

    # Parse file info and extract unique factors and anisotropy coefficients
    file_info = []
    factors = set()
    anisotropy_pairs = set()

    for filename in heterogeneity_files:
        parts = (
            filename.replace("errors_heterogeneity_", "").replace(".txt", "").split("_")
        )
        factor, mu_lam_index = float(parts[0]), int(parts[-1])
        file_info.append((factor, mu_lam_index, mu_lam_index, filename))
        factors.add(factor)
        anisotropy_pairs.add((mu_lam_index, mu_lam_index))

    # Sort the file info and factor/aniso lists
    file_info.sort()
    factors = sorted(factors)
    anisotropy_pairs = sorted(anisotropy_pairs)

    # Define custom styles for plotting
    custom_styles = {
        # Every other color in "displacement" and "traction" is black/grey. The others
        # are colors (pink, green, blue), where traction has a darker color and
        # displacement a lighter one.
        "displacement": ["#D8B5DE", "black", "#9FE6A2", "black", "#55A1FF", "black"],
        "traction": [
            "#A45892",
            "darkgray",
            "#01840C",
            "darkgray",
            "#01308E",
            "darkgray",
        ],
        "linestyles": ["-", ":", "-", ":", "-", ":"],
        "markers": ["o", "s", "o", "D", "o", "d"],
        "alpha": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    }

    # Map styles to each combination of factor and anisotropy
    style_map = {}
    for i, combo in enumerate(product(factors, anisotropy_pairs)):
        style_map[combo] = {
            "color_displacement": custom_styles["displacement"][
                i % len(custom_styles["displacement"])
            ],
            "color_traction": custom_styles["traction"][
                i % len(custom_styles["traction"])
            ],
            "linestyle": custom_styles["linestyles"][
                i % len(custom_styles["linestyles"])
            ],
            "marker": custom_styles["markers"][i % len(custom_styles["markers"])],
            "alpha": custom_styles["alpha"][i % len(custom_styles["alpha"])],
        }

    # Create figure for plotting
    fig, ax = plt.subplots(figsize=(10, 9))

    # Initialize lists for handles and labels
    handles_u, labels_u, handles_T, labels_T = [], [], [], []

    # Plot data for each file
    for factor, mu, lam, filename in file_info:
        filepath = os.path.join(output_dir, filename)
        num_cells, num_time_steps, displacement_errors, traction_errors = np.loadtxt(
            filepath, delimiter=",", skiprows=1, unpack=True, dtype=float
        )

        x_axis = (num_cells * num_time_steps) ** (1 / 4)
        style = style_map[(factor, (mu, lam))]
        label_name = label_dict.get(
            (factor, mu, lam), f"Unknown Case ({factor},{mu},{lam})"
        )

        # Plot displacement error
        (line_u,) = ax.loglog(
            x_axis,
            displacement_errors,
            linestyle=style["linestyle"],
            marker=style["marker"],
            color=style["color_displacement"],
            alpha=style["alpha"],
            markersize=8 if style["marker"] != "o" else 16,
            linewidth=5,
        )

        # Plot traction error
        (line_T,) = ax.loglog(
            x_axis,
            traction_errors,
            linestyle=style["linestyle"],
            marker=style["marker"],
            color=style["color_traction"],
            markersize=8 if style["marker"] != "o" else 16,
            linewidth=5,
        )

        handles_u.append(line_u)
        handles_T.append(line_T)

        # Draw convergence slope indicator for the last factor and anisotropy pair
        if len(x_axis) == 5:
            if factor == factors[-1] and (mu, lam) == anisotropy_pairs[-1]:
                draw_multiple_loglog_slopes(
                    fig,
                    ax,
                    origin=(1.05 * x_axis[-2], 1.05 * displacement_errors[-2]),
                    triangle_width=1.45,
                    slopes=[-2],
                    inverted=False,
                )
        # If less than five refinements is performed, the triangle is located elsewhere
        # to not overlap with the legend. Then the triangle is based on the first
        # heterogeneity factor and anisotropy pair.
        elif len(x_axis) < 5:
            if factor == factors[0] and (mu, lam) == anisotropy_pairs[0]:
                draw_multiple_loglog_slopes(
                    fig,
                    ax,
                    origin=(0.95 * x_axis[-1], 0.95 * displacement_errors[-1]),
                    triangle_width=1.45,
                    slopes=[-2],
                    inverted=True,
                )

    # Create 7 additional white lines and add them to the legend. This is done to ease
    # further customization with the legend.
    invisible_lines = [plt.Line2D([0], [0], color="white") for _ in range(7)]

    # These are the labels which are common for the displacement and traction errors for
    # each of the simulation runs:
    common_labels = [
        "",
        r"$r_h = 1$" + ",     " + r"$ r_{a} = 0$",
        r"$r_h = 1$" + ",     " + r"$ r_{a} = 10$",
        r"$r_h = 2^{-4}$" + ",  " + r"$r_{a} = 0$",
        r"$r_h = 2^{-4}$" + ",  " + r"$r_{a} = 10$",
        r"$r_h = 2^{-8}$" + ",  " + r"$ r_{a} = 0$",
        r"$r_h = 2^{-8}$" + ",  " + r"$ r_{a} = 10$",
    ]

    # Configure plot labels, title, grid and ticks
    ax.set_xlabel(r"$(N_x \cdot N_t)^{1/4}$", fontsize=16)
    ax.set_ylabel("Relative $L^2$ error", fontsize=16)
    ax.set_title("Convergence analysis results", fontsize=22)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    ax.xaxis.set_tick_params(which="both", labelsize=15)
    ax.yaxis.set_tick_params(which="both", labelsize=15)
    ax.set_ylim(top=2e1)

    # Create custom legend with headers for displacement and traction errors
    handles_u_header = plt.Line2D([0], [0], color="white", linewidth=0)
    handles_T_header = plt.Line2D([0], [0], color="white", linewidth=0)
    labels_u_header, labels_T_header = (
        r"$\mathcal{E}_u$",
        r"$\mathcal{E}_T$",
    )

    # Modify labels_u to have empty strings for the first column
    labels_u_empty = [""] * 6

    # Create the legend, adjust spacing and center the headers
    handles = (
        [handles_u_header]
        + handles_u
        + [handles_T_header]
        + handles_T
        + invisible_lines
    )
    labels = (
        [labels_u_header]
        + labels_u_empty
        + [labels_T_header]
        + labels_u_empty
        + common_labels
    )

    # Legend for the plot. The fontsize and handleheight is slightly smaller when we have
    # less than 5 data points/refinements.
    leg = ax.legend(
        handles,
        labels,
        fontsize=13 if len(x_axis) < 5 else 15,
        loc="upper right",
        bbox_to_anchor=(0.999, 0.999),
        ncol=3,
        frameon=True,
        handleheight=1.5 if len(x_axis) < 5 else 1.7,
        handlelength=2.2,
        columnspacing=0.8,
        labelspacing=0.275,
    )

    # Adjust alignment and further refine the positioning of legend entries
    for vpack in leg._legend_handle_box.get_children():
        for hpack in vpack.get_children()[:1]:
            hpack.get_children()[0].set_width(0)

    for vpack in leg._legend_handle_box.get_children()[:1]:
        for hpack in vpack.get_children():
            hpack.get_children()[0].set_width(0)

    for vpack in leg._legend_handle_box.get_children()[2:3]:
        for hpack in vpack.get_children():
            handle = hpack.get_children()[0]
            handle.set_visible(False)

    # Save the plot to a file
    figures_folder = "figures"
    figures_output_dir = os.path.join(script_dir, figures_folder)
    os.makedirs(figures_output_dir, exist_ok=True)

    plt.savefig(
        os.path.join(
            figures_output_dir,
            "space_time_convergence_analysis_absorbing_boundaries.png",
        ),
        bbox_inches="tight",
    )
