import sys
from copy import deepcopy

import numpy as np
import porepy as pp
from porepy.utils.txt_io import TxtData, export_data_to_txt

sys.path.append("../")

import run_models.run_linear_model as rlm


def run_analysis(self) -> list:
    """Run convergence analysis.

    See the original method (found in porepy.applications.convergence_analysis) for more
    detailed documentation. This function is a copy of the PorePy-one, with the
    extension to also append the number of cells to the convergence results.

    Returns:
        List of results (i.e., data classes containing the errors) for each
        refinement level.

    """
    convergence_results: list = []
    for level in range(self.levels):
        setup = self.model_class(deepcopy(self.model_params[level]))
        rlm.run_linear_model(setup, deepcopy(self.model_params[level]))

        setattr(setup.results[-1], "cell_diameter", setup.mdg.diameter())
        setattr(setup.results[-1], "dt", setup.time_manager.dt)
        setattr(
            setup.results[-1],
            "num_cells",
            setup.mdg.subdomains(dim=setup.nd)[0].num_cells,
        )

        convergence_results.append(setup.results[-1])
    return convergence_results


def export_errors_to_txt(
    self,
    list_of_results: list,
    variables_to_export=None,
    file_name="error_analysis.txt",
) -> None:
    """Write errors into a ``txt`` file.

    See the original method (found in porepy.applications.convergence_analysis) for more
    detailed documentation. This function is a copy of the PorePy-one, with the
    extension to also export number of cells in the domain.

    """
    # Filter variables from the list of results
    var_names: list[str] = self._filter_variables_from_list_of_results(
        list_of_results=list_of_results,
        variables=variables_to_export,
    )

    # Filter errors to be exported
    errors_to_export: dict[str, np.ndarray] = {}
    for name in var_names:
        # Loop over lists of results
        var_error: list[float] = []
        for result in list_of_results:
            var_error.append(getattr(result, name))
        # Append to the dictionary
        errors_to_export[name] = np.array(var_error)

    # Prepare to export
    list_of_txt_data: list[TxtData] = []
    # Append cell diameters
    cell_diameters = np.array([result.cell_diameter for result in list_of_results])
    list_of_txt_data.append(
        TxtData(
            header="cell_diameter",
            array=cell_diameters,
            format=self._set_column_data_format("cell_diameter"),
        )
    )

    # Append cell number
    cell_diameters = np.array([result.num_cells for result in list_of_results])
    list_of_txt_data.append(
        TxtData(
            header="num_cells",
            array=cell_diameters,
            # Want integer value for cell number
            format="%d",
        )
    )

    time_steps = np.array([result.dt for result in list_of_results])
    list_of_txt_data.append(
        TxtData(
            header="time_step",
            array=time_steps,
            format=self._set_column_data_format("time_step"),
        )
    )

    for key in errors_to_export.keys():
        list_of_txt_data.append(
            TxtData(
                header=key,
                array=errors_to_export[key],
                format=self._set_column_data_format(key),
            )
        )

    # Finally, call the function to write into the txt
    export_data_to_txt(list_of_txt_data, file_name)


def traction_error_volume_weight(
    sd: pp.GridLike,
    true_array: np.ndarray,
    approx_array: np.ndarray,
    is_scalar: bool,
    relative: bool = False,
) -> pp.number:
    """To be documented."""
    face_nodes = sd.face_nodes
    all_meas = np.zeros(sd.num_faces)
    for face_number in range(sd.num_faces):
        # The nodes which belong to a certain face
        face_node_indices = face_nodes.indices[
            face_nodes.indptr[face_number] : face_nodes.indptr[face_number + 1]
        ]

        node_coordinates = []
        for node in face_node_indices:
            node_coordinates.append(sd.nodes[:, node])

        neighboring_cells = np.where(sd.cell_faces[face_number].todense() != 0)[
            1
        ].tolist()

        def compute_volume(sd, cell, nodes):
            cc = sd.cell_centers[:, cell].reshape(-1, 1)
            polygon = np.stack(nodes, axis=1)
            distance = pp.geometry.distances.points_polygon(
                p=cc, poly=polygon, tol=1e-8
            )[0]
            area = 1 / 3 * distance * sd.face_areas[face_number]
            return area

        def compute_area(sd, cell, nodes):
            cc = sd.cell_centers[:, cell]
            starting_point = nodes[0]
            end_point = nodes[1]
            distance, _ = pp.geometry.distances.points_segments(
                p=cc, start=starting_point, end=end_point
            )
            area = 1 / 2 * distance * sd.face_areas[face_number]
            return area

        area_local = 0
        for cell in neighboring_cells:
            if sd.dim == 3:
                area_local += compute_volume(sd, cell, node_coordinates)
            elif sd.dim == 2:
                area_local += compute_area(sd, cell, node_coordinates)
        all_meas[face_number] = area_local
    meas = all_meas

    if not is_scalar:
        meas = meas.repeat(sd.dim)

    # Obtain numerator and denominator to determine the error.
    numerator = np.sqrt(np.sum(meas * np.abs(true_array - approx_array) ** 2))
    denominator = np.sqrt(np.sum(meas * np.abs(true_array) ** 2)) if relative else 1.0

    return numerator / denominator
