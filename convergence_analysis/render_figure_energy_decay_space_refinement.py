import sys

import numpy as np

sys.path.append("../")
import matplotlib.pyplot as plt
import plotting.plot_utils as pu

import runscript_ABC_energy_vary_dx

relative_path = "convergence_analysis/energy_values/"

# Tuple value in dictionary:
#   * Legend text
#   * Color
#   * Dashed/not dashed line
#   * Logarithmic y scale/not logarithmic y scale.
index_dx_dict = {
    12: ("$\Delta x = 1/512$", pu.RGB(0, 0, 0), False, True),
    11: ("$\Delta x = 1/256$", pu.RGB(216, 27, 96), False, True),
    10: ("$\Delta x = 1/128$", pu.RGB(30, 136, 229), True, True),
    9: ("$\Delta x = 1/64$", pu.RGB(255, 193, 7), True, True),
    8: ("$\Delta x = 1/32$", pu.RGB(0, 0, 0), True, True),
}

for key, value in index_dx_dict.items():
    filename = f"{relative_path}energy_values_{key}.txt"
    energy_values = (
        pu.read_float_values(filename=filename)
        / pu.read_float_values(filename=filename)[0]
    )
    final_time = 15
    time_values = np.linspace(0, final_time, len(energy_values))

    plt.yscale("log" if value[3] else "linear")
    plt.plot(
        time_values,
        energy_values,
        label=value[0],
        color=value[1],
        linestyle="-" if not value[2] else "--",
    )
    print(value[0], energy_values[-1] * 100)

plt.axvline(
    x=10 / np.sqrt(3),
    ymin=0,
    ymax=5,
    color=(0.65, 0.65, 0.65),
    linestyle="--",
    linewidth=1,
)
plt.axvline(
    x=10 * np.sqrt(6) / 3,
    ymin=0,
    ymax=5,
    color=(0.65, 0.65, 0.65),
    linestyle="--",
    linewidth=1,
)

plt.axhline(
    y=0,
    xmin=0,
    xmax=12,
    color=(0, 0, 0),
    linewidth=0.5,
)

plt.xlabel("Time [s]", fontsize=12)
plt.ylabel("$\\frac{E}{E_0}$", fontsize=16)
plt.title("Energy evolution with respect to time")
plt.legend()
plt.show()
