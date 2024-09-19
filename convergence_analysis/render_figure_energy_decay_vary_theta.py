import os
import sys

import numpy as np

sys.path.append("../")
import matplotlib.pyplot as plt
import plotting.plot_utils as pu


import runscript_ABC_energy_vary_theta

folder_name = "energy_values"
script_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(script_dir, folder_name)
os.makedirs(output_dir, exist_ok=True)

# Tuple value in dictionary:
#   * Legend text
#   * Color
#   * Dashed/not dashed line
#   * Logarithmic y scale/not logarithmic y scale.
index_angle_dict = {
    0: ("$\\theta = 0$", pu.RGB(216, 27, 96), False, False),
    1: ("$\\theta = \pi/6$", pu.RGB(30, 136, 229), False, False),
    2: ("$\\theta = \pi/3$", pu.RGB(255, 193, 7), True, False),
    3: ("$\\theta = \pi/4$", pu.RGB(0, 0, 0), True, False),
    4: ("$\\theta = \pi/8$", pu.RGB(216, 27, 96), True, False),
}

for key, value in index_angle_dict.items():
    filename = os.path.join(output_dir, f"energy_values_{key}.txt")
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
