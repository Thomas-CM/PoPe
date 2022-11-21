import os
import glob
import numpy as np
import compress_pickle
import matplotlib.pyplot as plt


"""
This script create a list of Poincare plots out of a list of simulation files
"""

list_of_simulation_files = glob.glob("../data/simulation/*.p.gz")
if not list_of_simulation_files:
    print("There are not simulation to study in ../data/simulation")

for simulation_file in list_of_simulation_files:

    print("Producing Poincare plot of file "+simulation_file)

    # load data
    data = compress_pickle.load(simulation_file, compression="gzip")

    X = data['X']
    parameters = data['parameters']
    numerical_scheme = data['numerical_scheme']
    N_iter_per_time_unit = data['N_iter_per_time_unit']
    N_iter_of_time_unit_transient = int(data['N_iter_of_time_unit']/2)

    plt.plot(np.mod(X[0, :, N_iter_of_time_unit_transient:] + .5, 1) - .5,
             X[1, :, N_iter_of_time_unit_transient:], '.k', markersize=1)
    plt.plot(np.mod(X[0, :, 0] + .5, 1) - .5,
             X[1, :, 0], '*r', markersize=3)

    name = f"Poincaré plot, sigma={parameters['sigma_chir']}, nu={parameters['nu']}, " + \
           f"scheme={numerical_scheme}, N={N_iter_per_time_unit}"

    plt.xlim([-.5, .5])

    if parameters["sigma_chir"] == 7 and parameters["nu"] == .2:

        if numerical_scheme == "rk2":
            if N_iter_per_time_unit <= 16:
                plt.ylim([-200, 200])
            elif N_iter_per_time_unit < 64:
                plt.ylim([-40, 40])
            else:
                plt.ylim([-10, 10])
        elif numerical_scheme == "rk4":
            if N_iter_per_time_unit < 12:
                plt.ylim([-70, 70])
            else:
                plt.ylim([-10, 10])

    elif parameters["sigma_chir"] == 2.3 and parameters["nu"] == .8:
        plt.ylim([-5, 5])

    plt.xlabel('$\Theta$')
    plt.ylabel('J')
    plt.title(name)
    os.makedirs('../figures/plot_poincare_plots', exist_ok=True)
    plt.savefig("../figures/plot_poincare_plots/"+data["name_of_simulation"]+".png")
    plt.close()