import os
import compress_pickle
import utils


"""
This script runs an ensemble of simulations that are defined by lists of parameters.
The output are a list of files containing the results of each simulation.
Those results can be post process to plot Poincaré plots or to verify simulations
"""


N_iter_of_time_unit = 2_000 # default = 2000
N_points = 30 # default = 30

# parameters of the system simulated
list_of_parameters = [
    {"sigma_chir": 7,
     "nu": .2},
    {"sigma_chir": 2.3,
     "nu": .8},
]

list_of_numerical_schemes = ["rk2", "rk4"]

list_of_N_iter_per_time_unit = [6, 8, 10, 12, 16, 20, 24, 28, 32, 40, 48, 56, 64, 96, 128, 192, 256, 512, 1024]

# initial perturbation used to compute the Lyapunov exponent,
# if = 0, this is not computed
d0 = 10**-8

save_PoPe_data = True

# loops over different lists

for parameters in list_of_parameters:
    parameters['B'] = 0.25 * (parameters['sigma_chir']**2)

    for numerical_scheme in list_of_numerical_schemes:
        if numerical_scheme == "rk2":
            integrator = utils.rk2
        elif numerical_scheme == "rk4":
            integrator = utils.rk4

        for N_iter_per_time_unit in list_of_N_iter_per_time_unit:
            parameters['dt'] = 1. / N_iter_per_time_unit

            name_of_simulation = f"sigma_{parameters['sigma_chir']}__nu_{parameters['nu']}__" + \
                                 f"scheme_{numerical_scheme}__N_{N_iter_per_time_unit}__d0__{d0}"
            simulation_filename = '../data/simulation/'+name_of_simulation+'.p.gz'

            need_to_compute = True
            if os.path.exists(simulation_filename):
                try:
                    data = compress_pickle.load(simulation_filename, compression="gzip")
                except:
                    print(
                        "Simulation file " + simulation_filename +
                        " already exists but the file is corrupted. The simulation will be redone.")
                else:
                    print("Simulation already exists, remove file "+simulation_filename+" if you want to recompute it.")
                    need_to_compute = False

            if need_to_compute:
                # compute trajectories
                t_stored, X_stored_trajectory, Lyapunov_data, PoPe_data = utils.compute_trajectory(
                    N_iter_of_time_unit,
                    N_iter_per_time_unit,
                    N_points,
                    integrator,
                    utils.RHS_strange_attractor,
                    parameters,
                    d0=d0,
                    save_PoPe_data=save_PoPe_data
                )

                # save data
                os.makedirs('../data/simulation', exist_ok=True)
                compress_pickle.dump(
                    {"t": t_stored,
                     "X": X_stored_trajectory,
                     "PoPe_data": PoPe_data,
                     "Lyapunov_data": Lyapunov_data,
                     "parameters": parameters,
                     "numerical_scheme": numerical_scheme,
                     "N_iter_per_time_unit": N_iter_per_time_unit,
                     "N_iter_of_time_unit": N_iter_of_time_unit,
                     "N_points": N_points,
                     "d0": d0,
                     "name_of_simulation": name_of_simulation
                     },
                    simulation_filename,
                    compression="gzip"
                )