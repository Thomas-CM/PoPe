import os
import glob
import numpy as np
import compress_pickle
import matplotlib.pyplot as plt
import pandas as pd


"""
This script compute lyapunov exponent and plot it
"""

list_of_simulation_files = glob.glob("../data/simulation/*.p.gz")
if not list_of_simulation_files:
    print("There are not simulation to study in ../data/simulation")

summary = pd.DataFrame(columns=['sigma_chir', 'nu', 'B', 'numerical_scheme', 'N_iter_per_time_unit', 'dt', 'd0',
                                'name_of_simulation', 'lyapunov_exponent'])

for simulation_file in list_of_simulation_files:

    print('Analysing Lyapunov exponent of file '+simulation_file)

    # load data
    data = compress_pickle.load(simulation_file, compression="gzip")

    Lyapunov_data = data['Lyapunov_data']
    parameters = data['parameters']
    N_iter_of_time_unit = data['N_iter_of_time_unit']
    d0 = data['d0']

    end_of_transient = int(N_iter_of_time_unit/2)

    dX1 = (Lyapunov_data[:,:,end_of_transient:,1] - Lyapunov_data[:,:,end_of_transient:,0])/d0
    dX2 = (Lyapunov_data[:,:,end_of_transient:,2] - Lyapunov_data[:,:,end_of_transient:,0])/d0

    d1 = np.sqrt(dX1[0,:]**2 + dX1[1,:]**2)
    d2 = np.sqrt(dX2[0,:]**2 + dX2[1,:]**2)

    d = np.maximum.reduce([d1,d2])

    lyapunov_exponent = np.mean(np.log10(d)/0.25)

    #d = np.maximum.reduce([d1[0,:], d2[0,:]])
    #lyapunov_exponent = np.mean(np.log10(d/d0)/0.25)

    summary = summary.append({
        'sigma_chir': parameters['sigma_chir'],
        'nu': parameters['nu'],
        'B': parameters['B'],
        'numerical_scheme': data['numerical_scheme'],
        'N_iter_per_time_unit': data['N_iter_per_time_unit'],
        'dt': parameters['dt'],
        'd0': data['d0'],
        'name_of_simulation': data['name_of_simulation'],
        'lyapunov_exponent': lyapunov_exponent
    }, ignore_index=True)


summary.sort_values(by=['sigma_chir','numerical_scheme','N_iter_per_time_unit'],inplace=True)

legends = []
plt.figure()
for B in summary.B.unique():
    for numerical_scheme in summary.numerical_scheme.unique():
        tmp = summary.loc[
            (summary.numerical_scheme == numerical_scheme) & (summary.B == B)
        ]
        plt.plot(tmp.N_iter_per_time_unit,tmp.lyapunov_exponent)
        legends.append(f"{numerical_scheme}, case {B}")
plt.legend(legends)
plt.xscale('log')
plt.xlabel('N = Number of time steps per unit of time (log scale)')
plt.ylabel('Lyapunov exponent')
plt.title("Lyapnuov exponent, Fig. 8 left")

os.makedirs('../figures', exist_ok=True)
plt.savefig("../figures/lyapunov_exponent.png")