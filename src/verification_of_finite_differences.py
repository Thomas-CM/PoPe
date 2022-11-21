import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import utils

"""
This script verifies the implementation of the finite differences functions by comparing the analytical derivative of 
sin'(x) = cos(x) to numerical estimations of this derivative
"""

list_of_N = [2**x for x in range(4,20)]

summary = pd.DataFrame(columns={'N', 'order', 'error'})

for N in list_of_N:
    x = np.linspace(0,2*np.pi,N)
    y = np.sin(x)
    y_prime_theoretical = np.cos(x)
    for order in [2,4,6,8]:
        y_prime_finite_differences = utils.first_derivative(y, x[1]-x[0], order)
        summary = summary.append(
            {'N': N,
             'order': order,
             'error': np.mean(np.abs(y_prime_finite_differences[4:-4] - y_prime_theoretical[4:-4]))
             },
            ignore_index=True
        )

legends = []
plt.figure()
for order in summary.order.unique():
    tmp = summary.loc[summary.order == order]
    plt.plot(tmp.N, tmp.error)
    legends.append("order "+str(order))
plt.xscale('log')
plt.xlabel('N = Number of time steps (log scale)')
plt.yscale('log')
plt.ylabel('Norm 1 of error')
plt.legend(legends)
plt.title('Verification of the functions computing finite differences, Fig. 9')
os.makedirs('../figures', exist_ok=True)
plt.savefig("../figures/verification_of_finite_differences.png")