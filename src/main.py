# computes a set of simulations and store the results in ./data/simulation
import compute_trajectory

# plots Poincare plots illustrating the system's chaotic behavior and store figures in ./figures/plot_poincare_plots
import plot_poincare_plots

# computes and plots the lyapunov exponents with respect to the numerical discretisation and stores it in ./figures
import plot_lyapunov_exponent

# verifies the accuracy of the finite differences functions written in utils.py and used in the PoPe analysis
import verification_of_finite_differences

# computes the PoPe and iPoPe analysis, stores results in ./data/PoPe, stores figures in ./figures/PoPe
import compute_PoPe

# uses results from compute_PoPe to plot trends in ./figures
import plot_PoPe
