# (i)PoPe example

This repository contains an example of how to use PoPe and iPoPe methods to verify a code / a solution by studying 
simulations outputs. As described in the paper `Verification and accuracy check of simulations with PoPe and iPoPe`. The
system chosen is chaotic (in chaotic regime) to stress the ability to verify 
simulations in non linear chaotic regimes, where other verification methods might fail.

# Installation

A few mainstream packages are needed. To install them: `python3 -m pip install -r requirements`

# Usage

In order to reproduce some plots of the paper XXX, run the script `src/main.py`. This script calls:
``` content of src/main.py
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
```

### Repository structure
```tree -A -I "__init__.py|venv|__pycache__|*.png|*.p.gz"
.
├── data                                                < where the data produced is stored
│   └── simulation                                      < where simulations are stored
├── figures                                             < where figures are stored
│   ├── lyapunov_exponent.png                           < Fig. 8 left
│   ├── verification_of_finite_differences.png          < Fig. 9 left
│   ├── PoPe                                            < where PoPe figures are stored, Fig. 10
│   └── plot_poincare_plots                             < where Poincare plots are stored, Fig. 2, Fig. 7
├── requirements.txt                                    < list of the required python packages
└── src                                                 < where the python sources are stored
    ├── compute_PoPe.py                                 < script computing PoPe analysis
    ├── compute_trajectory.py                           < script computing a set of simulations
    ├── main.py                                         < main script calling other scripts
    ├── plot_lyapunov_exponent.py                       < script plotting the lypunov exponent figure
    ├── plot_poincare_plots.py                          < script plotting Poincare plots
    ├── plot_PoPe.py                                    < script plotting PoPe results
    ├── README.md                                       < file you currently read
    ├── utils.py                                        < few utils functions used in different scripts
    └── verification_of_finite_differences.py           < script verifying the finite differences functions
```
