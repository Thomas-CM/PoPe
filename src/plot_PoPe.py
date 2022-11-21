import matplotlib.pyplot as plt
import pandas as pd


"""
This script plots PoPe and iPoPe trends
"""

summary = pd.read_csv("../data/PoPe_summary_0D.csv")

summary.sort_values(by=['N_iter_per_time_unit','B','numerical_scheme','reconstruction_order'], inplace=True)
color = {2:'r',4:'g',6:'b',8:'c'}
case = 1
for B in summary.B.unique():
    for numerical_scheme in summary.numerical_scheme.unique():
        for reconstruction_order in summary.reconstruction_order.unique():

            tmp = summary.loc[
                (summary.B==B) &
                (summary.numerical_scheme==numerical_scheme) &
                (summary.reconstruction_order==reconstruction_order)
            ]

            plt.figure(1,figsize=(10,10))
            plt.subplot(2,2,case)
            plt.plot(tmp.N_iter_per_time_unit, tmp.norm_2_x_delta,'--o'+color[reconstruction_order])

            plt.figure(2,figsize=(10,10))
            plt.subplot(2,2,case)
            plt.plot(tmp.N_iter_per_time_unit, tmp.norm_2_x_RHS_effective,'--o'+color[reconstruction_order])

            plt.figure(3,figsize=(10,10))
            plt.subplot(2,2,case)
            plt.plot(tmp.N_iter_per_time_unit, tmp.norm_2_x_RHS_theoretical,'--o'+color[reconstruction_order])

            plt.figure(4,figsize=(10,10))
            plt.subplot(2,2,case)
            plt.plot(tmp.N_iter_per_time_unit, tmp.norm_2_x_PoPe_residual,'--o'+color[reconstruction_order])

            plt.figure(5,figsize=(10,10))
            plt.subplot(2,2,case)
            plt.plot(tmp.N_iter_per_time_unit, tmp.norm_2_x_iPoPe_residual,'--o'+color[reconstruction_order])

            plt.figure(6,figsize=(10,10))
            plt.subplot(2,2,case)
            plt.plot(tmp.N_iter_per_time_unit, tmp.coefficient_x_full,'--x'+color[reconstruction_order])
            plt.figure(7,figsize=(10,10))
            plt.subplot(2,2,case)
            plt.plot(tmp.N_iter_per_time_unit, tmp.coefficient_x_delta,'--o'+color[reconstruction_order])
            plt.figure(8,figsize=(10,10))
            plt.subplot(2,2,case)
            plt.plot(tmp.N_iter_per_time_unit, tmp.coefficient_x_iPoPe,'--+'+color[reconstruction_order])

            plt.figure(101,figsize=(10,10))
            plt.subplot(2,2,case)
            plt.plot(tmp.N_iter_per_time_unit, tmp.norm_2_y_delta,'-*'+color[reconstruction_order])

            plt.figure(102, figsize=(10, 10))
            plt.subplot(2, 2, case)
            plt.plot(tmp.N_iter_per_time_unit, tmp.norm_2_y_RHS_effective, '--o' + color[reconstruction_order])

            plt.figure(103, figsize=(10, 10))
            plt.subplot(2, 2, case)
            plt.plot(tmp.N_iter_per_time_unit, tmp.norm_2_y_RHS_theoretical, '--o' + color[reconstruction_order])

            plt.figure(104, figsize=(10, 10))
            plt.subplot(2, 2, case)
            plt.plot(tmp.N_iter_per_time_unit, tmp.norm_2_y_PoPe_residual, '--o' + color[reconstruction_order])

            plt.figure(105, figsize=(10, 10))
            plt.subplot(2, 2, case)
            plt.plot(tmp.N_iter_per_time_unit, tmp.norm_2_y_iPoPe_residual, '--o' + color[reconstruction_order])

            plt.figure(106,figsize=(10,10))
            plt.subplot(2,2,case)
            plt.plot(tmp.N_iter_per_time_unit, tmp.coefficient1_y_full,'--x'+color[reconstruction_order])
            plt.plot(tmp.N_iter_per_time_unit, tmp.coefficient2_y_full,'-x'+color[reconstruction_order])
            plt.figure(107, figsize=(10, 10))
            plt.subplot(2, 2, case)
            plt.plot(tmp.N_iter_per_time_unit, tmp.coefficient1_y_delta,'--o'+color[reconstruction_order])
            plt.plot(tmp.N_iter_per_time_unit, tmp.coefficient2_y_delta,'-o'+color[reconstruction_order])
            plt.figure(108, figsize=(10, 10))
            plt.subplot(2, 2, case)
            plt.plot(tmp.N_iter_per_time_unit, tmp.coefficient1_y_iPoPe,'--+'+color[reconstruction_order])
            plt.plot(tmp.N_iter_per_time_unit, tmp.coefficient2_y_iPoPe,'-+'+color[reconstruction_order])

        for ifig in [1,2,3,4,5,6,7,8]:
            for offset in [0,100]:
                plt.figure(ifig+offset)
                plt.xscale('log')
                if ifig in [6,7,8]:
                    plt.yscale('linear')
                else:
                    plt.yscale('log')
                plt.title(f"sigma ={tmp.iloc[0].sigma_chir}, "
                          f"nu = {tmp.iloc[0].nu}, "
                          f"{numerical_scheme}")

        case = case + 1


def format_figure(id_fig, xlabel, ylabel, file_name):
    """
    function formating the figure and saving it
    :param id_fig: id of the figure concerned
    :param xlabel: xlabel to print
    :param ylabel:  ylabel to print
    :param file_name: file name used to save the figure
    """

    plt.figure(id_fig)
    for sub in [3, 4]:
        plt.subplot(2, 2, sub)
        plt.xlabel(xlabel)
    for sub in [1, 3]:
        plt.subplot(2, 2, sub)
        plt.ylabel(ylabel)

    plt.savefig(file_name)
    plt.close()

format_figure(id_fig=1,
              xlabel='N',
              ylabel='eq x : norm 2 of effective - theoretical RHS',
              file_name='../figures/summary_eq_x_norm_2_x_delta.png')

format_figure(id_fig=2,
              xlabel='N',
              ylabel='eq x : norm 2 of effective RHS',
              file_name='../figures/summary_eq_x_norm_2_x_RHS_effective.png')

format_figure(id_fig=3,
              xlabel='N',
              ylabel='eq x : norm 2 of theoretical RHS',
              file_name='../figures/summary_eq_x_norm_2_x_RHS_theoretical.png')

format_figure(id_fig=4,
              xlabel='N',
              ylabel='eq x : norm 2 of PoPe residual',
              file_name='../figures/summary_eq_x_norm_2_x_PoPe_residual.png')

format_figure(id_fig=5,
              xlabel='N',
              ylabel='eq x : norm 2 of iPoPe residual',
              file_name='../figures/summary_eq_x_norm_2_x_iPoPe_residual.png')

format_figure(id_fig=6,
              xlabel='N',
              ylabel='eq x : coefficients',
              file_name='../figures/summary_eq_x_coefficients_full.png')

format_figure(id_fig=7,
              xlabel='N',
              ylabel='eq x : coefficients',
              file_name='../figures/summary_eq_x_coefficients_delta.png')

format_figure(id_fig=8,
              xlabel='N',
              ylabel='eq x : coefficients',
              file_name='../figures/summary_eq_x_coefficients_iPoPe.png')

###

format_figure(id_fig=101,
              xlabel='N',
              ylabel='eq y : norm 2 of effective - theoretical RHS',
              file_name='../figures/summary_eq_y_norm_2_y_delta.png')

format_figure(id_fig=102,
              xlabel='N',
              ylabel='eq y : norm 2 of effective RHS',
              file_name='../figures/summary_eq_y_norm_2_y_RHS_effective.png')

format_figure(id_fig=103,
              xlabel='N',
              ylabel='eq y : norm 2 of theoretical RHS',
              file_name='../figures/summary_eq_y_norm_2_y_RHS_theoretical.png')

format_figure(id_fig=104,
              xlabel='N',
              ylabel='eq y : norm 2 of PoPe residual',
              file_name='../figures/summary_eq_y_norm_2_y_PoPe_residual.png')

format_figure(id_fig=105,
              xlabel='N',
              ylabel='eq y : norm 2 of iPoPe residual',
              file_name='../figures/summary_eq_y_norm_2_y_iPoPe_residual.png')

format_figure(id_fig=106,
              xlabel='N',
              ylabel='eq y : coefficients',
              file_name='../figures/summary_eq_y_coefficients_full.png')

format_figure(id_fig=107,
              xlabel='N',
              ylabel='eq y : coefficients',
              file_name='../figures/summary_eq_y_coefficients_delta.png')

format_figure(id_fig=108,
              xlabel='N',
              ylabel='eq y : coefficients',
              file_name='../figures/summary_eq_y_coefficients_iPoPe.png')