import os
import glob
import numpy as np
import compress_pickle
import matplotlib.pyplot as plt
import pandas as pd
import utils


"""
This script computes PoPe and iPoPe post treatment
It can plots some "high dimensional" results for each simulation, otherwise, the script plot_PoPe plots trends over
different simulations
"""

# variable used to store 0 dimensional data in one file for the entire set of simulation
summary_0D = pd.DataFrame()

plot_details_of_each_simulation = True

list_of_simulation_files = glob.glob("../data/simulation/*.p.gz")
if not list_of_simulation_files:
    print("There are not simulation to study in ../data/simulation")

else:
    for simulation_file in list_of_simulation_files:
        print("Studying file " + simulation_file)

        # load data
        data = compress_pickle.load(simulation_file, compression="gzip")

        t = data['t'].copy()
        X = data['X'].copy()
        PoPe_data = data['PoPe_data'].copy()
        parameters = data['parameters'].copy()
        N_points = data['N_points']
        name_of_simulation = data['name_of_simulation']
        N_iter_per_time_unit = data['N_iter_per_time_unit']
        dt = parameters['dt']
        nu = parameters['nu']
        B = parameters['B']
        sigma_chir = parameters['sigma_chir']
        numerical_scheme = data['numerical_scheme']

        del data

        # selection of the subset of iterations used to do the verification
        # the first and last iterations in time are not be used because centered time derivatives cannot be used for
        # those two points, also the first used iteration is chosen far enough into the simulation to be on the strange
        # attractor
        first_iter = int(np.shape(X)[2]/2)
        last_iter = -1

        # (re)compute RHS, see paper for definition

        # dX[0,:]/dt = X[1,:]
        x_operator = X[1,:,first_iter:last_iter].flatten()
        x_RHS_theoretical = x_operator
        # first and last time iteration are not verified because centered time derivatives cannot be computed,
        # thus first and last point of the array are omitted
        # in this equation, there is only one operator, which make the verification quite singular

        # dX[1,:]/dt = ... see paper
        y_operator_1_2 = - 2*np.pi*parameters['B']*(
                np.sin( 2*np.pi*X[0,:,first_iter:last_iter] )
                +np.sin( 2*np.pi*(X[0,:,first_iter:last_iter] - np.tile(t[first_iter:last_iter],(N_points,1))) )
        )
        y_operator_1_2 = y_operator_1_2.flatten()
        y_operator_3 = - parameters['nu']*X[1,:,first_iter:last_iter].flatten()
        y_RHS_theoretical = y_operator_1_2 + y_operator_3

        xs_for_dt = np.moveaxis(np.squeeze(PoPe_data[0,:,first_iter:last_iter,:]),[0,1,2],[1,2,0])
        ys_for_dt = np.moveaxis(np.squeeze(PoPe_data[1,:,first_iter:last_iter,:]),[0,1,2],[1,2,0])

        # check alignment of X and PoPe_data, error should be 0
        #error_of_alignment = np.max(np.abs(X[:,:,first_iter:last_iter]-PoPe_data[:,:,first_iter:last_iter,4]))

        plt.figure()
        legends = []
        for i,order in enumerate([2,4,6,8]):

            # compute time derivatives
            x_RHS_effective = utils.first_derivative(xs_for_dt,parameters['dt'],order)[4,:,:].flatten()
            y_RHS_effective = utils.first_derivative(ys_for_dt,parameters['dt'],order)[4,:,:].flatten()

            # 'A' matrix for the x equation
            A_x = np.zeros((np.size(x_RHS_effective),1))
            A_x[:, 0] = x_operator

            # projection of PoPe eq 1 on the full RHS
            b_x = x_RHS_effective
            coefficient_x_full = np.linalg.lstsq(A_x, b_x)[0]

            # projection of PoPe eq 1 on the error only (delta = RHS_effective - RHS_theoretical)
            b_x_delta = x_RHS_effective - x_RHS_theoretical
            coefficient_x_delta = 1+np.linalg.lstsq(A_x, b_x_delta)[0]

            norm_2_x_delta = np.sqrt(np.mean( (x_RHS_effective - x_RHS_theoretical) ** 2 ))
            norm_2_x_RHS_effective = np.sqrt(np.mean( x_RHS_effective ** 2 ))
            norm_2_x_RHS_theoretical = np.sqrt(np.mean( x_RHS_theoretical ** 2 ))
            x_PoPe_residual = x_RHS_theoretical-np.matmul(A_x,coefficient_x_full)
            norm_2_x_PoPe_residual = np.sqrt(np.mean( x_PoPe_residual ** 2 ))

            # projection of iPoPe eq 1 :: PoPe and iPoPe are equivalent as there is only one operator in this equation
            coefficient_x_iPoPe = 1+np.dot(b_x_delta,A_x[:,0]) / np.dot(A_x[:,0],A_x[:,0])
            x_iPoPe_residual = x_RHS_theoretical.T-((coefficient_x_iPoPe)*A_x).flatten()
            norm_2_x_iPoPe_residual = np.sqrt(np.mean( x_iPoPe_residual ** 2 ))

            # 'A' matrix for the y equation
            A_y = np.zeros((np.size(y_RHS_effective),2))
            A_y[:,0] = y_operator_1_2
            A_y[:,1] = y_operator_3

            # projection of PoPe eq 2 on the full RHS
            b_y = y_RHS_effective
            coefficient_y_full = np.linalg.lstsq(A_y, b_y)[0]
            coefficient1_y_full = coefficient_y_full[0]
            coefficient2_y_full = coefficient_y_full[1]

            # projection of PoPe eq 2 on the error only (delta = RHS_effective - RHS_theoretical)
            b_y_delta = y_RHS_effective - y_RHS_theoretical
            coefficient_y_delta = 1+np.linalg.lstsq(A_y, b_y_delta)[0]
            coefficient1_y_delta = coefficient_y_delta[0]
            coefficient2_y_delta = coefficient_y_delta[1]

            norm_2_y_delta = np.sqrt(np.mean( (y_RHS_effective - y_RHS_theoretical) ** 2 ))
            norm_2_y_RHS_effective = np.sqrt(np.mean( y_RHS_effective ** 2 ))
            norm_2_y_RHS_theoretical = np.sqrt(np.mean( y_RHS_theoretical ** 2 ))
            y_PoPe_residual = y_RHS_theoretical-np.matmul(A_y,coefficient_y_full)
            norm_2_y_PoPe_residual = np.sqrt(np.mean( y_PoPe_residual ** 2 ))

            # projection of iPoPe eq 2
            coefficient_y_iPoPe = np.zeros(coefficient_y_delta.shape)
            for operator in range(2):
                coefficient_y_iPoPe[operator] = 1 + np.dot(b_y_delta,A_y[:,operator]) \
                                                / np.dot(A_y[:,operator],A_y[:,operator])
            coefficient1_y_iPoPe = coefficient_y_iPoPe[0]
            coefficient2_y_iPoPe = coefficient_y_iPoPe[1]
            y_iPoPe_residual = y_RHS_theoretical-np.matmul(A_y,coefficient_y_iPoPe)
            norm_2_y_iPoPe_residual = np.sqrt(np.mean( y_iPoPe_residual ** 2 ))

            local_summary_0D = {
                    'N_iter_per_time_unit': N_iter_per_time_unit,
                    'dt': dt,
                    'nu': nu,
                    'B': B,
                    'numerical_scheme': numerical_scheme,
                    'sigma_chir': sigma_chir,
                    'name_of_simulation': name_of_simulation,
                    'reconstruction_order': order,
                    'coefficient_x_full':coefficient_x_full[0],
                    'coefficient_x_delta':coefficient_x_delta[0],
                    'norm_2_x_delta':norm_2_x_delta,
                    'norm_2_x_RHS_effective':norm_2_x_RHS_effective,
                    'norm_2_x_RHS_theoretical':norm_2_x_RHS_theoretical,
                    'norm_2_x_PoPe_residual':norm_2_x_PoPe_residual,
                    'coefficient_x_iPoPe':coefficient_x_iPoPe,
                    'norm_2_x_iPoPe_residual':norm_2_x_iPoPe_residual,
                    'coefficient1_y_full':coefficient1_y_full,
                    'coefficient2_y_full':coefficient2_y_full,
                    'coefficient1_y_delta':coefficient1_y_delta,
                    'coefficient2_y_delta':coefficient2_y_delta,
                    'norm_2_y_delta':norm_2_y_delta,
                    'norm_2_y_RHS_effective':norm_2_y_RHS_effective,
                    'norm_2_y_RHS_theoretical':norm_2_y_RHS_theoretical,
                    'norm_2_y_PoPe_residual':norm_2_y_PoPe_residual,
                    'coefficient1_y_iPoPe':coefficient1_y_iPoPe,
                    'coefficient2_y_iPoPe':coefficient2_y_iPoPe,
                    'norm_2_y_iPoPe_residual':norm_2_y_iPoPe_residual
                }

            summary_0D = summary_0D.append(local_summary_0D, ignore_index=True)
            summary_0D.to_csv('../data/PoPe_summary_0D.csv')

            if plot_details_of_each_simulation:

                os.makedirs('../figures/PoPe', exist_ok=True)
                PoPe_root_filename = simulation_file.replace('data','figures').replace('simulation','PoPe').replace('.p.gz','')

                plt.figure(figsize=(20,20))

                plt.subplot(3,3,1)
                plt.plot(x_RHS_theoretical, b_x_delta, '.')
                plt.title('x equation')
                plt.xlabel('theoretical RHS')
                plt.ylabel('effective - theoretical RHS')

                plt.subplot(3,3,2)
                plt.plot(x_RHS_theoretical, x_PoPe_residual, '.')
                plt.title('x equation')
                plt.xlabel('theoretical RHS')
                plt.ylabel('PoPe residual')

                plt.subplot(3,3,3)
                plt.plot(x_RHS_theoretical, x_iPoPe_residual, '.')
                plt.title('x equation')
                plt.xlabel('theoretical RHS')
                plt.ylabel('iPoPe residual')

                plt.subplot(3, 3, 5)
                plt.plot(b_x_delta, x_PoPe_residual, '.')
                plt.title('x equation')
                plt.xlabel('effective - theoretical RHS')
                plt.ylabel('PoPe residual')

                plt.subplot(3, 3, 6)
                plt.plot(b_x_delta, x_iPoPe_residual, '.')
                plt.title('x equation')
                plt.xlabel('effective - theoretical RHS')
                plt.ylabel('iPoPe residual')

                plt.subplot(3, 3, 9)
                plt.plot(x_PoPe_residual, x_iPoPe_residual, '.')
                plt.title('x equation')
                plt.xlabel('PoPe residual')
                plt.ylabel('iPoPe residual')

                plt.savefig(PoPe_root_filename+'__eq_x__order__'+str(order)+'.png')
                plt.close()

                plt.figure(figsize=(20,20))

                plt.subplot(3, 3, 1)
                plt.plot(y_RHS_theoretical, b_y_delta, '.')
                plt.title('y equation')
                plt.xlabel('theoretical RHS')
                plt.ylabel('effective - theoretical RHS')

                plt.subplot(3, 3, 2)
                plt.plot(y_RHS_theoretical, y_PoPe_residual, '.')
                plt.title('y equation')
                plt.xlabel('theoretical RHS')
                plt.ylabel('PoPe residual')

                plt.subplot(3, 3, 3)
                plt.plot(y_RHS_theoretical, y_iPoPe_residual, '.')
                plt.title('y equation')
                plt.xlabel('theoretical RHS')
                plt.ylabel('iPoPe residual')

                plt.subplot(3, 3, 5)
                plt.plot(b_y_delta, y_PoPe_residual, '.')
                plt.title('y equation')
                plt.xlabel('effective - theoretical RHS')
                plt.ylabel('PoPe residual')

                plt.subplot(3, 3, 6)
                plt.plot(b_y_delta, y_iPoPe_residual, '.')
                plt.title('y equation')
                plt.xlabel('effective - theoretical RHS')
                plt.ylabel('iPoPe residual')

                plt.subplot(3, 3, 9)
                plt.plot(y_PoPe_residual, y_iPoPe_residual, '.')
                plt.title('y equation')
                plt.xlabel('PoPe residual')
                plt.ylabel('iPoPe residual')

                plt.savefig(PoPe_root_filename + '__eq_y__order__' + str(order) + '.png')
                plt.close()
