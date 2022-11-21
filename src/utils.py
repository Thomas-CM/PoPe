import numpy as np


"""
List of utils functions called by different main functions
"""


def compute_trajectory(N_iter_of_time_unit: int,
                       N_iter_per_time_unit: int,
                       N_points: int,
                       integrator,
                       func,
                       parameters: dict,
                       d0: float = 0.0,
                       save_PoPe_data: bool = False):
    """
    Function computing the trajectory of a system defined by the function 'func' and the numerical method 'integrator'
    :param N_iter_of_time_unit: ~= duration of the simulation
    :param N_iter_per_time_unit: ~= discretisation of the simulation
    :param N_points: number of points simulated
    :param integrator: method chosen for the numerical integration
    :param func: function to integrate ~= equation of the system
    :param parameters: dictionary containing the parameter of the system
    :param d0: initial perturbation used to compute Lyapunov exponent
    :param save_PoPe_data: switch used to save or not PoPe data
    :return: t_stored, X_stored_trajectory, Lyapunov_data, PoPe_data
    """

    # initial conditions
    X = np.zeros((2, N_points))
    X[1, :] = np.linspace(-3, 3, N_points)
    t = 0

    # variable used to compute the Lyapunov exponent
    Nsmall = int(N_iter_per_time_unit / 4)  # number of integration time steps used to compute exponentiation
    if d0:
        Lyapunov_data = np.zeros((2,N_points,N_iter_of_time_unit, 3))
    else:
        Lyapunov_data = None

    # variable used to store data required for PoPe
    # dimensions:
    #   2 for the 2 dimensions of the problem,
    #   N_points for the N_points simulated,
    #   N_iter_of_time_unit for every point of the trajectory already saved in X_stored_trajectory
    #   9 for the 4 preceding and 4 following points in time with respect to the point saved in X_stored_trajectory
    # this approach is not optimal, an iterative on-line approach is used for regular production code, see example in
    # time_derivative_iterative_computation.py
    if save_PoPe_data:
        PoPe_data = np.zeros((2, N_points, N_iter_of_time_unit, 9))
    else:
        PoPe_data = None

    # main variable
    X_stored_trajectory = np.zeros((2, N_points, N_iter_of_time_unit))
    X_stored_trajectory[:, :, 0] = X
    t_stored = np.zeros((N_iter_of_time_unit,))

    # time integration "over time unit"
    for i_time_unit in range(N_iter_of_time_unit):

        X_stored_trajectory[:, :, i_time_unit] = X
        t_stored[i_time_unit] = t

        if np.mod(i_time_unit, N_iter_of_time_unit / 20) == 0:
            print(f"{parameters['sigma_chir']}, {parameters['nu']}, {integrator.__name__}, {N_iter_per_time_unit} : " + \
                  f"i_time_unit = {i_time_unit} out of {N_iter_of_time_unit}")

        # compute lyapunov using a few timestep
        if d0:

            # initialize 3 sets of points with a gap of d0 in each dimension
            X0 = X.copy()
            X1 = X.copy()
            X1[0, :] = X[0, :] + d0
            X2 = X.copy()
            X2[1, :] = X[1, :] + d0

            # integrating trajectories on Nsmall time steps
            for step in range(Nsmall):
                X0 = integrator(X0, t + parameters['dt'] * step, parameters, func)
                X1 = integrator(X1, t + parameters['dt'] * step, parameters, func)
                X2 = integrator(X2, t + parameters['dt'] * step, parameters, func)

            Lyapunov_data[:,:,i_time_unit,0] =X0
            Lyapunov_data[:,:,i_time_unit,1] =X1
            Lyapunov_data[:,:,i_time_unit,2] =X2

        # time integration "within time unit"
        for i_iter_in_per_time_unit in range(N_iter_per_time_unit):

            if save_PoPe_data and i_time_unit != 0 and i_time_unit != (N_iter_of_time_unit-1):
                # first (0) and last (N_iter_of_time_unit-1) are not stored as there are necessarily incomplete when
                # using centered finite differences in post treatment
                PoPe_data = build_PoPe_data(X,
                                            PoPe_data,
                                            i_iter_in_per_time_unit,
                                            N_iter_per_time_unit,
                                            i_time_unit)

            # compute trajectory
            X = integrator(X, t, parameters, func)

            t = t + parameters['dt']

    return t_stored, X_stored_trajectory, Lyapunov_data, PoPe_data

def RHS_strange_attractor(x, t: float, parameters: dict):
    """
    Computing the RHS of the strange attractor.
    This function is used multiple times to build rk2 and rk4 time integration
    eq 18a and 18b of paper
    :param x: coordinates
    :param t: time
    :param parameters of the equation
    :return: y = dx/dt
    """
    y = 0. * x
    y[0, :] = x[1, :]
    y[1, :] = - 2 * np.pi * parameters['B'] * (np.sin(2 * np.pi * x[0, :]) + np.sin(2 * np.pi * (x[0, :] - t))) \
              - parameters['nu'] * x[1, :]
    return y

def rk4(x, t: float, parameters: dict, func):
    """
    Computing rk4 time integration
    :param x: coordinates
    :param t: time
    :param parameters of the equation
    :param func: function to integrate
    :return: y: time integration of x coordinates
    """

    rhs1 = func(x, t, parameters)
    rhs2 = func(x + parameters['dt'] / 2 * rhs1, t + parameters['dt'] / 2, parameters)
    rhs3 = func(x + parameters['dt'] / 2 * rhs2, t + parameters['dt'] / 2, parameters)
    rhs4 = func(x + parameters['dt'] * rhs3, t + parameters['dt'], parameters)
    y = x + parameters['dt'] * (rhs1 + 2 * rhs2 + 2 * rhs3 + rhs4) / 6

    return y

def rk2(x, t: float, parameters: dict, func):
    """
    Computing rk2 time integration
    :param x: coordinates
    :param t: time
    :param parameters of the equation
    :param func: function to integrate
    :return: y: time integration of x coordinates
    """

    rhs1 = func(x, t, parameters)
    rhs2 = func(x + parameters['dt'] / 2 * rhs1, t + parameters['dt'] / 2, parameters)
    y = x + parameters['dt'] * rhs2

    return y

def build_PoPe_data(X,
                    PoPe_data,
                    i_iter_in_per_time_unit,
                    N_iter_per_time_unit,
                    i_time_unit):
    """
    Function selecting the necessary data for the "offline" PoPe study.
    This approach is not optimal, an iterative on-line approach is used for regular production code, see example in
    time_derivative_iterative_computation.py
    :param X: coordinate of the trajectory
    :param PoPe_data: array containing selected copies of X
    :param i_iter_in_per_time_unit: current iteration within a unit of time
    :param N_iter_per_time_unit: number of iteration within a time unit
    :param i_time_unit: current iteration of time
    :return: PoPe_data
    """

    # a list of "if" is used because several "if"s could be valid at once

    if i_iter_in_per_time_unit == N_iter_per_time_unit - 4:
        # 4 time steps before saving X the next i_time_unit
        PoPe_data[:,:,i_time_unit+1,0] = X

    if i_iter_in_per_time_unit == N_iter_per_time_unit - 3:
        # 3 time steps before saving X the next i_time_unit
        PoPe_data[:,:,i_time_unit+1,1] = X

    if i_iter_in_per_time_unit == N_iter_per_time_unit - 2:
        # 2 time steps before saving X the next i_time_unit
        PoPe_data[:,:,i_time_unit+1,2] = X

    if i_iter_in_per_time_unit == N_iter_per_time_unit - 1:
        # 1 time step before saving X the next i_time_unit
        PoPe_data[:,:,i_time_unit+1,3] = X

    if i_iter_in_per_time_unit == 0:
        # the X at i_time_unit has just been saved
        # this time step is already saved in stored_X_trajectory and it is not require to compute centered first
        # derivatives but it is still saved for pedagogic reasons and the "offline" verification of the synchronization
        # between X and PoPe_data
        PoPe_data[:,:,i_time_unit,4] = X

    if i_iter_in_per_time_unit == 1:
        # 1 time step after saving X
        PoPe_data[:,:,i_time_unit,5] = X

    if i_iter_in_per_time_unit == 2:
        # 2 time steps after saving X
        PoPe_data[:,:,i_time_unit,6] = X

    if i_iter_in_per_time_unit == 3:
        # 3 time steps after saving X
        PoPe_data[:,:,i_time_unit,7] = X

    if i_iter_in_per_time_unit == 4:
        # 4 time steps after saving X
        PoPe_data[:,:,i_time_unit,8] = X

    return PoPe_data

def first_derivative(x,dx,order):
    if order == 2:
        return first_derivative_order_2(x,dx)
    elif order == 4:
        return first_derivative_order_4(x,dx)
    elif order == 6:
        return first_derivative_order_6(x,dx)
    elif order == 8:
        return first_derivative_order_8(x,dx)
    else:
        raise Exception("first_derivative function called with an order not implemented, it should be 2, 4, 6 or 8")

def first_derivative_order_2(x,dx):
    y = np.zeros(x.shape)
    y[1:-1,] = (1/2*x[2:,]
               -1/2*x[:-2,])/dx
    return y

def first_derivative_order_4(x,dx):
    y = np.zeros(x.shape)
    y[2:-2,] = (-1/12*x[4:,] +2/3*x[3:-1,]
               -2/3*x[1:-3,] +1/12*x[:-4,])/dx
    return y

def first_derivative_order_6(x,dx):
    y = np.zeros(x.shape)
    y[3:-3,] = (1/60*x[6:,] -3/20*x[5:-1,] +3/4*x[4:-2,]
               -3/4*x[2:-4,] +3/20*x[1:-5,] -1/60*x[:-6,])/dx
    return y

def first_derivative_order_8(x,dx):
    y = np.zeros(x.shape)
    y[4:-4,] = (-1/280*x[8:,] +4/105*x[7:-1,] -1/5*x[6:-2,] +4/5*x[5:-3,]
               -4/5*x[3:-5,] +1/5*x[2:-6,] -4/105*x[1:-7,] +1/280*x[:-8,])/dx
    return y
