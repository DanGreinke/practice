import numpy as np
# Note: credit to 3blue1brown -> https://github.com/3b1b/videos/blob/1385f6d12cabf0c7cd61bdd9b78987e77b9edfbf/_2019/diffyq/solve_pendulum_ode_sample_code.py
# Physical constants
g = 9.8
L = 2
mu = 0.1

THETA_0 = np.pi / 3  # 60 degrees
THETA_DOT_0 = 0  # No initial angular velocity

# Definition of ODE
def get_theta_double_dot(theta, theta_dot):
    return -mu * theta_dot - (g / L) * np.sin(theta)


# Solution to the differential equation
def theta(t):
    # Initialize changing values
    theta = THETA_0
    theta_dot = THETA_DOT_0
    delta_t = 0.01  # Some time step
    for time in np.arange(0, t, delta_t):
        # Take many little time steps of size delta_t
        # until the total time is the function's input
        theta_double_dot = get_theta_double_dot(
            theta, theta_dot
        )
        theta += theta_dot * delta_t
        theta_dot += theta_double_dot * delta_t
    return theta


print(theta(2)) # Output: -0.5623420732009209






























