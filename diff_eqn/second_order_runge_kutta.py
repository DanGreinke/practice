import numpy as np
import pandas as pd
import time
"""
From Chapter 3.5 of "A First Course in Differential Equations" (1975) by Frank G. Hagin

This chapter applies the Runge-Kutta method to solve second order differential equations

    y'' = g(x,y,y')

Next, y'' is substituted for u' and y' is substituted for u

    y' = u
    u' = g(x,y,u)


The Runge-Kutta method requires these equations:
    Y_1 = Y_0 + (h/6)*(f_0 + 2*f_1 + 2*f_2 + f_3)
    U_1 = U_0 + (h/6)*(g_0 + 2*g_1 + 2*g_2 + g_3)

Where:
    f_0 = f(x_0, Y_0, U_0)                                g_0 = g(x_0, Y_0, U_0)
    f_1 = f(x_0 + h/2, Y_0 + (h/2)*f_0, U_0 + (h/2)*g_0)  g_1 = g(x_0 + h/2, Y_0 + (h/2)*f_0, U_0 + (h/2)*g_0)
    f_2 = f(x_0 + h/2, Y_0 + (h/2)*f_1, U_0 + (h/2)*g_1)  g_2 = g(x_0 + h/2, Y_0 + (h/2)*f_1, U_0 + (h/2)*g_1)
    f_3 = f(x_0 + h, Y_0 + h*f_2, U_0 + h*g_2)            g_3 = g(x_0 + h, Y_0 + h*f_2, U_0 + h*g_2)

    and h is the x increment

Note that I am borrowing the problem from 2.9, where Improved Euler is applied to solve the problem,
    y'' = 2y - (1/x)y' + x^2

Results should appear as follows:

X increment: 0.1

  x        U         Y
1.0 3.000000 -2.000000
1.1 2.478207 -1.727346
1.2 2.090329 -1.499919
1.3 1.810946 -1.305680
1.4 1.622761 -1.134702
1.5 1.514136 -0.978485
1.6 1.477556 -0.829480
1.7 1.508645 -0.680724
1.8 1.605527 -0.525563
1.9 1.768421 -0.357423
2.0 1.999384 -0.169614

Where Improved Euler gives a solution accurate to 3 decimal places with delta_x = 0.01,
Runge-Kutta gives a similar result with delta_x = 0.1

That's hella cool!
"""


# Initial conditions
x_0 = 1
y_0 = -2
y_dot_0 = 3
delta_x = 0.1


# Let y' = u, and y'' = u'
def u_dot(x, y, u):
    return 2*y - (1/x)*u + x**2


def Runge_Kutta(x, y, u, delta_x, func):
    """
    Runge-Kutta (RK4) implementation for second order differential equations
    """
    f_0 = u
    g_0 = func(x,y,u)
    
    f_1 = u + (delta_x/2)*g_0
    g_1 = func(x + delta_x/2, y + (delta_x/2)*f_0, u + (delta_x/2)*g_0)
    
    f_2 = u + (delta_x/2)*g_1
    g_2 = func(x + delta_x/2, y + (delta_x/2)*f_1, u + (delta_x/2)*g_1)
    
    f_3 = u + delta_x*g_2
    g_3 = func(x + delta_x, y + delta_x*f_2, u + delta_x*g_2)

    new_y = y + (delta_x/6)*(f_0 + 2*f_1 + 2*f_2 + f_3)
    new_u = u + (delta_x/6)*(g_0 + 2*g_1 + 2*g_2 + g_3)

    return (new_y, new_u)


def Y(x):
    y = y_0
    u = y_dot_0

    #Initialize output lists for dataframe
    X_list = []
    U_list = []
    Y_list = []
    
    # Iterate over range from x_0 to input x
    for x_i in np.arange(x_0, x + delta_x, delta_x):
        # Append values to output lists
        X_list.append(x_i)
        U_list.append(u)
        Y_list.append(y)

        # Use Runge-Kutta method to solve differential equation
        y, u = Runge_Kutta(x_i, y, u, delta_x, u_dot)
        
    return (X_list, Y_list, U_list)


def main(x):
    # Print result dataframe
    x_list, y_list, v_list = Y(x)
    data = {
        "x" : x_list,
        "U" : v_list,
        "Y" : y_list
    }
    df = pd.DataFrame(data)
    print("X increment: " + str(delta_x) + "\n")
    print(df.to_string(index=False))

if __name__ == "__main__":
    tic = time.time()
    main(2)
    toc = time.time()
    print("\nElapsed time: " + str((toc - tic)*1000) + " ms")