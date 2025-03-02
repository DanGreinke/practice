import numpy as np
import pandas as pd
import math
"""
From Chapter 2.9 of "A First Course in Differential Equations" (1975) by Frank G. Hagin

This chapter applies the improved Euler method to solving second order differential equations

    y'' = f(x,y,y')

Next, y'' is substituted for v' and y' is substituted for v

    y' = v
    v' = f(x,y,v)


The Improved Euler method requires these equations:
    y_bar = Y_0 + h*V_0
    v_bar = V_0 + h*f(x_0, Y_0, V_0)
    Y_1 = Y_0 + (h/2)*(V_0 + v_bar)
    V_1 = V_0 + (h/2)*(f(x_0, Y_0, V_0) + f(x_0 + h, y_bar, v_bar))

Where h is the x increment
    
Results should appear as follows, reproducing Table 9-1 (p. 121):
______________________
X increment: 0.01

  x        V         Y
1.0 3.000000 -2.000000
1.1 2.478228 -1.727364
1.2 2.090362 -1.499950
1.3 1.810985 -1.305721
1.4 1.622801 -1.134750
1.5 1.514173 -0.978540
1.6 1.477587 -0.829542
1.7 1.508666 -0.680793
1.8 1.605537 -0.525640
1.9 1.768417 -0.357509
2.0 1.999361 -0.169712
______________________
"""


# Initial conditions
x_0 = 1
y_0 = -2
y_dot_0 = 3
delta_x = 0.01
display_increment = 0.1


# Let y' = v, and y'' = v'
def v_dot(x, y, v):
    return 2*y - (1/x)*v + x**2


def is_display_increment(num):
    # Note: This breaks for delta_x < 0.01
    #
    # Checks that each row is a multiple of the increment, to abbreviate
    # result table for smaller increment sizes
    tmp = num/display_increment
    return math.isclose(tmp, int(tmp), abs_tol=1e-9)


def Y(x):
    y = y_0
    v = y_dot_0

    #Initialize output lists for dataframe
    X_list = []
    V_list = []
    Y_list = []
    
    # Iterate over range from x_0 to input x
    for x_i in np.arange(x_0, x + delta_x, delta_x):
        # Limit display dataframe so we can summarize results for smaller increments
        if is_display_increment(x_i) == True:
            X_list.append(x_i)
            V_list.append(v)
            Y_list.append(y)
        
        # Use Improved Euler method to solve differential equation
        v_bar = v + delta_x*v_dot(x_i, y, v)
        y_bar = y + delta_x*v
        new_v = v + (delta_x/2)*(v_dot(x_i, y, v) + v_dot(x_i + delta_x, y_bar, v_bar))
        new_y = y + (delta_x/2)*(v + v_bar)
        
        # Update v and y
        v = new_v
        y = new_y

    return (X_list, Y_list, V_list)


def main(x):
    # Print result dataframe
    x_list, y_list, v_list = Y(x)
    data = {
        "x" : x_list,
        "V" : v_list,
        "Y" : y_list
    }
    df = pd.DataFrame(data)
    print("X increment: " + str(delta_x) + "\n")
    print(df.to_string(index=False))

if __name__ == "__main__":
    main(2)