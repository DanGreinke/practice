import numpy as np
import pandas as pd

pd.options.display.float_format = '{:.7f}'.format

"""
From Chapter 3.2 of "A First Course in Differential Equations" (1975) by Frank G. Hagin

The chapter introduces the 4th order Runge-Kutta method for numerically solving differential equations. 
This enables highly accurate approximations with low computation al cost.

Solve:

y' = 2xy^2

with initial condition y(0) = 0.5

Results should appear as follows:

  x  1/(2-x^2)    Euler  Improved Euler  Runge-Kutta
0.0   0.500000 0.500000        0.500000     0.500000
0.1   0.502513 0.500000        0.502500     0.502513
0.2   0.510204 0.505000        0.510177     0.510204
0.3   0.523560 0.515201        0.523513     0.523560
0.4   0.543478 0.531127        0.543397     0.543478
0.5   0.571429 0.553695        0.571284     0.571429
0.6   0.609756 0.584352        0.609486     0.609756
0.7   0.662252 0.625328        0.661720     0.662252
0.8   0.735294 0.680073        0.734192     0.735295
0.9   0.840336 0.754073        0.837895     0.840336
1.0   1.000000 0.856426        0.994063     0.999996

"""

x_0 = 0
y_0 = 0.5
delta_x = 0.05

y_dot = lambda x, y: 2 * x * pow(y, 2)

def Y_Euler(x):
    Y = y_0
    result_list = [Y]
    for x_i in np.arange(x_0, x, delta_x):
        
        Y += delta_x * y_dot(x_i,Y)
        
        result_list.append(Y)
    return result_list

def Y_Euler_Improved(x):
    Y = y_0
    result_list = [Y]
    for x_i in np.arange(x_0, x, delta_x):
        
        Y_bar = Y + delta_x * y_dot(x_i,Y)
        
        Y += (1/2) * delta_x * (y_dot(x_i,Y) + y_dot(x_i + delta_x,Y_bar))
        
        result_list.append(Y)
    return result_list

def Y_Runge_Kutta(x):
    Y = y_0
    result_list = [Y]
    for x_i in np.arange(x_0, x, delta_x):
        
        f_0 = y_dot(x_i,Y)
        f_1 = y_dot(x_i + (delta_x/2), Y + (delta_x/2)*f_0)
        f_2 = y_dot(x_i + (delta_x/2), Y + (delta_x/2)*f_1)
        f_3 = y_dot(x_i + delta_x, Y + delta_x*f_2)
        
        Y += (delta_x / 6)*(f_0 + 2*f_1 + 2*f_2 + f_3)

        result_list.append(Y)
    return result_list

def Y_analytic_solution(x):
    result_list = []
    for x_i in np.arange(0, x + delta_x, delta_x):
        result_list.append(1 / (2 - pow(x_i,2)))
    return result_list

def main(x):
    data = {
        "x": np.arange(0,1 + delta_x,delta_x),
        "1/(2-x^2)":Y_analytic_solution(x),
        "Euler": Y_Euler(x),
        "Improved Euler": Y_Euler_Improved(x),
        "Runge-Kutta": Y_Runge_Kutta(x)
    }
    df = pd.DataFrame(data)
    print(df.to_string(index=False))

if __name__ == "__main__":
    main(1)