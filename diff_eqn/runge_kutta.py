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

        x  1/(2-x^2)     Euler  Improved Euler  Runge-Kutta
0.0000000  0.5000000 0.5000000       0.5000000    0.5000000
0.1000000  0.5025126 0.5000000       0.5025000    0.5025126
0.2000000  0.5102041 0.5050000       0.5101772    0.5102041
0.3000000  0.5235602 0.5152010       0.5235132    0.5235603
0.4000000  0.5434783 0.5311269       0.5433973    0.5434784
0.5000000  0.5714286 0.5536946       0.5712841    0.5714288
0.6000000  0.6097561 0.5843524       0.6094856    0.6097564
0.7000000  0.6622517 0.6253285       0.6617198    0.6622522
0.8000000  0.7352941 0.6800735       0.7341918    0.7352948
0.9000000  0.8403361 0.7540735       0.8378954    0.8403365
1.0000000  1.0000000 0.8564263       0.9940628    0.9999956

"""

x_0 = 0
y_0 = 0.5
delta_x = 0.1

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