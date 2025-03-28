import numpy as np
import pandas as pd

# From Chapter 1.8 of "A First Course in Differential Equations" (1975) by Frank G. Hagin
# 
# This chapter introduces numerical methods for solving differential equations and contrasts the Euler Method
# with an "Improved Euler Method"
#
#     Euler:          Y_1 = Y_0 + (x_1 - x_0)*f(x_0,Y_0)
#
#     Improved Euler: Y_1 = Y_0 + 0.5*(x_1 - x_0)*[f(x_0,Y_0) + f(x_1,Y_bar)]
#              where, Y_bar = Y_0 + (x_1 - x_0)*f(x_0,Y_0)
#
# The code below reproduces the results from tables 8-1 and 8-3 (p. 47 and 49), highlighting the contrast between 
# these two methods in solving:
#
#     y' = exp(-x^2)/y
#
# Results should appear as follows:
#
#   x     Euler  Euler Improved
# 0.0  1.000000        1.000000
# 0.1  1.100000        1.095002
# 0.2  1.190005        1.180735
# 0.3  1.270743        1.257628
# 0.4  1.342664        1.325992
# 0.5  1.406131        1.386133
# 0.6  1.461517        1.438412
# 0.7  1.509253        1.483264
# 0.8  1.549844        1.521208
# 0.9  1.583867        1.552836
# 1.0  1.611954        1.578791

# initial conditions
x_0 = 0
y_0 = 1.0
delta_x = 0.1

def y_dot(x,y):
    return np.exp(-x**2)/y

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
    for x_i in np.arange(0, x, delta_x):
        Y_bar = Y + delta_x * y_dot(x_i,Y)
        Y += (1/2) * delta_x * (y_dot(x_i,Y) + y_dot(x_i + delta_x,Y_bar))
        result_list.append(Y)
    return result_list

def main(x):
    data = {
        "x" : [x for x in np.arange(0, x + delta_x, delta_x)],
        "Euler" : Y_Euler(x),
        "Euler Improved" : Y_Euler_Improved(x)
    }
    df = pd.DataFrame(data)
    print(df.to_string(index=False))

if __name__ == "__main__":
    main(1.0)