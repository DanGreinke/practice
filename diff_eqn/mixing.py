import numpy as np
# https://sites.math.washington.edu//~conroy/m125-general/mixingTankExamples/mixingTankExamples01.pdf
# Constraints
S_initial = 10 # kg
dS_dt_initial = S_initial / 10
delta_t = 0.003


print("Example 1")
"""
A tank has pure water flowing into it at 10 l/min. The contents of the tank are kept
thoroughly mixed, and the contents flow out at 10 l/min. Initially, the tank contains 10 kg of
salt in 100 l of water.

How much salt will there be in the tank after 30 minutes?
"""
print(50*"-")
def get_S_dot_1(S):
    return -S / 10

def s(t):
    S = S_initial
    for time in np.arange(0, t, delta_t):
        S += get_S_dot_1(S) * delta_t
    return S

print(s(30))
print(50*"-")

print("Example 2")
"""
A tank has pure water flowing into it at 10 l/min. The contents of the tank are kept
thoroughly mixed, and the contents flow out at 10 l/min. Salt is added to the tank at the rate
of 0.1 kg/min. Initially, the tank contains 10 kg of salt in 100 l of water.

How much salt is in the tank after 30 minutes?
"""
print(50*"-")

def get_S_dot_2(S):
    return -S/10 + 0.1

def s_2(t):
    S = S_initial
    for time in np.arange(0, t, delta_t):
        S += get_S_dot_2(S) * delta_t
    return S

print(s_2(30))
print(50*"-")

print("Example 3")
"""
A tank has pure water flowing into it at 12 l/min. The contents of the tank are kept
thoroughly mixed, and the contents flow out at 10 l/min. Initially, the tank contains 10 kg of
salt in 100 l of water.

In this case, the inflow rate is greater than the outflow rate. As a result, the volume is not
constant.


"""
print(50*"-")

def get_S_dot_3(S, t):
    return -10*S/(100 + 2*t)

def s_3(t):
    S = S_initial
    for time in np.arange(0, t, delta_t):
        S += get_S_dot_3(S, time) * delta_t
    return S

print(s_3(30))
print(50*"-")
