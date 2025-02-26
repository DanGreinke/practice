import numpy as np
"""
For the game, Oxygen Not Included.

Simulates one radbolt generator dumping its heat into an infinite nuclear waste liquid storage, while 1kg/sec of nuclear waste is pumped in at 100 C (373 Kelvin). Initial nuclear waste is 1000kg at 373 Kelvin.
Neglects adjacent blocks such as tempshift plates, pipes, etc that would be present in-game.
"""
C_nuke_waste = 7.440 # kDTU/kg/C
Heat_radbolt = 5.000 # kDTU/sec

nuke_waste_in = 1 # kg/sec
nuke_waste_initial = 1000 # kg
T_initial = 373 # Kelvin, or 100 C

Heat_initial = T_initial*C_nuke_waste*nuke_waste_initial
HC_initial = C_nuke_waste*nuke_waste_initial

delta_t = 1
dT_dt_initial = (1/(C_nuke_waste*nuke_waste_initial)) * (Heat_radbolt)

def get_H_dot_HC_dot(t):
    """
    Change in total heat in system divided by change in total heat capacity
    """
    Heat_added = (Heat_radbolt + (C_nuke_waste*nuke_waste_in*T_initial))*t
    HC_added = C_nuke_waste*nuke_waste_in*t
    return Heat_added, HC_added

def T(t):
    T_dot = dT_dt_initial
    H = Heat_initial
    HC = HC_initial
    for time in np.arange(0, t, delta_t):
        H += get_H_dot_HC_dot(time)[0]*delta_t
        HC += get_H_dot_HC_dot(time)[1]*delta_t

    return H/HC

print(T(1000*600)) # 1000 cycles, 600 seconds each.