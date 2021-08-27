"""

Simulator of heater water tank

Finn Aakre Haugen
finn.haugen@usn.no

2020 10 03

"""

#%% Imports:

import numpy as np
import matplotlib.pyplot as plt

#%% Model params:

c = 4200 # [J/(kg*K)]
rho = 1000 # [kg/m3]
V = 0.2 # [m3]
U = 1000 # [W/K]
F = 0.25e-3  # [m3/s]
T_in = 20  # [deg C]
T_env = 20  # [deg C]
t_delay = 60  # [s]

T_min = 0
T_max = 100

#%% Sim time settings:

dt = 1
t_start = 0
t_stop = 6000
N_sim = int((t_stop - t_start)/dt) + 1

#%% Array for transport delay:
    
P_delayed_init = 0  # [W]
N_delay = int(round(t_delay/dt)) + 1
P_delay_array = np.zeros(N_delay) + P_delayed_init

#%% Preallocation of arrays for storing:

t_array = np.zeros(N_sim)
T_array = np.zeros(N_sim)
T_in_array = np.zeros(N_sim)
T_env_array = np.zeros(N_sim)
P_array = np.zeros(N_sim)

#%% Sim loop:

P0 = 0  # [W]
P1 = 1000  # [W]
T_k = T_init = 20  # [deg C]
t0 = 500  # [s] Time of step in P

for k in range(0, N_sim):

    t_k = k*dt
    
    if (0 <= t_k <= t0):
        P_k = P0
        T_in_k = T_in
        T_env_k = T_env
    else:
        P_k = P1
        T_in_k = T_in
        T_env_k = T_env
    # Moving array elements one step:
    P_delayed_k = P_delay_array[-1]
    P_delay_array[1:] = P_delay_array[0:-1]
    P_delay_array[0] = P_k

    
    dT_dt_k = ((1/(c*rho*V))
               *(P_delayed_k
                 + (c*rho*F)*(T_in-T_k) 
                 + U*(T_env-T_k)))
    T_kp1 = T_k + dt*dT_dt_k
    T_kp1 = np.clip(T_kp1, T_min, T_max)
    
    t_array[k] = t_k
    T_array[k] = T_k
    T_in_array[k] = T_in_k
    T_env_array[k] = T_env_k
    P_array[k] = P_k

    
    # Time index shift:
    T_k = T_kp1


# %% Plotting:

plt.close('all')
plt.figure(1)

plt.subplot(2, 1, 1)
plt.plot(t_array, T_array, 'r', label='T')
plt.plot(t_array, T_in_array, 'b', label='T_in')
plt.plot(t_array, T_env_array, 'g', label='T_env')
plt.legend()
plt.grid()
plt.xlabel('t [s]')
plt.ylabel('[C]')

plt.subplot(2, 1, 2)
plt.plot(t_array, P_array, 'm', label='P')
plt.legend()
plt.grid()
plt.xlabel('t [s]')
plt.ylabel('[W]')

# plt.savefig('plot_sim_heated_water_tank.pdf')
plt.show()

#%% Dynamics:

# From model:
gain_model = 1/(c*rho*F + U)  # [deg C/W]
time_const_model = c*rho*V/(c*rho*F + U)  # [s]
time_delay_model = t_delay  # [s]

# From simulation:
dT = T_array[-1] - T_init
dP = P1 - P0
gain_sim = dT/dP

T_at_time_const = T_init + (1-np.exp(-1))*dT
k0 = np.min(np.where(T_array > T_at_time_const))
time_const_sim = t_array[k0] - t_delay - t0

k1 = np.min(np.where(T_array > T_init))
time_delay_sim = t_array[k1] - t0

# Comparison:
print('gain_model =', f'{gain_model:.3e}')
print('gain_sim =', f'{gain_sim:.3e}')
print('time_const_model =', f'{time_const_model:.3e}')
print('time_const_sim =', f'{time_const_sim:.3e}')
print('time_delay_model =', f'{time_delay_model:.3e}')
print('time_delay_sim =', f'{time_delay_sim:.3e}')
