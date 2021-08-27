"""

Simulator of wood-chip tank with PI controller

Author: Finn Haugen, USN, finn.haugen@usn.no

Updated 2020 11 07

"""

# %% Import

import matplotlib.pyplot as plt
import numpy as np

# %% Def of PI controller function:

def fun_pi_contr(y_sp_k, y_k, u_i_km1, contr_params, dt):
    
    (Kc, Ti, u_man, u_min, u_max) = contr_params
    e_k = y_sp_k - y_k  # Control error
    u_p_k = Kc*(y_sp_k - y_k)  # P term
    u_i_k = u_i_km1 + (Kc*dt/Ti)*e_k  # PI term
    u_i_min = u_min - u_man
    u_i_max = u_max - u_man

    # Limiting ui (integrator anti windup):
    u_i_k = np.clip(u_i_k, u_i_min, u_i_max)

    u_k = u_man + u_p_k + u_i_k  # Total control signal

    # Limiting applied control signal:
    u_k = np.clip(u_k, u_min, u_max)
    
    return (u_k, u_i_k)

# %% Def of process model simulator

def process_sim(h_k,
                u_k,
                delay_array_k,
                process_params,
                dt):

    # Reading process params:
    (rho, A, t_delay, h_min, h_max) = process_params
    
    # Limiting the state:
    h_k = np.clip(h_k, h_min, h_max)

    # Time delay:
    u_delayed_k = delay_array_k[-1]
    delay_array_k[1:] = delay_array_k[0:-1]
    delay_array_k[0] = u_k
    delay_array_kp1 = delay_array_k

    # Euler-forward integration (Euler step):
    w_in_k = u_delayed_k  # kg/s
    dh_dt_k = (1/(rho*A))*(w_in_k - w_out_k)
    h_kp1 = h_k + dt*dh_dt_k

    return (h_kp1, delay_array_kp1)

# %% Time settings:

dt = 1  # Time-step [s]
t_start = 0.0  # [s]
t_stop = 5000.0  # [s]
N_sim = int((t_stop-t_start)/dt) + 1

# %% Arrays for plotting:

t_array = np.zeros(N_sim)
u_array = np.zeros(N_sim)
h_sp_array = np.zeros(N_sim)
h_array = np.zeros(N_sim)
w_out_array = np.zeros(N_sim)

# %% Process params:

rho = 145  # [kg/m^3]
A = 13.4  # [m^2]
t_delay = 250.0  # [s]
h_min = 0  # [m]
h_max = 15  # [m]
process_params = (rho, A, t_delay, h_min, h_max)

# %% Initialization of time delay:

u_delayed_init = 25  # [kg/s]
N_delay = int(round(t_delay/dt)) + 1
delay_array_k = np.zeros(N_delay) + u_delayed_init

# %% PI controller settings:

Kc = 4  # [(kg/s)/m]
Ti = 1000  # [s]
print('Kc = ', f'{Kc:.2e}')
print('Ti = ', f'{Ti:.2e}')
u_man = 25  # [kg/s]
u_max = 50  # [kg/s]
u_min = 0  # [kg/s]
contr_params = (Kc, Ti, u_man, u_min, u_max)

# %% Initial state:

h_k = 10.0  # m
u_i_km1 = 0  # [kg/s]

# %% Simulation for-loop:

for k in range(0, N_sim):

    t_k = k*dt

    if t_k <= 500:
        h_sp_k = 10  # [m]
        w_out_k = 25  # [kg/s]
    else:
        h_sp_k = 11  # [m]
        w_out_k = 25  # [kg/s]
    
    # PI controller:
    (u_k, u_i_k) = fun_pi_contr(h_sp_k,
                                h_k, 
                                u_i_km1, 
                                contr_params,
                                dt)

    # Process simulator:
    (h_kp1, delay_array_kp1) = process_sim(h_k,
                                       u_k,
                                       delay_array_k,
                                       process_params,
                                       dt)
    
    # Storage for plotting:
    t_array[k] = t_k
    u_array[k] = u_k
    h_sp_array[k] = h_sp_k
    h_array[k] = h_k
    w_out_array[k] = w_out_k

    # Time index shift:
    h_k = h_kp1
    u_i_km1 = u_i_k
    delay_array_k = delay_array_kp1

# %% Plotting:

plt.close('all')
plt.figure(1)

plt.subplot(3, 1, 1)
plt.plot(t_array, h_sp_array, 'r', label='h_sp')
plt.plot(t_array, h_array, 'b', label='h')
plt.legend()
plt.grid()
plt.xlim(t_start, t_stop)
# plt.ylim(20, 25)
plt.xlabel('t [s]')
plt.ylabel('[m]')

plt.subplot(3, 1, 2)
plt.plot(t_array, u_array, 'g', label='u')
plt.plot(t_array, u_array*0 + u_min, 'm', label='u_min')
plt.plot(t_array, u_array*0 + u_max, 'm', label='u_max')
plt.legend()
plt.grid()
plt.xlim(t_start, t_stop)
# plt.ylim(20, 25)
plt.xlabel('t [s]')
plt.ylabel('[kg/s]')

plt.subplot(3, 1, 3)
plt.plot(t_array, w_out_array, 'c', label='w_out')
plt.legend()
plt.grid()
plt.xlim(t_start, t_stop)
# plt.ylim(20, 25)
plt.xlabel('t [s]')
plt.ylabel('[kg/s]')

# plt.savefig('plot_sim_woodchiptank_pi_control.pdf')
# plt.savefig('plot_sim_woodchiptank_pi_control.png')
plt.show()