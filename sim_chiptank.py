"""

Simulator of wood-chip tank

Author: Finn Haugen, USN, finn.haugen@usn.no

Updated 2020 09 27

"""

# %% Import
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

# %% Time settings:

ts = 1  # Time-step [s]
t_start = 0.0  # [s]
t_stop = 2000.0  # [s]
N_sim = int((t_stop-t_start)/ts) + 1

# %% Process params:

rho = 145  # [kg/m^3]
A = 13.4  # [m^2]
t_delay = 250.0  # [s]
h_min = 0  # [m]
h_max = 15  # [m]

#P controller
z = np.zeros(N_sim)
e = np.zeros(N_sim)
sps=np.zeros(N_sim+1)
SP=10
u_man = 1
Kp=50


# %% Initialization of time delay:

u_delayed_init = 25  # [kg/s]
N_delay = int(round(t_delay/ts)) + 1
delay_array = np.zeros(N_delay) + u_delayed_init

# %% Arrays for plotting:

t_array = np.zeros(N_sim+1)
h_array = np.zeros(N_sim+1)
u_array = np.zeros(N_sim+1)

# %% Initial state:

h_k = 10.0  # m

# %% Simulation for-loop:

for k in range(0, N_sim):
    
    t_k = k*ts
    
    error = SP-h_k
    h_kp1= u_man + Kp*error
    
    if t_k <= 500:
        u_k = 25  # kg/s
        w_out_k = 25  # kg/s
    else:
        u_k = 40  # kg/s
        w_out_k = 25  # kg/s
    
    # Time delay:
    u_delayed_k = delay_array[-1]
    delay_array[1:] = delay_array[0:-1]
    delay_array[0] = u_k

    # Euler-forward integration (Euler step):
    w_in_k = u_delayed_k  # kg/s
    dh_dt_k = (1/(rho*A))*(w_in_k - w_out_k)
    h_kp1 = h_k + ts*dh_dt_k
    h_kp1 = np.clip(h_kp1, h_min, h_max)
    
    # Storage for plotting:
    t_array[k] = t_k
    u_array[k] = u_k
    h_array[k] = h_k
    
    # Time shift:
    h_k = h_kp1
    
    sps[k+1] = SP
    
    
# %% Plotting:

plt.close('all')
plt.figure(1)

plt.subplot(2, 1, 1)
plt.plot(t_array,sps,'k-',linewidth=2)
plt.plot(t_array,h_array,'b-',linewidth=3)
plt.legend()
plt.grid()
plt.xlim(t_start, t_stop)
# plt.ylim(20, 25)



