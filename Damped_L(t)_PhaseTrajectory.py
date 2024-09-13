import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

g = 9.81
L_0 = 3
w = 1
delta = 0.08
c_c = 911.47
c_u = 400
c_o = 1500
m = 28
k = 1

# Underdamped
def Standing_Pumping_New(phi, omega, t):
    L = L_0 + delta * np.cos(2 * w * t)
    L_dot = - 2 * w * delta * np.sin(2 * w * t)
    d_phi_dt = omega
    d_omega_dt = (-g * np.sin(phi))/L - (2 * L_dot * d_phi_dt)/L - (k * phi)/(m * L**2) - (c_u * d_phi_dt)/(m * L**2)
    return d_phi_dt, d_omega_dt

def euler_method(phi0, omega0, t0, dt, steps):
    phi_values = [phi0]
    omega_values = [omega0]
    t_values = [t0]

    for _ in range(steps):
        d_phi_dt, d_omega_dt = Standing_Pumping_New(phi0, omega0, t0)
        phi0 += dt * d_phi_dt
        omega0 += dt * d_omega_dt
        t0 += dt
     
    
        phi_values.append(phi0)
        omega_values.append(omega0)
        t_values.append(t0)
        
    return phi_values, omega_values, t_values

phi0 = np.pi / 4  
omega0 = 0.0
t0 = 0
dt = 0.008
steps = 6000

phi_values, omega_values, t_values = euler_method(phi0, omega0, t0, dt, steps)

# Overdamped
def Standing_Pumping_New2(phi2, omega2, t):
    L = L_0 + delta * np.cos(2 * w * t)
    L_dot2 = - 2 * w * delta * np.sin(2 * w * t)
    d_phi_dt2 = omega2
    d_omega_dt2 = (-g * np.sin(phi2))/L - (2 * L_dot2 * d_phi_dt2)/L - (k * phi2)/(m * L**2) - (c_o * d_phi_dt2)/(m * L**2)
    return d_phi_dt2, d_omega_dt2

def euler_method2(phi_2, omega_2, t_2, dt, steps):
    phi_values2 = [phi_2]
    omega_values2 = [omega_2]
    t_values2 = [t_2]

    for _ in range(steps):
        d_phi_dt2, d_omega_dt2 = Standing_Pumping_New2(phi_2, omega_2, t_2)
        phi_2 += dt * d_phi_dt2
        omega_2 += dt * d_omega_dt2
        t_2 += dt
     
    
        phi_values2.append(phi_2)
        omega_values2.append(omega_2)
        t_values2.append(t_2)
        
    return phi_values2, omega_values2, t_values2

phi_2 = np.pi / 4  
omega_2 = 0.0
t_2 = 0
dt = 0.008
steps = 6000

phi_values2, omega_values2, t_values2 = euler_method2(phi_2, omega_2, t_2, dt, steps)

# Critically damped
def Standing_Pumping_New1(phi1, omega1, t):
    L = L_0 + delta * np.cos(2 * w * t)
    L_dot1 = - 2 * w * delta * np.sin(2 * w * t)
    d_phi_dt1 = omega1
    d_omega_dt1 = (-g * np.sin(phi1))/L - (2 * L_dot1 * d_phi_dt1)/L - (k * phi1)/(m * L**2) - (c_c * d_phi_dt1)/(m * L**2)
    return d_phi_dt1, d_omega_dt1

def euler_method1(phi_1, omega_1, t_1, dt, steps):
    phi_values1 = [phi_1]
    omega_values1 = [omega_1]
    t_values1 = [t_1]

    for _ in range(steps):
        d_phi_dt1, d_omega_dt1 = Standing_Pumping_New1(phi_1, omega_1, t_1)
        phi_1 += dt * d_phi_dt1
        omega_1 += dt * d_omega_dt1
        t_1 += dt
     
    
        phi_values1.append(phi_1)
        omega_values1.append(omega_1)
        t_values1.append(t_1)
        
    return phi_values1, omega_values1, t_values1

phi_1 = np.pi / 4  
omega_1 = 0.0
t_1 = 0
dt = 0.008
steps = 6000

phi_values1, omega_values1, t_values1 = euler_method1(phi_1, omega_1, t_1, dt, steps)

# Set up custom color palette
colors = plt.cm.viridis(np.linspace(0, 1, 2))


# Plot phase trajectory
plt.figure(figsize=(20, 10))
plt.plot(phi_values, omega_values, label='$c_u$ = 400', color='blue', linewidth=2, linestyle=':', marker='s', markersize=1)
plt.plot(phi_values1, omega_values1, label='$c_c$ = 911.47', color='green', linewidth=2, linestyle='-.', marker='s', markersize=1)
plt.plot(phi_values2, omega_values2, label='$c_o$ = 1500', color='red', linewidth=2, linestyle='--', marker='s', markersize=1)
plt.xlabel('$\phi$', fontsize=40, fontweight='bold', color='green')
plt.ylabel('$\dot{\phi}$', fontsize=40, fontweight='bold', color='green')
plt.grid(True, linestyle='--', color='gray', alpha=0.5)
plt.legend(fontsize=20, loc='upper right')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.gca().set_facecolor('#f0f0f0')  # Background color
plt.show()