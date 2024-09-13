import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

g = 9.81
L_0_1 = 3
L_0_2 = 3
L_0_3 = 3
w = 1
delta = 0.08
c_c = 911.47
c_u = 400
c_o = 1500
m = 28
k = 1


def Standing_Pumping_New_1(phi, omega, t):
    L = L_0_1 + delta * np.cos(2 * w * t)
    L_dot = - 2 * w * delta * np.sin(2 * w * t)
    d_phi_dt_1 = omega
    d_omega_dt_1 = (-g * np.sin(phi))/L - (2 * L_dot * d_phi_dt_1)/L - (k * phi)/(m * L**2) - (c_c * d_phi_dt_1)/(m * L**2)
    return d_phi_dt_1, d_omega_dt_1

def euler_method_1(phi0_1, omega0_1, t0_1, dt, steps):
    phi_values_1 = [phi0_1]
    omega_values_1 = [omega0_1]
    t_values_1 = [t0_1]

    for _ in range(steps):
        d_phi_dt_1, d_omega_dt_1 = Standing_Pumping_New_1(phi0_1, omega0_1, t0_1)
        phi0_1 += dt * d_phi_dt_1
        omega0_1 += dt * d_omega_dt_1
        t0_1 += dt
     
    
        phi_values_1.append(phi0_1)
        omega_values_1.append(omega0_1)
        t_values_1.append(t0_1)
        
    return phi_values_1, omega_values_1, t_values_1

def Standing_Pumping_New_2(phi, omega, t):
    L = L_0_2 + delta * np.cos(2 * w * t)
    L_dot = - 2 * w * delta * np.sin(2 * w * t)
    d_phi_dt_2 = omega
    d_omega_dt_2 = (-g * np.sin(phi))/L - (2 * L_dot * d_phi_dt_2)/L - (k * phi)/(m * L**2) - (c_o * d_phi_dt_2)/(m * L**2)
    return d_phi_dt_2, d_omega_dt_2

def euler_method_2(phi0_2, omega0_2, t0_2, dt, steps):
    phi_values_2 = [phi0_2]
    omega_values_2 = [omega0_2]
    t_values_2 = [t0_2]

    for _ in range(steps):
        d_phi_dt_2, d_omega_dt_2 = Standing_Pumping_New_2(phi0_2, omega0_2, t0_2)
        phi0_2 += dt * d_phi_dt_2
        omega0_2 += dt * d_omega_dt_2
        t0_2 += dt
     
    
        phi_values_2.append(phi0_2)
        omega_values_2.append(omega0_2)
        t_values_2.append(t0_2)
        
    return phi_values_2, omega_values_2, t_values_2

def Standing_Pumping_New_3(phi, omega, t):
    L = L_0_3 + delta * np.cos(2 * w * t)
    L_dot = - 2 * w * delta * np.sin(2 * w * t)
    d_phi_dt_3 = omega
    d_omega_dt_3 = (-g * np.sin(phi))/L - (2 * L_dot * d_phi_dt_3)/L - (k * phi)/(m * L**2) - (c_u * d_phi_dt_3)/(m * L**2)
    return d_phi_dt_3, d_omega_dt_3

def euler_method_3(phi0_3, omega0_3, t0_3, dt, steps):
    phi_values_3 = [phi0_3]
    omega_values_3 = [omega0_3]
    t_values_3 = [t0_3]

    for _ in range(steps):
        d_phi_dt_3, d_omega_dt_3 = Standing_Pumping_New_3(phi0_3, omega0_3, t0_3)
        phi0_3 += dt * d_phi_dt_3
        omega0_3 += dt * d_omega_dt_3
        t0_3 += dt
     
    
        phi_values_3.append(phi0_3)
        omega_values_3.append(omega0_3)
        t_values_3.append(t0_3)
        
    return phi_values_3, omega_values_3, t_values_3

t0_1 = t0_2 = t0_3 = 0
phi0_2 = phi0_3 = phi0_1 = np.pi / 4  
omega0_2 = omega0_1 = omega0_3 = 0.0
dt = 0.003
steps = 6000


phi_values_1, omega_values_1, t_values_1 = euler_method_1(phi0_1, omega0_1, t0_1, dt, steps)
phi_values_2, omega_values_2, t_values_2 = euler_method_2(phi0_2, omega0_2, t0_2, dt, steps)
phi_values_3, omega_values_3, t_values_3 = euler_method_3(phi0_3, omega0_3, t0_3, dt, steps)


# Plot phase trajectory
plt.figure(figsize=(20, 10))
plt.plot(t_values_3, phi_values_3, label='$c_u$ = 400', color="blue", linewidth=2, linestyle=':', markersize=1)
plt.plot(t_values_1, phi_values_1, label='$c_c$ = 911.47', color="green", linewidth=2, linestyle='-.', markersize=1)
plt.plot(t_values_2, phi_values_2, label='$c_o$ = 1500', color="red", linewidth=2, linestyle='--', markersize=1)
plt.xlabel('t', fontsize=40, fontweight='bold', color='green')
plt.ylabel('$\phi$', fontsize=40, fontweight='bold', color='green')
plt.grid(True, linestyle='--', color='gray', alpha=0.5)
plt.legend(fontsize=20, loc='upper right')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.gca().set_facecolor('#f0f0f0')  # Background color
plt.show()