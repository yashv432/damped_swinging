import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

g = 9.81
L_0_1 = 3
L_0_2 = 3
L_0_3 = 3
c_c = 911.47
c_u = 400
c_o = 1500
m = 28
k = 1

#Underdamping
def L(theta):
    return L_0_1 - L_0_1 * np.abs(theta) / np.pi

def Standing_Pumping(theta, omega):
    dtheta_dt = omega
    dL_dt = (-L_0_1 * theta * omega) / (np.pi * np.abs(theta))
    domega_dt = (- c_u / (m * L_0_1**2))*dtheta_dt - (g / L_0_1 + k / (m * L_0_1**2))*theta
    return dtheta_dt, domega_dt

def euler_method(theta0, omega0, dt, steps):
    theta_values = [theta0]
    omega_values = [omega0]

    for _ in range(steps):
        dtheta_dt, domega_dt = Standing_Pumping(theta0, omega0)
        theta0 += dt * dtheta_dt
        omega0 += dt * domega_dt

        theta_values.append(theta0)
        omega_values.append(omega0)

    return theta_values, omega_values

theta0 = np.pi / 4  
omega0 = 0.0
dt = 0.01
steps = 1000

time = np.linspace(0, 5, 1001)

theta_values, omega_values = euler_method(theta0, omega0, dt, steps)

#Critial Damping
def L(theta1):
    return L_0_2 - L_0_2 * np.abs(theta1) / np.pi

def Standing_Pumping(theta1, omega1):
    dtheta_dt1 = omega1
    dL_dt = (-L_0_2 * theta1 * omega1) / (np.pi * np.abs(theta1))
    domega_dt1 = (- c_c / (m * L_0_2**2))*dtheta_dt1 - (g / L_0_2 + k / (m * L_0_2**2))*theta1
    return dtheta_dt1, domega_dt1

def euler_method(theta_1, omega_1, dt, steps):
    theta_values1 = [theta_1]
    omega_values1 = [omega_1]

    for _ in range(steps):
        dtheta_dt1, domega_dt1 = Standing_Pumping(theta_1, omega_1)
        theta_1 += dt * dtheta_dt1
        omega_1 += dt * domega_dt1

        theta_values1.append(theta_1)
        omega_values1.append(omega_1)

    return theta_values1, omega_values1

theta_1 = np.pi / 4  
omega_1 = 0.0
dt = 0.01
steps = 1000

time = np.linspace(0, 5, 1001)

theta_values1, omega_values1 = euler_method(theta_1, omega_1, dt, steps)

#Overdamping
def L(theta2):
    return L_0_3 - L_0_3 * np.abs(theta2) / np.pi

def Standing_Pumping(theta2, omega2):
    dtheta_dt2 = omega2
    dL_dt2 = (-L_0_3 * theta2 * omega2) / (np.pi * np.abs(theta2))
    domega_dt2 = (- c_o / (m * L_0_3**2))*dtheta_dt2 - (g / L_0_3 + k / (m * L_0_3**2))*theta2
    return dtheta_dt2, domega_dt2

def euler_method(theta_2, omega_2, dt, steps):
    theta_values2 = [theta_2]
    omega_values2 = [omega_2]

    for _ in range(steps):
        dtheta_dt2, domega_dt2 = Standing_Pumping(theta_2, omega_2)
        theta_2 += dt * dtheta_dt2
        omega_2 += dt * domega_dt2

        theta_values2.append(theta_2)
        omega_values2.append(omega_2)

    return theta_values2, omega_values2

theta_2 = np.pi / 4  
omega_2 = 0.0
dt = 0.01
steps = 1000

time = np.linspace(0, 5, 1001)

theta_values2, omega_values2 = euler_method(theta_2, omega_2, dt, steps)

# Set up custom color palette
colors = plt.cm.viridis(np.linspace(0, 1, 2))


# Plot phase trajectory (Omega vs Theta)
plt.figure(figsize=(20, 10))
plt.plot(theta_values, omega_values, label='$c_u$ = 400', color='blue', linewidth=2, linestyle=':', marker='s', markersize=1)
plt.plot(theta_values1, omega_values1, label='$c_c$ = 911.47', color='green', linewidth=2, linestyle='-.', marker='s', markersize=1)
plt.plot(theta_values2, omega_values2, label='$c_o$ = 1500', color='red', linewidth=2, linestyle='--', marker='s', markersize=1)
plt.xlabel('$\phi$', fontsize=40, fontweight='bold', color='green')
plt.ylabel('$\dot{\phi}$', fontsize=40, fontweight='bold', color='green')
plt.grid(True, linestyle='--', color='gray', alpha=0.5)
plt.legend(fontsize=20, loc='upper right')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.gca().set_facecolor('#f0f0f0')  # Background color
plt.show()