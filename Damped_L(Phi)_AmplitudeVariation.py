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

def L(theta):
    return L_0_1 * (1 - (np.abs(theta) / np.pi))

def Standing_Pumping(theta, omega):
    dtheta_dt = omega
    dL_dt = (-L_0_1 * theta * omega) / (np.pi * np.abs(theta))
    domega_dt = (- c_c / (m * L_0_1**2))*dtheta_dt - (g / L_0_1 + k / (m * L_0_1**2))*theta
    return dtheta_dt, domega_dt

def euler_method(theta0, omega0, dt, steps):
    theta_values = [theta0]  # Include the initial value
    omega_values = [omega0]

    for _ in range(steps):  # Run for steps iterations
        dtheta_dt, domega_dt = Standing_Pumping(theta0, omega0)
        theta0 += dt * dtheta_dt
        omega0 += dt * domega_dt
        
        theta_values.append(theta0)
        omega_values.append(omega0)

    return theta_values, omega_values

def L2(theta2):
    return L_0_2 * (1 - (np.abs(theta2) / np.pi))

def Standing_Pumping2(theta2, omega2):
    dtheta_dt2 = omega2
    dL_dt2 = (-L_0_2 * theta2 * omega2) / (np.pi * np.abs(theta2))
    domega_dt2 = (- c_o / (m * L_0_1**2))*dtheta_dt2 - (g / L_0_1 + k / (m * L_0_1**2))*theta2
    return dtheta_dt2, domega_dt2

def euler_method2(theta02, omega02, dt, steps):
    theta_values2 = [theta02]  # Include the initial value
    omega_values2 = [omega02]

    for _ in range(steps):  # Run for steps iterations
        dtheta_dt2, domega_dt2 = Standing_Pumping2(theta02, omega02)
        theta02 += dt * dtheta_dt2
        omega02 += dt * domega_dt2

        theta_values2.append(theta02)
        omega_values2.append(omega02)

    return theta_values2, omega_values2

def L3(theta3):
    return L_0_3 * (1 - (np.abs(theta3) / np.pi))

def Standing_Pumping3(theta3, omega3):
    dtheta_dt3 = omega3
    dL_dt3 = (-L_0_3 * theta3 * omega3) / (np.pi * np.abs(theta3))
    domega_dt3 = (- c_u / (m * L_0_1**2))*dtheta_dt3 - (g / L_0_1 + k / (m * L_0_1**2))*theta3
    return dtheta_dt3, domega_dt3

def euler_method3(theta03, omega03, dt, steps):
    theta_values3 = [theta03]  # Include the initial value
    omega_values3 = [omega03]

    for _ in range(steps):  # Run for steps iterations
        dtheta_dt3, domega_dt3 = Standing_Pumping3(theta03, omega03)
        theta03 += dt * dtheta_dt3
        omega03 += dt * domega_dt3

        theta_values3.append(theta03)
        omega_values3.append(omega03)

    return theta_values3, omega_values3

theta02 = np.pi / 4  
omega02 = 0.0
theta03 = np.pi / 4  
omega03 = 0.0
dt = 0.01
steps = 1000  # Adjusted to match the length of time

theta0 = np.pi / 4  # Reusing the same initial conditions as in euler_method2
omega0 = 0.0  # Reusing the same initial conditions as in euler_method2

time = np.linspace(0, 5, 1001)  # Adjusted to have 4000 points

theta_values, omega_values = euler_method(theta0, omega0, dt, steps)
theta_values2, omega_values2 = euler_method2(theta02, omega02, dt, steps)
theta_values3, omega_values3 = euler_method3(theta03, omega03, dt, steps)


# Plot phase trajectory (Theta vs Time)
plt.figure(figsize=(20, 10))
plt.plot(time, theta_values3, label='$c_u$ = 400', color="blue", linewidth=2, linestyle=':', markersize=1)
plt.plot(time, theta_values, label='$c_c$ = 911.47', color="green", linewidth=2, linestyle='-.', markersize=1)
plt.plot(time, theta_values2, label='$c_o$ = 1500', color="red", linewidth=2, linestyle='--', markersize=1)
plt.xlabel('t', fontsize=40, fontweight='bold', color='green')
plt.ylabel('$\phi$', fontsize=40, fontweight='bold', color='green')
plt.grid(True, linestyle='--', color='gray', alpha=0.5)
plt.legend(fontsize=20, loc='upper right')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.gca().set_facecolor('#f0f0f0')  # Background color
plt.show()