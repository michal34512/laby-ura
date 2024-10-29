from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def zad1():
    R1 = 2
    R2 = 5
    C1 = 0.5
    L1 = 2
    L2 = 0.5
    A = np.array([[-R1/L1, 0, -1/L1],[0, -R2/L2, 1/L2],[1/C1, -1/C1, 0]])
    B = np.array([[1/L1],[0],[0]])
    C = np.array([[0, 1, 0]])
    D = 0
    u = 1
    def model(x, t):
        x = x.reshape(3,1)
        dxdt = A @ x + B * u
        return dxdt.reshape(3)
    t = np.linspace(0, 5, 100)
    x0 = np.array([0.0, 0.0, 0.0])
    x = odeint(model, x0, t)
    y = (C @ x.T).flatten() + D * u
    plt.figure(figsize=(10, 5))
    plt.plot(t, x[:, 0], label='x1(t)')
    plt.plot(t, x[:, 1], label='x2(t)')
    plt.plot(t, y, label='y(t)', linestyle='--')
    plt.xlabel('Czas t')
    plt.ylabel('Wartości')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.show()
def zad2():
    R1 = 2
    R2 = 5
    C1 = 0.5
    L1 = 2
    L2 = 0.5
    yd = 3
    kp = 1
    ki = 4
    kd = 0.1

    def model(xe, t): # e_dot = -y_dot
        x1, x2, x3, e_int = xe
        # PID
        e_dot = -(x2*(-R2/L2) + x3 * (-1/L2))
        e = yd - x2
        u = kp * e + ki * e_int + kd * e_dot

        x1_dot = x1*(-R1/L1) + x3 * (-1/L1) + u * (1/L1)
        x2_dot = x2*(-R2/L2) + x3 * (1/L2)
        x3_dot = x1*(1/C1) + x2 * (-1/C1)

        return [x1_dot, x2_dot, x3_dot, e]
    t = np.linspace(0, 25, 500)
    x0 = np.array([0.0, 0.0, 0.0, 0.0])
    sol = odeint(model, x0, t)
    x1, x2, x3, e_int = sol.T
    y = x2
    plt.figure(figsize=(10, 5))
    plt.plot(t, x1, label='x1(t)')
    plt.plot(t, x2, label='x2(t)')
    plt.plot(t, y, label='y(t)', linestyle='--')
    plt.xlabel('Czas t')
    plt.ylabel('Wartości')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.show()
zad2()