import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def calc_P(kp, T):
    return (kp/2, 0, 0)
def calc_PI(kp, T):
    return(0.45*kp, 0.54*kp/T, 0)
def calc_PD(kp, T):
    return(0.8*kp, 0, 0.1*kp*T)
def calc_PID(kp, T):
    return(0.6*kp, 1.2*kp/T, 0.075*kp*T)

def zad1():
    # Niekończące sie oscylacje dla kp = 75.5
    # Okres 100/62 = 1.61
    # Ziegler-Nichols
    g_kp = 75.5
    g_T = 1.61
    P = calc_P(g_kp, g_T)
    PI = calc_PI(g_kp, g_T)
    PD = calc_PD(g_kp, g_T)
    PID = calc_PID(g_kp, g_T)
    print(f"P: {P}\nPI: {PI}\nPD: {PD}\nPID: {PID}\n")

def plot_PID(nastawy):
    R1 = 2
    R2 = 5
    C1 = 0.5
    L1 = 2
    L2 = 0.5
    yd = 1
    (kp, ki, kd) = nastawy

    def model(xe, t): # e_dot = -y_dot
        x1, x2, x3, e_int = xe
        # PID
        e_dot = -(x2*(-R2/L2) + x3 * (1/L2))
        e = yd - x2
        u = kp * e + ki * e_int + kd * e_dot

        x1_dot = x1*(-R1/L1) + x3 * (-1/L1) + u * (1/L1)
        x2_dot = x2*(-R2/L2) + x3 * (1/L2)
        x3_dot = x1*(1/C1) + x2 * (-1/C1)

        return [x1_dot, x2_dot, x3_dot, e]
    t = np.linspace(0, 100, 50000)
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
def zad2():
    g_kp = 75.5
    g_T = 1.61
    plot_PID(calc_P(g_kp, g_T))
    plot_PID(calc_PI(g_kp, g_T))
    plot_PID(calc_PD(g_kp, g_T))
    plot_PID(calc_PID(g_kp, g_T))

zad2()