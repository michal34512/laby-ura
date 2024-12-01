import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def calc_P(ku, T):
    return (ku/2, 0, 0)
def calc_PI(ku, T):
    return(0.45*ku, 0.54*ku/T, 0)
def calc_PD(ku, T):
    return(0.8*ku, 0, 0.1*ku*T)
def calc_PID(ku, T):
    return(0.6*ku, 1.2*ku/T, 0.075*ku*T)

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
    yd = 3
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
    t = np.linspace(0, 20, 50000)
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
    e = y - yd
    plt.figure(figsize=(10, 5))
    plt.plot(t, e, label='e(t)', linestyle='--')
    plt.xlabel('Czas t')
    plt.ylabel('Uchyb')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.show()

def zad2():
    g_kp = 75.5
    g_T = 1.61

    print(calc_P(g_kp, g_T)) # (37.75, 0, 0)
    print(calc_PI(g_kp, g_T)) # (33.975, 25.32298136645963, 0)
    print(calc_PD(g_kp, g_T)) # (60.400000000000006, 0, 12.155500000000002)
    print(calc_PID(g_kp, g_T)) # (45.3, 56.273291925465834, 9.116625)
    
    plot_PID(calc_P(g_kp, g_T))
    plot_PID(calc_PI(g_kp, g_T))
    plot_PID(calc_PD(g_kp, g_T))
    plot_PID(calc_PID(g_kp, g_T))

#zad2()
plot_PID((1000,66.8,158))