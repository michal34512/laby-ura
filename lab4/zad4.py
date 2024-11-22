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

def plot_PID(nastawy):
    R1 = 2
    R2 = 5
    C1 = 0.5
    L1 = 2
    L2 = 0.5
    yd = 1
    (kp, ki, kd) = nastawy
    kp = 1000
    ki = 34.7
    kd = 194.4

    def model(xe, t): # e_dot = -y_dot
        x1, x2, x3, e_int, I_ISE, I_ITSE, I_IAE, I_ITAE  = xe
        # PID
        e_dot = -(x2*(-R2/L2) + x3 * (1/L2))
        e = yd - x2
        u = kp * e + ki * e_int + kd * e_dot
        ISE = e**2
        ITSE = t*e**2
        IAE = abs(e)
        ITAE = t*abs(e)


        x1_dot = x1*(-R1/L1) + x3 * (-1/L1) + u * (1/L1)
        x2_dot = x2*(-R2/L2) + x3 * (1/L2)
        x3_dot = x1*(1/C1) + x2 * (-1/C1)

        return [x1_dot, x2_dot, x3_dot, e, ISE, ITSE, IAE, ITAE]
    t = np.linspace(0, 10, 50000)
    x0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    sol = odeint(model, x0, t)
    x1, x2, x3, e_int, i_ise, i_itse, i_iae, i_itae = sol.T
    print(f"{i_ise[-1]}\n{i_itse[-1]}\n{i_iae[-1]}\n{i_itae[-1]}\n")
    y = x2
    plt.figure(figsize=(10, 5))
    # plt.plot(t, i_ise, label='ise(t)')
    # plt.plot(t, x1, label='x1(t)')
    # plt.plot(t, x2, label='x2(t)')
    plt.plot(t, y, label='y(t)', linestyle='--')
    plt.xlabel('Czas t')
    plt.ylabel('Wartości')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.show()

def simulate_PID(kp, ki, kd):
    # Parametry układu
    R1 = 2
    R2 = 5
    C1 = 0.5
    L1 = 2
    L2 = 0.5
    yd = 1  # wartość zadana

    def model(xe, t): # e_dot = -y_dot
        x1, x2, x3, e_int, I_ISE, I_ITSE, I_IAE, I_ITAE  = xe
        # PID
        e_dot = -(x2*(-R2/L2) + x3 * (1/L2))
        e = yd - x2
        u = kp * e + ki * e_int + kd * e_dot
        ISE = e**2
        ITSE = t*e**2
        IAE = abs(e)
        ITAE = t*abs(e)


        x1_dot = x1*(-R1/L1) + x3 * (-1/L1) + u * (1/L1)
        x2_dot = x2*(-R2/L2) + x3 * (1/L2)
        x3_dot = x1*(1/C1) + x2 * (-1/C1)

        return [x1_dot, x2_dot, x3_dot, e, ISE, ITSE, IAE, ITAE]

    t = np.linspace(0, 100, 50000)
    x0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    sol = odeint(model, x0, t)
    i_ise = sol[:, 7]  # I_ISE

    return i_ise[-1]

def optimize_PID():
    kp_values = np.linspace(10000, 20000, 10)  # Możliwe wartości kp
    ki_values = np.linspace(100, 500, 10)  # Możliwe wartości ki
    kd_values = np.linspace(100, 800, 10)  # Możliwe wartości kd
    
    best_kp, best_ki, best_kd = 0, 0, 0
    best_ISE = float('inf')

    for kp in kp_values:
        for ki in ki_values:
            for kd in kd_values:
                ISE = simulate_PID(kp, ki, kd)
                if ISE < best_ISE:
                    best_ISE = ISE
                    best_kp, best_ki, best_kd = kp, ki, kd
                    print(f"Nowe najlepsze nastawy: kp={kp}, ki={ki}, kd={kd} z ISE={ISE}")

    print("\nNajlepsze znalezione nastawy PID:")
    print(f"kp = {best_kp}, ki = {best_ki}, kd = {best_kd}")
    print(f"Minimalny wskaźnik ISE = {best_ISE}")

    return best_kp, best_ki, best_kd, best_ISE

def zad1():
    g_kp = 75.5
    g_T = 1.61
    
    # ISE = 0.63
    # ITSE = 0.49
    # IAE = 1.3
    # ITAE = 2.11
    plot_PID(calc_PID(g_kp, g_T))
optimize_PID()