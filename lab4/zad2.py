from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def czas_ustalania(y, t, y_d, tolerance=0.02):
    min_value = y_d * (1 - tolerance)
    max_value = y_d * (1 + tolerance)
    
    for i in range(len(y)):
        if np.all(y[i:] >= min_value) and np.all(y[i:] <= max_value):
            return t[i]
    return np.inf  

def model(xe, t, kp, ki, kd, yd, R1, R2, C1, L1, L2):
    x1, x2, x3, e_int = xe
    e_dot = -(x2 * (-R2 / L2) + x3 * (1 / L2))
    e = yd - x2
    u = kp * e + ki * e_int + kd * e_dot
    x1_dot = x1 * (-R1 / L1) + x3 * (-1 / L1) + u * (1 / L1)
    x2_dot = x2 * (-R2 / L2) + x3 * (1 / L2)
    x3_dot = x1 * (1 / C1) + x2 * (-1 / C1)
    return [x1_dot, x2_dot, x3_dot, e]

def grid_search_pid(yd=3):
    R1, R2, C1, L1, L2 = 2, 5, 0.5, 2, 0.5
    t = np.linspace(0, 20, 5000)
    x0 = np.array([0.0, 0.0, 0.0, 0.0])
    
    kp_values = np.linspace(120, 130, 20)  # np.linspace(start, stop, liczba_prób)
    ki_values = np.linspace(22, 27, 10)
    kd_values = np.linspace(35, 40, 10)

    najlepszy_czas = np.inf
    najlepsze_parametry = (0, 0, 0)

    for kp in kp_values:
        for ki in ki_values:
            for kd in kd_values:
                sol = odeint(model, x0, t, args=(kp, ki, kd, yd, R1, R2, C1, L1, L2))
                y = sol[:, 1]
                ts = czas_ustalania(y, t, yd)
                
                if ts < najlepszy_czas:
                    najlepszy_czas = ts
                    najlepsze_parametry = (kp, ki, kd)
                    print(f"Znaleziono lepsze: kp={kp}, ki={ki}, kd={kd}, czas ustalania={ts:.2f}s")

    return najlepsze_parametry, najlepszy_czas
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
    yd = 1
    kp = 45.3
    ki = 56.273291925465834
    kd = 9.116625

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

    ts = czas_ustalania(y, t, yd)
    print(f"Czas ustalania: {ts:.2f} sekund" if ts is not None else "Czas ustalania nie został osiągnięty.")
    

    plt.figure(figsize=(10, 5))
    # plt.plot(t, x1, label='x1(t)')
    # plt.plot(t, x2, label='x2(t)')
    # plt.plot(t, x3, label='x3(t)')
    plt.plot(t, y, label='y(t)', linestyle='--')
    plt.xlabel('Czas t')
    plt.ylabel('Wartości')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.show()

# Pkt. 3
    # - pochodna = układ bardziej dynamiczny
    #   całka = sprowadzanie uchybu do zera / przeregulowania
    #   proporcja = początkowy "kop"
    # - zerowy uchyb nie jest możliwy w przypadku regulatora typu P (wzmocnienei kp musiało by być nieskończone)
    # - 
zad2()
#grid_search_pid()