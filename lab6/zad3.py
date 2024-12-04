import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def zad1():
    kp = 2
    omega = 4
    ksi = 0.25
    u = 1
    def model(x, t):
        dydt = x[1]
        dy2dt2 = -(2*ksi/omega) * x[1] - np.sqrt(x[0])/omega + (kp/omega**2) * u
        return [dydt, dy2dt2]
    t = np.linspace(0, 50, 500)
    yn0 = [0, 0]
    yn = odeint(model, yn0, t)

    plt.plot(t, yn[:, 0], label=r'y(t)')
    plt.plot(t, yn[:, 1], label=r'$\dot{y}(t)$')
    plt.xlabel('Czas t [s]')
    plt.ylabel('y(t)')
    plt.title(r'Rozwiązanie równania różniczkowego')
    plt.legend()
    plt.grid(True)
    plt.show()
zad1()