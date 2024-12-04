import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def zad1():
    def model(y, t):
        dydt = t**2
        return dydt

    t = np.linspace(0, 10, 100)
    y0 = 0
    y = odeint(model, y0, t)

    plt.plot(t, y, label='Rozwiązanie y(t)')
    plt.xlabel('Czas t [s]')
    plt.ylabel('y(t)')
    plt.title(r'Rozwiązanie równania różniczkowego $\dot{y}(t) = t^2$')
    plt.legend()
    plt.grid(True)
    plt.show()

def zad2():
    def model(y, t):
        dydt = t**2
        return dydt
    t = np.linspace(0, 10, 100)
    ya = 1/3 * t**3
    yn0 = 0
    yn = odeint(model, yn0, t)

    plt.plot(t, ya, label='Rozwiązanie analityczne y(t)')
    plt.plot(t, yn, label='Rozwiązanie numeryczne y(t)')
    plt.xlabel('Czas t [s]')
    plt.ylabel('y(t)')
    plt.title(r'Rozwiązanie równania różniczkowego $\dot{y}(t) = t^2$')
    plt.legend()
    plt.grid(True)
    plt.show()

# zad1()
zad2()