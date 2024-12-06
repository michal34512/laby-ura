import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def zad1():
    u = 1
    def model(y, t):
        z1 = y[0]
        z2 = y[1]
        z1dot = z1 * np.log(z2)
        z2dot = -z2 * np.log(z1) + z2 * u
        return [z1dot, z2dot]

    t = np.linspace(0, 10, 100)
    y0 = [1 ,1]
    y = odeint(model, y0, t)

    plt.plot(t, y[:, 0], label='x1(t)')
    plt.plot(t, y[:, 1], label='x2(t)')
    plt.xlabel('Czas t [s]')
    plt.ylabel('y(t)')
    plt.legend()
    plt.grid(True)
    plt.show()

def zad3():
    u = 1
    def model(y, t):
        x1 = y[0]
        x2 = y[1]
        x1dot = x2
        x2dot = -x1 + u
        return [x1dot, x2dot]

    t = np.linspace(0, 10, 100)
    y0 = [np.log(1), np.log(1)]
    y = odeint(model, y0, t)

    plt.plot(t, y[:, 0], label='x1(t)')
    plt.plot(t, y[:, 1], label='x2(t)')
    plt.xlabel('Czas t [s]')
    plt.ylabel('y(t)')
    plt.legend()
    plt.grid(True)
    plt.show()

    z1 = np.exp(y[:, 0])
    z2 = np.exp(y[:, 1])

    plt.plot(t, z1, label='x1(t)')
    plt.plot(t, z2, label='x2(t)')
    plt.xlabel('Czas t [s]')
    plt.ylabel('y(t)')
    plt.legend()
    plt.grid(True)
    plt.show()
zad2()