import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def zad1():
    L = 0.65
    m = 0.01187
    km = 1.16e-4
    R = 27.7
    g = 9.81

    
    h0 = 0.1
    x0 = [h0, 0 , np.sqrt(g*m/km)*h0]
    u = x0[2] * (R - (2*x0[1]*km)/(x0[0]**2))
    
    def model(x, t):
        x1, x2, x3 = x
        x1dot = x2
        x2dot = g - (km/m) * (x3**2)/(x1**2)
        x3dot = (2*km/L) * (x2*x3)/(x1**2) - R*x3/L + u/L
        return [x1dot, x2dot, x3dot]
    t = np.linspace(0, 2, 1000)
    x = odeint(model, x0, t)
    plt.figure(figsize=(10, 5))
    plt.plot(t, x[:, 0], label='x1(t)')
    plt.xlabel('Czas t')
    plt.ylabel('theta(t)')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.grid(True)
    plt.show()

def zad2():
    L = 0.65
    m = 0.01187
    km = 1.16e-4
    R = 27.7
    g = 9.81

    
    h0 = 0.1
    x0 = [h0, 0 , h0*np.sqrt(g*m/km)]
    u0 = x0[2] * (R - (2*x0[1]*km)/(x0[0]**2))

    A0 = np.array([[0, 1, 0], 
          [(2*km*x0[2]**2)/(m*x0[0]**3), 0, (-2*km*x0[2])/(m*x0[0]**2)], 
          [(-4*km*x0[1]*x0[2])/(L*x0[0]**3), (2*km*x0[2])/(L*x0[0]**2), (2*km*x0[1])/(L*x0[0]**2)-R/L]])
    B0 = np.array([[0],[0],[1/L]])

    ufala = 0

    def model(x, t):
        x1, x2, x3 = x
        x1dot = x2
        x2dot = g - (km/m) * (x3**2)/(x1**2)
        x3dot = (2*km/L) * (x2*x3)/(x1**2) - R*x3/L + u0/L
        return [x1dot, x2dot, x3dot]
    
    def lin_model(x, t):
        xdot = A0 @ x + B0.flatten() * ufala
        return xdot
    t = np.linspace(0, 3, 1000)
    xfala = odeint(lin_model, [0, 0, 0], t)
    x_lin = xfala + x0
    x = odeint(model, x0, t)
    plt.figure(figsize=(10, 5))
    plt.plot(t, x[:, 0], label='x1(t)')
    plt.plot(t, x_lin[:, 0], label='x1_lin(t)')
    plt.xlabel('Czas t')
    plt.ylabel('theta(t)')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.grid(True)
    plt.show()

zad2()