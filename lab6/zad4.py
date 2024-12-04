import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def zad1():
    kp = 2
    T = 2
    kob = 4
    xd = 1
    def feedback(y, t):
        e = xd - y[0]
        u = e * kp
        uc = np.clip(u, -0.1, 0.1)
        ydot = 1/T * (uc * kob - y[0])
        return [ydot]
    y0 = [0]
    t = np.linspace(0, 10, 100)
    y = odeint(feedback, y0, t)

    plt.plot(t, y, label=r'y(t)')
    plt.xlabel('Czas t [s]')
    plt.ylabel('y(t)')
    plt.title(r'Rozwiązanie równania różniczkowego')
    plt.legend()
    plt.grid(True)
    plt.show()
def zad2():
    kp = 2
    T = 2
    kob = 4
    xd_values = [1, 2, 3]

    t = np.linspace(0, 10, 100)
    plt.figure(figsize=(10, 6))
    
    for xd in xd_values:
        def feedback(y, t):
            e = xd - y[0]
            u = e * kp
            uc = np.clip(u, -0.1, 0.1)
            ydot = 1 / T * (uc * kob - y[0])
            return [ydot]
        y0 = [0]
        y = odeint(feedback, y0, t)
        plt.plot(t, y, label=f'$x_d = {xd}$')
    
    plt.xlabel('Czas t [s]')
    plt.ylabel('y(t)')
    plt.title('Rozwiązanie równania różniczkowego dla różnych $x_d$')
    plt.legend()
    plt.grid(True)
    plt.show()
    # zasada superpozycji  iskalowania nie jesat zachowana - mamy kilkukrotnie wiekszy 
    # sygnal zadany a nic sie nie zmienia na wyjscu
    # 
    # uklad nie ma charakteru liniowego, bo nie jest zachowana zasada superpozycji 
def zad3():
    kp = 2
    T = 2
    kob = 4
    xd_values = [1, 2, 3]

    t = np.linspace(0, 10, 100)
    plt.figure(figsize=(10, 6))
    
    for xd in xd_values:
        def feedback(y, t):
            e = xd - y[0]
            u = e * kp
            ydot = 1 / T * (u * kob - y[0])
            return [ydot]
        y0 = [0]
        y = odeint(feedback, y0, t)
        plt.plot(t, y, label=f'$x_d = {xd}$')
    
    plt.xlabel('Czas t [s]')
    plt.ylabel('y(t)')
    plt.title('Rozwiązanie równania różniczkowego dla różnych $x_d$')
    plt.legend()
    plt.grid(True)
    plt.show()
    # uklad jest liniowy

zad3()