import numpy as np
import matplotlib . pyplot as plt
from scipy import signal
from scipy.integrate import odeint

def zad21():
    kp = 2
    T = 2

    ##### Transfer function #####
    G = signal.TransferFunction([kp],[T, 1])
    t, y = signal.step(G)
    plt.plot(t,y, label="Transfer function", color='r')
    # Tak, odpowiedz skokowa odpowiada zalozeniom. Ma charakter odpowiedzi układu pierwszego rzędu

def zad23():
    kp = 2
    T = 2

    ##### State-equasions #####
    A = np.array([-1/T])
    B = np.array([1])
    C = np.array([[kp/T]])
    D = 0
    SS = signal.StateSpace(A,B,C,D)
    t2, y2 = signal.step(SS)
    plt.plot(t2,y2, label="State Space", color='b')
    # Spostrzeżenie - zamiana B z C nic nie zmienia

def zad24():
    kp = 2
    T = 2

    t = np.linspace(0, 15, 100)
    def model(y, t):
        dydt = (-1/T) * y + (kp/T) * (1 if t >= 0 else 0)
        return dydt
    y0 = 0
    y = odeint(model, y0, t) 
    plt.plot(t ,y, label="odeint", color='g')
    

if __name__ == "__main__":
    zad21()
    zad23()
    zad24()
    plt.legend()
    plt.show()
    # odpowiedzi skokowe się pokrywają (odeint może delikatnie odstawać w zależności od rozdzielczości czasu)
    # zastosowane przekształcenia nie zmianiają charakteru wyjścia systemu
