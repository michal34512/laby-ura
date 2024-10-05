import numpy as np
import matplotlib . pyplot as plt
from scipy import signal
from scipy.integrate import odeint

L = 12
R = 1
Ca = 100e-6

def zad31():
    ##### Transfer function #####
    G = signal.TransferFunction([1, 0],[L, R, 1/Ca])
    t, y = signal.step(G)
    t2, y2 = signal.impulse(G)
    plt.plot(t, y, label = "TF step")
    plt.plot(t2, y2, label = "TF impulse")

def zad32():
    ##### State-equasions #####
    A = np.array([[0, 1],[-1/(L*Ca), -(R/L)]])
    B = np.array([[0],[1/L]])
    C = np.array([[0, 1]])
    D = 0
    SS = signal.StateSpace(A,B,C,D)
    t, y = signal.step(SS)
    t2, y2 = signal.impulse(SS)
    plt.plot(t, y, label = "SS step")
    plt.plot(t2, y2, label = "SS impulse")
    # Wykresy sie pokrywają, oscylacyjny charakter odpowiedzi wynika z oscylacyjnego 
    # charakteru symulowanego układu (RLC). Już po samej transmitancji widać że ona również 
    # ma postać obiektu oscylacyjnego

    # Odpowiedzi impulsowe przyjmują niezerowe wartości początkowe co wynika z twierdzenia o
    # wartości początkowej. Jeśli pomnożymy transmitancję razy 's' i obliczymy granicę przy 
    # 's' dążącym do nieskończoności to uzyskamw wartość poczatkową 1/12 

def zad33():
    ##### SS -> TF #####
    A = np.array([[0, 1],[-1/(L*Ca), -(R/L)]])
    B = np.array([[0],[1/L]])
    C = np.array([[0, 1]])
    D = 0
    G = signal.ss2tf(A,B,C,D)
    print(G)

    ##### TF -> SS #####
    SS = signal.tf2ss([1, 0],[L, R, 1/Ca])
    print(SS)

    # Zamiana SS -> TF daje taką samą transmitancję (może być jedynie mianownik i licznik pomnożony
    # razy jakąś stałą). Zamiana TF -> SS daje inne równania stanu, bo równania stanu nie są jednoznaczne. 
    # W zaleności co przyjmiemy za zmienne stanu to równania stanu mogą wyglądać inaczej. Najlepiej jest za 
    # zmienne stanu przyjąć wielkości które mają odzwierciedlenie w rzeczywistości (np 'q' i 'I')

    # Da się to zrobić, trzebaby zrobić transpozycję macierzy A, zamienić miejscami wartości B i D (odpowiada to zamianie miejscami zmiennych stanu), 
    # 
if __name__ == "__main__":
    zad31()
    zad32()
    plt.legend()
    plt.show()
    zad33()
    