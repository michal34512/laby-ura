from matplotlib.pylab import det
import numpy as np
import matplotlib . pyplot as plt
from scipy import signal
from scipy.integrate import odeint
from sympy import Matrix, simplify, symbols

def characteristic_equation(A, B, K):
    s = symbols('s')
    A_closed = A - B @ K
    char_matrix = s * Matrix.eye(A.shape[0]) - A_closed
    characteristic_eq = char_matrix.det()
    return characteristic_eq

# Dla sterowalnej postaci układu z Rys. 2 wyznaczyć równanie charakterystyczne układu
# zamkniętego ze sprzężeniem od stanu (11). W tym celu można skorzystać z relacji (13).
def zad31():
    C1 = 1
    R1 = 1
    C2 = 2
    R2 = 1
    C3 = 3
    R3 = 1

    ##### State-equasions #####
    k1, k2, k3 = symbols("k1 k2 k3")

    A = np.array([[-1/(R1*C1), 0, 0],[0, -1/(R2*C2), 0],[0, 0, -1/(R3*C3)]])
    B = np.array([[1/(R1*C1)],[1/(R2*C2)],[1/(R3*C3)]])
    K = np.array([[k1, k2, k3]])
    print(simplify(characteristic_equation(A, B, K)))
    # k1*s**2 + 0.83333*k1*s + 0.16666*k1 + 0.5*k2*s**2 + 0.6666*k2*s + 0.1666*k2 + 0.3333*k3*s**2 + 0.5*k3*s + 0.16666*k3 + s**3 +1.83333*s**2 + s + 0.1666

# Wyznaczyć takie wartości macierzy K, aby bieguny układu zamkniętego wynosiły odpowiednio
# s1 = −1, s1 = −2, s3 = −5. Uwaga: idealne równanie charakterystyczne dla przyjętych 
# biegunów wynosi s^3 + 8s^2 + 17s + 10
def zad32():
    desired_poles = [-1, -2, -5]
    C1 = 1
    R1 = 1
    C2 = 2
    R2 = 1
    C3 = 3
    R3 = 1
    
    A = np.array([[-1/(R1*C1), 0, 0],[0, -1/(R2*C2), 0],[0, 0, -1/(R3*C3)]])
    B = np.array([[1/(R1*C1)],[1/(R2*C2)],[1/(R3*C3)]])

    K = signal.place_poles(A, B, desired_poles).gain_matrix

    print(f'K vector: {K}')

# Przeprowadzić symulację odpowiedzi obiektu na wymuszenie (11).
# • Jaki jest charakter odpowiedzi układu zamkniętego?
# • Jak poszczególne wartości wzmocnień K wpływają na pozycje biegunów układu zamkniętego?
def zad33():
    desired_poles = [-1, -2, -5]
    C1 = 1
    R1 = 1
    C2 = 2
    R2 = 1
    C3 = 3
    R3 = 1
    
    A = np.array([[-1/(R1*C1), 0, 0],[0, -1/(R2*C2), 0],[0, 0, -1/(R3*C3)]])
    B = np.array([[1/(R1*C1)],[1/(R2*C2)],[1/(R3*C3)]])
    C = np.array([1, 2, 3]) # nie było podane jaki jest sygnał wyjściowy
    D = 0
    
    K = np.array(signal.place_poles(A, B, desired_poles).gain_matrix)
    #K = np.array([[10, 0, 0]])
    
    current_A = (A - B @ K)
    print(f"Bieguny nowej macierzy A: {np.linalg.eigvals(current_A)}") # Bieguny teraz są dokładnie w tym miejscu co chcieliśmy żeby były
    current_C = (C - D * K)

    SS_before = signal.StateSpace(A, B, C, D)
    SS_after = signal.StateSpace(current_A, B, current_C, 0)

    t,y = signal.step(SS_before)
    plt.plot(t,y)
    plt.show()
    t,y = signal.step(SS_after)
    plt.plot(t,y)
    plt.show()
    # • Odpowiedź układu się stabilizuje na 0.1. Dynamika jest większa, ale pojawia się przeregulowanie (coś za coś).
    # • Wpływają zgodnie z równaniem wyzmaczonym w zadaniu 3.1 (miejsca zerowe wielomianu to miejsca biegunów)

if __name__ == "__main__":
    zad31()    
   
