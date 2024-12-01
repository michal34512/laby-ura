import numpy as np
from scipy.linalg import solve_continuous_are
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def get_K_mat(_Q, _R):
    L = 0.2
    R = 0.5
    C = 0.5
    A = np.array([[0, 1],
                [-1/(L*C), -R/L]])
    B = np.array([[0],
                [1/L]])
    P = solve_continuous_are(A, B, _Q, _R)
    K = np.linalg.inv(_R) @ B.T @ P
    return K

def zad1():
    L = 0.2
    R = 0.5
    C = 0.5
    A = np.array([[0, 1],
                [-1/(L*C), -R/L]])
    B = np.array([[0],
                [1/L]])

    Q = np.eye(2)
    R = np.array([[1]])

    P = solve_continuous_are(A, B, Q, R)
    # Funkcja tworzy macierz Hamiltona i używa dekompozycji Shura do znalezienia rozwiązań

    K = np.linalg.inv(R) @ B.T @ P

    # Wyniki
    print("Macierz P (rozwiązanie równania Riccatiego):\n", P)
    print("Wzmocnienia K:\n", K)

def zad2():
    L = 0.2
    R = 0.5
    C = 0.5
    A = np.array([[0, 1],
                [-1/(L*C), -R/L]])
    B = np.array([[0],
                [1/L]])


    u = 1

    def model(x, t):
        dxdt = A @ x + B.flatten() * u
        return dxdt
    x0 = [0, 0]
    t = np.linspace(0, 5, 500)

    x = odeint(model, x0, t)

    plt.figure(figsize=(10, 6))
    plt.plot(t, x[:, 0], label="x1(t)")
    plt.plot(t, x[:, 1], label="x2(t)")
    plt.title("Odpowiedź skokowa układu otwartego")
    plt.xlabel("Czas [s]")
    plt.ylabel("Stan")
    plt.legend()
    plt.grid()
    plt.show()

def zad3(_Q, _R):
    L = 0.2
    R = 0.5
    C = 0.5
    A = np.array([[0, 1],
                [-1/(L*C), -R/L]])
    B = np.array([[0],
                [1/L]])
    K = get_K_mat(_Q, _R)
    
    saved_u = []
    saved_t = []
    def model(x, t):
        u = - K @ x
        dxdt = A @ x + B.flatten() * u
        saved_u.append(u)
        saved_t.append(t)
        return dxdt
    x0 = [1, 1]
    t = np.linspace(0, 5, 500)

    x = odeint(model, x0, t)

    plt.figure(figsize=(10, 6))
    plt.plot(t, x[:, 0], label="x1(t)")
    plt.plot(t, x[:, 1], label="x2(t)")
    #plt.plot(saved_t, saved_u, label="u(t)")
    plt.title("Odpowiedź skokowa układu zamkniętego")
    plt.xlabel("Czas [s]")
    plt.ylabel("Stan")
    plt.legend()
    plt.grid()
    plt.show()

def zad4(_Q, _R):
    L = 0.2
    R = 0.5
    C = 0.5
    A = np.array([[0, 1],
                [-1/(L*C), -R/L]])
    B = np.array([[0],
                [1/L]])

    K = get_K_mat(_Q, _R)

    def model(x, t):
        u = -K @ x[:2]
        dxdt = A @ x[:2] + B.flatten() * u
        dJdt = x[:2].T @ _Q @ x[:2] + u.T @ _R @ u 
        
        return np.concatenate([dxdt, [dJdt]])

    x0 = [1, 1, 0]
    t = np.linspace(0, 5, 500) 

    x = odeint(model, x0, t)

    J_value = x[-1, 2]

    print(f"Wartość wskaźnika jakości J: {J_value}")

    # Wykres wyników
    plt.figure(figsize=(10, 6))
    plt.plot(t, x[:, 0], label="x1 (stan 1)")
    plt.plot(t, x[:, 1], label="x2 (stan 2)")
    plt.title("Odpowiedź układu z wskaźnikiem jakości")
    plt.xlabel("Czas [s]")
    plt.ylabel("Stan")
    plt.legend()
    plt.grid()
    plt.show()

#zad1()
# zad2()
zad3([[1, 0], [0, 1]], [[1]])
# zad4([[1, 0], [0, 0.001]], [[1]])
