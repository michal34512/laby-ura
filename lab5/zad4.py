import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.linalg import solve_continuous_are
from scipy.interpolate import interp1d

def zad1(Q, R, qd):
    L = 0.2
    _R = 0.5
    C = 0.5
    A = np.array([[0, 1],
                [-1/(L*C), -_R/L]])
    B = np.array([[0],
                [1/L]])
    P = solve_continuous_are(A, B,Q, R)
    K = np.linalg.inv(R) @ B.T @ P
    xd = [qd, 0]
    def model(x, t):
        u = (qd/C) - K @ (x - xd)
        dxdt = A @ x + B.flatten() * u
        return dxdt
    x0 = [1, 1]
    t = np.linspace(0, 5, 500)

    x = odeint(model, x0, t)

    plt.figure(figsize=(10, 6))
    plt.plot(t, x[:, 0], label="x1(t)")
    plt.plot(t, x[:, 1], label="x2(t)")
    plt.title("Odpowiedź skokowa układu zamkniętego")
    plt.xlabel("Czas [s]")
    plt.ylabel("Stan")
    plt.legend()
    plt.grid()
    plt.show()

def zad2(S, qd):
    L = 0.2
    R = 0.5
    C = 0.5
    A = np.array([[0, 1],
                [-1/(L*C), -R/L]])
    B = np.array([[0],
                [1/L]])

    Q = np.eye(1)  # Przykładowa macierz Q
    R = np.array([[1]])
    
    

    def riccati(p, t):
        P = np.reshape(p, (2, 2))
        P_dot = -(P @ A  - P @ B @ np.linalg.inv(R) @ B.T @ P + A.T @ P + Q)
        return P_dot.reshape(-1)

    P0_vec = S.reshape(-1)

    t_riccati = np.linspace(1, 0, 500)
    P_vec = odeint(riccati, P0_vec, t_riccati)

    P_matrices = np.array([np.reshape(p, (2, 2)) for p in P_vec])
    P11_interp = interp1d(t_riccati, P_matrices[:, 0, 0], fill_value="extrapolate")
    P12_interp = interp1d(t_riccati, P_matrices[:, 0, 1], fill_value="extrapolate")
    P22_interp = interp1d(t_riccati, P_matrices[:, 1, 1], fill_value="extrapolate")

    saved_u = []
    saved_t = []

    xd = [qd, 0]
    def model(x, t):
        P_t = np.array([[P11_interp(t), P12_interp(t)],
                        [P12_interp(t), P22_interp(t)]])
        
        K_t = np.linalg.inv(R) @ (B.T @ P_t)

        u = (qd/C) - K_t @ (x - xd)
        dxdt = A @ x + B.flatten() * u
        if t <= 2:
            saved_u.append(u)
            saved_t.append(t)
        return dxdt

    # Warunki początkowe i symulacja
    x0 = [1, 1]
    t_sim = np.linspace(0, 2, 500)
    x = odeint(model, x0, t_sim)
    
    # Wykres odpowiedzi układu
    plt.figure(figsize=(12, 6))
    plt.plot(t_sim, x[:, 0], label="x1(t)")
    plt.plot(t_sim, x[:, 1], label="x2(t)")
    plt.plot(saved_t, saved_u, label="u(t)")
    #plt.plot(np.linspace(0, 5, len(u_values)), u_values, label="u")

    plt.title("Odpowiedź układu z P(t) w czasie")
    plt.xlabel("Czas [s]")
    plt.ylabel("Stan")
    plt.legend()
    plt.grid()
    plt.show()

zad2(np.eye(2)*1000, 20)

