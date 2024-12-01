import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.linalg import solve_continuous_are
from scipy.interpolate import interp1d
# horyzont 1s, S na diagonali 100 -> jesli sie nei wyreguluje to  cos jest nie tak

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

    def riccati(p, t):
        P = np.reshape(p, (2, 2))
        P_dot = A.T @ P + P @ A - P @ B @ np.linalg.inv(R) @ B.T @ P + Q
        return P_dot.reshape(-1)

    P0 = np.eye(2)
    P0_vec = P0.reshape(-1)
    t = np.linspace(0, 5, 500)

    P_vec = odeint(riccati, P0_vec, t)
    P_matrices = [np.reshape(p, (2, 2)) for p in P_vec]

    print("Macierz P(t) na końcu symulacji:")
    print(P_matrices[-1])

def zad2(S):
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

    t = np.linspace(5, 0, 500)
    P_vec = odeint(riccati, P0_vec, t)
    P_matrices = [np.reshape(p, (2, 2)) for p in P_vec]

    P_matrices = np.array([np.reshape(p, (2, 2)) for p in P_vec])

    P11 = [P[0, 0] for P in P_matrices]
    P12 = [P[0, 1] for P in P_matrices]
    P21 = [P[1, 0] for P in P_matrices]
    P22 = [P[1, 1] for P in P_matrices]

    plt.figure(figsize=(12, 6))
    plt.plot(t, P11, label="P11 (t)")
    plt.plot(t, P12, label="P12 (t)")
    plt.plot(t, P21, label="P21 (t)")
    plt.plot(t, P22, label="P22 (t)")
    plt.title("Przebieg elementów macierzy P(t) w czasie")
    plt.xlabel("Czas [s]")
    plt.ylabel("Wartość elementów macierzy P(t)")
    plt.legend()
    plt.grid()
    plt.show()

    P_t0 = P_matrices[-1]
    print("Macierz P(t) na końcu symulacji (t = 0):")
    print(P_t0)


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


# nie mozna uzyc care -> uzyj odeint(riccatti rownanie 10). Nie wiemy jaki jest warunek poczatkow ale wiemy warunek koncowy (S), dlatego odwracamy czas
# P flip zeby bylo chronologiczne
# uzyc interpolacje jesli bedzie blad systemu
def zad3(S):
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
    def model(x, t):
        P_t = np.array([[P11_interp(t), P12_interp(t)],
                        [P12_interp(t), P22_interp(t)]])
        
        K_t = np.linalg.inv(R) @ (B.T @ P_t)
        
        u = -K_t @ x
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


def zad4(t_hor):
    L = 0.2
    R = 0.5
    C = 0.5
    A = np.array([[0, 1],
                [-1/(L*C), -R/L]])
    B = np.array([[0],
                [1/L]])

    Q = np.eye(2)  # Przykładowa macierz Q
    R = np.array([[1]])
    
    def riccati(p, t):
        P = np.reshape(p, (2, 2))
        P_dot = -(P @ A  - P @ B @ np.linalg.inv(R) @ B.T @ P + A.T @ P + Q)
        return P_dot.reshape(-1)

    P0_vec = 100*np.eye(2).reshape(-1)

    t_riccati = np.linspace(t_hor, 0, 1000)
    P_vec = odeint(riccati, P0_vec, t_riccati)

    P_matrices = np.array([np.reshape(p, (2, 2)) for p in P_vec])
    P11_interp = interp1d(t_riccati, P_matrices[:, 0, 0], fill_value="extrapolate")
    P12_interp = interp1d(t_riccati, P_matrices[:, 0, 1], fill_value="extrapolate")
    P22_interp = interp1d(t_riccati, P_matrices[:, 1, 1], fill_value="extrapolate")

    saved_u = []
    saved_t = []
    def model(x, t):
        P_t = np.array([[P11_interp(t), P12_interp(t)],
                        [P12_interp(t), P22_interp(t)]])
        
        K_t = np.linalg.inv(R) @ (B.T @ P_t)
        
        u = -K_t @ x
        dxdt = A @ x + B.flatten() * u
        if t <= t_hor + 1:
            saved_u.append(u)
            saved_t.append(t)
        return dxdt

    # Warunki początkowe i symulacja
    x0 = [1, 1]
    t_sim = np.linspace(0, t_hor+1, 1000)
    x = odeint(model, x0, t_sim)
    
    # Wykres odpowiedzi układu
    plt.figure(figsize=(12, 6))
    plt.plot(t_sim, x[:, 0], label="x1(t)", color="lightgreen", linestyle="--")
    plt.plot(t_sim, x[:, 1], label="x2(t)", color="lightgreen", linestyle="--")
    plt.plot(saved_t, saved_u, label="u(t)", color="forestgreen")
    #plt.plot(np.linspace(0, 5, len(u_values)), u_values, label="u")
    def get_K_mat(_Q, _R):
        P = solve_continuous_are(A, B, _Q, _R)
        K = np.linalg.inv(_R) @ B.T @ P
        return K
    K = get_K_mat(Q, R)
    
    saved_u = []
    saved_t = []
    def model(x, t):
        u = - K @ x
        dxdt = A @ x + B.flatten() * u
        if t <= t_hor + 1:
            saved_u.append(u)
            saved_t.append(t)
        return dxdt
    x0 = [1, 1]
    t = np.linspace(0, t_hor+1, 500)

    x = odeint(model, x0, t)

    plt.plot(t, x[:, 0], label="x1(t)", color="lightblue", linestyle="--")
    plt.plot(t, x[:, 1], label="x2(t)", color="lightblue", linestyle="--")
    plt.plot(saved_t, saved_u, label="u(t)", color="dodgerblue")
    
    plt.title("Odpowiedź układu z P(t) w czasie")
    plt.xlabel("Czas [s]")
    plt.ylabel("Stan")
    plt.legend()
    plt.grid()
    plt.show()



zad4(1)


### DLA CHETNYCH - PLOT u -> dodac do sprawozdania
