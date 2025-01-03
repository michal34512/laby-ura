import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.linalg import solve_continuous_are
from scipy.interpolate import interp1d

def zad1(u):
    k = 3
    J1 = 0.04
    J2 = 0.3
    m = 0.5
    l = 0.5
    g = 9.81

    x0 = [0, 0, 0, 0]
    t = np.linspace(0, 10, 1000)    
    def model(x, t):
            x1, x2, x3, x4 = x
            x1dot = x2
            x3dot = x4
            x2dot = -m*g*l*np.sin(x1)/J1 - k*(x1-x3)/J1
            x4dot = k*(x1- x3)/J2 + u/J2
            return [x1dot, x2dot, x3dot, x4dot]
    x = odeint(model, x0, t)
    plt.figure(figsize=(10, 5))
    plt.plot(t, x[:, 0], label='kat_nieliniowy(t)')
    plt.xlabel('Czas t')
    plt.ylabel('theta(t)')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.grid(True)
    plt.show()


def zad2(u):
    k = 3
    J1 = 0.04
    J2 = 0.3
    m = 0.5
    l = 0.5
    g = 9.81

    x_approx = [0, 0, 0, 0]
    t = np.linspace(0, 10, 1000)
    A = np.array([[0, 1, 0, 0], 
                  [-(k + m*g*l*np.cos(x_approx[0]))/J1, 0, k/J1, 0], 
                  [0, 0, 0, 1], 
                  [k/J2, 0, -k/J2, 0]])
    B = np.array([[0],[0],[0],[1/J2]])

    
    def model(x, t):
            xdot = A @ x + B.flatten() * u
            return xdot
    x = odeint(model, [np.pi/4, 0, 0, 0], t)
    plt.figure(figsize=(10, 5))
    plt.plot(t, x[:, 0], label='kat_liniowy(t)')
    plt.xlabel('Czas t')
    plt.ylabel('theta(t)')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.grid(True)
    plt.show()
#zad2(0)

def get_P_mat(A, B, t_end, Q, R, S):
    def riccati(p, t):
        P = np.reshape(p, A.shape)
        P_dot = -(P @ A  - P @ B @ np.linalg.inv(R) @ B.T @ P + A.T @ P + Q)
        return P_dot.reshape(-1)
    P0_vec = S.reshape(-1)
    t_riccati = np.linspace(t_end, 0, 1000)
    P_vec = odeint(riccati, P0_vec, t_riccati)
    P_matrices = np.array([np.reshape(p, A.shape) for p in P_vec])
    interpolators = {}
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            key = (i, j)
            interpolators[key] = interp1d(
                t_riccati,
                P_matrices[:, i, j],
                fill_value="extrapolate"
            )
    def P_t(t):
        return np.array([
            [interpolators[(i, j)](t) for j in range(A.shape[1])]
            for i in range(A.shape[0])
        ])
    return P_t

def get_K_mat(A, B, Q, R):
    P = solve_continuous_are(A, B, Q, R)
    K = np.linalg.inv(R) @ B.T @ P
    return K
def zad3_inf(Q, R):
    k = 3
    J1 = 0.04
    J2 = 0.3
    m = 0.5
    l = 0.5
    g = 9.81

    x_approx = [0, 0, 0, 0]
    t = np.linspace(0, 10, 1000)
    A = np.array([[0, 1, 0, 0], 
                  [-(k + m*g*l*np.cos(x_approx[0]))/J1, 0, k/J1, 0], 
                  [0, 0, 0, 1], 
                  [k/J2, 0, -k/J2, 0]])
    B = np.array([[0],[0],[0],[1/J2]])
    K = get_K_mat(A, B, Q, R)
    def model(x, t):
            u = -K @ x
            xdot = A @ x + B.flatten() * u
            return xdot
    x = odeint(model, [np.pi/4, 0, 0, 0], t)
    plt.figure(figsize=(10, 5))
    plt.plot(t, x[:, 0], label='kat_liniowy(t)')
    plt.xlabel('Czas t')
    plt.ylabel('theta(t)')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.grid(True)
    plt.show()

def zad3_fin(Q, R, S):
    k = 3
    J1 = 0.04
    J2 = 0.3
    m = 0.5
    l = 0.5
    g = 9.81

    x_approx = [0, 0, 0, 0]
    t = np.linspace(0, 10, 1000)
    A = np.array([[0, 1, 0, 0], 
                  [-(k + m*g*l*np.cos(x_approx[0]))/J1, 0, k/J1, 0], 
                  [0, 0, 0, 1], 
                  [k/J2, 0, -k/J2, 0]])
    B = np.array([[0],[0],[0],[1/J2]])

    P_t = get_P_mat(A, B, 10, Q, R, S)
    
    def model(x, t):
            K_t = np.linalg.inv(R) @ (B.T @ P_t(t))
            u = -K_t @ x
            xdot = A @ x + B.flatten() * u
            return xdot
    x = odeint(model, [np.pi/4, 0, 0, 0], t)
    plt.figure(figsize=(10, 5))
    plt.plot(t, x[:, 0], label='kat_liniowy(t)')
    plt.xlabel('Czas t')
    plt.ylabel('theta(t)')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.grid(True)
    plt.show()
# zad3_fin(np.eye(4) , np.array([[1]]), np.eye(4))
# zad3_inf(np.eye(4) , np.array([[1]]))

def zad4(Q, R, S, t_end, x0, u0=0):
    k = 3
    J1 = 0.04
    J2 = 0.3
    m = 0.5
    l = 0.5
    g = 9.81

    x_approx = [0, 0, 0, 0]

    t = np.linspace(0, 10, 1000)
    A = np.array([[0, 1, 0, 0], 
                  [-(k + m*g*l*np.cos(x_approx[0]))/J1, 0, k/J1, 0], 
                  [0, 0, 0, 1], 
                  [k/J2, 0, -k/J2, 0]])
    B = np.array([[0],[0],[0],[1/J2]])

    finite_P_t = get_P_mat(A, B, t_end, Q, R, S) # finite
    infinite_K = get_K_mat(A, B, Q, R) # infinite
    
    def finite_model(x, t):
        K_t = np.linalg.inv(R) @ (B.T @ finite_P_t(t))
        u = -K_t @ (x - x0) + u0
        xdot = A @ x + B.flatten() * u
        return xdot
    
    def infinite_model(x, t):
        u = -infinite_K @ (x - x0) + u0
        xdot = A @ x + B.flatten() * u
        return xdot
    
    finite_x = odeint(finite_model, x0, t) + x0
    infinite_x = odeint(infinite_model, x0, t) + x0
    plt.figure(figsize=(10, 5))
    plt.plot(t, finite_x[:, 0], label='horyzont skończony')
    plt.plot(t, infinite_x[:, 0], label='horyzont nieskończony')
    plt.xlabel('Czas t')
    plt.ylabel('theta(t)')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.grid(True)
    plt.show()

zad4(np.eye(4) , np.array([[1]]), np.eye(4), 5, np.array([np.pi, 0, np.pi/2, 0]))