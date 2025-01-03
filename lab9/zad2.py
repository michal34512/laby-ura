import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.linalg import solve_continuous_are

def zad1(u): # nieliniowe wahadlo
    m = 9 # kg
    l = 1 # [m] dlugosc wahadla
    d = 0.5 # [Nms^2/rad^2]
    g = 9.81
    J = 1
    
    # common
    t = np.linspace(0, 10, 1000)
    x0 = [np.pi/4, 0]

    # nieliniowy
    def model(x, t):
            x1, x2 = x
            x1dot = x2
            x2dot = (u - d * x2)/J + (m * g * l * np.sin(x1))
            return [x1dot, x2dot]
    
    res = odeint(model, x0, t)
    theta = res[: , 0]
    thetadot = res[:, 1]

    plt.figure(figsize=(10, 5))
    plt.plot(t, theta, label='kat_nielinowy')
    plt.xlabel('Czas t')
    plt.ylabel('theta(t)')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.grid(True)
    plt.show()
    
# zad1(0)

def zad2(u): # wahadło SDC
    m = 9 # kg
    l = 1 # [m] dlugosc wahadla
    d = 0.5 # [Nms^2/rad^2]
    g = 9.81
    J = 1
    
    # common
    t = np.linspace(0, 10, 1000)
    x0 = [np.pi/4, 0]

    # liniowy
    def A(x0):
        return np.array([[0, 1],[m*g*l*np.sin(x0[0])/x0[0], -d/J]]) 
    B = np.array([[0], [1/J]])
    
    def model(x, t):
            xdot = A(x) @ x + B.flatten() * u
            return xdot
    
    res = odeint(model, x0, t)
    theta = res[: , 0]
    thetadot = res[:, 1]

    plt.figure(figsize=(10, 5))
    plt.plot(t, theta, label='kat_nielinowy')
    plt.xlabel('Czas t')
    plt.ylabel('theta(t)')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.grid(True)
    plt.show()


# zad2(0)
def get_K_mat(A, B, Q, R):
    P = solve_continuous_are(A, B, Q, R)
    K = np.linalg.inv(R) @ B.T @ P
    return K

def zad3(Q, R): # wahadło SDRE
    m = 9 # kg
    l = 1 # [m] dlugosc wahadla
    d = 0.5 # [Nms^2/rad^2]
    g = 9.81
    J = 1
    
    # common
    t = np.linspace(0, 5, 1000)
    x0 = [2*np.pi, 0]

    # liniowy
    saved_u = []
    saved_t = []
    def A(x0):
        return np.array([[0, 1],[m*g*l*np.sin(x0[0])/x0[0], -d/J]]) 
    B = np.array([[0], [1/J]])
    
    def model(x, t):
            K = get_K_mat(A(x), B, Q, R)
            u = -K @ x
            xdot = A(x) @ x + B.flatten() * u
            saved_u.append(u)
            saved_t.append(t)
            return xdot
    
    res = odeint(model, x0, t)
    theta = res[: , 0]
    thetadot = res[:, 1]

    plt.figure(figsize=(10, 5))
    plt.plot(t, theta, label='sqre')
    plt.plot(saved_t, saved_u, label="u(t)")
    plt.xlabel('Czas t')
    plt.ylabel('theta(t)')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.grid(True)
    plt.show()

zad3(np.eye(2), np.array([[1]]))

def zad4(R): # wahadło SDRE + zmienne Q
    m = 9 # kg
    l = 1 # [m] dlugosc wahadla
    d = 0.5 # [Nms^2/rad^2]
    g = 9.81
    J = 1
    
    # common
    t = np.linspace(0, 5, 1000)
    x0 = [2*np.pi, 0]

    # liniowy
    saved_u = []
    saved_t = []
    def A(x0):
        return np.array([[0, 1],[m*g*l*np.sin(x0[0])/x0[0], -d/J]]) 
    B = np.array([[0], [1/J]])
    
    def model(x, t):
            Q = np.array([[x[0]**2, 0],[0, x[1]**2]]) 
            K = get_K_mat(A(x), B, Q, R)
            u = -K @ x
            xdot = A(x) @ x + B.flatten() * u
            saved_u.append(u)
            saved_t.append(t)
            return xdot
    
    res = odeint(model, x0, t)
    theta = res[: , 0]
    thetadot = res[:, 1]

    plt.figure(figsize=(10, 5))
    plt.plot(t, theta, label='sqre')
    plt.plot(saved_t, saved_u, label="u(t)")
    plt.xlabel('Czas t')
    plt.ylabel('theta(t)')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.grid(True)
    plt.show()

zad4(np.array([[1]]))