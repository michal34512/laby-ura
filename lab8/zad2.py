import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.linalg import solve_continuous_are

# Przetestować globalność
# Zad 2.2  LQR + niel, LQR + lin

def zad1(u):
    m = 9 # kg
    l = 1 # [m] dlugosc wahadla
    d = 0.5 # [Nms^2/rad^2]
    g = 9.81
    J = 1
    
    # common
    t = np.linspace(0, 10, 1000)
    x0 = [np.pi/4, 0]
    x0aprox = [np.pi, 0]

    # nieliniowy
    def model(x, t):
            x1, x2 = x
            x1dot = x2
            x2dot = (1/J) * (u - d * x2 - (m * g * l * np.sin(x1)))
            return [x1dot, x2dot]
    
    res = odeint(model, x0, t)
    theta = res[: , 0]
    thetadot = res[:, 1]


    # liniowy
    u0 = 0
    ufala = u - u0

    A = np.array([[0, 1],[-(m*g*l*np.cos(x0aprox[0]))/J, -d/J]])
    B = np.array([[0],[1/J]])
    print(np.linalg.eig(A))
    print(B)
    def model(xfala, t):
        xfaladot = A @ xfala + B.flatten() * ufala
        return xfaladot
    xfala = odeint(model, x0, t)
    x = xfala + x0

    plt.figure(figsize=(10, 5))
    plt.plot(t, x[:, 0], label='kat_liniowy(t)')
    plt.plot(t, theta, label='kat_nielinowy')
    plt.xlabel('Czas t')
    plt.ylabel('theta(t)')
    plt.title('Symulacja układu równań stanu')
    plt.legend()
    plt.grid(True)
    plt.show()

def get_K_mat(A, B, _Q, _R):
    P = solve_continuous_are(A, B, _Q, _R)
    K = np.linalg.inv(_R) @ B.T @ P
    return K

def zad2(x0, Q, R):
    m = 9 # kg
    l = 1 # [m] dlugosc wahadla
    d = 0.5 # [Nms^2/rad^2]
    g = 9.81
    J = 1
    
    # common
    t = np.linspace(0, 2, 1000)
    u = 0
    x0aprox = np.array([np.pi, 0])

    # liniowy
    u0 = 0
    A = np.array([[0, 1],[-(m*g*l*np.cos(x0aprox[0]))/J, -d/J]])
    B = np.array([[0],[1/J]])
    K = get_K_mat(A, B, Q, R)
    print(K)

    # def model(xfala, t):
    #     ufala = -K @ xfala
    #     xfaladot = A @ xfala + B.flatten() * ufala
    #     return xfaladot
    def model(x, t):
            x1, x2 = x
            ufala = -K[0, 0] * x1 - K[0, 1] * x2
            x1dot = x2
            x2dot = (1/J) * (ufala - d * x2 - (m * g * l * np.sin(x1)))
            return [x1dot, x2dot]
    
    x0 = np.array([0, 0])
    plt.figure(figsize=(10, 5))
    x = odeint(model, x0 -  x0aprox, t) + x0aprox
    plt.plot(t, x[:, 0], label=r'$x(0) = (0, 0)$')
    x0 = np.array([-np.pi/2, 0])
    x = odeint(model, x0 -  x0aprox, t) + x0aprox
    plt.plot(t, x[:, 0], label=r'$x(0) = (-\pi/2, 0)$')
    x0 = np.array([np.pi*3/4, 0])
    x = odeint(model, x0 -  x0aprox, t) + x0aprox
    plt.plot(t, x[:, 0], label=r'$x(0) = (3/4\pi, 0)$')
    x0 = np.array([np.pi*3/2, 0])
    x = odeint(model, x0 -  x0aprox, t) + x0aprox
    plt.plot(t, x[:, 0], label=r'$x(0) = (3/2\pi, 0)$')
    plt.xlabel('Czas t')
    plt.ylabel('x_0(t)')
    plt.title('Symulacja nieliniowego układu z regulatorem LQR')
    plt.legend()
    plt.grid(True)
    plt.show()

# zad1(0)
zad2(np.array([0, 0]), np.eye(2), np.eye(1))