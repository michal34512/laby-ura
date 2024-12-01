import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.linalg import solve_continuous_are

def zad3(Q, R, qd):
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
zad3(np.eye(2)*1000, [[1]], 20)



# def zad3(_Q, _R, qd):
#     L = 0.2
#     R = 0.5
#     C = 0.5
#     A = np.array([[0, 1],
#                 [-1/(L*C), -R/L]])
#     B = np.array([[0],
#                 [1/L]])
#     P = solve_continuous_are(A, B, _Q, _R)
#     K = np.linalg.inv(_R) @ B.T @ P
#     xd = [qd, 0]
#     saved_u = []
#     saved_t = []
#     def model(x, t):
#         u = (qd/C) - K @ (x - xd)
#         dxdt = A @ x + B.flatten() * u
#         saved_u.append(u)
#         saved_t.append(t)
#         return dxdt
#     x0 = [1, 1]
#     t = np.linspace(0, 5, 500)

#     x = odeint(model, x0, t)

#     plt.figure(figsize=(10, 6))
#     plt.plot(t, x[:, 0], label="x1(t)")
#     plt.plot(t, x[:, 1], label="x2(t)")
#     plt.plot(saved_t, saved_u, label="u(t)")
#     plt.title("Odpowiedź skokowa układu zamkniętego")
#     plt.xlabel("Czas [s]")
#     plt.ylabel("Stan")
#     plt.legend()
#     plt.grid()
#     plt.show()
# zad3(np.eye(2)*1000, [[1]], 20)