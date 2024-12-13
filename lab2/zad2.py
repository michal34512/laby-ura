import numpy as np
import matplotlib . pyplot as plt
from scipy import signal
from scipy.integrate import odeint

def calman_matrix(A, B):
    columns = [np.linalg.matrix_power(A, i) @ B for i in range(A.shape[1])]
    result_matrix = np.column_stack(columns)
    return result_matrix

def getAs(A):
    eigenvalues = np.linalg.eigvals(A)
    polynomial_coeffs = np.poly(eigenvalues)
    rows = []
    for i in range(0, polynomial_coeffs.size-2):
        new_vec = np.zeros(polynomial_coeffs.size-1)
        new_vec[i+1] = 1
        rows.append(new_vec)
    rows.append(-polynomial_coeffs[1:][::-1])
    return np.vstack(rows)

def getBs(B):
    Bs = np.zeros_like(B)
    Bs[-1, 0] = 1
    return Bs

def getPinv(As, Bs, A, B):
    return calman_matrix(A, B) @ np.linalg.inv(calman_matrix(As, Bs))

def plot_all_state_vars_Pinv(A,B, Pinv, utype = 0):
    D = 0
    for i in range(A.shape[1]):
        C = np.zeros(A.shape[1])
        C[i] = 1
        if Pinv is not None:
            C = C @ Pinv
        
        SS = signal.StateSpace(A,B,C,D)
        if utype == 0:
            t, y = signal.step(SS)
        elif utype == 1:
            time = np.linspace(0, 40, 100)
            u = np.ones_like(time) * 2
            t, y, x_out = signal.lsim(SS, u, time)
        else:
            time = np.linspace(0, 40, 100)
            u = [ (np.sin(x) if x > 0 else 0) - 1/2 for x in time]
            t, y, x_out = signal.lsim(SS, u, time)
        plt.plot(t,y, label=f"x{i+1}")
    plt.legend()
    plt.show()
def plot_all_state_vars(A,B,C, u = 1):
    C = np.array(C)  
    def model(x, t):
        xdot = A @ x + B.flatten() * u
        return xdot
    t = np.linspace(0, 22, 1000)
    x = odeint(model, [0, 0, 0], t)
    y = x @ C.T 
    plt.plot(t, x[:, 0], label='x1(t)')
    plt.plot(t, x[:, 1], label='x2(t)')
    plt.plot(t, x[:, 2], label='x3(t)')
    plt.plot(t, y, label='y(t)')
    plt.legend()
    plt.show()

def zad21_rys2():
    C1 = 1
    R1 = 1
    C2 = 2
    R2 = 1
    C3 = 3
    R3 = 1

    ##### State-equasions #####
    A = np.array([[-1/(R1*C1), 0, 0],[0, -1/(R2*C2), 0],[0, 0, -1/(R3*C3)]])
    B = np.array([[1/(R1*C1)],[1/(R2*C2)],[1/(R3*C3)]])

    As = getAs(A)
    Bs = getBs(B)
    Pinv = getPinv(As, Bs, A, B)
    C = [[0.5, 1, 1.5]]
    Cs = C @ Pinv
    plot_all_state_vars(As, Bs, Cs)
    print(As)
    plot_all_state_vars(A, B, C)
    # jest sterowalny

def zad21_rys4():
    R1 = 2
    L1 = 0.5
    R2 = 1
    L3 = 1
    C4 = 2

    ##### State-equasions #####
    A = np.array([[-R1/L1, 0, -1/L1],[0, 0, 1/L3],[1/C4, -1/C4, -1/(R2*C4)]])
    B = np.array([[1/L1],[0],[0]])
    
    As = getAs(A)
    Bs = getBs(B)
    Pinv = getPinv(As, Bs, A, B)
    plot_all_state_vars_Pinv(As, Bs, Pinv , utype = 0)
    print(As)
    plot_all_state_vars_Pinv(A, B, None, utype = 0)
    # jest sterowalny

def calman_theorem_rank(A, B):
    columns = [np.linalg.matrix_power(A, i) @ B for i in range(A.shape[1])]
    result_matrix = np.column_stack(columns)
    print(f"Calman matrix: {result_matrix}")
    return np.linalg.matrix_rank(result_matrix)

def zad21_rys1():
    C1 = 1
    R1 = 2
    C2 = 0.5
    R2 = 4

    ##### State-equasions #####
    A = np.array([[-1/(R1*C1), 0],[0, -1/(R2*C2)]])
    B = np.array([[1/(R1*C1)],[1/(R2*C2)]])
    As = getAs(A)
    Bs = getBs(B)
    Pinv = getPinv(As, Bs, A, B)
    plot_all_state_vars_Pinv(As, Bs, Pinv , utype = 0)
    print(As)
    plot_all_state_vars_Pinv(A, B, None, utype = 0)
    # jest niesterowalny

def zad21_rys3():
    R = 1
    C1 = 0.1
    L1 = 0.1
    C2 = 0.1
    L2 = 0.1

    ##### State-equasions #####
    A = np.array([[0, 1/L1, 0, 0],[-1/C1, -1/(R*C1), 0, -1/(R*C1)],[0, 0, 0, 1/L2], [0, -1/(R*C2), -1/C2, -1/(R*C2)]])
    B = np.array([[0],[-1/(R*C1)],[0],[-1/(R*C2)]])
    As = getAs(A)
    Bs = getBs(B)
    Pinv = getPinv(As, Bs, A, B)
    plot_all_state_vars_Pinv(As, Bs, Pinv , utype = 0)
    plot_all_state_vars_Pinv(A, B, None, utype = 0)
    print(f"Rank before: {calman_theorem_rank(A, B)}")
    print(As)
    print(f"Rank after: {calman_theorem_rank(As, Bs)}")
    # nie jest sterowalny

if __name__ == "__main__":
    zad21_rys2()
    # Opisują obiekty w różny sposób ale efekt jest ten sam
    # Przebiegi są takie same, to oznacza że do projektowania układu regulacji możemy użyć postaci sterowalnej równań stanu
    
   
