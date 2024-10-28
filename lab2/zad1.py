import numpy as np
import matplotlib . pyplot as plt
from scipy import signal
from scipy.integrate import odeint

def calman_theorem_rank(A, B):
    columns = [np.linalg.matrix_power(A, i) @ B for i in range(A.shape[1])]
    result_matrix = np.column_stack(columns)
    print(f"Calman matrix: {result_matrix}")
    return np.linalg.matrix_rank(result_matrix)

def plot_all_state_vars(A,B, utype = 0):
    D = 0
    for i in range(A.shape[1]):
        C = np.zeros(A.shape[1])
        C[i] = 1
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

def zad11_1():
    # U = uR + Uc
    # U = RC(duC)/(dt) + uC
    # uC' = U/(RC) - uC/(RC)
    C1 = 1
    R1 = 2
    C2 = 0.5
    R2 = 4

    ##### State-equasions #####
    A = np.array([[-1/(R1*C1), 0],[0, -1/(R2*C2)]])
    B = np.array([[1/(R1*C1)],[1/(R2*C2)]])
    print(f"Rank: {calman_theorem_rank(A, B)}")
    plot_all_state_vars(A, B, 2)
    # jest niesterowalny

def zad11_2():
    C1 = 1
    R1 = 1
    C2 = 2
    R2 = 1
    C3 = 3
    R3 = 1

    ##### State-equasions #####
    A = np.array([[-1/(R1*C1), 0, 0],[0, -1/(R2*C2), 0],[0, 0, -1/(R3*C3)]])
    B = np.array([[1/(R1*C1)],[1/(R2*C2)],[1/(R3*C3)]])
    print(f"Rank: {calman_theorem_rank(A, B)}")
    plot_all_state_vars(A,B)
    # jest sterowalny

def zad11_3():
    R = 1
    C1 = 0.1
    L1 = 0.1
    C2 = 0.1
    L2 = 0.1

    ##### State-equasions #####
    A = np.array([[0, 1/L1, 0, 0],[-1/C1, -1/(R*C1), 0, -1/(R*C1)],[0, 0, 0, 1/L2], [0, -1/(R*C2), -1/C2, -1/(R*C2)]])
    B = np.array([[0],[1/(R*C1)],[0],[1/(R*C2)]])
    
    print(f"Rank: {calman_theorem_rank(A, B)}")
    plot_all_state_vars(A,B, 2)
    # nie jest sterowalny

def zad11_4():
    R1 = 2
    L1 = 0.5
    R2 = 1
    L3 = 1
    C4 = 2

    ##### State-equasions #####
    A = np.array([[-R1/L1, 0, -1/L1],[0, 0, 1/L3],[1/C4, -1/C4, -1/(R2*C4)]])
    B = np.array([[1/L1],[0],[0]])
    
    print(f"Rank: {calman_theorem_rank(A, B)}")
    plot_all_state_vars(A,B, 2)
    # jest sterowalny

if __name__ == "__main__":
    zad11_2()
   
