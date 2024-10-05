
import numpy as np
import matplotlib . pyplot as plt
from scipy import signal
from scipy.integrate import odeint

m = 1
L = 0.5
d = 0.1
J = 1/3 * (m*L**2)

def zad41():
    ##### State-equasions #####
    A = np.array([[0, 1],[0, -(d/J)]])
    B = np.array([[0],[1/J]])
    C = np.array([[1, 0]])
    D = 0
    SS = signal.StateSpace(A,B,C,D)
    t, y = signal.step(SS)
    plt.plot(t, y, label = "SS step")
    # odpowiedź skokowa ma charakter narastający stabilizujący się do liniowego

def zad42():
    ##### State-equasions #####
    A = np.array([[0, 1],[0, -(d/J)]])
    B = np.array([[0],[1/J]])
    C = np.array([[1, 0]])
    D = 0
    SS = signal.StateSpace(A,B,C,D)
    t = np.linspace(0, 10, 100)
    tau = t # Liniowo narastajacy sygnał zadany y = x
    t_out, y_out, x_out = signal.lsim(SS, tau, t) # x_out zawiera stany dla każdej chwili
    plt.plot(t_out, y_out, label = "SS linear")

def zad43():
    ##### State-equasions #####
    A = np.array([[0, 1],[0, -(d/J)]])
    B = np.array([[0],[1/J]])
    C = np.array([[1, 0]])
    D = 0
    SS = signal.StateSpace(A,B,C,D)
    w, mag, phase = signal.bode(SS)
    plt.subplot(2, 1, 1)
    plt.semilogx(w, mag) 
    plt.title('Bode')
    plt.ylabel('Gain [dB]')
    plt.grid()
    plt.subplot(2, 1, 2)
    plt.semilogx(w, phase)
    plt.ylabel('Phase [deg]')
    plt.xlabel('Freq [rad/s]')
    plt.grid()
    # W tak zamodelowanym układzie mamy do czynienia z obiektem inercyjnym z astatyzmem. 
    # Wykresy bodgo to potwierdzają (faza w nieskończoności jest równa -180 co świadczy o 
    # astatyźmie, widać też lekkie przełamanie na wykiesie modułu, co świadczy o inercji. 
    # Przełamanie występuje dla częstotliwości d/J = 1.2)

if __name__ == "__main__":
    zad41()
    plt.legend()
    plt.show()
    zad42()
    plt.legend()
    plt.show()
    zad43()
    plt.legend()
    plt.show()
