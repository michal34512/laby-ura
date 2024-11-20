from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt


def main(): 
# Krok (a): Ustawienie trybu optymalizacji dynamicznej
    m = GEKKO(remote=False)
    m.options.IMODE = 6  # Tryb optymalizacji dynamicznej

    m.time = np.linspace(0, 1.01, 101)  # Dyskretyzacja przedziału [0, 1]
                                        # Dodajemy jedną próbkę dla dokładnego
                                        # rozwiązania numerycznego
    x = m.Var(value=1)
    m.fix(x, pos=0, val=1)
    m.fix(x, pos=len(m.time)-1, val=3)

    J = m.Var(value=0)
    t = m.Var(value=0)
    m.Equation(t.dt() == 1)
    m.Equation(J.dt() == 24*x*t+2*(x.dt()**2)-4*t)  

    Jf = m.FV()
    Jf.STATUS = 1 
    m.Connection(Jf ,J , pos2 = 'end')
    m.Obj(Jf)

    m.solve(disp=False)

    #analityczne
    def J_a(t):
        return t**3+t+1
    y_a = J_a(m.time)

    plt.figure(figsize=(10, 5))
    # Ucinamy odtatnią próbkę
    plt.plot(m.time[0:100], x.value[0:100], label="Rozwiązanie numeryczne x(t)")
    plt.plot(m.time[0:100], y_a[0:100], label="Analitycznie x(t)")
    plt.xlabel('Czas t')
    plt.ylabel('Wartość')
    plt.legend()
    plt.title('Rozwiązanie problemu optymalizacji dynamicznej')
    plt.show()

    # Dzwina próbka na końcu wynikała z -1 przy m.fix(x, pos=len(m.time)-1, val=3)

main()
