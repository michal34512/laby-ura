import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def calman_matrix(A, B):
    columns = [np.linalg.matrix_power(A, i) @ B for i in range(A.shape[1])]
    result_matrix = np.column_stack(columns)
    return result_matrix

def zad1(u):
    m = 9 # kg
    l = 1 # [m] dlugosc wahadla
    d = 0.5 # [Nms^2/rad^2]
    g = 10
    J = 1
    
    # common
    t = np.linspace(0, 20, 1000)
    y0 = [0, 0]
    
    # model nieliniowy
    def lin_model(x, t):
            x1, x2 = x
            x1dot = x2
            x2dot = (1/J) * (u - d * x2 - (m * g * l * np.sin(x1)))
            return [x1dot, x2dot]
    
    # model liniowy
    A = np.array([[0, 1],[-(m*g*l)/J, -d/J]])
    B = np.array([[0],[1/J]])
    def model(x, t):
        xdot = A @ x + B.flatten() * u
        return xdot
    plt.figure(figsize=(10, 5))
    u=5
    y = odeint(model, y0, t)
    y_lin = odeint(lin_model, y0, t)
    plt.plot(t, y_lin[:, 0], label=r'$układ\ zlinearyzowany,\ u=5$', color='blue')
    plt.plot(t, y[:, 0], label=r'$układ\ nieliniowy,\ u=5$', color=(0.5, 0.7, 1))
    u=20
    y = odeint(model, y0, t)
    y_lin = odeint(lin_model, y0, t)
    plt.plot(t, y_lin[:, 0], label=r'$układ\ zlinearyzowany,\ u=20$', color='red')
    plt.plot(t, y[:, 0], label=r'$układ\ nieliniowy,\ u=20$', color=(1.0, 0.6, 0.6))
    plt.xlabel('Czas t')
    plt.ylabel('x1(t)')
    plt.title(r'$u(t) = 5, u(t) = 20$')
    plt.legend()
    plt.grid(True)
    plt.show()

    print(f"Rank ukladu:{np.linalg.matrix_rank(calman_matrix(A, B))}")


# wartosci u (nieliniowy):
# 0 = wahadlo sie nie rusza
# 5, 20, 45*sqrt(2) = buja sie i stabilizuje
# 70 = kreci sie w nieskonczonosc
# wartosci u (nieliniowy):
# wahadlo zawsze sie na czyms stabilizuje dla wymuszenia do 20 nie widac duzej roznicy

# - Czy model zlinearyzowany wiernie odzwierciedla zachowanie nieliniowego obiektu
# dla wszystkich wymuszeń?
# Nie, tylko dla malych wzmocnien
# 
# - Jaka jest interpretacja wartości wymuszenia, w odniesieniu do punktu równowagi?
# im mocniejsze wymuszenie tym dalej wahadlo wychyla sie od punktu rownowagi (punktu pracy)
# 
# Czy jest sterowalny?
# Tak, mozna uzyc macierz kalmana

def zad2(u):
    m = 9 # kg
    l = 1 # [m] dlugosc wahadla
    d = 0.5 # [Nms^2/rad^2]
    g = 10
    J = 1
    
    # common
    t = np.linspace(0, 20, 1000)
    x0 = [np.pi/4, 0]

    # nieliniowy
    def model(x, t):
            x1, x2 = x
            x1dot = x2
            x2dot = (1/J) * (u - d * x2 - (m * g * l * np.sin(x1)))
            return [x1dot, x2dot]
    
    # y = odeint(model, x0, t)


    # liniowy
    u0 = 45 * np.sqrt(2)
    ufala = u - u0

    A = np.array([[0, 1],[-(m*g*l*np.sqrt(2))/(2*J), -d/J]])
    B = np.array([[0],[1/J]])
    def lin_model(xfala, t):
        xfaladot = A @ xfala + B.flatten() * ufala
        return xfaladot
    # y_lin = odeint(model, [0, 0], t) 

    plt.figure(figsize=(10, 5))
    u = 45*np.sqrt(2)
    ufala = u - u0
    y = odeint(model, x0, t)
    y_lin = odeint(lin_model, [0, 0], t) + x0
    plt.plot(t, y_lin[:, 0], label=r'$układ\ zlinearyzowany\ u(t)=45\sqrt{2}$', color='green')  # Zielony
    plt.plot(t, y[:, 0], label=r'$układ\ nieliniowy\ u(t)=45\sqrt{2}$', color=(0.3, 0.8, 0.3))  # Ciemniejszy jasnozielony
    u = 45*np.sqrt(2) + 10
    ufala = u - u0
    y = odeint(model, x0, t)
    y_lin = odeint(lin_model, [0, 0], t) + x0
    plt.plot(t, y_lin[:, 0], label=r'$układ\ zlinearyzowany\ u(t)=45\sqrt{2}+10$', color='blue')  # Niebieski
    plt.plot(t, y[:, 0], label=r'$układ\ nieliniowy\ u(t)=45\sqrt{2}+10$', color=(0.5, 0.7, 1))  # Jasnoniebieski
    plt.xlabel('Czas t')
    plt.ylabel('x1(t)')
    plt.title(r'$Symulacja\ układu\ zlinearyzowanego\ w\ punkcie\ x_0(\pi/4,0), u_0=45\sqrt{2}$')
    plt.legend()
    plt.grid(True)
    plt.show()

def zad3(u):
    m = 9  # kg
    l = 1  # m (długość wahadła)
    d = 0.5  # Nms^2/rad^2 (damping)
    g = 10  # m/s^2 (przyspieszenie grawitacyjne)
    J = 1  # moment bezwładności
    t = np.linspace(0, 20, 1000)
    x0 = [np.pi / 4, 0]  # Początkowy kąt i prędkość
    def model_nieliniowy(x, t):
        x1, x2 = x
        x1dot = x2
        x2dot = (1 / J) * (u - d * x2 - (m * g * l * np.sin(x1)))
        return [x1dot, x2dot]
    res_nieliniowy = odeint(model_nieliniowy, x0, t)
    theta_nieliniowy = res_nieliniowy[:, 0]
    
    def model_SDC(x, t):
        x1, x2 = x
        A = np.array([[0, 1], 
                      [-(m * g * l * np.sin(x1)) / x1, -d / J]])
        B = np.array([0, 1 / J])
        dxdt = A @ x + B * u
        return dxdt
    res_SDC = odeint(model_SDC, x0, t)
    theta_SDC = res_SDC[:, 0]
    
    # Wykresy
    plt.figure(figsize=(10, 6))
    plt.plot(t, theta_nieliniowy, label="Oryginalny model nieliniowy", linewidth=2)
    plt.plot(t, theta_SDC, '--', label="Model SDC", linewidth=2)
    plt.xlabel("Czas [s]")
    plt.ylabel("Kąt [rad]")
    plt.title("Porównanie odpowiedzi układu nieliniowego i SDC")
    plt.legend()
    plt.grid()
    plt.show()

    # Analiza dla warunku początkowego x1(0)=0, x2(0)=0
    x0_zero = [0, 0]
    res_nieliniowy_zero = odeint(model_nieliniowy, x0_zero, t)
    res_SDC_zero = odeint(model_SDC, x0_zero, t)
    
    plt.figure(figsize=(10, 6))
    plt.plot(t, res_nieliniowy_zero[:, 0], label="Nieliniowy (x1(0)=0, x2(0)=0)", linewidth=2)
    plt.plot(t, res_SDC_zero[:, 0], '--', label="SDC (x1(0)=0, x2(0)=0)", linewidth=2)
    plt.xlabel("Czas [s]")
    plt.ylabel("Kąt [rad]")
    plt.title("Porównanie przy warunkach początkowych x1(0)=0, x2(0)=0")
    plt.legend()
    plt.grid()
    plt.show()

zad1(5)



