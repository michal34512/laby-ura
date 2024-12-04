import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def wahadlo_amin(res, R, t):
    theta = res[:, 0]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(-1.5 * R, 1.5 * R)
    ax.set_ylim(-1.5 * R, 1.5 * R)
    ax.set_aspect('equal')
    ax.grid()

    line, = ax.plot([], [], 'o-', lw=2) 
    time_template = 'Czas: {:.1f} s'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    delta_t = t[1] - t[0]
    interval = delta_t * 1000  

    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text

    def update(frame):
        x = R * np.sin(theta[frame]) 
        y = -R * np.cos(theta[frame])  
        line.set_data([0, x], [0, y])
        time_text.set_text(time_template.format(t[frame]))
        return line, time_text

    ani = FuncAnimation(fig, update, frames=len(t), init_func=init, blit=True, interval=interval)
    plt.show()

def zad2():
    m = 0.5 # kg
    R = 1 # [m] dlugosc wahadla
    d = 0.5 # [Nms^2/rad^2]
    g = 10
    J = R**2 * m
    
    # czas
    t = np.linspace(0, 20, 1000)
    
    # sygnal zadany
    A = 1.5
    omega= 0.65 # [rad / s] 
    def tau(time):
        return A * np.cos(omega * time)

    def model(x, t):
            x1, x2 = x
            x1dot = x2
            x2dot = (1/J) * (tau(t) - d * x2 - (m * g * R * np.sin(x1)))
            return [x1dot, x2dot]
    y0 = [0, 0]
    res = odeint(model, y0, t)
    theta = res[: , 0]
    thetadot = res[:, 1]
    plt.plot(t, theta, label='Rozwiązanie y(t)')
    plt.xlabel('Czas t [s]')
    plt.ylabel('y(t)')
    plt.title('Rozwiązanie')
    plt.legend()
    plt.grid(True)
    plt.show()
    # - charakter oscylacyjny (logiczne bo bujamy wahadlo)
    #  - interpretacja fizyczna to wahadlo wprawiane w ruch wymuszeniem oscylacyjnym
    wahadlo_amin(res, R, t)
zad2()