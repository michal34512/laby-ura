import numpy as np
from scipy.optimize import minimize
from scipy.optimize import dual_annealing

def method1(): # local
    def objective_function(x):
        return x**4 - 4*x**3 - 2*x**2 + 12*x + 9
    bounds = [(0, None)]
    x0 = np.array([0.5])
    result = minimize(objective_function, x0, bounds=bounds, method='L-BFGS-B')
    print("Punkt, w którym osiągnięto minimum:", result.x)
    print("Wartość minimalna funkcji:", result.fun)

def method2(): # global
    def objective_function(x):
        return x[0]**4 - 4*x[0]**3 - 2*x[0]**2 + 12*x[0] + 9
    bounds = [(0, 10000)]
    result = dual_annealing(objective_function, bounds)
    print("Punkt, w którym osiągnięto minimum globalne:", result.x)
    print("Wartość minimalna funkcji:", result.fun)
method1()
method2()