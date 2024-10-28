import numpy as np
from scipy.optimize import minimize

def zad21():
    def objective_function(z):
        x, y = z
        return y
    def constraint1(z):
        x, y = z
        return 4 - (2 * x - y)
    def constraint2(z):
        x, y = z
        return y + x - 3
    def constraint3(z):
        x, y = z
        return y + 4 * x + 2
    x0 = np.array([0, 0])
    constraints = [
        {'type': 'ineq', 'fun': constraint1},  # 2x - y <= 4
        {'type': 'ineq', 'fun': constraint2},  # y + x > 3
        {'type': 'ineq', 'fun': constraint3}   # y + 4x >= -2
    ]
    result = minimize(objective_function, x0, constraints=constraints)
    print("Punkt, w którym osiągnięto minimum:", result.x)
    print("Wartość funkcji celu (maksymalizacja y):", result.fun)
    #minimum lokalne - algorytm zaczyna w x0

def main():
    zad21()
    
if __name__ == "__main__":
    main()
