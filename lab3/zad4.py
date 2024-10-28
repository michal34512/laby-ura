import numpy as np
from scipy.optimize import LinearConstraint
from scipy.optimize import minimize

def model(t, a0, a1, a2, a3):
    """Zwraca wartość funkcji x(t) w danym czasie t."""
    return a0 + a1 * t + a2 * t**2 + a3 * t**3

def problem_dyn(a):
    """Oblicza całkę funkcji model w horyzoncie t ∈ [0, 1] z zadanymi parametrami a."""
    a0, a1, a2, a3 = a
    t = np.linspace(0, 1, 100) 
    x_values = model(t, a0, a1, a2, a3)
    integral = np.trapz(x_values, t)  
    return integral

# Ograniczenia liniowe
A = np.array([[1, 1, 1, 1]])  # a0 + a1 + a2 + a3
b = np.array([3])  # Musi być równe 3
bounds = [(None, None), (None, None), (None, None), (None, None)]  # Nieograniczone
constraint_eq = LinearConstraint(A, b, b)  # Ograniczenie równości

a_initial = [1, 0, 0, 0]  # Zakładamy, że a0=1, a1=a2=a3=0

# Optymalizacja
result = minimize(problem_dyn, a_initial, method='trust-constr', constraints=constraint_eq)

# Wynik
print("Optymalne wartości parametrów a0, a1, a2, a3:", result.x)
print("Minimalna wartość funkcji celu:", result.fun)
