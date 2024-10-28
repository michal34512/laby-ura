from gekko import GEKKO

def method1():  # lokalne minimum
    m = GEKKO(remote=False)
    
    # Definiowanie zmiennej
    x = m.Var(value=0.5, lb=0)  # Wartość początkowa 0.5, ograniczenie dolne x >= 0
    
    # Funkcja celu do minimalizacji
    m.Obj(x**4 - 4*x**3 - 2*x**2 + 12*x + 9)
    
    # Uruchomienie solvera
    m.solve(disp=False)
    
    # Wyświetlanie wyników
    print("Punkt, w którym osiągnięto lokalne minimum:", x.value[0])
    print("Wartość minimalna funkcji:", m.options.objfcnval)

def method2():  # globalne minimum
    m = GEKKO(remote=False)
    
    # Definiowanie zmiennej z szerszymi granicami poszukiwań
    x = m.Var(value=5, lb=0, ub=10)  # wartość początkowa 5, granice poszukiwań 0 <= x <= 10000
    
    # Funkcja celu do minimalizacji
    m.Obj(x**4 - 4*x**3 - 2*x**2 + 12*x + 9)
    
    # Ustawienie solvera na globalny (opcjonalne)
    m.options.SOLVER = 1  # APOPT, który wspiera optymalizację globalną w niektórych przypadkach
    
    # Uruchomienie solvera
    m.solve(disp=False)
    
    # Wyświetlanie wyników
    print("Punkt, w którym osiągnięto globalne minimum:", x.value[0])
    print("Wartość minimalna funkcji:", m.options.objfcnval)

# Uruchomienie metod
method1()
method2()
