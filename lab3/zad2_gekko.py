from gekko import GEKKO

def zad21():
    # Tworzenie modelu Gekko
    m = GEKKO(remote=False)

    # Definiowanie zmiennych
    x = m.Var(value=0)  # Wartość początkowa x
    y = m.Var(value=0)  # Wartość początkowa y

    # Funkcja celu: minimalizujemy -y, aby zmaksymalizować y
    m.Obj(y)

    # Definicje ograniczeń
    m.Equation(2 * x - y <= 4)  # 2x - y <= 4
    m.Equation(y + x >= 3)      # y + x > 3
    m.Equation(y + 4 * x >= -2) # y + 4x >= -2

    # Wybór solvera i uruchomienie optymalizacji
    m.solve(disp=False)

    # Wyniki
    print("Punkt, w którym osiągnięto maksimum y:", x.value[0], y.value[0])
    print("Wartość funkcji celu (maksymalizacja -y):", y.value[0])

def main():
    zad21()

if __name__ == "__main__":
    main()
