import numpy as np
import matplotlib.pyplot as plt

def CoefDF(k, xbar, x):
    x = np.array(x)
    n = len(x)
    A = np.zeros((n, n))
    B = np.zeros(n)
    h = min(x[i] - x[i - 1] for i in range(1, n))
    #h = min(x[1:] - x[:-1])
    h2 = min(np.abs(x - xbar))
    if h2 > 0:
        h = min(h, h2)
    p = n - k
    for i in range(n):
        for j in range(n):
            A[i, j] = (x[j] - xbar)**(i) / np.math.factorial(i)

    B[k] = 1

    coef = np.linalg.solve(A, B)
    coef = coef * h**k

    return coef
# Exemple de données
h = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
y = np.array([1.0, 4.0, 9.0, 16.0, 25.0])  # Par exemple, y = x^2

# Point autour duquel nous voulons calculer la dérivée
xbar = 3.0

# Degré de la dérivée que nous voulons calculer (1 pour la première dérivée)
k = 1

# Calcul des coefficients pour la première dérivée
coefficients = CoefDF(k, xbar, h)

# Affichage des coefficients
print("Coefficients de la première dérivée :", coefficients)

# Calcul de la première dérivée approximative en utilisant les coefficients
# Cette dérivée approximative est un polynôme de degré k
first_derivative = np.polyval(coefficients, xbar)

# Affichage de la première dérivée
print("Dérivée approximative à xbar :", first_derivative)

# Trace de la fonction et de la tangente à xbar
plt.plot(h, y, label="y = x^2")
plt.plot(xbar, xbar ** 2, "ro", label=f"xbar = {xbar}")
plt.plot(h, np.polyval(coefficients, h), label="Tangente à xbar")
plt.legend()
plt.xlabel("h")
plt.ylabel("y")
plt.title("Approximation de la première dérivée")
plt.show()
