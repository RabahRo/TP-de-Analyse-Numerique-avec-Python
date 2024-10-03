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

def fct_u(x, y):
    return np.cos(x**2 + y**2)

def fct_u_Lap(x, y):
    return -4 * ((x**2 + y**2) * np.cos(x**2 + y**2) + np.sin(x**2 + y**2))

def Laplacien_2D(cas):
    if cas == 1:
        x0, y0 = 1, 1
        odd = 4
    else:
        x0, y0 = 0, 0
        odd = 6

    h = [5e-1, 1e-1, 5e-2, 1e-2, 5e-3]
    nn = len(h)
    err1 = np.zeros(nn)
    ordre1 = np.zeros(nn - 1)
    err2 = np.zeros(nn)
    ordre2 = np.zeros(nn - 1)

    k = 2
    coef1 = CoefDF(k, 0, [-1, 0, 1])
    nk1 = (len(coef1) + 1) // 2

    coef2 = CoefDF(k, 0, [-2, -1, 0, 1, 2])
    nk2 = (len(coef2) + 1) // 2

    L = fct_u_Lap(x0, y0)

    for i in range(nn):
        Lh1 = 0
        for j in range(len(coef1)):
            Lh1 += coef1[j] * (fct_u(x0 + (j - nk1) * h[i], y0) + fct_u(x0, y0 + (j - nk1) * h[i]))
        Lh1 /= h[i] ** k

        Lh2 = 0
        for j in range(len(coef2)):
            Lh2 += coef2[j] * (fct_u(x0 + (j - nk2) * h[i], y0) + fct_u(x0, y0 + (j - nk2) * h[i]))
        Lh2 /= h[i] ** k

        err1[i] = Lh1 - L
        err2[i] = Lh2 - L

        if i > 0:
            ordre1[i - 1] = np.log(err1[i - 1] / err1[i]) / np.log(h[i - 1] / h[i])
            ordre2[i - 1] = np.log(err2[i - 1] / err2[i]) / np.log(h[i - 1] / h[i])

    fig, ax = plt.subplots()
    ax.loglog(h, [h**2 for h in h], 'k-*', h, np.abs(err1), 'b-d', h, [h**odd for h in h], 'k-*', h, np.abs(err2), 'g-d')
    ax.grid(True)
    ax.set(xlabel='h', ylabel='erreur', title='Graphe log-log de la précision de la méthode')
    ax.legend(['h^2', 'Lh1', f'h^{odd}', 'Lh2'], loc='upper left')
    plt.show()

    coef_err1 = np.polyfit(np.log(h), np.log(np.abs(err1)), 1)
    coef_err2 = np.polyfit(np.log(h), np.log(np.abs(err2)), 1)

    pente_moyenne_err1 = coef_err1[0]
    pente_moyenne_err2 = coef_err2[0]

    print("Pente moyenne pour Lh1:", pente_moyenne_err1)
    print("Pente moyenne pour Lh2:", pente_moyenne_err2)
    print("Ordre pour Lh1:", ordre1)
    print("Ordre pour Lh2:", ordre2)

# Exemple d'appel
Laplacien_2D(0)  # Vous pouvez appeler la fonction avec l'argument de votre choix
