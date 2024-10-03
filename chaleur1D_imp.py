import numpy as np
import matplotlib.pyplot as plt

def chaleur1D_imp(cfl):
    # Paramètres du problème
    L = 10
    nu = 1
    T = 10

    Nx = int(input("Nombre de mailles : "))
    h = L / Nx

    if cfl == 0:
        Nt = int(input("Nombre de pas de temps : "))
        deltat = T / Nt
        Nt2 = Nt
    else:
        deltat = cfl * h ** 2 / nu
        Nt = int(T / deltat)
        if Nt * deltat != T:
            Nt2 = Nt + 1
        else:
            Nt2 = Nt
    print("Le nombre de pas de temps est :",Nt2)

    cas = int(input("Quel type de donnée initiale (1, 2 ou 3) : "))

    lambda_ = nu * deltat / h ** 2

    x = np.zeros(Nx + 1)
    t = np.zeros(Nt2 + 1)

    for i in range(Nx + 1):
        x[i] = i * h
    for j in range(Nt + 1):
        t[j] = j * deltat
    if Nt2 != Nt:
        t[Nt2] = T - t[Nt2]

    u = np.zeros((Nx + 1, Nt2 + 1))
    
    for i in range(Nx + 1):
        if cas == 1:
            u[i, 0] = np.exp(-5 * (x[i] - L / 2) ** 2)
        elif cas == 2:
            if L / 2 - 1 <= x[i] <= L / 2 + 1:
                u[i, 0] = 1
        else:
            u[i, 0] = np.sin(np.pi * x[i] / L) + np.sin(10 * np.pi * x[i] / L)

    A = np.zeros((Nx - 1, Nx - 1))
    for i in range(Nx - 1):
        A[i, i] = 1 + 2 * lambda_
        if i > 0:
            A[i, i - 1] = -lambda_
        if i < Nx - 2:
            A[i, i + 1] = -lambda_

    A_inv = np.linalg.inv(A)

    for n in range(1, Nt2 + 1):
        u[1:Nx, n] = np.dot(A_inv, u[1:Nx, n - 1])
        
    print('x = ',x)
    print('t = ',t)
    print('A = ',A)
    print('A^(-1) = ',A_inv)
    print('u = ',u)
    
    plt.figure()
    nn2 = 0
    for nn in range(10, 0, -1):
        nn2 = nn2 + int(np.ceil((Nt2 + 1) / 2**nn))
        if nn2 > Nt2 + 1:
            break
        plt.plot(x, u[:, nn2], 'b')
    if nn2 != Nt2 + 1:
        plt.plot(x, u[:, Nt2], 'r')
    plt.plot(x, u[:, 0], 'k')

    plt.show()

chaleur1D_imp(0)  # Vous pouvez appeler la fonction avec les arguments de votre choix
