import numpy as np
import matplotlib.pyplot as plt
def chaleur1D_imp():
    L = 10
    nu = 1
    T = 10
    Nx = int(input("nombre de pas d'espaces est :"))
    h = L/Nx
    Nt = int(input("nombre de pas de temps est :"))
    dt = T/Nt
    lamda = nu * dt/h**2
    x = np.zeros(Nx+1)
    t = np.zeros(Nt+1)
    for i in range(Nx+1):
        x[i] = i*h
    for j in range(Nt+1):
        t[j] = j*dt

    A = np.zeros((Nx-1, Nx-1))
    for i in range(Nx-1):
        A[i,i] = 1 + 2 * lamda
        if i > 0:
            A[i, i-1] = -lamda
        if i < Nx-2:    
            A[i, i+1] = -lamda 
    B = np.linalg.inv(A)               
    u = np.zeros((Nx+1, Nt+1))
    cas = int(input("choisir le type de la condition initiale(1, 2 ou 3) :"))
    for i in range(Nx+1):
        if cas == 1 :
            u[i,0] = np.exp(-5*(x[i]-L/2)**2)
        elif cas == 2 :
            if L/2 - 1 <= x[i] <= L/2 + 1 :
                u[i,0] = 1
        else :
            u[i,0] = np.sin(np.pi * x[i]/L) + np.sin(10 * np.pi * x[i]/L)
    for j in range(1, Nt+1):
        u[1 : Nx,j] = np.dot(B,u[1 : Nx, j-1])
    print("x = ",x)
    print("t = ",t)   
    print("A = ",A)
    print("u = ",u)
    #for i in range(Nt+1):
    plt.plot(x, u, 'r')
    plt.xlabel('x')
    plt.ylabel('u')
    plt.title('Graphe de la solution u')
    plt.legend(['u'])
    plt.show() 
chaleur1D_imp()    