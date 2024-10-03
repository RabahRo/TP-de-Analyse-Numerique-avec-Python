import numpy as np
import matplotlib.pyplot as plt
def chaleur1D_exp():
    L = 10
    v = 1
    T = 10
    Nx = int(input("donner le nombre de pas d'espces :"))
    h = L/Nx
    CFL = float(input("donner le Coefficient de CFL :"))
    dt = CFL*(h**2)/v
    Nt = int(T/dt)
    if (Nt-1)*dt < T :
        dt2 = T- (Nt-1)*dt
        Nt2 = Nt +1
    else :
        Nt2 = Nt    
    print("nombre de pas de temps Nt =", Nt2)
    x = np.zeros(Nx +1)
    for i in range(Nx + 1):
        x[i] = i*h    
    t = np.zeros(Nt2 +1)
    for j in range(Nt +1):
        t[j] = j*dt 
    if Nt2 < Nt :
        t[Nt2] = T 
    print('t = ',t)       
    u = np.zeros((Nx+1, Nt2+1))
    cas = int(input("choisir le type de la condition initiale (1, 2 ou 3) :"))
    for i in range(Nx+1):
        if cas == 1 :
            u[i,0] = np.exp(-5*(x[i]-L/2)**2)
        elif cas == 2 :
            if L/2-1 <= x[i] <= L/2+1 : 
                u[i,0] = 1 
        else :
            u[i,0] = np.sin(np.pi*x[i]/L) + np.sin(10*np.pi*x[i]/L)   
    for j in range(1,Nt+1):
        for i in range(1,Nx):
            u[i,j] = u[i,j-1] + (v*dt/h**2)*(u[i+1,j-1] -2*u[i,j-1] + u[i-1,j-1])
    for i in range(1,Nx):
        u[i,Nt2] = u[i,Nt] + (v*dt2/h**2)*(u[i+1,Nt] -2*u[i,Nt] + u[i-1,Nt])                          
    print("u = ",u) 
    plt.plot(x, u[:,Nt2], 'b')
    plt.xlabel(x)
    plt.ylabel(u)
    plt.title('')
    plt.legend(['u'])
    plt.show()  
chaleur1D_exp()        