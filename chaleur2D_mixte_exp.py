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
    for i in range(n):
        for j in range(n):
            A[i, j] = (x[j] - xbar)**(i) / np.math.factorial(i)
    B[k] = 1
    coef = np.linalg.solve(A, B)
    coef = coef * h**k
    return coef
a = CoefDF(1, 0, [0,3,6])
print('a =',a)

def chaleur_mixte_exp(Nx,Ny,cfl):
    cas = int(input('choisir la place de radiateur (1,2,3 ou 4):'))
    L = 1
    v = 1
    T = 3
    T_in = 20
    T_ext = 5
    dx = (L-0)/Nx
    dy = (L-0)/Ny
    dt = cfl*(dx**2*dy**2)/(v*(dx**2+dy**2))
    Nt = int(T/dt)
    if (Nt-1)*dt<T:
        dt2 = T-((Nt-1)*dt)
        Nt2 = Nt+1
    else:
        Nt2 = Nt
    print('le nombre de pas de temps est :',Nt2) 

    lamda_x = (v*dt)/dx**2
    lamda_y = (v*dt)/dy**2
    lamda_x2 = (v*dt2)/dx**2
    lamda_y2 = (v*dt2)/dy**2
    
    x = np.zeros(Nx+1)
    for i in range(Nx+1):
        x[i] = i*dx
    y = np.zeros(Ny+1)
    for j in range(Ny+1):
        y[j] = j*dy
    t = np.zeros(Nt2+1)
    for k in range(Nt+1):
        t[k] = k*dt 
    if Nt2>Nt:
        t[Nt2] = t[Nt] + dt2
    
    #les murs ne permettent pas de flux de chaleur, i.e. α(x, t) = 0
    u = np.zeros((Nx+1,Ny+1,Nt2+1))    
    # condition initiale :
    for j in range(Ny+1):
        for i in range(Nx+1):
            u[i,j,0] = T_in 
            if (0.4<=y[j]<=0.6):
                u[0,j,0] = T_ext      
    for k in range(1,Nt2+1):
        for i in range(1,Nx):
            for j in range(1,Ny):
                u[i,j,k] = (1 - 2*(lamda_x + lamda_y))*u[i,j,k-1] + lamda_x*(u[i+1,j,k-1] + u[i-1,j,k-1]) + lamda_y*(u[i,j+1,k-1] + u[i,j-1,k-1]) + phi_u(x[i],y[j],u[i,j,k-1], cas)*dt
                if Nt2>Nt :
                    u[i,j,Nt2] = (1 - 2*(lamda_x2 + lamda_y2))*u[i,j,Nt2-1] + lamda_x2*(u[i+1,j,Nt2-1] + u[i-1,j,Nt2-1]) + lamda_y2*(u[i,j+1,Nt2-1] + u[i,j-1,Nt2-1]) + phi_u(x[i],y[j],u[i,j,Nt2-1], cas)*dt2      
        for i in range(1,Nx):
            u[i,0,k] = (4/3)*u[i,1,k] -(1/3)*u[i,2,k]
            u[i,Ny,k] = (4/3)*u[i,Ny-1,k] -(1/3)*u[i,Ny-2,k]
            
        for j in range(1,Ny):
            if (0.4<=y[j]<=0.6):
                u[0,j,k] = T_ext
            else :         
                u[0,j,k] = (4/3)*u[1,j,k] -(1/3)*u[2,j,k]
            u[Nx,j,k] = (4/3)*u[Nx-1,j,k] -(1/3)*u[Nx-2,j,k]
        u[0,0,k] = (2/3)*(u[1,0,k] + u[0,1,k]) -(1/6)*(u[2,0,k] + u[0,2,k])
        u[Nx,0,k] = (2/3)*(u[Nx-1,0,k] + u[Nx,1,k]) -(1/6)*(u[Nx-2,0,k] + u[Nx,2,k])
        u[0,Ny,k] = (2/3)*(u[0,Ny-1,k] + u[1,Ny,k]) -(1/6)*(u[0,Ny-2,k] + u[2,Ny,k])
        u[Nx,Ny,k] = (2/3)*(u[Nx-1,Ny,k] + u[Nx,Ny-1,k]) -(1/6)*(u[Nx-2,Ny,k] + u[Nx,Ny-2,k])    
    
   
    print('u(x,y,t) =',u)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(x,y)
    #for k in range(Nt2+1):
        #ax.plot_surface(X, Y, u[:, :, k].transpose())
    ax.plot_surface(X, Y, u[:, :, Nt2].transpose())
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('u')
    plt.show()
    return u
def phi_u(x, y,u,cas):
    T_rad = 40
    f = 0
    if cas ==1 :
        x1 = 10
        x2 = 20  
    elif cas == 2:
        x1 = 0.9
        x2 = 1    
    elif cas ==3:
        x1 = 0.45
        x2 = 0.55
    else :
        x1 = 0
        x2 = 0.1
    if x1<=x<=x2 :
        if 0.4<=y<= 0.6:
            f = (T_rad- u)**3
    return f        
     
#chaleur_mixte_exp(5,5,10)
u = chaleur_mixte_exp(20,20,0.45)
print('la température maximale et minimale dans la chambre est:',[u.max(),u.min()])   
 

