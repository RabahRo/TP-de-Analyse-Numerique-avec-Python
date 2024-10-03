import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
def chaleur2D(Nx,Ny,CFL):
    #coefficient de diffusion termique v :
    v = 1
    T = 0.01
    #Ny = Nx
    dx = (1-0)/Nx
    dy = (1-0)/Ny
    dt = CFL*(dx**2*dy**2)/(v*(dx**2 + dy**2))
    Nt = int(T/dt)
    if (Nt*dt) < T:
        dt2 = T - Nt*dt
        Nt2 = Nt +1
    else :
        Nt2 = Nt    
    print(' le nombre de pas de temps est :',Nt2)
    #dt = (T-0)/Nt
    lamda_x = v*dt/dx**2 
    lamda_y = v*dt/dy**2 
    lamda_x2 = v*dt2/dx**2 
    lamda_y2 = v*dt2/dy**2
    x = np.zeros(Nx+1)
    for i in range(Nx+1):
        x[i] = i*dx     
    y = np.zeros(Ny+1)
    for i in range(Ny+1):
        y[i] = i*dy      
    t = np.zeros(Nt2+1)
    for i in range(Nt+1):
        t[i] = i*dt
    if Nt2 > Nt:
        t[Nt2] = t[Nt] + dt2
        #t[Nt2] = T
    print('le vecteur t =',t) 
        
    u = np.zeros((Nx+1, Ny+1, Nt2+1))
    for i in range(Nx+1):
        for j in range(Ny+1):
            if (0.4 <= y[j]<= 0.6 and 0.4 <= x[i]<= 0.6):
                u[i,j,0] = 1
    for k in range(1,Nt+1):
        for i in range(1,Nx):
            for j in range(1,Ny):
                u[i,j,k] = (1 - 2*(lamda_x + lamda_y))*u[i,j,k-1] + lamda_x*(u[i+1,j,k-1] + u[i-1,j,k-1]) + lamda_y*(u[i,j+1,k-1] + u[i,j-1,k-1])
    if Nt2>Nt:
        for i in range(1,Nx):
            for j in range(1,Ny): 
                u[i,j,Nt2] = (1 - 2*(lamda_x2 + lamda_y2))*u[i,j,Nt] + lamda_x2*(u[i+1,j,Nt] + u[i-1,j,Nt]) + lamda_y2*(u[i,j+1,Nt] + u[i,j-1,Nt])            
    #la stabilité l_infini :   dt <= 1/2(dx**2*dy**2/(dx**2+dy**2))
    #comme dx = dy alors dt<=dx**2/4
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(x, y)
    # le graphe pour chaque t :
    for k in range(Nt2+1):
        ax.plot_surface(X, Y, u[:,:,k].transpose(), cmap='viridis')
    # le graphe sur un seul t_j :    
    #ax.plot_surface(X, Y, u[:,:,:].transpose(), cmap='viridis')    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('u')
    ax.set_title('Évolution de u dans le temps')

    plt.show()
            

chaleur2D(40, 40, 0.45)  
