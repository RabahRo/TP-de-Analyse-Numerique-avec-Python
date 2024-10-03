import numpy as np
import matplotlib.pyplot as plt
def chaleur2D_imp(M, P,cfl):
    v = 1
    T = 0.01
    dx = (1-0)/M
    dy = (1-0)/P
    #dt = (T-0)/Nt
    dt = cfl*(dx**2*dy**2)/(v*(dx**2 + dy**2))
    Nt = int(T/dt)
    if (Nt*dt<T):
        dt2 = T-Nt*dt
        Nt2 = Nt+1
    print('le nombre de pas de temps est :',Nt2)   
    x = np.zeros(M+1)
    for i in range(M+1):
        x[i] = i*dx
    print('x=',x)    
    y = np.zeros(P+1)
    for i in range(P+1):
        y[i] = i*dy  
    print('y=',y)    
    t = np.zeros(Nt2+1)
    for j in range(Nt+1):
        t[j] = j*dt
    if Nt2>Nt:
        t[Nt2] =t[Nt] + dt2      
    print('t=',t) 
    lamda_x  = (v*dt)/dx**2
    lamda_y  = (v*dt)/dy**2
    lamda_x1  = (v*dt2)/dx**2
    lamda_y1  = (v*dt2)/dy**2
    m = M-1
    p = P-1
    q = m*p 
    u = np.zeros((M+1,P+1,Nt2+1))
    for i in range(M+1):
        for j in range(P+1):
            if (0.4 <= y[j]<= 0.6 and 0.4 <= x[i]<= 0.6):
                u[i,j,0] = 1
    u1 = np.zeros((q,Nt2+1))
    for i in range(q):
        u1[i,0] = 1
    A = np.zeros((q,q))
    for i in range(1,m+1):
        for j in range(1,p+1):
            k = (j-1)*m +i-1

            A[k,k] = 1+2*(lamda_x + lamda_y)     
            if i < m:
                A[k,k+1] = -lamda_x
            if i > 1:
                A[k,k-1] = -lamda_x
            if j > 1 :
                A[k,k-m] = -lamda_y
            if (j < p) :
                A[k,k+m] = -lamda_y
    print('A=',A) 
    B = np.linalg.inv(A)
    print('inverse de A est :',B) 
    A2 = np.zeros((q,q))
    for i in range(1,m+1):
        for j in range(1,p+1):
            k = (j-1)*m +i-1
            A2[k,k] = 1+2*(lamda_x1 + lamda_y1)
            if i < m:
                A2[k,k+1] = -lamda_x1
            if i > 1:
                A2[k,k-1] = -lamda_x1
            if j > 1 :
                A2[k,k-m] = -lamda_y1
            if (j < p) :
                A2[k,k+m] = -lamda_y1 
    B2 = np.linalg.inv(A2) 
    for n in range(1, Nt+1):
                u1[:, n] = np.dot(B, u1[:, n - 1])
    if Nt2>Nt:                                
        u1[:,Nt2] =np.dot(B2,u1[:,Nt2-1])   
    print('U=',u1)
    for n in range(1,Nt2+1):
        for i in range(1,m+1):
            for j in range(1,p+1):
                k = (j-1)*m+i-1
                u[i,j,n] = u1[k,n] 
    print('u=',u) 

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(x, y)
    #for n in range(Nt2+1):
        #ax.plot_surface(X, Y, u[:,:,n])
    ax.plot_surface(X, Y, u[:,:,Nt2])
    plt.show()
    return u[:,:,Nt2]               
chaleur2D_imp(10,10,1)
u_fin = chaleur2D_imp(10,10,0.1)
print('le max de la solution final est:',u_fin.max() )

