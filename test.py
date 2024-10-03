import numpy as np
import matplotlib.pyplot as plt
def AppDf_exp(m, cfl):
    alfa = 0.05*(4*np.pi**2)
    beta = 0.05
    T = 40
    dx = (3-1)/m
    dt = cfl * (dx**2)/alfa
    Nt = int(T/dt)
    if (Nt-1)*dt < T :
        dt2 = T - (Nt-1)*dt
        Nt2 = Nt+1
    else:
        Nt2 = Nt
    print("le nombre de pas de temps est :",Nt2)     
    x = np.zeros(m+1)
    for i in range(m+1):
        x[i] = i*dx
    t = np.zeros(Nt2+1)
    for j in range(Nt+1):
        t[j] = j*dt
    if Nt2> Nt:
        t[Nt2] = t[Nt]+dt2 
    # la solution analytique :
    u_ex = np.zeros((m+1,Nt2+1))
    for i in range(m+1):
        for j in range(Nt2+1):
            u_ex[i,j] = np.sin(2*np.pi*x[i])*np.exp(-1*(4*np.pi**2*alfa+beta)*t[j])  
    u_ap = np.zeros((m+1,Nt2+1)) 
    lamda_1 = (alfa*dt)/dx**2
    lamda_2 = beta*dt 
    lamda_1_1 = (alfa*dt2)/dx**2
    lamda_2_2 = beta*dt2 
    for j in range(Nt2+1):
        u_ap[0,j] = u_ap[m,j]
        u_ap[m, j] = u_ap[0, j]
    for i in range(m+1):
        u_ap[i,0] = np.sin(2*np.pi*x[i]) 
    for i in range(1,m):     
        for j in range(1,Nt+1):
            u_ap[i,j] = u_ap[i,j-1]*(1-2*lamda_1-lamda_2) +lamda_1*(u_ap[i-1,j-1]+u_ap[i+1,j-1])
        if Nt2>Nt:
            u_ap[i,Nt2] = u_ap[i,Nt2-1]*(1-2*lamda_1_1-lamda_2_2) +lamda_1_1*(u_ap[i-1,Nt2-1]+u_ap[i+1,Nt2-1])  
    print('la solution analytique est :',u_ex)
    print('la solution approché est :',u_ap)
    plt.plot(x,u_ap[:,Nt],'r',x,u_ex[:,Nt2],'b')
    plt.xlabel('x')
    plt.ylabel('u')
    plt.title('graphe de comparison')
    plt.legend(['u_approché','u_exact'])
    plt.show()  
AppDf_exp(10,0.49)                       




