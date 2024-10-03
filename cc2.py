import numpy as np
import matplotlib.pyplot as plt
def appDF_exp(n_p,cfl):
    c =1
    v =0.1
    T = 0.5
    L = 2
    #n_p = 50
    dx = (L-0)/(n_p)
    dt = cfl*(dx**2)/(c*dx+2*v)
    Nt = int(T/dt)
    if (Nt-1)*dt< T:
        dt2 = T-((Nt-1)*dt)
        Nt2 = Nt+1
    else:
        Nt2 =Nt
    print("le nombre de pas de temps est :",Nt2)    
    x = np.zeros(n_p+1)
    for i in range(n_p+1):
        x[i] = i*dx
    t = np.zeros(Nt2+1) 
    for j in range(Nt+1):
        t[j] = j*dt
    if Nt2>Nt:
        t[Nt2] = t[Nt]+ dt2
    u = np.zeros((n_p+1,Nt2+1))
    for i in range(n_p+1):
        u[i,0] = np.sin(2*np.pi*x[i]) 
    lamda1 = (c*dt)/dx
    lamda2 = (v*dt)/dx**2    
    for j in range(1,Nt+1):
        for i in range(1,n_p):
            u[i,j]= u[i,j-1] -lamda1*(u[i+1,j-1]-u[i,j-1]) +lamda2*(u[i-1,j-1]-2*u[i,j-1]+u[i+1,j-1])
    if Nt2>Nt:
        for i in range(1,n_p):
            u[i,Nt2]= u[i,Nt] -lamda1*(u[i+1,Nt]-u[i,Nt]) +lamda2*(u[i-1,Nt]-2*u[i,Nt]+u[i+1,Nt])        
    u_ex = np.zeros((n_p+1,Nt2+1))
    for j in range(Nt2+1):
        for i in range(n_p+1): 
            u_ex[i,j] = np.exp(-4*np.pi**2*v*t[j])*np.sin(2*np.pi*(x[i]-c*t[j])) 
    #calcul de l'erreur l_infini,l2, l1 à t = 0.5:         
    #err_l_infini = np.zeros(n_p+1)
    #err_l_2 = np.zeros(n_p+1)
    #err_l_1 = np.zeros(n_p+1)
    #for i in range(n_p+1) :
        #err_l_infini = np.max(np.abs(u[i,Nt2]-u_ex[i,Nt2]))   
        #err_l_2 = np.sqrt(u[i,Nt2]**2 + u_ex[i,Nt2]**2) 
        #err_l_1 = np.abs(u[i,Nt2]-u_ex[i,Nt2])

    err_l_infini = np.max(np.abs(u[:,Nt2] - u_ex[:,Nt2]))   
    #err_l_2 = np.sqrt(u[:,Nt2]**2 + u_ex[:,Nt2]**2) 
    #err_l_1 = np.abs(u[:,Nt2] - u_ex[:,Nt2])
    err_l_2 = np.sqrt(np.sum((u[:, Nt2] - u_ex[:, Nt2])**2) * dx)
    err_l_1 = np.sum(np.abs(u[:, Nt2] - u_ex[:, Nt2])) * dx
    print("l'erreur infini est :",err_l_infini)
    print("l'erreur L2 est :",err_l_2)
    print("l'erreur L1 est :",err_l_1)
    #for j in range(Nt2+1):
        #for i in range(n_p):
    plt.plot(x, u[:,Nt2],'r',x,u_ex[:,Nt2],'b')
    plt.xlabel('x')
    plt.ylabel('u')
    plt.title('graphe de comparison')
    plt.legend(['u_approché','u_exacte'])
    plt.show() 
    print("err = ",err_l_infini)
    return err_l_infini

#appDF_exp(20,0.9)
#appDF_exp(50,1)
vect_n_p = [25,50,100,200,400]
vect_err =[]
for n_p in vect_n_p:
    err = appDF_exp(n_p,0.9)
    vect_err.append(err)
vect_dx = [(2-0)/(n_p) for n_p in vect_n_p]
print("levecteur_dx =",vect_dx) 
print("le vecteur erreur =",vect_err)
coefficients = np.polyfit(np.log(vect_dx), np.log(vect_err),1)
pente = coefficients[0]

plt.loglog(vect_dx,vect_err,'ro-') 
plt.xlabel('dx')
plt.ylabel('erreur')
plt.title("graphe loglog de convergence à t = 0.5 ")
plt.legend(['erreur L_infini'])
plt.show()
print("l'ordre de convergence ou la pente de graphe loglog est :",pente)

      


                