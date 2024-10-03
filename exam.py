import numpy as np
import matplotlib.pyplot as plt
def u_ex(x,t):
    alpha =0.05/(4*np.pi**2)
    beta =0.05
    return np.sin(2*np.pi*x)*np.exp(-((4*np.pi**2*alpha+beta)*t))

def appDF_exp(m,cfl,cas):
    
    alpha =0.05/(4*np.pi**2)
    beta =0.05
    T = 40
    a = 1
    b = 3
    h = (b-a)/(m)
    if cas == 1 :
        dt = cfl*(h**2)/(alpha)
    else :
        dt = cfl*(1/(3*(1+(beta/4*np.pi**2*alpha))**2))*(h**2/alpha)
    Nt = int(T/dt)
    if (Nt-1)*dt< T:
        dt2 = T-((Nt-1)*dt)
        Nt2 = Nt+1
    else:
        Nt2 =Nt
    print("le nombre de pas de temps est :",Nt2)    
    x = np.zeros(m+1)
    for i in range(m+1):
        x[i] = a + i*h
    t = np.zeros(Nt2+1) 
    for j in range(Nt+1):
        t[j] = j*dt
    if Nt2>Nt:
        t[Nt2] = t[Nt]+ dt2
    u = np.zeros((m+1,Nt2+1))
    for i in range(m+1):
        u[i,0] = np.sin(2*np.pi*x[i]) 
    for j in range(Nt2+1):
        u[0,j] = u[m,j] 
              
    u_exact = np.zeros((m+1,Nt2+1))
    for i in range(m+1):
        for j in range(Nt2+1):
            u_exact[i,j] = u_ex(x[i],t[j])
    
    lamda = (alpha*dt)/(h**2)
    lamda2 = (alpha*dt2)/(h**2) 
    for i in range(1,m):   
        for j in range(1,Nt+1):
            u[i,j]= u[i,j-1]*(1-2*lamda - beta*dt) +lamda*(u[i+1,j-1]+u[i-1,j-1]) 
    if Nt2>Nt:        
        for i in range(1,m):
            u[i,Nt2]= u[i,Nt2-1]*(1-2*lamda2 - beta*dt2) +lamda2*(u[i+1,Nt2-1]+u[i-1,Nt2-1])
    err_t = np.zeros(m+1)
    #err_l_infini = np.zeros(m+1)
    err_l_2 = np.zeros(m+1)
    err_l_1 = np.zeros(m+1)
    for i in range(m+1) :
        err_t[i] = np.abs(u[i,Nt2]-u_exact[i,Nt2])   
        #err_l_2[i] = np.sqrt((u[i,Nt2] - u_exact[i,Nt2])**2) 
        #err_l_1[i] = np.abs(u[i,Nt2]-u_exact[i,Nt2])
    plt.plot(x,u[:,0],'r+',x,u_exact[:,0],'b')
    #plt.plot(x,u[:,Nt2],'r',x,u_exact[:,Nt2],'b')
    plt.xlabel('x')
    plt.ylabel('u')
    plt.title('graphe de comparison')
    plt.legend(['u_approché','u_exacte'])
    plt.show() 
    err_infini = np.max(np.abs(err_t))
    return err_infini

appDF_exp(100,0.45,1)
appDF_exp(100,0.55,1)
vect_m = [20,40,80,160,320]
vect_err =[]
for m in vect_m:
    err = appDF_exp(m,0.45,1)
    vect_err.append(err)
vect_h = [(3-1)/(m) for m in vect_m]
print("levecteur_dx =",vect_h) 
print("le vecteur erreur =",vect_err)
coefficients = np.polyfit(np.log(vect_h), np.log(vect_err),1)
pente = coefficients[0]
plt.loglog(vect_h,vect_err,'ro-') 
plt.xlabel('dx')
plt.ylabel('erreur')
plt.title("graphe loglog de convergence à t = 0.5 ")
plt.legend(['erreur L_infini'])
plt.show()
print("l'ordre de convergence ou la pente de graphe loglog est :",pente)