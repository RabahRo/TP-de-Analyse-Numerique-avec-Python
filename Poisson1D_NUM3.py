import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

def F(x):
    return x**2
    #return (np.sin(np.pi*x))
    #return np.exp(x)

# probleme de poisson1D avec condition de dérichlet :
# {-u"(x) =f(x) sur ]0,10[ et u(0)=u(10) =0} 
def poisson(n) :
       
    #n = le nombre de pas
    h = (10-0)/n

    x = np.zeros(n + 1)
    for i in range(n+1):
        x[i] = i*h
        
    A = 2*np.eye(n-1)-np.diag(np.ones(n-2),1)-np.diag(np.ones(n-2),-1)
    A = A*(1/h)
    print("A =",A)

    X = np.linspace(0,10,1000)
    U_ex = (-1/12)*X**4 +(1000/12)*X
    #U_ex= (1/(np.pi**2))*np.sin(np.pi*X) 
    #U_ex = -np.exp(X) + ((np.exp(10)-1)/10)*X+1  

    u_ex = np.zeros(n+1)
    for i in range(n+1):
        u_ex[i] = (-1/12)*x[i]**4 +(1000/12)*x[i]
        #u_ex[i]= (1/(np.pi**2))*np.sin(np.pi*x[i]) 
        #u_ex[i] = -np.exp(x[i]) + ((np.exp(10)-1)/10)*x[i]+1 

    B1 = np.zeros(n-1)
    for i in range(1,n):
        B1[i-1] = h*F(x[i])
    
    U = np.linalg.solve(A, B1)
    u_ap1 = np.zeros(n+1)
    u_ap1[1:n] = U
    #avec une quadrature de Python :
    B2 = np.zeros(n-1)
    for i in range(1,n):
        a = x[i-1]
        b = x[i]
        c = x[i+1]
        def phi_1(x):
            return (x-a)/h
        def phi_2(x):
            return (c-x)/h
        def g_1(x):
            return F(x)*phi_1(x)
        def g_2(x):
            return F(x)*phi_2(x)

        B2[i-1]= integrate.quad(g_1, a, b)[0] + integrate.quad(g_2 ,b,c)[0]

    
    U2 = np.linalg.solve(A, B2)
    u_ap2 = np.zeros(n+1)
    u_ap2[1:n] = U2
    #print('u_exact',u_ex)
    #print('u_approche1',u_ap1)
    #print('u_approche2',u_ap2)
    
    plt.plot(x,u_ap1,'b',X,U_ex,'r',x,u_ap2,'k')
    plt.legend(['sol_approche1', 'sol_exact','sol_approche2'])
    plt .title('graphe de comparaison')
    plt.show()
    err_t = np.zeros(n + 1)
    for i in range(n+1):
        
        err_t[i] = np.abs(u_ap2[i]-u_ex[i])
    err_l2 = h*(np.sum(np.sqrt(err_t**2)))
    err_inf = np.max(err_t) 
    print('erreur infini est :',err_inf)
    print('erreur L2 est :',err_l2)
    #return err_l2
    return err_inf 
#poisson(20)  
vect_n = [20, 30,50,100]
vect_err = []
for n in vect_n:
    err = poisson(n)
    vect_err.append(err)
    
h_vect = [(10-0)/n for n in vect_n]
print("le vecteur h =",h_vect)
print("le vecteur erreur =",vect_err)
plt.loglog(h_vect, vect_err, 'bo-', h_vect,h_vect,'ko-',h_vect,[h_vect**2 for h_vect in h_vect],'ro-',h_vect,[h_vect**3 for h_vect in h_vect],'go-') 
plt.xlabel('h')
plt.ylabel('erreur')
plt.title('graphe de convergence')
plt.legend(['err','h','h^2','h^3'])
coef = np.polyfit(np.log(h_vect), np.log(vect_err), 1)
pente = coef[0]
print("Ordre de la méthode :", coef[0])
plt.show()







        
        
   
