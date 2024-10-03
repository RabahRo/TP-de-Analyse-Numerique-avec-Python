import numpy as np
import matplotlib.pyplot as plt
def CoefDF(k, xbar, x):
    x = np.array(x)
    n = len(x)
    A = np.zeros((n, n))
    B = np.zeros(n)
    h = min(x[i] - x[i - 1] for i in range(1, n))
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
    
def u(x,y):
    return np.cos(x**2 + y**2)
def lap(x,y):
    return -4*(np.sin(x**2 + y**2) + np.cos(x**2 + y**2)*(x**2 + y**2))
def laplacien2D(x0, y0):
    h = [5*10**(-1), 10**(-1), 5*10**(-2), 10**(-2), 5*10**(-3)]
    #h = [9*10**(-1), 8*10**(-1), 7*10**(-1), 6*10**(-1), 5*10**(-1)]
    n = len(h)
    k = 2
    err_t = np.zeros(n)
    err_t2 = np.zeros(n)
    ordre = np.zeros(n-1)
    ordre2 = np.zeros(n-1)
    Lu = lap(x0, y0)   
    for i in range(n):
        Lh1 = 0
        Lh1 = Lh1 + (u(x0 + h[i],y0)+u(x0, y0 + h[i])-4*u(x0,y0)+u(x0 - h[i],y0)+u(x0, y0 - h[i])) 
        Lh1 = Lh1/h[i]**k
        Lh2 = 0
        Lh2 = Lh2 + (-u(x0 + 2*h[i],y0) -u(x0, y0 + 2*h[i]) + 16*u(x0 + h[i], y0) + 16*u(x0, y0 + h[i]) -60*u(x0,y0) -u(x0 - 2*h[i],y0) -u(x0, y0 - 2*h[i]) + 16*u(x0 - h[i], y0) + 16*u(x0, y0 - h[i]))
        Lh2 = Lh2/(12*h[i]**k)     
        err_t[i] = Lh1 - Lu
        err_t2[i] = Lh2 - Lu
        ordre[i-1] = np.log(err_t[i-1]/err_t[i]) / np.log(h[i-1]/h[i])
        ordre2[i-1] = np.log(err_t2[i-1]/err_t2[i]) / np.log(h[i-1]/h[i])    
    plt.loglog(h, [h**2 for h in h] , 'b-',h,np.abs(err_t), 'r-', h, np.abs(err_t2), 'g-', h, [h**4 for h in h],'y-', h, [h**6 for h in h],'k-') 
    plt.xlabel('h') 
    plt.ylabel('erreur')
    plt.title("le graphe loglog de la précision de la méthode ")  
    plt.legend(['h^2', 'Lh', 'Lh2', 'h^4', 'h^6'])
    plt.show()
    print("l'ordre pour Lh est :", ordre)
    print("l'ordre pour Lh2 est : ", ordre2)
laplacien2D(0,0)    
#print('coef1 =',CoefDF(2, 1, [-2, 1, 4]))
#print('coef2 =',CoefDF(2, 1, [-9,-4, 1, 6, 11]))
print('coef3 =',CoefDF(2, 0, [-10,-5, 0, 5, 10]))    
#laplacien2D(0,0) 
