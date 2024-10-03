import numpy as np
import matplotlib.pyplot as plt

def poisson(nint):
    a = 0
    b = 3
    alfa = -5
    beta = 3
    
    h = (b-a)/nint
    #INITIALISATION :
    x = np.zeros(nint + 1)
    A = np.zeros((nint, nint))
    F = np.zeros(nint)
    Uex = np.zeros(nint + 1)
    errt = np.zeros(nint + 1)
    # resolution de susteme :
    #print('coef1 =',CoefDF(1, 0, [0, 3, 6]))
    A[0, 0] = 2
    A[0, 1] = 1
    x[0] = a
    F[0] = alfa/h
    for i in range(1,nint-1):
        #print('coef1 =',CoefDF(2, 0, [-3, 0, 3]))
        A[i, i] = -2
        A[i, i-1] = 1    
        A[i, i+1] = 1
        x[i] = a + i*h
        F[i] = np.exp(x[i])
    A[nint-1, nint-1] = -2
    A[nint-1, nint-2] = 1
    x[nint-1] = a + (nint-1)*h
    F[nint-1] = np.exp(x[nint-1]) - beta/h**2
    x[nint] = b
    A = A/h**2 
    print(A)
    U = np.linalg.solve(A, F)
    V = np.zeros(nint + 1)
    V[:nint]= U
    V[nint] = beta
    #print('La silution aprocher U_h= ',V)

    for i in range(nint + 1):
        Uex[i] = np.exp(x[i]) - np.exp(b) + (alfa - np.exp(a)) * (x[i] - b) + beta
       
    for i in range(nint + 1):
         errt[i] = V [i] - Uex[i]
    err = np.max(np.abs(errt))
    #print('erreur maximale est : err=', err)
    plt.plot(x, V, 'b*-', x, Uex, 'r')
    plt.xlabel('x')
    plt.ylabel('U')
    plt.legend(['Solution exacte', 'Solution numérique'])
    plt.show()
    return err
poisson(5)
