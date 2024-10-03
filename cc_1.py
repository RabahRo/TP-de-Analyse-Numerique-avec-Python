import numpy as np
import matplotlib.pyplot as plt
def u(x,y):
    u1 = np.cos(x**4 + y**2)
    u2 = np.sin(x**4 + y**2)
    return np.array([[u1],[u2]])
def grad_u(x,y):
    d_u1 = [-4*x**3*np.sin(x**4 + y**2), -2*y*np.sin(x**4 + y**2)]
    d_u2 = [4*x**3*np.cos(x**4 + y**2), 2*y*np.cos(x**4 + y**2)]
    return np.array([[d_u1],[d_u2]])
def v(x,y):
    return np.exp(x**4+y**4)
def grad_v(x,y):
    return ([[4*x**3*v(x,y)],[4*y**3*v(x,y)]])
def DF(x0,y0):
    h = [30*10**(-2), 32*10**(-2), 34*10**(-2), 35*10**(-2), 36*10**(-2)]
    n = len(h)
    err1_inf = np.zeros(n)
    err2_1 = np.zeros(n)
    err2_2 = np.zeros(n)
    err2_inf = np.zeros(n)
    ordre1 = np.zeros(n-1)
    ordre2 = np.zeros(n-1)
    L1 = grad_u(x0,y0)
    L2 = grad_v(x0,y0)
    Lh1 = np.zeros((1,1))
    Lh2 = np.array([[0],[0]])
    for i in range(n):
        Lh1 = (1/(2*h[i]))*np.array([[u(x0 +h[i],y0)[0,0]-u(x0-h[i],y0)[0,0],u(x0,y0+h[i])[0,0]-u(x0,y0-h[i])[0,0]],
                                    [u(x0 +h[i],y0)[1,0]-u(x0-h[i],y0)[1,0],u(x0,y0+h[i])[1,0]-u(x0,y0-h[i])[1,0]]])
        Lh2 = (1/(2*h[i]))*np.array([[v(x0 +h[i],y0)-v(x0-h[i],y0)],[v(x0,y0 +h[i])-v(x0,y0-h[i])]])
        err1_inf[i] = np.max(np.abs(Lh1-L1))
        ordre1[i-1] = np.log(err1_inf[i-1]/err1_inf[i]) / np.log(h[i-1]/h[i])
        ordre2[i-1] = np.log(err2_2[i-1]/err2_2[i]) / np.log(h[i-1]/h[i])
        err2_1[i] = np.abs(Lh2[0]-L2[0]) + np.abs(Lh2[1]-L2[1])
        err2_2[i] = np.sqrt((Lh2[0]-L2[0])**2+(Lh2[1]-L2[1])**2)
        err2_inf[i] = np.max(np.abs(Lh2-L2))
        print("errer infini =",err1_inf)
        print('ordre1 =',ordre1)
        print('ordre2 =',ordre2)
    plt.loglog(h,[h**2 for h in h], 'r', h, err1_inf,'b',h,err2_1,'k', h,err2_2,'y',h,err2_inf,'g')
    plt.xlabel('h')
    plt.ylabel('erreur')  
    plt.title('graphe loglog de comparison') 
    plt.legend(['h^2', 'err1_inf','err2_1','err2_2','err2_inf'])
    plt.show() 
DF(1,1)
