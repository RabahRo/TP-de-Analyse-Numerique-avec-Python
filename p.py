import numpy as np
import matplotlib.pyplot as plt

def F(x):
    return np.sin(np.pi*x)

    
n = int(input("donner le nombre de pas n= "))
h = (10-0)/n
#INITIALISATION :
u_ex = np.zeros(n + 1)
x = np.zeros(n + 1)
for i in range(n+1):
    x[i] = i*h
        
A = 2*np.eye(n-1)-np.diag(np.ones(n-2),1)-np.diag(np.ones(n-2),-1)
    
A = A/h**2
#print("A =",A)
for i in range(n+1):
    u_ex[i] = 1/(np.pi**2)*np.sin(np.pi*x[i])           
B = np.zeros(n-1)
U = np.zeros(n-1)
for i in range(n-1):
    B[i] = h*F(x[i])
U = np.linalg.solve(A, B)
u_ap = np.zeros(n+1)
u_ap[1:n] = U
plt.plot(x, u_ap, 'b*-', x, u_ex, 'r')
plt.show()