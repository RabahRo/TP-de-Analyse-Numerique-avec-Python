import numpy as np

def point_milieu(a, b, m, f):
    h = (b-a)/m
    x = [a + i * h for i in range(m)]
    integral = 0
    for i in range(m-1):
        integral += h*f((x[i]+ x[i+1])/2)
    return integral
def fonction(x):
    return x**2
# Déclaration de l'intervalle [a,b] et nombre de petit intervalle :
a = 0
b = 1
n = 100
# Appel la fonction point_milieu pour l'integration :
resultat = point_milieu(a, b, n, fonction)
print("le resulta de l'approximation est :",resultat)

    

