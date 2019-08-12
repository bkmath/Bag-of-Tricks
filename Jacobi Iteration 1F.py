import numpy as np

A = np.matrix([[2, -1, 0],
               [-1, 2, -1],
               [0, -1, 2]])

D = np.diag(1/np.diag(A))
I = np.identity(3)
b = np.matrix([[1],
               [1],
               [1]])

def Jacobi(A,m):
    J = I - D*A
    x0 = np.matrix([[0],
                    [0],
                    [0]])

    xold = J*x0 + D*b
    for i in np.arange(0,m):
        xnew = J*xold +D*b
        xold = xnew
    return xnew

print(Jacobi(A,50))