
import numpy as np

n=11
nx = 11
ny = 11
pi = np.pi
dx = 1/nx
dy = 1/ny

#Stencil Construction

def Laplacian(nx,ny):

    A = np.diag(np.ones(ny)) * (-2/(dx**2)-(2/(dy**2)))
    A = np.asmatrix(A)

    temp1y = np.asmatrix((np.diag(np.ones(ny-1),1)*(1/dy**2)))
    temp2y = np.asmatrix((np.diag(np.ones(ny-1),-1)*(1/dy**2)))
    A = A + temp1y +temp2y

    K = np.kron(np.identity(ny), A)
    K = np.asmatrix(K)
    temp1x = np.asmatrix((np.diag(np.ones((nx-1)*(ny)),ny))*(1/dx**2))
    temp2x = np.asmatrix((np.diag(np.ones((nx-1)*(ny)),-ny))*(1/dx**2))
    K = K + temp1x + temp2x
    return K

def boundary_cond(x,y):
    return np.sin(2*pi*(x*dx+y*dy))

def boundary_impose(u):
    A = u.reshape(nx,ny)
    
    for i in np.arange(0,nx):
        A[i,0]=boundary_cond(i,0)
        
    for i in np.arange(0,nx):
        A[i,-1] = boundary_cond(i,nx)
        
    for j in np.arange(0,ny):
        A[0,j] = boundary_cond(0,j)
        
    for j in np.arange(0,ny):
        A[-1,j] = boundary_cond(ny,j)
        
    return np.asmatrix(A.flatten()).T

A = np.diag(np.zeros(n))
for i in np.arange(0,nx):
    A[i,0]=boundary_cond(i,0)
        
for i in np.arange(0,nx):
    A[i,-1] = boundary_cond(i,nx)
        
for j in np.arange(0,ny):
    A[0,j] = boundary_cond(0,j)
        
for j in np.arange(0,ny):
    A[-1,j] = boundary_cond(ny,j)
    

def Jacobi(A,m): 
    D = np.diag((1/np.diag(Laplacian(nx,ny))))
    J = np.identity(nx*ny)-D*Laplacian(nx,ny)
    Aold = np.asmatrix(A.flatten()).T
    for i in np.arange(0,m):
        Anew = boundary_impose(J*Aold)
        Aold = Anew
    B = Anew.reshape(nx,ny)
    return B

