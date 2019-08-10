#                                   Construction.
#OVERALL: We want to develop a Laplacian Operator for some 2d (uniform) mesh and employ Jacobi's Simultaneous Relaxation.
#1. We need to develop the Laplacian Operator. This will be done in 2 stages:
#1. (cont.) First, we will develop the tridiagonal in the X direction. Then we will add the outriggers after expanding the matrix into an n^2 by n^2 matrix.
#2. After we have our operator, we are now free to develop Jacobi. But here's a problem: The boundary will change after every iteration unless we impose the boundary conditions at each step.
#2. (cont.) Hence, we will construct it in two parts: The boundary value function (bound_cond) and the imposing function (boundary_impose).
#3. Finally, after having an boundary imposing routine, we can start employing a for-loop that iteratively employs the Jacobi proceedure for m steps.
#Note to self: Add a error measurement function.
#########################################################################################
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
#Here we apply constants:
#Let nx and ny are the number of points in the x and y direction, if one chooses a non-uniform mesh.
#Let dx and dy denote the partition length in the x and y direction.
nx = 11
ny = 11
pi = np.pi
dx = 1/nx
dy = 1/ny

#The Operator

def Laplacian(nx,ny):   #The construction of the Laplacian operator is divided into 2 parts

    TRI = np.diag(np.ones(ny)) * (-2/(dx**2)-(2/(dy**2)))       #We construct the main diagonal which is -2u[i,j]/dx^2 -2u[i,j]dy^2
    TRI = np.asmatrix(TRI)

    temp1x = np.asmatrix((np.diag(np.ones(nx-1),1)*(1/dx**2)))
    temp2x = np.asmatrix((np.diag(np.ones(nx-1),-1)*(1/dx**2))) #We assemble the off-diagonals, to make a tridiagonal! Now we need to make our difference equation for the y direction.
    TRI = TRI + temp1x +temp2x

    PENTA = np.kron(np.identity(ny), TRI)                             #We make the full sized matrix. We employ a kronecker product.
    PENTA = np.asmatrix(PENTA)
    temp1y = np.asmatrix((np.diag(np.ones((ny-1)*(nx)),nx))*(1/dy**2))
    temp2y = np.asmatrix((np.diag(np.ones((ny-1)*(nx)),-nx))*(1/dy**2)) #We make the outriggers like we made the off-diagonal, but seperated by distance nx.
    PENTA = PENTA + temp1y + temp2y
    return PENTA

#Boundary Routines
def boundary_cond(x,y):
    return np.sin(2*pi*(x*dx+y*dy))

def boundary_impose(u):     #We're replacing each boundary's entry (first row, last row, first column, last column) with the boundary conditions. Then we flatten again (which is useful later).
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

#Inital A with 0's in the boundary.
A = np.diag(np.zeros(nx))
for i in np.arange(0,nx):
    A[i,0]=boundary_cond(i,0)
        
for i in np.arange(0,nx):
    A[i,-1] = boundary_cond(i,nx)
        
for j in np.arange(0,ny):
    A[0,j] = boundary_cond(0,j)
        
for j in np.arange(0,ny):
    A[-1,j] = boundary_cond(ny,j)
    
#The Jacobi Iteration for this particular problem:

def Jacobi(A,m): 
    D = np.diag((1/np.diag(Laplacian(nx,ny))))  #this is diag(laplacian)^-1
    J = np.identity(nx*ny)-D*Laplacian(nx,ny)   #J = I - D^(-1)*Lap
    Aold = np.asmatrix(A.flatten()).T           #Flatten A and apply over and over again.
    for i in np.arange(0,m):
        Anew = boundary_impose(J*Aold)
        Aold = Anew
    B = Anew.reshape(nx,ny)
    return B

Answer = Jacobi(A,11) #We store our answers in a matrix (to graph) by doing 11 iterations of the Jacobi method.

#Now we plot

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = y = np.linspace(0,1,nx)
X, Y = np.meshgrid(x, y)
zs = np.array([Answer[i,j] for i in np.arange(0,nx) for j in np.arange(0,ny)])
Z = zs.reshape(X.shape)

ax.plot_surface(X, Y, Z)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Energy')

plt.show()