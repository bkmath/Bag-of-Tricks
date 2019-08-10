    #Basic Operations
import numpy as np
import scipy.linalg as linalg

v = np.matrix([-1.0,3.0]).T
A = np.matrix([[1.0,2.0],[3.0,4.0]])
B = np.matrix([[3.0,-2.0],[2.0,1.0]])
C = np.matrix([[1.0,1.5,-2.0],[2.0,1.0,-1.0],[3.0,-1.0,2.0]])

#Not necessary because we can simply do A+B, but this is to gain insight.
def m_sum(A,B):
    return np.matrix([[A[i,j]+B[i,j] for j in np.arange(len(B))] for i in np.arange(len(A))])  # the matrix A+B

sum(v[i]*v[i] for i in range(len(v)))  # the dot product v.v

#Not necessary because we can simply write A*v or A*v.T, but this is the inner mechanisms of it.
def T(A,v):
    return np.array([sum(A[i,j]*v[j] for j in np.arange(len(v))) for i in np.arange(len(A))])  # the vector A*v

#Not necessary because we can simply type A*B, but this is to gain insight.
def m_mult(A,B):
    return np.matrix([[sum(A[i,k]*B[k,j] for k in np.arange(len(B))) for j in np.arange(len(B[0]))] for i in range(len(A))])  # the matrix A*B


#Before going any further, we write a pair of short routines to display matrices and vectors in the more familiar format:
import string
def vectfmt(thevect,fmt):
    return string.join(['[']+[fmt.format(x) for x in thevect]+[']'],'')

def matfmt(themat,fmt):
    return string.join([vectfmt(x,fmt)+'\n' for x in themat],'')

#       So we can now format our matrices as:
#       >>> print vectfmt(v,'{:10.4f}')
#       [   -1.0000    3.0000]
#       >>> print matfmt(B,'{:8.4f}')
#       [  3.0000 -2.0000]
#       [  2.0000  1.0000]

#The fatal flaw in this gaussian elimination function is that if there's a zero in a pivot, then this doesn't work. Hence, we need to develop a pivot function.
def gausselim1(A):      #inputs A as an np.matrix
    for i in np.arange(len(A)):
        for j in np.arange(i+1,len(A)):
         m = A[j,i]/A[i,i]   # Ratio of (i,j) elt by (i,i) (diagonal) elt
         A[j] = [A[j,k]-m*A[i,k] for k in np.arange(len(A))]
         return A

#This is our improved gauss elimination method.
def gausselim2(A): # Basic row pivoting
   m = A.shape[0]; n = A.shape[1]
   for j in np.arange(min(n,m)):  # for each column on the main diag
      if(A[j,j]==0): # Find a non-zero pivot and swap rows
         thecolumn = [A[k,j] for k in np.arange(j,m)]
         ipivot = thecolumn.index(max(thecolumn))
         temp = A[j]; A[j] = A[ipivot]; A[ipivot] = temp
      for i in np.arange(j+1,m):
         c = A[i,j]/A[j,j]   # Ratio of (i,j) elt by (j,j) (diagonal) elt
         A[i] = [A[i,k]-c*A[j,k] for k in np.arange(n)]
   return A

#We need to make an LU decomposition method.















