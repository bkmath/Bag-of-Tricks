#import linalg package of the SciPy module for the LU decomp 
import scipy.linalg as linalg 
import numpy as np


#define A same as before 
A = np.array([[2., 1., 1.], [1., 3., 2.], [1., 0., 0.]]) 

#define B 
B = np.array([4., 5., 6.]) 

#call the lu_factor function 
LU = linalg.lu_factor(A) 

#solve given LU and B 
x = linalg.lu_solve(LU, B) 
x 

#we get the same solution as before 
#Solutions: [ 6. 15. -23.] 

#now we want to see how A has been factorized, P is the so called Permutation matrix 
P, L, U = linalg.lu(A) 

#print P 
#[[ 1. 0. 0.] [ 0. 1. 0.] [ 0. 0. 1.]] 
#
#print L 
#[[ 1. 0. 0. ] [ 0.5 1. 0. ] [ 0.5 -0.2 1. ]] 
#
#print U 
#[[ 2. 1. 1. ] [ 0. 2.5 1.5] [ 0. 0. -0.2]]