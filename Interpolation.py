#To do:
#https://docs.scipy.org/doc/scipy/reference/interpolate.html
#https://docs.scipy.org/doc/numpy/reference/routines.polynomials.polynomial.html#polynomial-class
from scipy.interpolate import lagrange
import numpy as np
from numpy.polynomial.polynomial import Polynomial

x = np.linspace(-5,5,100)
y = np.array([1-(k**2) for k in x])
f = Polynomial(np.arange(0,10))
p = lagrange(x,y)