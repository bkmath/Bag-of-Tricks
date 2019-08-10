#The idea: We want to develop a function so that it gives us the first n amount of nonzero terms as a polynomial estimation.
#So, we will define a routine that gives a non-zero coefficient (fn),
#Second, we will define a routine that constructs the estimated polynoimal (meta_coef) using the coefficients, (fn)
#Finally, we will evaluate each function (meta_graph) at each point x in some interval I.
#########################################################################################

import numpy as np
import matplotlib.pyplot as plt



#Let n denote the amount of terms in the taylor series
#Let k denote the number of partitions in the interval, I
#Let I = [-pi,pi]
n = 200
k = 1000
pi = np.pi
I = np.linspace((-1)*pi,pi,k)

#for convenience
def sin(x):
    y = np.sin(x)
    return y 

#nth coefficient of the Taylor series expansion of sine.
def a(m):
    y = ((-1)**(m))/(np.math.factorial(2*m+1))
    return y


#Prototype: coef = [coef(i) for i in np.arange(1,k)]
#The meta coefficients are the non-zero coefficients for the Taylor seires of sine.
def meta_coef(m):
    coef = np.empty(0)
    for i in np.arange(0,m+1):
        coef = np.append(coef,np.array([0, a(i)]))
    return coef
#test = [meta_coef(i) for i in np.arange(1,5)]

#Prototype: p = np.polynomial.polynomial.Polynomial(coef) 

#Outputs y=fn(x)  for m non-zero terms of fn.
def fn(x,m):
    y = np.polynomial.polynomial.polyval(x,meta_coef(m))
    return y

#Outputs graph of sine.
graph_Sin = np.array([sin(x) for x in I])

#Outputs graph of fn.
def meta_graph(m):
    y = np.array([fn(x,m) for x in I])
    return y

#Prototype: Set_graphs = np.array([meta_graph(i) for i in np.arange(1,n+1)])

fig = plt.figure()
plt.plot(I,graph_Sin,label='Sine Function')
plt.plot(I,meta_graph(n),label='Estimate')
plt.title('Approximation of with $k terms')
plt.grid(True)
axes = plt.gca()
axes.set_xlim([I[0],I[-1]])
axes.set_ylim([-2,2])

plt.show()