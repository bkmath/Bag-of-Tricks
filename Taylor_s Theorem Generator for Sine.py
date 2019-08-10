
#                                   Construction.
#OVERALL: We want to develop a function so that it gives us the first n amount of nonzero terms as a polynomial estimation.
#1. We will define a routine that gives us the values of the m-th non-zero coefficient (a(m)),
#2. We will define a routine that constructs the estimated polynoimal (meta_coef) using the coefficients. Because sine only has odd polynomials, we will alternate between 0 and the function a(m)
#3. Finally, we will define a routine that will spit out the graph of the polynomial defined in 2.
#4. We will graph all the approximations between 1 and 16.
#########################################################################################
#                                   Debug Space
#test = [meta_coef(i) for i in np.arange(1,5)]
#Prototype: p = np.polynomial.polynomial.Polynomial(coef) 
#Prototype: coef = [coef(i) for i in np.arange(1,k)]
#Prototype: Set_graphs = np.array([meta_graph(i) for i in np.arange(1,n+1)])
#plt.plot(I,meta_graph(n),label='Estimate')
#plt.title('Approximation of with {k} terms')
#########################################################################################


import numpy as np
import matplotlib.pyplot as plt


#                                   Variables
#Let n denote the amount of non-zero terms in the Taylor series for sine.
#Let I denote the closed interval, [-pi,pi]
#Let k denote the number of partitions in the interval I.
n = 16
k = 1000
pi = np.pi
I = np.linspace((-1)*pi,pi,k)

def sin(x):             #we simply define this for convenience
    y = np.sin(x)
    return y 

def a(m):               #this is the mth non-zero term for sine's Taylor series.
    y = ((-1)**(m))/(np.math.factorial(2*m+1))
    return y



def meta_coef(m):       #I wanted a function that can give me the coordinate vectors for the polynomial.
    coef = np.empty(0)
    for i in np.arange(0,m+1):
        coef = np.append(coef,np.array([0, a(i)]))
    return coef

def fn(x,m):            #Outputs y=fn(x)  for m non-zero terms of fn. Esentially, a dot product between the standard basis for polynomials and the coordinates for sine.
    y = np.polynomial.polynomial.polyval(x,meta_coef(m))
    return y

graph_Sin = np.array([sin(x) for x in I])    #Outputs graph of sine.


def meta_graph(m):                      #Outputs the graph of estimated sine.
    y = np.array([fn(x,m) for x in I])
    return y

#                       #Graph
fig = plt.figure()
plt.plot(I,graph_Sin, 'g--', label='Sine Function')
plt.grid(True)
axes = plt.gca()
axes.set_xlim([I[0],I[-1]])
for k in np.arange(1,n+1):
    plt.plot(I, meta_graph(k), label='{k} Term Approx.'.format(k=k))
#plt.legend(loc='best')             #I commented out the legend because it looked bad.
plt.show()