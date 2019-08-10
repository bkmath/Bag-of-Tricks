
#                                   Construction
#OVERALL: For the first order ODE y'=f(t,y(t)), we wish to use Euler's method to reconstuct y(t).
# 1. We need to define f, which is easy since y'=f(t,y(t))=y (given by the question) NOTE: it is autonomous so we can kick out our t(n).
# 2. We need to define an iterative routine for Euler's method: y(n+1)=y0+dt*f(t(n),y(n)). (e_method)
# 2. (cont.) Further, because our f is autonomous, let's define this routine to start constructing itself depending on initial conditions (y0)
# 3. Graph for various y0's {-3,-2,-1,0,1,2,3}. Cool idea for future: let's start implementing random initial conditions. :) 
#########################################################################################
#                                   Debug Space
#plt.plot(I,forward_euler(y0,n),'g-')
#plt.plot(I,solution_graph,'y--', label='Exact Solution')
#Let solution_graph be the graph of the solution for y'=y(t)
#def solution(t):    #The exactly solution (for comparison)
#    z = np.exp(t)
#    return z
#solution_graph = np.array([solution(t) for t in I])
#########################################################################################


import numpy as np
import matplotlib.pyplot as plt

#                                   Variables
#Let I denote the closed interval, [a,b].
#Let n denote the amount of partitions in the interval I.
#Let dt be the length of each partition for I.

a=0
b=1
n=100
dt = (b-a)/n
I = np.linspace(a,b,n)
initials = np.arange(-3,4)


def f(t,y):         #Integral Curve function y'=f(t,y). In this example, y'=y.
    y
    return y

def e_method(t,y):         #Iteration Function (to make new estimate of y's). This is part 1 for the function below.
    z = y+dt*f(t,y)
    return z

def forward_euler(y0,n):    #This function spits out the graph of e_method with starting value y0 and n partitions.
    ys = np.empty(0)                    # ys = {}
    ys = np.append(ys,y0)               # np.append(ys,y0) -> ys={y0}
    for i in np.arange(1,n):
        ys = np.append(ys,e_method(I[i-1],ys[i-1]))     
    return ys               


plt.grid(True)
for k in initials:
    plt.plot(I, forward_euler(k,n), label='y0 = {k}'.format(k=k))
    

plt.show()