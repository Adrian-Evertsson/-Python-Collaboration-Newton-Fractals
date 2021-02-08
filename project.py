import scipy as sp
import numpy as np
import pylab as pl
from sympy import symbols, diff
from scipy.linalg import solve
from sympy import *

def f(x):
    return (x[0]**3)-(3*x[0]*x[1]**2)-1
def g(x):
    return (3*x[0]**2*x[1]-(x[1]**3))

#def f(x):
#    return x[0]**8-28*x[0]**6*x[1]**2+70*x[0]**4*x[1]**4+15*x[0]**4-28*x[0]**2*x[1]**6-90*x[0]**2*x[1]**2+x[1]**8+15*x[1]**4-16
#def g(x):
#    return 8*x[0]**7*x[1]-56*x[0]**5*x[1]**3+56*x[0]**3*x[1]**5+60*x[0]**3*x[1]-8*x[0]*x[1]**7-60*x[0]*x[1]**3

class fractal2D(object):
    def __init__ (self, f, g):
        self.f = f
        self.g = g
        self.tol= 1.e-9
        self.listempty=True
        self.xz=[]
#        self.xinitial=xinitial
    def __call__(self, x):
        return f(x), g(x)
    
    def __repr__(self):
        return ("({},{})".format(self.f,self.g))
    
    def __Jacobian__ (self, x):
        f, g = self.f, self.g 
        fvec = np.array([f,g])
        xold=x
        
        xl = symbols('x0 x1')
        #Mw=np.array([[diff(f(xl), x0),diff(f(xl), x1)],
        #     [diff(g(xl), x0),diff(g(xl), x1)]])
        J = np.zeros([len(fvec),len(x)])

        for i in range(2):
            for j in range(2):
                Jsym = diff(fvec[i](xl),xl[j])
                J[i,j] = Jsym.subs([(xl[0],x[0]),(xl[1],x[1])])
        return J
    
    def __Newton__(self,x0):
        x = np.transpose(x0)
        xt = np.transpose(x)
        f,g,tol= self.f,self.g,self.tol

        for i in range(10000):
            prev=xt #row vec
            fvec = np.array([f(xt),g(xt)])
            J = self.__Jacobian__(prev)
            xt= xt - np.linalg.solve(J,fvec) # col vec
            if abs(xt-prev).all() < tol:
                return xt 
        else:
            return "No conv detected"
            
    def __getzeroes__(self, xintial):
        if self.listempty==True:
            self.xz.append(self.__Newton__(xintial[0]))
            self.listempty=False
        for i in range(0,len(xintial)):
            N = self.__Newton__(xintial[i])
            if type(N)== str:
                self.xz.append(N)
            else:
                C = abs(self.xz-N) > self.tol
                for i in range(len(C)):
                    if C[i].any() != True:
                        break
                else:
                    self.xz.append(N)
        return self.xz
    def __plot__(self, N, a,b,c,d, S):
        xvalues = np.linspace(a,b,N)
        yvalues = np.linspace(c,d,N)
        X, Y = np.meshgrid(xvalues, yvalues)
        P=np.zeros((N,N), dtype='float,float')
        for i in range(N):
            for j in range(N):
                    P[i,j]=(X[i,j],Y[i,j])

#        P=np.fliplr(P)
        P=np.transpose(P)
#        A=np.zeros(np.shape(P), dtype='float,float')
        A=np.zeros(np.shape(P))
        B=np.zeros(np.shape(P))
        for j in range(N):   
            for k in range(N):
                if S != True:
                    current=(self.__Newton__(tuple(P[j,k])))
                else:
                    current=(self.__SimplifiedNewton__(tuple(P[j,k])))
#                    current=(self.__SNM__(tuple(P[j,k])))
                if isinstance(current, (str)):
                    raise TypeError(current)
                A[j,k]=current[0]
                B[j,k]=current[1]
                print(k)
            print(j)
        pl.pcolor(A)
        pl.pcolor(B)
        return A,B
    
    def __SimplifiedNewton__(self,x0):
        x = np.transpose(x0)
        xt = np.transpose(x)
        f,g,tol= self.f,self.g,self.tol
        J = self.__Jacobian__(xt)
        for i in range(1000000):
            prev=xt #row vec
            fvec = np.array([f(xt),g(xt)])
            xt= xt - np.linalg.solve(J,fvec) # col vec
            if abs(xt-prev).all() < tol:
                return xt 
        else:
            return "No conv detected"
    
    def __SNM__(self,x0):
        x = np.transpose(x0)
        xt = np.transpose(x)
        f,g,tol= self.f,self.g,self.tol
        J = self.__Jacobian__(xt)
        for i in range(1000000):
            prev=xt #row vec
            fvec = np.array([f(xt),g(xt)])
            xt= xt - np.linalg.solve(J,fvec) # col vec
            if abs(xt-prev).all() < tol:
                return xt
        else:
            return "No conv detected"
        
    def __itPlot__(self,N,a,b):
        xvalues = np.linspace(a,b,N)
        values= []
        for i in range(len(xvalues)):
            values.append(np.array([[xvalues[i]],[xvalues[i]]]))
        A=[]
        for i in range(len(xvalues)):
            A.append(self.__Newton__(values[i])[1])
#            B.append(self.__inpoint__(xvalues[i],True)[1])
        pl.plot(xvalues,A,'.')
        pl.xlabel("Starting Value")
        pl.ylabel("Iteration needed")
        pl.show()
        return A



      
I=fractal2D(f, g)
#print(I.__Jacobian__(np.array([2,3]))) #WORKS
#define x0 as list of row vectors containing initial points
A=tuple((5.1,6.5))

T = I.__Newton__(A)
#T = I.__SNM__([[5],[6]])
#
#for i in range(-1000,1000):
#    if i==0:
#        i+=1

As=I.__itPlot__(250,-50,50,-50,50, True)

#print(T)
