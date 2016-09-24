from __future__ import division
from sympy import symbols,Eq,init_printing,pprint,solve,simplify,Function,diff,cancel,apart,powsimp,collect,factor,fraction,integrate,Sum,oo
from sympy.functions import exp
from sympy.core.containers import Tuple


x1 = symbols(r'x1')
x2 = symbols(r'x2')
x3 = symbols(r'x3')

r1 = Function(r'r1')
r2 = Function(r'r2')

i = symbols(r'i',integer=True)
t = symbols(r't')
n = symbols(r'n',integer=True)
M = symbols(r'M',integer=True)
f = Function(r'f')
w = Function(r'w')
c = Function(r'c')


#p = Function(r'p')
#P = Function(r'P')
def P(x1,x2):
    return Sum(w(i)*exp(-((r1(i)-x1)**2 + (r2(i)-x2)**2)/c(i)**2),(i,1,M))
def p(x1):
    return Sum(w(i)*exp(-(r1(i)-x1)**2/c(i)**2),(i,1,M))



Eq1 = Eq(diff(p(x1),t),diff(diff(p(x1),x1) + (n-1)*integrate(P(x1,x2)*f(x1,x2),x2), x1))
#Eq2 = Eq(diff(P(x1,x2),t),diff(diff(P(x1,x2),x1) + f(x1,x2)*P(x1,x2) +  (n-2)*(P(x1,x2)/(p(x1)*p(x2)))*integrate((P(x1,x3)*P(x2,x3)/p(x3))*f(x1,x3),x3), x1) + diff(diff(P(x1,x2),x2) + f(x2,x1)*P(x1,x2) +  (n-2)*(P(x1,x2)/(p(x1)*p(x2)))*integrate((P(x1,x3)*P(x2,x3)/p(x3))*f(x2,x3),x3), x2))


init_printing()

print '------------------------------------------------'
print '       EQUATIONS'
print '------------------------------------------------'
pprint(Eq1)
#pprint(Eq2)


