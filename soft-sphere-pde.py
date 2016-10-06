from __future__ import division
from sympy import symbols,Eq,init_printing,pprint,solve,simplify,Function,diff,cancel,apart,powsimp,collect,factor,fraction,integrate,Sum,oo,MatrixSymbol,Inverse,FunctionMatrix,Lambda,Transpose,cse,Integral,Derivative
from sympy.functions import exp
from sympy.core.containers import Tuple


x1 = symbols(r'x1')
x2 = symbols(r'x2')
x3 = symbols(r'x3')

r1 = Function(r'r1')
r2 = Function(r'r2')

i = symbols(r'i',integer=True)
j = symbols(r'j',integer=True)
t = symbols(r't')
n = symbols(r'n',integer=True)
M = symbols(r'M',integer=True)
f = Function(r'f')
w = Function(r'w')
c = Function(r'c')
g = Function(r'g')


p = Function(r'p')
P = Function(r'P')
#def P(x1,x2):
#    return Sum(w(i)*exp(-((r1(i)-x1)**2 + (r2(i)-x2)**2)/c(i)**2),(i,1,M))
#def p(x1):
#    return Sum(w(i)*exp(-(r1(i)-x1)**2/c(i)**2),(i,1,M))

#k = Function(r'k')
#w_p = MatrixSymbol(r'w_p',M,1)
#w_P = MatrixSymbol(r'w_P',M,1)
#
#def K1(x1):
#    return FunctionMatrix(M,M,Lambda((i,j), k(i,j,x1))))
#def K2(x1,x2):
#    return FunctionMatrix(M,M,Lambda((i,j), k(i,j,x1,x2)))
#def K3(x1,x2,x3):
#    return FunctionMatrix(M,M,Lambda((i,j), k(i,j,x1,x2,x3)))
#def F(x1,x2):
#    return FunctionMatrix(M,M,Lambda((i,j), f(i,j)))
#
#def P(x1,x2):
#    return K2(x1,x2)*w_P
#def p(x1):
#    return K1(x1)*w_p
#def diff_disc(expr,by,K):
#    return diff(K,by)*Inverse(K)*expr
#def integrate_disc(expr,by,K):
#    return integrate(K,by)*Inverse(K)*expr
#

Eq1 = Eq(diff(p(x1,t),t),Derivative(diff(p(x1),x1) + (n-1)*Integral(P(x1,x2)*f(x1,x2),x2), x1))
#Eq2 = Eq(diff(P(x1,x2),t),diff(diff(P(x1,x2),x1) + f(x1,x2)*P(x1,x2) +  (n-2)*(P(x1,x2)/(p(x1)*p(x2)))*integrate((P(x1,x3)*P(x2,x3)/p(x3))*f(x1,x3),x3), x1) + diff(diff(P(x1,x2),x2) + f(x2,x1)*P(x1,x2) +  (n-2)*(P(x1,x2)/(p(x1)*p(x2)))*integrate((P(x1,x3)*P(x2,x3)/p(x3))*f(x2,x3),x3), x2))
Eq2 = Eq(diff(P(x1,x2,t),t),Derivative(diff(P(x1,x2),x1) + f(x1,x2)*P(x1,x2) +  (n-2)*(P(x1,x2)/(p(x1)*p(x2)))*integrate(g(x1,x2,x3)*f(x1,x3),x3), x1) + Derivative(diff(P(x1,x2),x2) + f(x2,x1)*P(x1,x2) +  (n-2)*(P(x1,x2)/(p(x1)*p(x2)))*integrate(g(x1,x2,x3)*f(x2,x3),x3), x2))

Eq3 = Eq(g(x1,x2,x3),P(x1,x3)*P(x2,x3)/p(x3))
#Eq1 = diff_disc(diff_disc(p(x1),x1,K1(x1)) + (n-1)*integrate_disc(F(x1,x2)*P(x1,x2),x2,K2(x1,x2)), x1,K1(x1))
#
#Eq2 = diff_disc(diff_disc(P(x1,x2),x1,K2(x1,x2)) + F(x1,x2)*P(x1,x2) +  (n-2)*(P(x1,x2)/(p(x1).dot(p(x2))))*integrate_disc(F(x1,x3)*(P(x1,x3)*Transpose(P(x2,x3))/p(x3)),x3,K3(x1,x2,x3)),x1,K2(x1,x2)) + diff_disc(diff_disc(P(x1,x2),x2,K2(x1,x2)) + F(x2,x1)*P(x1,x2) +  (n-2)*(P(x1,x2)/(p(x1)*Transpose(p(x2))))*integrate_disc(F(x2,x3)*(P(x1,x3)*P(x2,x3)/p(x3)),x3,K2(x1,x2)), x2,K2(x1,x2))




init_printing()

print '------------------------------------------------'
print '       EQUATIONS'
print '------------------------------------------------'
pprint(Eq1)
pprint(Eq2)
pprint(Eq3)


