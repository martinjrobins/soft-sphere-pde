from __future__ import division
from sympy import *
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
g1_i = Function(r'g^1_i')
g2_i = Function(r'g^2_i')
g3_i = Function(r'g^3_i')


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
Eq2 = Eq(diff(P(x1,x2,t),t),Derivative(diff(P(x1,x2),x1) + f(x1,x2)*P(x1,x2) +  (n-2)*(P(x1,x2)/(p(x1)*p(x2)))*integrate(g(x1,x2,x3)*f(x1,x3),x3), x1) + Derivative(diff(P(x1,x2),x2) + f(x2,x1)*P(x1,x2) +  (n-2)*(P(x1,x2)/(p(x1)*p(x2)))*integrate(g(x1,x2,x3)*f(x2,x3),x3), x2))
Eq3 = Eq(g(x1,x2,x3),P(x1,x3)*P(x2,x3)/p(x3))



init_printing()

print '------------------------------------------------'
print '       EQUATIONS'
print '------------------------------------------------'
pprint(Eq1)
pprint(Eq2)
pprint(Eq3)
with open('equations.tex','w') as fil:
    fil.write('------------------------------------------------\n')
    fil.write('   EQUATIONS \n')
    fil.write('------------------------------------------------\n')
    fil.write('\n')
    fil.write(latex(Eq1)+'\n')
    fil.write(latex(Eq2)+'\n')
    fil.write(latex(Eq3)+'\n')
    fil.write('\n')

Eq1 = Eq(diff(p(x1,t),t),diff(diff(p(x1),x1) + (n-1)*Integral(f(x1,x2)*P(x1,x2),x2), x1))
Eq2 = Eq(diff(P(x1,x2,t),t),diff(diff(P(x1,x2),x1) + f(x1,x2)*P(x1,x2) +  (n-2)*(P(x1,x2)/(p(x1)*p(x2)))*integrate(f(x1,x3)*g(x1,x2,x3),x3), x1) + diff(diff(P(x1,x2),x2) + f(x2,x1)*P(x1,x2) +  (n-2)*(P(x1,x2)/(p(x1)*p(x2)))*integrate(f(x2,x3)*g(x1,x2,x3),x3), x2))
Eq3 = Eq(g(x1,x2,x3),P(x1,x3)*P(x2,x3)/p(x3))


print '------------------------------------------------'
print '       EQUATIONS EXPANDED'
print '------------------------------------------------'
#Eq1 = Eq1.expand()
#Eq2 = Eq2.expand()
#Eq3 = Eq3.expand()
pprint(Eq1)
pprint(Eq2)
pprint(Eq3)
with open('equations.tex','a') as fil:
    fil.write('------------------------------------------------\n')
    fil.write('   EQUATIONS EXPANDED \n')
    fil.write('------------------------------------------------\n')
    fil.write('\n')
    fil.write(latex(Eq1)+'\n')
    fil.write(latex(Eq2)+'\n')
    fil.write(latex(Eq3)+'\n')



k = Function('k')
w_P = MatrixSymbol('w_P',M,1)
w_p = MatrixSymbol('w_p',M,1)
w_g = MatrixSymbol('w_g',M,1)

w_P_0 = MatrixSymbol('w_P_0',M,1)
w_p_0 = MatrixSymbol('w_p_0',M,1)
w_g_0 = MatrixSymbol('w_g_0',M,1)

xi = Indexed('x',i)
yi = Indexed('y',i)
zi = Indexed('z',i)
Pi = Indexed('P',i)
pi = Indexed('p',i)

dt = Symbol('dt')

def as_mat(expr):
    return FunctionMatrix(M,M,Lambda((i,j),expr))

def mat(name):
    return MatrixSymbol(name,M,M)

def my_sum():
    return Sum

K1 = mat('K1')
K2 = mat('K2')
K3 = mat('K3')
K4 = mat('K4')
K5 = mat('K5')
K6 = mat('K6')
K7 = mat('K7')
K8 = mat('K8')

EqK1 = Eq(K1,FunctionMatrix(M,M,Lambda((i,j),k(x1))))
EqK2 = Eq(K2,FunctionMatrix(M,M,Lambda((i,j),k(x1,x2))))
EqK3 = Eq(K3,FunctionMatrix(M,M,Lambda((i,j),k(x1,x2,x3))))

EqK4 = Eq(K4,FunctionMatrix(M,M,Lambda((i,j),diff(diff(k(x1),x1,x1)))))
EqK5 = Eq(K5,FunctionMatrix(M,M,Lambda((i,j),diff(f(x1,x2)*k(x1,x2),x1))))

EqK7 = Eq(K7,FunctionMatrix(M,M,Lambda((i,j),diff(f(x1,x3)*k(x1,x2,x3),x1))))


Eq1 = dt*(n-1)*K5*w_P + dt*K4*w_p - K1*(w_p - w_p_0)

#Eq2 = (n-1)*(K2*w_P).multiply_elementwise(K7*w_g)

#Eq1 = dt*(n-1)*Sum((integrate(diff(f(x1,x2),x1)*k(x1,x2),x2) + integrate(f(x1,x2)*diff(k(x1,x2),x1),x2))*Pi,(i,0,M)) + Sum(diff(diff(k(x1),x1),x1)*pi,(i,0,M))



#replacements = []
#
#for a in [f(x1,x3),f(x2,x3),diff(f(x1,x3),x1),diff(f(x2,x3),x2)]:
#    replacements.append((Integral(a*g(x1,x2,x3),x3),FunctionMatrix(M,M,Lambda((i,j), integrate(a*k(x1,x2,x3),x3)))*w_g))
#    replacements.append((Integral(a*diff(g(x1,x2,x3),x1),x3),FunctionMatrix(M,M,Lambda((i,j), integrate(a*diff(k(x1,x2,x3),x1),x3)))*w_g))
#    replacements.append((Integral(a*diff(g(x1,x2,x3),x2),x3),FunctionMatrix(M,M,Lambda((i,j), integrate(a*diff(k(x1,x2,x3),x2),x3)))*w_g))
#    replacements.append((Integral(g(x1,x2,x3),x3)*a,FunctionMatrix(M,M,Lambda((i,j), integrate(a*k(x1,x2,x3),x3)))*w_g))
#
#
#for replace in replacements:
#    Eq1 = Eq1.xreplace({replace[0]:replace[1]})
#    Eq2 = Eq2.xreplace({replace[0]:replace[1]})
#    Eq3 = Eq3.xreplace({replace[0]:replace[1]})
#    #Eq1 = Eq1.xreplace({replace[0]:replace[1]})
#    #Eq2 = Eq2.xreplace(replace[0],replace[1])
#    #Eq3 = Eq3.xreplace(replace[0],replace[1])
#
#replacements = []
#
#for a in [diff(f(x1,x2),x2),diff(f(x1,x2),x1),diff(f(x2,x1),x2),diff(f(x2,x1),x1)]:
#    replacements.append((a*P(x1,x2),FunctionMatrix(M,M,Lambda((i,j), a*k(x1,x2)))*w_P))
#
#for a in [f(x1,x2),f(x2,x1)]:
#    replacements.append((a*diff(P(x1,x2),x1),FunctionMatrix(M,M,Lambda((i,j), diff(k(x1,x2),x1)))*w_P))
#    replacements.append((a*diff(P(x1,x2),x2),FunctionMatrix(M,M,Lambda((i,j), diff(k(x1,x2),x2)))*w_P))
#
#for replace in replacements:
#    Eq1 = Eq1.replace(replace[0],replace[1])
#    Eq2 = Eq2.replace(replace[0],replace[1])
#    Eq3 = Eq3.replace(replace[0],replace[1])
#
#
#for a in [x1,x2,x3]:
#    replacements.append((diff(p(a),a),FunctionMatrix(M,M,Lambda((i,j), diff(k(a),a)))*w_p))
#for a in [x1,x2]:
#    replacements.append((diff(P(x1,x2),a),FunctionMatrix(M,M,Lambda((i,j), diff(k(x1,x2),a)))*w_P))
#
#for a in [x1,x2,x3]:
#    replacements.append((p(a),FunctionMatrix(M,M,Lambda((i,j), k(a)))*w_p))
#
#replacements.append((P(x1,x2),FunctionMatrix(M,M,Lambda((i,j), k(x1,x2)))*w_P))



print '------------------------------------------------'
print '       EQUATIONS DISCRETISED'
print '------------------------------------------------'

#pprint(cse(Eq1))
#pprint(Eq2)
#pprint(Eq3)

#Eq1 = Eq(diff(p(x1,t),t),diff(diff(p(x1),x1) + (n-1)*g1_i(x1), x1))
#Eq2 = Eq(diff(P(x1,x2,t),t),diff(diff(P(x1,x2),x1) + f(x1,x2)*P(x1,x2) +  (n-2)*(P(x1,x2)/(p(x1)*p(x2)))*g2_i(x1,x2), x1) + diff(diff(P(x1,x2),x2) + f(x2,x1)*P(x1,x2) +  (n-2)*(P(x1,x2)/(p(x1)*p(x2)))*g3_i(x1,x2), x2))
#
#Eq3 = Eq(g(x1,x2,x3),P(x1,x3)*P(x2,x3)/p(x3))
#Eq4 = Eq(g1_i(x1),Integral(P(x1,x2)*f(x1,x2),x2))
#Eq5 = Eq(g2_i(x1,x2),integrate(g(x1,x2,x3)*f(x1,x3),x3))
#Eq6 = Eq(g3_i(x1,x2),integrate(g(x1,x2,x3)*f(x2,x3),x3))
#
#print '------------------------------------------------'
#print '       EQUATIONS'
#print '------------------------------------------------'
#pprint(Eq1)
#pprint(Eq2)
#pprint(Eq3)
#pprint(Eq4)
#pprint(Eq5)
#pprint(Eq6)


