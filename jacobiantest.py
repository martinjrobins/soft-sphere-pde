from __future__ import division
from sympy import *


n = 2
i = Symbol('i',integer=True)
j = Symbol('j',integer=True)
ks = symbols('k:'+str(n)+'(:'+str(n)+')')
ls = symbols('l:'+str(n)+'(:'+str(n)+')')
ms = symbols('m:'+str(n)+'(:'+str(n)+')')
ss = symbols('s:'+str(n)+'(:'+str(n)+')')
xs = symbols('x:'+str(n))
ys = symbols('y:'+str(n))
zs = symbols('z:'+str(n))
ws = symbols('w:'+str(n))

K = Matrix(n,n,ks)
L = Matrix(n,n,ls)
M = Matrix(n,n,ms)
S = Matrix(n,n,ss)

x = Matrix(n,1,xs)
y = Matrix(n,1,ys)
z = Matrix(n,1,zs)
w = Matrix(n,1,ws)

inv = lambda i: 1/i
inv2 = lambda i: 1/(i**2)

print '-------------------------------------'

term = (K*x).multiply_elementwise(L*y).multiply_elementwise((M*z).applyfunc(inv)).multiply_elementwise((S*z).applyfunc(inv))

jacobian = term.jacobian(z)
jacobian_action = term.jacobian(z)*w
jacobian_eq_sub = -(K*x).multiply_elementwise(L*y).multiply_elementwise((M*z).applyfunc(inv2)).multiply_elementwise((S*z).applyfunc(inv2))
jacobian_action_eq = jacobian_eq_sub.multiply_elementwise(
        (S*z).multiply_elementwise(M*w) + (M*z).multiply_elementwise(S*w))

pprint(simplify(jacobian_action-jacobian_action_eq))


print '-------------------------------------'
print '     Kz.Lz.(Mx)-1'
print '-------------------------------------'

term = (K*z).multiply_elementwise(L*z).multiply_elementwise((M*x).applyfunc(inv))

jacobian = term.jacobian(z)
jacobian_action = term.jacobian(z)*w
jacobian_action_eq = ((K*w).multiply_elementwise(L*z) + (L*w).multiply_elementwise(K*z)).multiply_elementwise((M*x).applyfunc(inv))

pprint(simplify(jacobian_action-jacobian_action_eq))


jacobian = term.jacobian(x)
jacobian_action = term.jacobian(x)*w

jacobian_action_eq = -((M*w).multiply_elementwise(K*z).multiply_elementwise(L*z)).multiply_elementwise((M*x).applyfunc(inv2))

pprint(jacobian_action)
pprint(simplify(jacobian_action))
pprint(simplify(jacobian_action-jacobian_action_eq))

