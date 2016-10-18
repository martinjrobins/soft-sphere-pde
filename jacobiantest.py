from __future__ import division
from sympy import *


n = 3
i = Symbol('i',integer=True)
j = Symbol('j',integer=True)
ks = symbols('k:'+str(n)+'(:'+str(n)+')')
ls = symbols('l:'+str(n)+'(:'+str(n)+')')
ms = symbols('m:'+str(n)+'(:'+str(n)+')')
ss = symbols('s:'+str(n)+'(:'+str(n)+')')
xs = symbols('x:'+str(n))
ys = symbols('y:'+str(n))
zs = symbols('z:'+str(n))

K = Matrix(n,n,ks)
L = Matrix(n,n,ls)
M = Matrix(n,n,ms)
S = Matrix(n,n,ss)

x = Matrix(n,1,xs)
y = Matrix(n,1,ys)
z = Matrix(n,1,zs)

inv = lambda i: 1/i
inv2 = lambda i: 1/(i**2)

term = (K*x).multiply_elementwise(L*y).multiply_elementwise((M*z).applyfunc(inv)).multiply_elementwise((S*z).applyfunc(inv))

jacobian = term.jacobian(z)
#jacobian = Matrix(n,n,lambda i,j: diff(term[i],zs[j]))
jacobian_eq_sub = -(K*x).multiply_elementwise(L*y).multiply_elementwise((M*z).applyfunc(inv2)).multiply_elementwise((S*z).applyfunc(inv2))
jacobian_eq_sub_Sz = jacobian_eq_sub.multiply_elementwise(S*z)
jacobian_eq_sub_Mz = jacobian_eq_sub.multiply_elementwise(M*z)
jacobian_eq_sub_Sz_full = jacobian_eq_sub_Sz
jacobian_eq_sub_Mz_full = jacobian_eq_sub_Mz
for i in range(n-1):
    jacobian_eq_sub_Sz_full = jacobian_eq_sub_Sz_full.row_join(jacobian_eq_sub_Sz)
    jacobian_eq_sub_Mz_full = jacobian_eq_sub_Mz_full.row_join(jacobian_eq_sub_Mz)

jacobian_eq = M.multiply_elementwise(jacobian_eq_sub_Sz_full) + S.multiply_elementwise(jacobian_eq_sub_Mz_full)
#pprint(simplify(jacobian))
#pprint(simplify(jacobian_eq))
pprint(jacobian==jacobian_eq)

