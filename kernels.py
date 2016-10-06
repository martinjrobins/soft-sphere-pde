from __future__ import division
from sympy import *
from sympy.functions import exp
from sympy.core.containers import Tuple


xa = symbols(r'x_a')
ya = symbols(r'y_a')
za = symbols(r'z_a')

xb = symbols(r'x_b')
yb = symbols(r'y_b')
zb = symbols(r'z_b')
cb = symbols(r'c_b',positive=true)
k = Function('k')(xa,ya,za)

r_a = Matrix([xa,ya,za])
r_b = Matrix([xb,yb,zb])
dx = r_a

kernel = exp(-dx.dot(dx)/cb**2)
#kernel = 1/sqrt(dx.dot(dx) + cb**2)

pprint(Eq(k,kernel))
for var1 in [xa,ya,za]:
    pprint(Eq(diff(k,var1),diff(kernel,var1)))

for var1 in [xa,ya,za]:
    pprint(Eq(integrate(k,var1),simplify(integrate(kernel,(var1,-oo,oo)))))



#for var1 in [xa,ya,za]:
#    for var2 in [xa,ya,za]:
#        pprint(Eq(diff(diff(k,var1),var2),diff(diff(kernel,var1),var2)))

