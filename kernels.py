from __future__ import division
from sympy import *
from sympy.functions import exp
from sympy.core.containers import Tuple


xi = symbols(r'x_i',real=True)
yi = symbols(r'y_i',real=True)
zi = symbols(r'z_i',real=True)

xj = symbols(r'x_j',real=True)
yj = symbols(r'y_j',real=True)
zj = symbols(r'z_j',real=True)

cj = symbols(r'c_j',positive=true)
e = symbols(r'e',positive=true)

ri = Matrix([xi,yi,zi])
rj = Matrix([xj,yj,zj])
dx = ri-rj

k = exp(-dx.dot(dx)*cj)
#kernel = 1/sqrt(dx.dot(dx) + cb**2)
f = diff(exp(-abs(xi-yi)*e),xi)

print '----------------------------------------------'
print '               K1'
print '----------------------------------------------'
#pprint(integrate(diff(f*(k.subs(zi-zj,0)),xi),yi))

print '----------------------------------------------'
print '               K2'
print '----------------------------------------------'
pprint(simplify(diff(diff((k.subs(zi-zj,0).subs(yi-yj,0)),xi),xi)))

print '----------------------------------------------'
print '               K3'
print '----------------------------------------------'
pprint(simplify(k.subs(zi-zj,0)))

print '----------------------------------------------'
print '               K4'
print '----------------------------------------------'
pprint(simplify(k.subs(zi-zj,0).subs(yi-yj,0)))

print '----------------------------------------------'
print '               K5'
print '----------------------------------------------'
pprint(simplify(k.subs(zi-zj,0).subs(yi-yj,0).subs(xi-xj,yi-xj)))

print '----------------------------------------------'
print '               K6'
print '----------------------------------------------'
#pprint(integrate(diff((f.subs(yi,zi))*k,xi),zi))

print '----------------------------------------------'
print '               K7'
print '----------------------------------------------'
#pprint(integrate(diff((f.subs(yi,zi).subs(xi,yi))*k,xi),zi))

print '----------------------------------------------'
print '               K8'
print '----------------------------------------------'
pprint(simplify(diff(k.subs(zi-zj,0).subs(yi-yj,0).subs(xi-xj,yi-xj),yi)))

print '----------------------------------------------'
print '               K9'
print '----------------------------------------------'
#pprint(integrate((f.subs(yi,zi).subs(xi,yi))*k,zi))

print '----------------------------------------------'
print '               K10'
print '----------------------------------------------'
pprint(integrate((f.subs(yi,zi))*k,zi))

print '----------------------------------------------'
print '               K11'
print '----------------------------------------------'
#

#pprint(Eq(k,kernel))
#for var1 in [xa,ya,za]:
#    pprint(Eq(diff(k,var1),diff(kernel,var1)))
#
#for var1 in [xa,ya,za]:
#    pprint(Eq(integrate(k,var1),simplify(integrate(kernel,(var1,-oo,oo)))))
#
#for var1 in [xa,ya,za]:
#    pprint(Eq(diff(diff(k,var1),var1),simplify(diff(diff(kernel,var1),var1))))




#for var1 in [xa,ya,za]:
#    for var2 in [xa,ya,za]:
#        pprint(Eq(diff(diff(k,var1),var2),diff(diff(kernel,var1),var2)))

