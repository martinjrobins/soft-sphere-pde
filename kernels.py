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
e = symbols(r'\epsilon',positive=true)

ri = Matrix([xi,yi,zi])
rj = Matrix([xj,yj,zj])
dx = ri-rj

k = exp(-dx.dot(dx)*cj)
#k = 1/sqrt(dx.dot(dx) + cj)
#k = sqrt(dx.dot(dx) + cj)
f = diff(exp(-abs(xi-yi)*e),xi)
fxigtyi = diff(exp(-(xi-yi)*e),xi)
fxiltyi = diff(exp((xi-yi)*e),xi)

def printeq(name,eq):
    eq = simplify(eq)
    pprint(eq)
    with open('K_equations.tex','a') as fil:
        fil.write('------------------------------------------------\n')
        fil.write('   %s\n'%name)
        fil.write('------------------------------------------------\n')
        fil.write('\n')
        fil.write(latex(eq)+'\n')
        fil.write('\n')



with open('K_equations.tex','w') as fil:
    fil.write('equations for matrix calculations\n')
    fil.write('------------------------------------------------\n')
    fil.write('   k\n')
    fil.write('------------------------------------------------\n')
    fil.write('\n')
    fil.write(latex(k)+'\n')
    fil.write('\n')
    fil.write('------------------------------------------------\n')
    fil.write('   f\n')
    fil.write('------------------------------------------------\n')
    fil.write('\n')
    fil.write(latex(f)+'\n')
    fil.write('\n')

print '----------------------------------------------'
print '               K1'
print '----------------------------------------------'
printeq('K1',integrate(diff(fxigtyi*(k.subs(zi-zj,0)),xi),(yi,-oo,xi)) + integrate(diff(fxiltyi*(k.subs(zi-zj,0)),xi),(yi,xi,oo)))

print '----------------------------------------------'
print '               K2'
print '----------------------------------------------'
printeq('K2',diff(diff((k.subs(zi-zj,0).subs(yi-yj,0)),xi),xi))

print '----------------------------------------------'
print '               K3'
print '----------------------------------------------'
printeq('K3',simplify(k.subs(zi-zj,0)))

print '----------------------------------------------'
print '               K4'
print '----------------------------------------------'
printeq('K4',simplify(k.subs(zi-zj,0).subs(yi-yj,0)))

print '----------------------------------------------'
print '               K5'
print '----------------------------------------------'
printeq('K5',simplify(k.subs(zi-zj,0).subs(yi-yj,0).subs(xi-xj,yi-xj)))

print '----------------------------------------------'
print '               K6'
print '----------------------------------------------'
printeq('K6',integrate(diff((fxigtyi.subs(yi,zi))*k,xi),(zi,-oo,xi))+integrate(diff((fxiltyi.subs(yi,zi))*k,xi),(zi,xi,oo)))

print '----------------------------------------------'
print '               K7'
print '----------------------------------------------'
printeq('K7',integrate(diff((f.subs(yi,zi).subs(xi,yi))*k,xi),(zi,-oo,oo)))

print '----------------------------------------------'
print '               K8'
print '----------------------------------------------'
printeq('K8',simplify(diff(k.subs(zi-zj,0).subs(yi-yj,0).subs(xi,yi),yi)))

print '----------------------------------------------'
print '               K18'
print '----------------------------------------------'
printeq('K18',simplify(diff(k.subs(zi-zj,0).subs(yi-yj,0),xi)))

print '----------------------------------------------'
print '               K9'
print '----------------------------------------------'
printeq('K9',integrate((f.subs(yi,zi).subs(xi,yi))*k,(zi,-oo,oo)))

print '----------------------------------------------'
print '               K10'
print '----------------------------------------------'
printeq('K10',integrate((f.subs(yi,zi))*k,(zi,-oo,oo)))

print '----------------------------------------------'
print '               K11'
print '----------------------------------------------'
printeq('K11',diff(k.subs(zi-zj,0),xi))

print '----------------------------------------------'
print '               K12'
print '----------------------------------------------'
printeq('K12',diff(k.subs(zi-zj,0),yi))

print '----------------------------------------------'
print '               K13'
print '----------------------------------------------'
printeq('K13',diff(f*k.subs(zi-zj,0),yi) + diff((f.subs(xi,zi).subs(yi,xi).subs(zi,yi))*k.subs(zi-zj,0),yi) +  diff(diff(k.subs(zi-zj,0),xi),xi) + diff(diff(k.subs(zi-zj,0),yi),yi))

print '----------------------------------------------'
print '               K14'
print '----------------------------------------------'
printeq('K14',k)

print '----------------------------------------------'
print '               K15'
print '----------------------------------------------'
printeq('K15',k.subs(zi-zj,0).subs(yi,zi))

print '----------------------------------------------'
print '               K16'
print '----------------------------------------------'
printeq('K16',k.subs(zi-zj,0).subs(yi,zi).subs(xi,yi))

print '----------------------------------------------'
print '               K17'
print '----------------------------------------------'
printeq('K17',k.subs(zi-zj,0).subs(yi-yj,0).subs(xi,zi))



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

