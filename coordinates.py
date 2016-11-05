from __future__ import division
from sympy import *
from sympy.functions import exp
from sympy.core.containers import Tuple


x1 = symbols(r'x1')
x2 = symbols(r'x2')
x3 = symbols(r'x3')

x = symbols(r'x')
t = symbols(r't')
r = symbols(r'r')


xv = Matrix([x1,x2,x3])
ix = Matrix([1,1,1])/sqrt(3)
i2 = Matrix([-1,1,0])/sqrt(2)
i3 = ix.cross(i2)
#i3 = Matrix([0,0,1])
xv = x*ix + r*cos(t)*i2 + r*sin(t)*i3


pprint(Eq(x1,simplify(xv[0])))
pprint(Eq(x2,simplify(xv[1])))
pprint(Eq(x3,simplify(xv[2])))
