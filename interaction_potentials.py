from __future__ import division
from sympy import *

x, y, e = symbols('x y e',real=True)

r = abs(x-y)
u = exp(-r/e)
f = diff(u,x)
pprint(f)
pprint(diff(f,x))
