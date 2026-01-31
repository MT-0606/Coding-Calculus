# CHAPTER 2: LIMITS

# STEP 1: Declare variables.
import sympy as sym
x = sym.symbols('x')

class Limits:
    def __init__(self,a,L):
        # a: value x approaches
        # L: value f(x) approaches
        self.a = a
        self.L = L

# STEP 2: Define the properties of limits.
class LimitLaws(Limits):
    def sum(self,f,g): # lim(f+g) = lim(f)+lim(g)
        return sym.limit(f,x,self.a) + sym.limit(g,x,self.a)
    def difference(self,f,g): # lim(f-g) = lim(f)-lim(g)
        return sym.limit(f,x,self.a) - sym.limit(g,x,self.a)
    def constant(self,c,f): # lim(c*f) = c*lim(f)
        return c*sym.limit(f,x,self.a)
    def product(self,f,g): # lim(f*g) = lim(f)*lim(g)
        return sym.limit(f,x,self.a) * sym.limit(g,x,self.a)
    def quotient(self,f,g): # lim(f/g) = lim(f)/lim(g)
        return sym.limit(f,x,self.a) / sym.limit(g,x,self.a)

    def power(self,f,n): # lim(f^n) = {lim(f)}^n
        return (sym.limit(f,x,self.a))**n
    def root(self,f,n): # lim(f^(1/n)) = {lim(f)}^(1/n)
        return (sym.limit(f,x,self.a))**(1/n)

    def sandwich(self,f,g,h): # i.e. the Sandwich Theorem
        # If f <= g <= h when x is near a and lim(f) = lim(g) = L.
        # then lim(g) = L.
        flag1, flag2 = False, False
        if (sym.limit(f,x,self.a) <= (sym.limit(g,x,self.a)) and (sym.limit(g,x,self.a) <= (sym.limit(h,x,self.a)):
            flag1 = True
        if sym.limit(f,x,self.a) == sym.limit(h,x,self.a):
            flag2 = True
        if (flag1 == True) and (flag2 == True)
            return sym.limit(g,x,self.a)

# STEP 3: Graph limits to determine continuity.
import matplotlib.pyplot as plt
import numpy as np

# STEP 4: Check whether f(x) is continuous at certain points.
class Continuity(Limits):
    # If lim{x->a+}(f) = f(a), f is continuous from the right.
    def rightHand(self,f):
        flag = False
        if sym.limit(f,x,self.a,dir='+') == f(a):
            flag = True
        return flag
    # If lim{x->a-}(f) = f(a), f is continuous from the left.
    def leftHand(self,f):
        flag = False
        if sym.limit(f,x,self.a,dir='-') == f(a):
            flag = True
        return flag        

# STEP 5: Check for asymptotes.
class LimitsAtInfinity(Limits):
    # lim{x->+/-oo}(f) = L means f(x) has a horizontal asymptote.
    def horizontal(f):
        flag = False
        if (sym.limit(f,x,-sym.oo) == self.L) or (sym.limit(f,x,sym.oo) == self.L):
            flag = True
        return flag
    # lim{x->a}(f) = +/-oo means f(x) has a vertical asymptote.
    def vertical(f):
        flag = False
        if (sym.limit(f,x,self.a) == -sym.oo) or (sym.limit(f,x,self.a) == sym.oo):
            flag = True
        return flag
