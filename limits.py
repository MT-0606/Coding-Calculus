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
        return sym.limit(f,x,sym.a) + sym.limit(g,x,sym.a)
    def difference(self,f,g): # lim(f-g) = lim(f)-lim(g)
        return sym.limit(f-g,x,sym.a) - sym.limit(g,x,sym.a)
    def constant(self,c,f): # lim(c*f) = c*lim(f)
        return c*sym.limit(f,x,sym.a)
    def product(self,f,g): # lim(f*g) = lim(f)*lim(g)
        return sym.limit(f,x,sym.a) * sym.limit(g,x,sym.a)
    def quotient(self,f,g): # lim(f/g) = lim(f)/lim(g)
        return sym.limit(f,x,sym.a) / sym.limit(g,x,sym.a)

    def squeeze(self,f,g):
        # If f <= g <= h when x is near a and lim(f) = lim(g) = L.
        # then lim(g) = L.

# STEP 3: Check whether f(x) is continuous at certain points.
class Continuity(Limits):
    # If lim{x->a+}(f) = f(a), f is continuous from the right.

    # If lim{x->a-}(f) = f(a), f is continuous from the left.

# STEP 4: Check for asymptotes.
class LimitsAtInfinity(Limits):
    # lim{x->+/-oo}(f) = L means f(x) has a horizontal asymptote.

    # lim{x->a}(f) = +/-oo means f(x) has a vertival asymptote.
