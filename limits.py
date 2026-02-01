##############################
#### CHAPTER 2: LIMITS #######
##############################

# STEP 1: Declare variables.
import sympy as sym
x = sym.symbols('x')

class Limits:
    def __init__(self,a): # a: value x approaches
        self.a = a

# STEP 2: Define the properties of limits.
class LimitLaws(Limits):
##    def sum(self,f,g): # lim(f+g) = lim(f)+lim(g)
##        return sym.limit(f,x,self.a) + sym.limit(g,x,self.a)
##    def difference(self,f,g): # lim(f-g) = lim(f)-lim(g)
##        return sym.limit(f,x,self.a) - sym.limit(g,x,self.a)
##    def constant(self,c,f): # lim(c*f) = c*lim(f)
##        return c*sym.limit(f,x,self.a)
##    def product(self,f,g): # lim(f*g) = lim(f)*lim(g)
##        return sym.limit(f,x,self.a) * sym.limit(g,x,self.a)
##    def quotient(self,f,g): # lim(f/g) = lim(f)/lim(g)
##        return sym.limit(f,x,self.a) / sym.limit(g,x,self.a)

##    def power(self,f,n): # lim(f^n) = {lim(f)}^n
##        return (sym.limit(f,x,self.a))**n
##    def root(self,f,n): # lim(f^(1/n)) = {lim(f)}^(1/n)
##        return (sym.limit(f,x,self.a))**(1/n)

    def sandwich(self,f,g,h): # i.e. the Sandwich Theorem
        # If f <= g <= h when x is near a and lim(f) = lim(g) = L.
        # then lim(g) = L.
        flag1, flag2 = False, False
        if (sym.limit(f,x,self.a)<=sym.limit(g,x,self.a)) and (sym.limit(g,x,self.a)<=sym.limit(h,x,self.a)):
            flag1 = True
        if sym.limit(f,x,self.a)==sym.limit(h,x,self.a):
            flag2 = True
        if (flag1==True) and (flag2==True):
            return sym.limit(g,x,self.a)

    def evaluate(self, expr):
        return sym.limit(expr, x, self.a)
    
# STEP 3: Graph limits to determine continuity.
##import matplotlib.pyplot as plt
##import numpy as np

# STEP 4: Check whether f(x) is continuous at certain points.
class Continuity(Limits):
    # If lim{x->a}(f) = f(a), f is continuous at a.
    def check(self,f):
##        flag = False
        lim_f = sym.limit(f,x,self.a)
        f_a = f.subs(x,self.a)
        if lim_f == f_a:
            return f_a
        if value_at_a in (sym.nan, sym.zoo, sym.oo, -sym.oo):
            return None
    # If lim{x->a+}(f) = f(a), f is continuous at the right of a.
    def rightHand(self,f):
##        flag = False
        lim_f = sym.limit(f,x,self.a,dir='+')
        f_a = f.subs(x,self.a)
        if lim_f == f_a:
            return f_a
        if value_at_a in (sym.nan,sym.zoo,sym.oo,-sym.oo):
            return None
    # If lim{x->a-}(f) = f(a), f is continuous at the left of a.
    def leftHand(self,f):
##        flag = False
        lim_f = sym.limit(f,x,self.a,dir='-')
        f_a = f.subs(x,self.a)
        if lim_f == f_a:
            return f_a
        if value_at_a in (sym.nan,sym.zoo,sym.oo,-sym.oo):
            return None
    def piecewise(self,f):
##        if not isinstance(f, sym.Piecewise):
##            left_limit = sym.limit(f,x,self.a,dir='-')
##            right_limit = sym.limit(f,x,self.a,dir='+')
##            return left_limit,right_limit
        left_expr,right_expr = None,None    # for functions with conditions x<a and x>a
        exp_at_a,exp_not_at_a = None,None   # for functions with conditions x!=a and x==a
        for expr,cond in f.args:
            if (cond == (x<self.a)) or (cond == sym.Lt(x,self.a)): 
                left_expr = expr
            if (cond == (x>self.a)) or (cond == sym.Gt(x,self.a)):
                right_expr = expr
            if (cond == (x!=self.a)) or (cond == sym.Ne(x, self.a)):
                expr_at_a = expr if expr_at_a is None else expr_at_a
                expr_not_at_a = expr if expr_not_at_a is None else expr_not_at_a
        left_limit,right_limit = None,None
        limit_at_a,limit_not_at_a = None,None
        if left_expr is not None:
            left_limit = self.leftHand(left_expr)
            return left_limit
        if right_expr is not None:
            right_limit = self.rightHand(right_expr)
            return right_limit
        if expr_at_a is not None:
            limit_at_a = sym.limit(expr_at_a,x,self.a)
            return limit_at_a
        if expr_not_at_a is not None:
            limit_not_at_a = sym.limit(expr_at_a,x,self.a)
            return limit_not_at_a
##        return left_limit,right_limit    

# STEP 5: Check for asymptotes.
class LimitsAtInfinity(Limits):
    # lim{x->+/-oo}(f) = L means f(x) has a horizontal asymptote.
    def horizontal(self,f):
        if self.L is None: # i.e. if L does not exist
            return {
                "x->-oo": sym.limit(f,x,-sym.oo),
                "x->oo": sym.limit(f,x,sym.oo),
            }
        return (sym.limit(f,x,-sym.oo) == self.L) or (sym.limit(f,x,sym.oo) == self.L)
    # lim{x->a}(f) = +/-oo means f(x) has a vertical asymptote.
    def vertical(self,f):
        flag = False
        if (sym.limit(f,x,self.a) == -sym.oo) or (sym.limit(f,x,self.a) == sym.oo):
            flag = True
        return flag

### lim(x->-2){(x**3+2*x**2-1)/(5-3*x)}
##f1 = LimitLaws(-2)
##f1_num = x**3 + 2*x**2 - 1
##f1_denom = 5 - 3*x
##f1_exp = f1_num / f1_denom
##lim_f1 = f1.evaluate(sym.simplify(f1_exp))
##print(lim_f1) # Answer: -1/11
##
### lim(x->5){2*x**2-3*x+4}
##f2 = LimitLaws(5)
##f2_exp = 2*x**2 - 3*x + 4
##lim_f2 = f2.evaluate(sym.simplify(f2_exp))
##print(lim_f2)
##
### lim(x->1){(x**2-9)/(x-3)}
##f3 = LimitLaws(1)
##f3_num = x**2 - 1
##f3_denom = x - 1
##f3_exp = f3_num / f3_denom
##lim_f3 = f3.evaluate(sym.simplify(f3_exp))
##print(lim_f3) # Answer: 6

# A piecewise function
a = 4
g = Continuity(a)
g_eval = LimitLaws(a)
g_expr = sym.Piecewise(
    (8-2*x, x<a),
    ((x-4)**0.5, x>a)
    )
result = g.piecewise(g_expr)
print(result)

# Another piecewise function
a = 2
g = Continuity(a)
g_eval = LimitLaws(a)
g_expr = sym.Piecewise(
    ((x**2-x-2)/(x-2), x!=a),
    (1, x==a)
    )
result = g.piecewise(g_expr)
print(result)


