"""Reduce the eight solutions of the original draft (Mamuscript-Tusar.docx, Eqs. 12-19) to closed form with SymPy.
E = exp(-b*theta). Every case is written as a rational function of E and simplified."""
import sympy as sp
a, b, r, w, c, s1, i1, r0, r1, e1, th = sp.symbols('a b r w c sigma1 iota1 rho0 rho1 epsilon1 theta', positive=True)
u = b*th; ch, sh = sp.cosh(u), sp.sinh(u)
E = sp.Symbol('E', positive=True)            # E = exp(-b theta)
def inE(expr):
    return sp.factor(sp.simplify(expr.rewrite(sp.exp).subs(sp.exp(-u), E).subs(sp.exp(u), 1/E)))
cases = {
 'Case 1': 4*a*b*r*w**2*s1*e1*(sh-ch)/(c*s1**2 + 4*a**2*r**2*w**2*e1**2*(1+2*ch*sh-2*ch**2)),
 'Case 2': b*w*(a**2*e1-b**2*r1-b**2*r**2*e1*(1+2*ch*sh-2*ch**2))/(sp.sqrt(c)*(b**2*r1-a**2*e1-b**2*r**2*e1*(1+2*ch*sh-2*ch**2))),
 'Case 3': b*w*(a*e1-b*r0-b*r*e1*(ch-sh))/(2*sp.sqrt(c)*(b*r0-a*e1+b*r*e1*(ch-sh))),
 'Case 4': b*s1*(b**2*r**2-a**2+2*a*b*r*sh)/(2*a*r0*(a**2+b**2*r**2-2*a*b*r*ch)),
 'Case 5': -i1*(a*(i1*r1-s1*e1)+b*r*(i1*r1+s1*e1)*(sh-ch))/(e1*(a*(i1*r1-s1*e1)+b*r*(i1*r1+s1*e1)*(ch-sh))),
 'Case 6 (Mamuscript-Tusar: -6abr cosh)': 32*a**2*i1*(a**2-9*b**2*r**2+6*a*b*r*sh)/(3*b**2*r1*(a**2+9*b**2*r**2-6*a*b*r*ch)),
 'Case 6 (bare-maths file: +6abr cosh)': 32*a**2*i1*(a**2-9*b**2*r**2+6*a*b*r*sh)/(3*b**2*r1*(a**2+9*b**2*r**2+6*a*b*r*ch)),
 'Case 7': -b*w*(a*r0-b*r1+b*r*r0*(ch-sh))/(2*sp.sqrt(c)*(a*r0-b*r1+b*r*r0*(sh-ch))),
 'Case 8': b*w*(a**2*e1-b**2*r1-b**2*r**2*e1+2*b**2*r**2*e1*ch*(ch-sh))/(sp.sqrt(c)*(b**2*r1-a**2*e1-b**2*r**2*e1+2*b**2*r**2*e1*ch*(ch-sh))),
}
print('identities: sinh-cosh =', sp.simplify((sh-ch).rewrite(sp.exp)), ';  1+2 ch sh - 2 ch^2 =', sp.simplify((1+2*ch*sh-2*ch**2).rewrite(sp.exp)))
red = {}
for k, v in cases.items():
    red[k] = inE(v); print(f'{k}:  {red[k]}')
print('Case 8 - Case 2 =', sp.simplify(red['Case 8'] - red['Case 2']))
print('Case 3 is constant:', sp.diff(red['Case 3'], E) == 0)
