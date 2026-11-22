"""Numerical check: each original solution of the equation file (Manuscript.docx) equals its final short form."""
import numpy as np
rng = np.random.default_rng(1)
th = np.linspace(-3, 3, 41)
sech = lambda u: 1/np.cosh(u); coth = lambda u: 1/np.tanh(u)
def rel(x, y): return np.max(np.abs(x - y))/np.max(np.abs(y))
ok = True
def report(name, orig, final):
    global ok
    e = rel(orig, final); ok &= e < 1e-10; print('%-34s max rel. difference %.1e' % (name, e))
for trial in range(3):
    a, b, r, w, s1, e1, i1, r0, r1 = rng.uniform(0.3, 1.5, 9)
    ch, sh = np.cosh(b*th), np.sinh(b*th)
    # L1, c<0 -> sech
    c = -rng.uniform(0.3, 2)
    L1 = 4*a*b*r*w**2*s1*e1*(sh-ch)/(c*s1**2 + 4*a**2*r**2*w**2*e1**2*(1+2*ch*sh-2*ch**2))
    th0 = np.log(2*a*r*w*e1/(s1*np.sqrt(-c)))/b
    report('L1 (c<0) = bw/sqrt(-c) sech', L1, b*w/np.sqrt(-c)*sech(b*(th-th0)))
    # L1, c>0 -> -csch
    cp = rng.uniform(0.3, 2)
    L1p = 4*a*b*r*w**2*s1*e1*(sh-ch)/(cp*s1**2 + 4*a**2*r**2*w**2*e1**2*(1+2*ch*sh-2*ch**2))
    th0p = np.log(2*a*r*w*e1/(s1*np.sqrt(cp)))/b
    m = np.abs(th-th0p) > 0.05
    report('L1 (c>0) = -bw/sqrt(c) csch', L1p[m], -b*w/np.sqrt(cp)/np.sinh(b*(th[m]-th0p)))
    # L2, c>0: A = a^2 e1 - b^2 rho1 ; A/e1 < 0 -> -tanh ; A/e1 > 0 -> -coth
    for rho1, kind in ((a**2*e1/b**2 + 0.7, 'tanh'), (a**2*e1/b**2 - 0.2*a**2*e1/b**2, 'coth')):
        A = a**2*e1 - b**2*rho1
        L2 = b*w*(a**2*e1-b**2*rho1-b**2*r**2*e1*(1+2*ch*sh-2*ch**2))/(np.sqrt(cp)*(b**2*rho1-a**2*e1-b**2*r**2*e1*(1+2*ch*sh-2*ch**2)))
        t0 = np.log(b**2*r**2*e1/abs(A))/(2*b); u = b*(th-t0); m = np.abs(u) > 0.05
        f = -np.tanh(u) if kind == 'tanh' else -coth(u)
        report('L2 (A/e1 %s 0) = -bw/sqrt(c) %s' % ('<' if kind == 'tanh' else '>', kind), (L2)[m], (b*w/np.sqrt(cp)*f)[m])
    # L3 = constant
    L3 = b*w*(a*e1-b*r0-b*r*e1*(ch-sh))/(2*np.sqrt(cp)*(b*r0-a*e1+b*r*e1*(ch-sh)))
    report('L3 = -bw/(2 sqrt(c)) (constant)', L3, -b*w/(2*np.sqrt(cp))*np.ones_like(th))
    # L4, br/a > 0 -> -coth(u/2); amplitude b s1/(2 a r0)
    L4 = b*s1*(b**2*r**2-a**2+2*a*b*r*sh)/(2*a*r0*(a**2+b**2*r**2-2*a*b*r*ch))
    t0 = np.log(b*r/a)/b; u = b*(th-t0)/2; m = np.abs(u) > 0.05
    report('L4 = -(b s1/(2a r0)) coth(u/2)', L4[m], (-b*s1/(2*a*r0)*coth(u))[m])
    # L6: A = a(i1 r1 - s1 e1), B = br(i1 r1 + s1 e1); B/A>0 -> -(i1/e1) tanh(u/2)
    r1b = s1*e1/i1 + 0.5                        # makes A > 0, so B/A > 0
    A = a*(i1*r1b-s1*e1); B = b*r*(i1*r1b+s1*e1)
    L6 = -i1*(A + B*(sh-ch))/(e1*(A + B*(ch-sh)))
    t0 = np.log(B/A)/b
    report('L6 (B/A>0) = -(i1/e1) tanh(u/2)', L6, -i1/e1*np.tanh(b*(th-t0)/2))
    r1c = 0.5*s1*e1/i1                          # makes A < 0, so B/A < 0
    A = a*(i1*r1c-s1*e1); B = b*r*(i1*r1c+s1*e1)
    L6c = -i1*(A + B*(sh-ch))/(e1*(A + B*(ch-sh)))
    t0 = np.log(abs(B/A))/b; u = b*(th-t0)/2; m = np.abs(u) > 0.05
    report('L6 (B/A<0) = -(i1/e1) coth(u/2)', L6c[m], (-i1/e1*coth(u))[m])
    # L7 (equation file, +6abr cosh): (32 a^2 i1/(3 b^2 rho1)) tanh(u/2), th0 = ln(3br/a)/b
    L7 = 32*a**2*i1*(a**2-9*b**2*r**2+6*a*b*r*sh)/(3*b**2*r1*(a**2+9*b**2*r**2+6*a*b*r*ch))
    t0 = np.log(3*b*r/a)/b
    report('L7 = (32a^2 i1/(3b^2 r1)) tanh(u/2)', L7, 32*a**2*i1/(3*b**2*r1)*np.tanh(b*(th-t0)/2))
    # L8: P = a r0 - b r1 ; brr0/P < 0 -> -(bw/2sqrt c) tanh(u/2)
    r1d = a*r0/b + 0.6                          # P < 0
    P = a*r0 - b*r1d
    L8 = -b*w*(P + b*r*r0*(ch-sh))/(2*np.sqrt(cp)*(P + b*r*r0*(sh-ch)))
    t0 = np.log(abs(b*r*r0/P))/b
    report('L8 (P<0) = -(bw/2sqrt c) tanh(u/2)', L8, -b*w/(2*np.sqrt(cp))*np.tanh(b*(th-t0)/2))
    # L9 = L2
    rho1 = a**2*e1/b**2 + 0.7
    L2 = b*w*(a**2*e1-b**2*rho1-b**2*r**2*e1*(1+2*ch*sh-2*ch**2))/(np.sqrt(cp)*(b**2*rho1-a**2*e1-b**2*r**2*e1*(1+2*ch*sh-2*ch**2)))
    L9 = b*w*(a**2*e1-b**2*rho1-b**2*r**2*e1+2*b**2*r**2*e1*ch*(ch-sh))/(np.sqrt(cp)*(b**2*rho1-a**2*e1-b**2*r**2*e1+2*b**2*r**2*e1*ch*(ch-sh)))
    report('L9 = L2', L9, L2)
print('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')
