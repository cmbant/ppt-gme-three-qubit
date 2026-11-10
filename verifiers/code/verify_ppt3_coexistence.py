#!/usr/bin/env python3
"""Exact PG-2 and PG-4 checks for phase_coexistence.tex (paper Section "One support, both phases").

PG-2: the 13-term cyclic/parity biseparable decomposition sigma_BS on the rank-five support R_1, its strict PPT on all
three cuts by exact principal minors, and the same-support witness values.
PG-4: the exposed biseparable mixture with tr(H_5 X) = 0, its strict PPT by shifted Sylvester minors, and the exact
phase line -1/200 <= lambda <= 1 (Hilbert-Schmidt bound).

Uses SymPy only. The research-note checks PG-1, PG-5, PG-6, PG-7 are in verify_ppt3_support_geometry.py.
"""
from __future__ import annotations
import sympy as sp

I, sqrt3, R = sp.I, sp.sqrt(3), sp.Rational
w = (-1 + I*sqrt3)/2
s, z = sp.symbols('s z', positive=True)

def q(z0, s0=1):
    return sp.Matrix([2*z0, -z0*z0*(z0*z0+3*s0), z0*z0-s0,
                      -z0*(z0*z0+s0), z0*z0+s0])


B=sp.zeros(8,5)
for j,val in ((0,1),(3,-1),(5,-1),(6,-1)): B[j,0]=val
B[7,1]=1
for j in (1,2,4): B[j,2]=1
for j,val in ((3,1),(5,w),(6,w**2)): B[j,3]=val
for j,val in ((4,1),(2,w),(1,w**2)): B[j,4]=val
T=sp.diag(4,1,3,3,3)
assert (B.conjugate().T*B-T).applyfunc(sp.simplify)==sp.zeros(5)
# Cyclic/parity twirls preserve entries with matching character labels.
cyclic=(0,0,0,2,2)
parity=(1,-1,-1,1,-1)
def orbit(z0):
    v=q(z0)
    norm=sp.simplify((v.conjugate().T*T*v)[0])
    M=v*v.conjugate().T/norm
    return sp.Matrix(5,5,lambda j,k:
        sp.simplify(M[j,k]) if cyclic[j]==cyclic[k] and parity[j]==parity[k] else 0)

assert (q(2).conjugate().T*T*q(2))[0]==1250
assert (q(2*I).conjugate().T*T*q(2*I))[0]==290
Ev=sp.zeros(5);Ev[1,1]=1
X=(Ev+(orbit(2)+orbit(-2))/2+(orbit(2*I)+orbit(-2*I))/2)/3
M=sp.Matrix([[1232,0,0,0,0], [0,30493,32,0,0],
             [0,32,1693,0,0], [0,0,0,3700,0], [0,0,0,0,925]])
assert (54375*X-M).applyfunc(sp.simplify)==sp.zeros(5)
assert all(M[:k,:k].det()>0 for k in range(1,6))
stateN=(B*M*B.conjugate().T).applyfunc(sp.simplify)
assert sp.trace(stateN)==54375
print('PG-2: PASS exact 13-term cyclic/parity biseparable decomposition and rank-five coefficient matrix')

def pt(A,cut):
    ans=sp.zeros(8)
    for x in range(8):
        for y in range(8):
            i=[(x>>j)&1 for j in (2,1,0)]
            j=[(y>>j)&1 for j in (2,1,0)]
            i[cut],j[cut]=j[cut],i[cut]
            ans[4*i[0]+2*i[1]+i[2], 4*j[0]+2*j[1]+j[2]]=A[x,y]
    return ans

expected=[
 [1232,4558400,11751097500,13766553000000],
 [2618,4698075,8087598750,83202897843750]
]
for indices,exp in zip(([0,3,5,6],[1,2,4,7]),expected):
    block=pt(stateN,0).extract(indices,indices)
    minors=[sp.simplify(block[:k,:k].det()) for k in range(1,5)]
    assert minors==exp,(minors,exp)
assert sp.simplify(pt(stateN,0).det())>0
# The cyclic map |abc> -> |cab> yields D=diag(1,1,1,w²,w²).
P=sp.zeros(8)
for i in range(8):
    a0,b0,c0=(i>>2)&1,(i>>1)&1,i&1
    P[4*c0+2*a0+b0,i]=1
D=sp.diag(1,1,1,w**2,w**2)
assert (P*B-B*D).applyfunc(sp.simplify)==sp.zeros(8,5)
Z=sp.diag(1,-1,-1,1,-1,1,1,-1)
assert Z*B==B*sp.diag(1,-1,-1,1,-1)
assert (P*stateN*P.conjugate().T-stateN).applyfunc(sp.simplify)==sp.zeros(8)
print('PG-2: PASS strict PPT on all three cuts by exact rational principal minors')

H5=sp.diag(-12,1,11,-6,14)
Xg=sp.diag(10,30,9,9,3)/133
assert sp.trace(H5*Xg)==-R(3,133)
assert sp.trace(H5*M)/54375==R(25082,54375)
lam=sp.simplify((R(25082,54375))/(R(25082,54375)+R(3,133)))
assert lam==R(3335906,3499031)
print('PG-2: PASS BS and GME same-support witnesses and exact transition upper bound')
# PG-4: an exposed BS state with H5=0 and strict PPT, including a
# rational-width interval crossing from BS to GME without changing ranks.
norms=[56,32,576]
assert [(q(z0).conjugate().T*T*q(z0))[0].simplify()
        for z0 in (1,I,I*sp.sqrt(5))] == norms
M0=sp.Matrix([[293,0,0,0,0],[0,577,-105,0,0],
              [0,-105,273,0,0],[0,0,0,188,0],[0,0,0,0,76]])
Q=lambda z0: (orbit(z0)+orbit(-z0))/2
Xpartial=Q(1)/5+Q(I)/2+R(3,10)*Q(I*sp.sqrt(5))
assert (3360*Xpartial-M0).applyfunc(sp.simplify)==sp.zeros(5)
assert sp.trace(H5*M0)==0
assert sp.trace(M0*T)==3360
assert all(sp.simplify((M0/R(3360)-sp.eye(5)/100)[:k,:k].det())>0
           for k in range(1,6))
print('PG-4: PASS exact exposed BS mixture with H5=0, full rank on R1')

PTpartial=pt((B*M0*B.conjugate().T).applyfunc(sp.simplify),0)
raw=[[293,55084,9077712,1338006528],
     [349,62244,2149584,150147072]]
shifted=[[1409,1163516,654702384,109183407616],
         [1689,1363796,97132144,7881970176]]
for indices,mins,shift in zip(([0,3,5,6],[1,2,4,7]),raw,shifted):
    block=PTpartial.extract(indices,indices)
    assert [sp.simplify(block[:k,:k].det()) for k in range(1,5)]==mins
    block2=(5*PTpartial-56*sp.eye(8)).extract(indices,indices)
    assert [sp.simplify(block2[:k,:k].det()) for k in range(1,5)]==shift
print('PG-4: PASS strict PPT with min PT eigenvalue > 1/300 by shifted Sylvester minors')

Delta=X-M0/R(3360) # X is the PG-2 BS coefficient matrix
HS2=sp.simplify(sp.trace(Delta*T*Delta*T))
assert HS2==R(575589900887,2318006250000) and HS2<R(1,4)
assert sp.trace(H5*Delta)==R(25082,54375)
# For -1/200 <= lambda <= 0, the operator and every partial transpose
# stay strictly positive by the HS bound; for 0 <= lambda <= 1
# positivity is preserved by convex combination of the two endpoints.
assert R(1,300)-R(1,400)==R(1,1200)>0
assert R(1,100)-R(1,400)==R(3,400)>0
print('PG-4: PASS exact phase line -1/200 <= lambda <= 1: BS iff lambda >= 0; GME otherwise')

print('ALL EXACT CHECKS PASS')
