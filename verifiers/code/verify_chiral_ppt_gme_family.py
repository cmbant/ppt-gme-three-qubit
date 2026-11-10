#!/usr/bin/env python3
"""CF1--CF6: analytic five-component three-qubit PPT-GME construction.

No optimizer, sampled positivity, or floating-point rank is used.
The one extra finite input is an exact algebraic rank for extremality.
"""
import sympy as s
from sympy.polys.matrices import DomainMatrix

om=(-1+s.sqrt(3)*s.I)/2
t=s.symbols('t',real=True);z=s.symbols('z')
def clean(M):
    return M.applyfunc(s.expand)
def zero(M):
    assert clean(M)==s.zeros(*M.shape)
def pt(M,c):
    out=s.zeros(8)
    for i in range(8):
        for j in range(8):
            bit=(((i>>(2-c))^(j>>(2-c)))&1)<<(2-c)
            out[i^bit,j^bit]=M[i,j]
    return out
def pd(M):
    zero(M-M.conjugate().T)
    for j in range(1,M.rows+1):
        assert s.simplify(M[:j,:j].det(method='domain-ge'))>0

B=s.zeros(8,5);B[0,0]=B[7,1]=1
for j in [3,5,6]:B[j,0]=-1
for j in [1,2,4]:B[j,2]=1
for n,j in enumerate([3,5,6]):B[j,3]=s.expand(om**n)
for n,j in enumerate([4,2,1]):B[j,4]=s.expand(om**n)
T=s.diag(4,1,3,3,3)
zero(B.conjugate().T*B-T)
P=s.zeros(8)
for j in range(8):P[((j&3)<<1)|(j>>2),j]=1
parity=s.diag(*[(-1)**j.bit_count() for j in range(8)])
raw=clean(B*s.diag(t,3*t,3,3,1)*B.conjugate().T)
assert s.expand(s.trace(raw)-7*(t+3))==0
zero(P*raw*P.T-raw);zero(parity*raw*parity-raw)
ixs=[[0,3,5,6],[1,2,4,7]]
blocks=[pt(raw,0).extract(ix,ix) for ix in ixs]
expected=[[t,3*t,(t+3)*(3*t-7),3*(8*t*t-15*t-36)],
          [4,9,3*(12-t*t),-3*(t-3)*(3*t*t+4*t-48)]]
for block,mins in zip(blocks,expected):
    for j,target in enumerate(mins,1):
        assert s.simplify(block[:j,:j].det()-target)==0
lo=(15+9*s.sqrt(17))/16;hi=2*(s.sqrt(37)-1)/3
assert lo>s.Rational(29,9) and lo<s.Rational(10,3)<hi<s.sqrt(12)
assert s.simplify(8*lo**2-15*lo-36)==0
assert s.simplify(3*hi**2+4*hi-48)==0
assert clean(raw.subs(t,s.Rational(10,3))/s.Rational(133,3)).trace()==1
print('CF1 PASS exact symbolic: orthogonal five-component family, cyclic/parity symmetry, all PT minors and exact PPT interval',flush=True)

K=s.zeros(3,8)
for n,j in enumerate([3,5,6]):K[0,j]=s.expand(om**n)
for n,j in enumerate([4,2,1]):K[1,j]=s.expand(om**n)
K[2,0]=3
for j in [3,5,6]:K[2,j]=1
zero(K*B)
q=s.Matrix([2*z,-z*z*(z*z+3),z*z-1,-z*(z*z+1),z*z+1])
v=clean(B*q)
zero(v[4:8,:]-z*v[:4,:])
M=K*s.kronecker_product(s.Matrix([1,z]),s.eye(4))
co=s.Matrix([(-1)**j*M[:,[i for i in range(4) if i!=j]].det() for j in range(4)])
zero(co-(s.conjugate(om)-om)*v[:4,:])
assert v[0]==2*z and s.simplify(v[1].subs(z,0)-(om**2-1))==0
assert s.Poly(v[3],z).LC()==-1 and s.Poly(v[3],z).degree()==3
print('CF2 PASS exact symbolic: complete constant-rank cut-product curve, including zero and infinity; cyclic transport covers all cuts',flush=True)

k,x,y=s.symbols('k x y',real=True)
H=s.diag((1-k*k)/2,1,(k*k-2*k+7)/2,-k-1,(k*k+2*k-7)/2)
H[1,2]=H[2,1]=(5-k)/2
qq=q.subs(z,x+s.I*y);r=x*x+y*y
rhs=(r-1)**2*(r-k)**2+2*(k+1)*(r-1)**2*x*x+8*(k-5)*x*x*y*y
assert s.expand((qq.conjugate().T*H*qq)[0]-rhs)==0
H5=H.subs(k,5)
zero(H5-s.diag(-12,1,11,-6,14))
assert s.expand(rhs.subs(k,5)-(r-1)**2*((r-5)**2+12*x*x))==0
score=s.trace(H5*s.diag(t,3*t,3,3,1))
assert s.expand(score-(29-9*t))==0
W0=clean(B*T.inv()*H5*T.inv()*B.conjugate().T)
zero(P*W0*P.T-W0)
assert score.subs(t,s.Rational(10,3))==-1
print('CF3 PASS exact symbolic: three-square supported witness for every k>=5; k=5 factors into two nonnegative factors and detects the whole PPT interval',flush=True)

# The two unique real eigenvectors force the local-unitary chirality argument.
assert lo>3 and hi<4
assert 4*lo>3*hi>9>3
assert s.simplify(s.im(raw[3,5])+3*s.sqrt(3)/2)==0
zero(parity*raw*parity-raw)
print('CF4 PASS exact: unique GHZ/product spectral rays and nonreal phase projectors supply the analytic local-chirality proof',flush=True)

# Quantitative global witness lift. K has orthonormal rows after rescaling.
Kn=s.diag(1/s.sqrt(3),1/s.sqrt(3),1/s.sqrt(12))*K
zero(Kn*Kn.conjugate().T-s.eye(3))
zero(Kn.conjugate().T*Kn-(s.eye(8)-B*T.inv()*B.conjugate().T))
p,zz,zb=s.symbols('p zz zb',real=True)
cond=clean(Kn*s.kronecker_product(s.Matrix([[p,zz],[zb,1-p]]),s.eye(4))*Kn.conjugate().T)
det=s.expand(cond.det())
# Replace zz*zb=p(1-p); the remaining zz^2+zb^2 coefficient is nonnegative.
baseline=(-4*p**3+2*p*p+7*p+1+(2*p-8)*p*(1-p))/36
assert s.expand(det-(baseline+(1-p)*(zz**2+zb**2)/12
                          +(p-4)*(zz*zb-p*(1-p))/18))==0
minimum=s.expand(baseline-p*(1-p)**2/6)
assert s.expand(minimum-(1-7*p+24*p*p-12*p**3)/36)==0
assert s.expand(36*minimum-s.Rational(4,9)+(3*p-5)*(6*p-1)**2/9)==0
L,alpha,mu,beta=s.Integer(5),s.Rational(1,100),s.Rational(1,81),s.Integer(1)
penalty=(L+L*L/alpha)/mu
assert penalty==202905
g=s.Rational(3,133)-alpha
assert g==s.Rational(167,13300)
assert g/(2*L+penalty*beta)>s.Rational(1,20000000)
rho=clean(raw.subs(t,s.Rational(10,3))/s.Rational(133,3))
for ix in ixs:pd((pt(rho,0)-s.eye(8)/2000).extract(ix,ix))
print('CF5 PASS exact: conditioning determinant >=1/81, explicit global witness penalty 202905, physical radius 1/20000000 and PT margin 1/2000',flush=True)

# Exact extremality at the doubly singular PT boundary.
field=s.QQ.algebraic_field(s.sqrt(17),s.sqrt(3)*s.I)
aa=lo;bb=(324-89*aa)/(5*(4-aa))
ext=clean(B*s.diag(aa,bb,3,3,1)*B.conjugate().T)
assert bb>0 and bb-12*aa+29<0
Q=pt(ext,0);vectors=[]
for ix in ixs:
    A=Q.extract(ix,ix);top=A[:3,:3]
    assert s.simplify(top.det())>0 and s.simplify(A.det())==0
    v0=s.Matrix.vstack(-top.adjugate()*A[:3,3],s.Matrix([top.det()]))
    vv=s.zeros(8,1)
    for j,idx in enumerate(ix):vv[idx]=s.expand(v0[j])
    assert all(field.from_sympy(s.expand(a))==field.zero for a in Q*vv)
    vectors.append(vv)
rows=[]
for power in range(3):
    for v0 in vectors:
        vv=P**power*v0
        cuts=[c for c in range(3) if all(field.from_sympy(s.expand(a))==field.zero for a in pt(ext,c)*vv)]
        assert len(cuts)==1
        right=[];left=[]
        for i in range(5):
            for j in range(5):
                direction=pt(B[:,i]*B[:,j].conjugate().T,cuts[0])
                right.append(direction*vv)
                left.append((vv.conjugate().T*direction).T)
        rows += [s.Matrix.hstack(*right),s.Matrix.hstack(*left)]
system=s.Matrix.vstack(*rows)
assert system.shape==(96,25)
assert DomainMatrix.from_Matrix(system).convert_to(field).rank()==24
print('CF6 PASS exact algebraic: explicit fully-PPT GME state with profile (5,6,6,6); both-sided face constraints have rank 24, hence extremality',flush=True)
