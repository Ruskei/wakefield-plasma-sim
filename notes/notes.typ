#import "@preview/physica:0.9.8": *

We want to derive the reduced Maxwell equations. The full fields are:

#grid(
  columns: (1fr,) * 2,
  $
    r E_theta (r,z,t) &= sum_(i,j) e_(i j)^theta hat(T)_(i j)^(theta,1) \
    bold(E)^p (r,z,t) &= sum_(i,j) (e_(i j)^(p,r) hat(bold(T))_(i j)^(p,r) + e_(i j)^(p,z) hat(bold(T))_(i j)^(p,z))
  $,
  $
    B_theta (r,z,t) &= sum_(i,j) b_(i j)^theta hat(T)_(i j)^(p,2) \
    r bold(B)^theta (r,z,t) &= sum_(i,j) (b_(i j)^(theta,r) hat(bold(T))_(i j)^(theta,r) + b_(i,j)^(theta,z) hat(bold(T))_(i j)^(theta,z)) \
  $
)

for the reduced version, first note that the reduced coefficients are defined as
$
  bold(e)^theta
    = R_theta^1 macron(bold(e))^theta
$
thus
$
  r E_theta
    = sum_(i,j) (R_theta^1 macron(bold(e))^theta)_(i j)
    hat(T)_(i j)^(theta,1) \
    = Lambda_theta^1 R_theta^1 macron(bold(e))^theta
$
where $Lambda_(theta,"full")^1$ is the full row vector of basis functions evaluated at $(r,z)$. from this we see
$
  macron(Lambda)^cal(l) = Lambda^cal(l) R^cal(l)
$
for any of our relevant bases. Recalling Gauss law in the full space
$
  cal(G)^T cal(M)^(p,1) bold(e)^p = -bold(rho)^0
$
remember however that $rho$ is a functional. If we have
$
  phi = R_p^0 macron(phi) \
  rho (phi)
    = rho(R_p^0 macron(phi)) \
    = macron(phi)^T (R_p^0)^T bold(rho)
$
thus
$
  macron(bold(rho))^0
    = (R_p^0)^T bold(rho)^0
$

now substitute for $bold(e)$ into Gauss law to get
$
  cal(G)^T cal(M)^(p,1) R_p^1 macron(bold(e))^p = -bold(rho)^0
$
we saw that testing against a functional with in reduced space is same as left-multiplying by transpose $R$ operator, thus
$
  (R_p^0)^T cal(G)^T cal(M)^(p,1) R_p^1 macron(bold(e))^p
    = -macron(bold(rho))^0
$
use the commutative property
$
  cal(G) R_p^0
    = R_p^1 macron(cal(G))
$

in the definition

== Mass Matrices

Using index sets $mu = (i,j), nu = (k,l)$
#grid(
  columns: (1fr,) * 2,
  $
    cal(M)_(mu,nu)^(p,0) &= 2 pi integral_hat(Omega) hat(T)_mu^(p,0) hat(T)_nu^(p,0) med r dif r dif z \
    cal(M)_(mu,nu)^(p,1) &= 2 pi integral_hat(Omega) hat(bold(T))_mu^(p,1) dot hat(bold(T))_nu^(p,1) med r dif r dif z \
    cal(M)_(mu,nu)^(p,2) &= 2 pi integral_hat(Omega) hat(T)_mu^(p,2) hat(T)_nu^(p,2) med r dif r dif z \
  $,
  $
    cal(M)_(mu,nu)^(theta,1) &= 2 pi integral_hat(Omega) hat(T)_mu^(theta,1) hat(T)_nu^(theta,1) med 1/r dif r dif z \
    cal(M)_(mu,nu)^(theta,2) &= 2 pi integral_hat(Omega) hat(bold(T))_mu^(theta,2) dot hat(bold(T))_nu^(theta,2) med 1/r dif r dif z \
    cal(M)_(mu,nu)^(theta,3) &= 2 pi integral_hat(Omega) hat(T)_mu^(theta,3) hat(T)_nu^(theta,3) med 1/r dif r dif z \
  $
)

our reduced versions are:
$
  macron(cal(M))^cal(l) = 2 pi integral_hat(Omega) (Lambda^cal(l) R_theta^cal(l))^T (Lambda^cal(l) R_theta^cal(l)) 1/r dif r dif z
$

notice that $Lambda^(p,1)$ is a matrix with 2 rows, and is not just a row vector. 
$
  bold(E)^p (r,z,t) &= sum_(i,j) (e_(i j)^(p,r) hat(bold(T))_(i j)^(p,r) + e_(i j)^(p,z) hat(bold(T))_(i j)^(p,z))
    = sum_(i,j)[mat(
      hat(bold(T))_(i,j)^(p,r),
      hat(bold(T))_(i,j)^(p,z),
    ) vec(
      e_(i,j)^(p,r),
      e_(i,j)^(p,z),
    )]
$
so we really have something more like
$
  bold(E)^p = 
    mat(
      Lambda^(p,1,r), 0;
      0, Lambda^(p,1,z)
    )
    vec(
      bold(e)^(p,r),
      bold(e)^(p,z),
    )
$
we just want $macron(cal(M))_(mu,nu)^(p,1)$, which is
$
  macron(cal(M))_(mu,nu)^cal(l)
    = 2 pi integral_hat(Omega)
    ((Lambda^cal(l) R_theta^cal(l))^T)_mu (Lambda^cal(l) R_theta^cal(l))_nu
    1/r dif r dif z \
  (Lambda^cal(l) R_theta^cal(l))_nu
    = sum_m Lambda_m^cal(l) 
$
