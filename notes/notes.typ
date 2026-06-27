#import "@preview/physica:0.9.8": *
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

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

but if we have potential then
$
  cal(G)^T cal(M)^(p,1) cal(G) phi = -bold(rho)^0
$
substitute the equation from a reduced fuctional
$
  cal(G)^T cal(M)^(p,1) cal(G) R_p^0 macron(phi) = -bold(rho)^0
$
use commutative property to get
$
  cal(G)^T cal(M)^(p,1) R_p^1 macron(cal(G)) macron(phi) = -bold(rho)^0
$
left multiply by transpose of reduction operator
$
  (R_p^0)^T cal(G)^T cal(M)^(p,1) R_p^1 macron(cal(G)) macron(phi) = -macron(bold(rho))^0
$
use transpose and commutation property to reach
$
  macron(cal(G))^T macron(cal(M))^(p,1) macron(cal(G)) macron(phi) = -macron(bold(rho))^0
$

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

=== Section 4 rz.pdf follow-along

(3.3):
$
  dot(bold(Xi)) &= NN (bold(Xi))^TT bold("V") \
  dot(bold(V)) &= MM_q MM_m^(-1) NN(bold(Xi))
    (tilde(bold("E")) (bold(Xi), t)
    + (NN(Xi)^TT bold("V")) times tilde(bold("B"))
    (bold(Xi), t))
$

(3.8a):
$
  tilde(sans("M"))_1 dot(tilde("e"))
    = sans("C")^TT tilde(sans("M"))_2
    tilde(bold("b")) - MM_q tilde(bb(L))^1 (bold(Xi))^TT
    NN(bold(Xi))^TT bold("V")
$
(3.9a):
$
  dot(tilde(bold("b"))) = - sans("C") tilde(bold("e"))
$

the Hamiltonian $cal(H)$ for the plasma system is
$
  cal(H) = sum_s m_s / 2
    integral abs(bold("v"))^2 f_s (bold(x), bold(v))
    dif bold("x") dif bold("v")
    + 1/2 integral abs(bold("E") (bold("x")))^2
    + abs(bold("B") (bold("x")))^2 dif bold("x")
$

then it's rewritten in curvilinear coordinates which can be written in matrix notation.

then the derivative of the discrete Hamiltonian 

=== Derivation for Evolution Method

a vector is written
$
  bold(v)
    = v^r "e"_r
    + v^theta "e"_theta
    + v^z "e"_z
$
and note that in Cartesian coordinates
$
  e_r &= (cos theta) hat(bold("x"))
    + (sin theta) hat(bold("y")) \
  e_theta &= -(r sin theta) hat(bold("x"))
    + (r cos theta) hat(bold("y")) \
$
let's start with a position vector and find time derivative of that for velocity
$
  bold(x) = r e_r + z e_z \
  dot(bold(x)) = bold(v)
    = dot(r) e_r + r dot(e)_r
    + dot(z) e_z \
  = dot(r) e_r + r (dot(theta) / r e_theta) + dot(z) e_z \
  = dot(r) e_r + dot(theta) e_theta + dot(z) e_z \
$
take time derivative of velocity
$
  dot(bold(v))
    = dot(v)^r e_r + v^r dot(e)_r
    + dot(v)^theta e_theta + v^theta dot(e)_theta
    + dot(v)^z e_z
$
$
  dot(e)_r = dv(e_r,t)
    = pdv(e_r,theta) pdv(theta,t)
    = dot(theta) / r e_theta \
  dot(e)_theta = dv(e_theta,t)
    = pdv(e_theta,r) pdv(r,t)
    + pdv(e_theta,theta) pdv(theta,t)
    = dot(r) / r e_theta - dot(theta) r e_r \
$
$
  bold(v)
    = dot(v)^r e_r + v^r (dot(theta) / r e_theta)
    + dot(v)^theta e_theta + v^theta (dot(r) / r e_theta - dot(theta) r e_r)
    + dot(v)^z e_z \
    = (dot(v)^r - (v^theta)^2 r) e_r
    + (dot(v)^theta + 2 (v^r v^theta)/r) e_theta
    + dot(v)^z e_z
$
this is with the covariant basis, but we want it with respect to a unit vector basis. this is basically the same as having $u^theta = r v^theta$. also in our formula we have "time derivative = acceleration + terms" so lets solve for that
$
  a^r = dot(v)^r - (v^theta)^2 r
    = dot(u)^r - (u^theta)^2 / r \
    dot(u)^r = a^r + (u^theta)^2 / r \
  a^theta = dot(v)^theta + 2 (v^r v^theta) / r
    = dv(,t) (r^(-1) u^theta) + 2 (u^r u^theta) / r^2 \
    = dot(u)^theta / r + (u^r u^theta) / r^2 \
    dot(u)^theta = r a^theta - (u^r u^theta) / r \
    qed
$

also let's rewrite the position time derivative with unit vectors
$
  bold(x)
    = r hat(e)_r + z hat(e)_z \
  dot(bold(x))
    = dot(r) hat(e)_r
    + r dot(theta) hat(e)_theta
    + dot(z) hat(e)_z
$
so we see that in total our equation are
$
  u^r = dot(r), u^theta = r dot(theta), u^z = dot(z)
$
$
  dot(bold(x))
    = vec(u^r, u^z) \
  dot(bold(u))
    = vec(a^r, r a^theta, a^z)
    + vec((u^theta)^2 \/ r, -u^r u^theta \/ r, 0)
$
awesome sauce. so let's fully summarize our equations of motion. for $N_p$ particles define $bold(Xi) = (bold(xi)_1, ..., bold(xi)_N_p)^TT$ where $bold(xi)_i in RR^2$ representing $r, z$ components. then $bold("V") = (bold(v)_1, ..., bold(v)_N_p)$ where $bold(v) in RR^3$ is the physical velocity. Where $omega_p, m_p, q_p$ is the number of particles, mass, and charge represented by each macroparticle, $MM_m := diag(omega_p m_p) times.o II_3, MM_q := diag(omega_p q_p) times.o II_3$.
$
  dot(bold(xi))_i &= bold("P") bold(v)_i,
    quad bold("P") = mat(
      1, 0, 0;
      0, 0, 1;
    ) \
  dot(bold(v))_i &= (omega_i q_i) / m_i
    (bold("E") (bold(xi_i), t) + bold(v)_i times
    bold("B") (bold(xi_i), t))
    + vec(v_(i, theta)^2 \/ xi_(i,r),
    - v_(i,r) v_(i, theta) \/ xi_(i,r), 0) \
$
Gauss law:
$
  macron(cal(G))^T macron(cal(M))^(p,1) macron(cal(G)) macron(phi) = -macron(bold(rho))^0
$
Ampere law:
$
  macron(cal(M))^(p,1) dot(macron(bold(e)))^p &= (macron(cal(C))^p)^T macron(cal(M))^(p,2) macron(bold(b))_theta - macron(bold(j))^p \
  macron(cal(M))^(theta,1) dot(macron(bold(e)))^theta &= (macron(cal(C))^theta)^T macron(cal(M))^(theta,2) macron(bold(b))^p - macron(bold(j))^theta \
$
Faraday and magnetic Gauss law:
$
  macron(cal(C))^p macron(bold(e))^p = - dot(macron(bold(b)))_theta, quad macron(cal(C))^theta macron(bold(e))^theta = - macron(dot(bold(b)))^p \
  macron(cal(D)) macron(bold(b))^p = 0
$

derivation of reduced Ampere's law:
$
  cal(M)^(p,1) dot(bold(e))^p = (cal(C)^p)^T cal(M)^(p,2) bold(b)^theta - bold(j)^p \
  bold(e)^p = R_p^1 macron(bold(e))^p, quad
  bold(b)^p = R_p^2 macron(bold(b))^p \
  cal(M)^(p,1) R_p^1 macron(bold(e))^p
    = (cal(C)^p)^T cal(M)^(p,2) R_p^2 macron(bold(b))^p - bold(j)^p \
$
again, testing against a functional is same as left multiplying by transpose $R$ operator, thus
$
  macron(cal(M))^(p,1) macron(bold(e))^p
    = (macron(cal(C))^p)^TT macron(cal(M))^(p,2) macron(bold(b))^p - macron(bold(j))^p \
$
the rest of the cases are similar.

These are our governing equations, and the goal now is to take an initial state and propagate it in time. Notice the dependencies for our equations $U_xi, U_v, U_e^p, U_e^theta, U_b^p, U_b^theta$

$
  U_xi : cases(
    v_r --> dot(xi)_r,
    v_z --> dot(xi)_z,
  )
$
$
  U_v : cases(
    bold(e)_r^p\, bold(xi)\, v_theta\, v_z\, bold(b)^theta\, bold(b)_z^p --> dot(v)_r,
    bold(e)^theta\, bold(v)\, bold(xi)\, bold(b)^p --> dot(v)_theta,
    bold(e)_z^p\, bold(xi)\, v_r\, v_theta\, bold(b)_r^p\, bold(b)^theta --> dot(v)_z,
  )
$
$
  U_e^p : bold(b)^theta, bold(j)^p --> dot(bold(e))^p \
  U_e^theta : bold(b)^p, bold(j)^theta --> dot(bold(e))^theta \
$
$
  U_b^p : bold(e)^theta --> dot(bold(b))^p \
  U_b^theta : bold(e)^p --> dot(bold(b))^theta \
$
