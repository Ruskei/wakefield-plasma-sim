#import "@preview/physica:0.9.8": *

== Introduction & De Rham Sequences

Our space is a reduced 2D3V axisymmetric cylindrical space, where we still have the map:
$
  vec(r, theta, z) |-> vec(r cos theta, r sin theta, z), quad det J = r
$

then we can write (where hatted index means physical component)
$
  bold(E) = E_hat(r) bold(e)_r + E_hat(theta) bold(e)_theta + E_hat(z) bold(e)_z = E_hat(r) partial_r + E_hat(theta) partial_theta / r + E_hat(z) partial_z \
  E^r = E_hat(r), quad E^theta = E_hat(theta) / r, quad E^z = E_hat(z) \
$
use, where $g$ is the metric tensor
$
  E_i = g_(i j) E^j, quad E^flat = E_i dif x^i
$
for cylindrical coordinates
$
  g = diagonalmatrix(1, r^2, 1) \
  E_r = E^r, quad E_theta = r^2 E^theta, quad E_z = E^z \
  bold(E)^flat = E^r dif r + r^2 E^theta dif theta + E^z dif z = E_hat(r) dif r + r E_hat(theta) dif theta + E_hat(z) dif z
$
which can be checked to be the right object for line integrals. Then we care about magnetic flux
$
  integral_S beta = integral_S bold(E) dot bold(n) dif S
$
this is a surface integral so it requires a 2-form, so we take
$
  beta = star#h(0em)bold(E)^flat = (r E_hat(r)) dif theta and dif z - (E_hat(theta)) dif r and dif z + (r E_hat(z)) dif r and dif theta
$
which clearly contains the correct coefficients. Then we have the density object, which we get by differentiating
$
  tau = dif beta = r (nabla dot bold(E)) dif r and dif theta and dif z
$

if we set $partial_theta = 0$
$
  tau = [partial_r (r E_hat(r)) + partial_z (r E_hat(z))] dif r and dif theta and dif z
$

Our regular de Rham chain would be
$
  V^0 stretch(-->)^G V^1 stretch(-->)^C V^2 stretch(-->)^D V^3
$
if we have a 0-form $phi$ then
$
  nabla phi = partial_r phi + partial_z phi
$
which corresponds to the poloidal $r$ and $z$ components of $bold(E)^flat$. Then if we diferentiate again
$
  dif bold(E) &= dif (E_hat(r) dif r + r E_hat(theta) dif theta + E_hat(z) dif z) \
  &= -partial_z (r E_hat(theta)) dif theta and dif z + (partial_r E_hat(z) - partial_z E_hat(r)) dif r and dif z + partial_r (r E_hat(theta)) dif r and dif theta \
$
then match terms with $beta$, renaming the earlier $bold(E)$
$
  B_hat(theta) = (partial_z E_hat(r) - partial_r E_hat(z)) dif r and dif z
$
this concludes the first chain:
$
  phi stretch(-->)^display(hat(nabla)) E_hat(r) dif r + E_hat(z) dif z stretch(-->)^display(hat(nabla)_p times) -B_hat(theta) dif r and dif z
$
then for the other chain we match and solve
$
  B_hat(r) &= - partial_z E_hat(theta) \
  B_hat(z) &= 1/r partial_r (r E_hat(theta)) \
$
since these are only dependent on toroidal $E$, and $tau$ is only dependent on poloidal $B$, we get the second chain
$
  (r E_hat(theta)) dif theta stretch(-->)^display(hat(nabla)_theta times) ((r B_hat(r)) dif theta and dif z + (r B_hat(z)) dif r and dif theta) stretch(-->)^display(hat(nabla) dot) tau
$

we call these the poloidal and toroidal chains

== Constraints

In our reduced space, our mapping is really just
$
  F_(r z) : vec(r, z) -> vec(r, z)
$

so some of the singularities of the paper don't apply on the poloidal chain. We still want all our vector fields to be smooth and physical, so first, any $f(x, y, z)$ is invariant under rotations around the z-axis by axisymmetry. If we have $f(r, z)$ we can have functions like $f(r) = r$ which are't smooth. So instead we have $f(r^2, z)$.

For transverse terms, we note:
$
  hat(e)_r = (cos theta, sin theta, 0) = (x/r, y/r, 0), quad hat(e)_theta = (- sin theta, cos theta, 0) = (-y/r, x/r, 0)
$
then a pure radial or toroidal vector where $arrow(V)(0, z) != 0$ leads to inconsistent values from different angles, breaking axisymmetry, thus
$
  E_theta (r, z) = r a (r^2, z) \
  V_perp = a(x^2 + y^2, z) (x, y) + b(x^2 + y^2, z) (-y, x) = r a (r^2, z) hat(e)_r + r b(r^2, z) hat(e)_theta
$
thus transverse(radial/toroidal) components are odd in $r$
$
  V_r, V_theta = r V_1 (z) + r^3 V_3 (z) + r^3 V_5 (z) + ...
$
while axial/scalar components are even in $r$
$
  V_z, phi, rho = V_0 (z) + r^2 V_2 (z) + r^4 V_4 (z) + ...
$
for quantities
$
  f = phi, E_z, B_z, rho
$
we require
$
  partial_r f (0, z) = 0
$
because otherwise $partial_r$ is dependeont on angle, and otherwise but these values are perfectly well defined on the axis. For transverse
$
  v = E_r, E_theta, B_r, B_theta => v(0, z) = 0
$
since if $v(0,z) != 0$, axisymmetry is clearly broken. Now if we have
$
  u = r E_theta, quad p = r B_r, quad q = r B_z
$
which are the stored components of the toroidal chain, we need
$
  u(0,z) = partial_r u (0,z) = 0 \
  p(0,z) = partial_r p(0,z) = 0 \
  q(0,z) = 0
$
the $partial_r$ constraint comes from the fact that if
$
  partial_r u(0,z) = a(z) != 0
$
then
$
  u(r,z) = a(z) r + cal(O)(r^2)
$
which would imply $E_theta (0, z) != 0$ which fails. Thus $q$ doesn't have this restriction as it can be nonzero at $r=0$.

== Splines

We consider the logical poloidal and toroidal de Rham sequences of tensor product splines:
$
  "poloidal" : &
  hat(W)_(p,h)^0 := SS_(p,p) (hat(Omega)) stretch(-->)^display(hat(nabla))
  hat(W)_(p,h)^1 := vec(SS_(p-1,p) (hat(Omega)), SS_(p,p-1) (hat(Omega))) stretch(-->)^display(hat(nabla)_p times)
  hat(W)_(p,h)^2 := SS_(p-1,p-1) (hat(Omega)) \
  "toroidal" : &
  hat(W)_(theta,h)^1 := SS_(p,p) (hat(Omega)) stretch(-->)^display(hat(nabla)_theta times)
  hat(W)_(theta,h)^2 := vec(SS_(p,p-1) (hat(Omega)), SS_(p-1,p) (hat(Omega))) stretch(-->)^display(hat(nabla) dot)
  hat(W)_(theta,h)^3 := SS_(p-1,p-1) (hat(Omega))
$

then our bases are
$
  "poloidal"
  cases(
    hat(W)_(p,h)^0 & : hat(T)_(i j)^(p,0) (r, z) = B_i (r) B_j (z),
    hat(W)_(p,h)^1 & : hat(bold(T))_(i j)^(p,r) = display(vec(M_i (r) B_j (z), 0)) "and" hat(bold(T))_(i,j)^(p,z) = display(vec(0, B_i (r) M_j (z))),
    hat(W)_(p,h)^2 & : hat(T)_(i j)^(p,2) (r, z) = M_i (r) M_j (z),
  ) \
  "toroidal"
  cases(
    hat(W)_(theta,h)^1 & : hat(T)_(i j)^(theta,1) (r, z) = B_i (r) B_j (z),
    hat(W)_(theta, h)^2 & : hat(bold(T))_(i,j)^(theta,r) = display(vec(B_i (r) M_j (z), 0)) "and" hat(bold(T))_(i,j)^(theta,z) = display(vec(0, M_i (r) B_j (z))),
    hat(W)_(theta,h)^3 & : hat(T)_(i j)^(theta,3) (r, z) = M_i (r) M_j (z),
  ) \
$

without constraints, these spaces clearly allow for non-physical objects in axisymmetric space. We define
$
  cal(A)^cal(l) := {"RZ coefficient forms whose axisymmetric extension is smooth at" r = 0} \
  hat(V)_(dot,h)^cal(l) := hat(W)_(dot,h)^cal(l) inter cal(A)^cal(l), quad
  dif hat(V)_(dot,h)^cal(l) subset.eq hat(V)_(dot,h)^cal(l + 1)
$
these spaces are
$
  "poloidal" &
  cases(
    hat(V)_(p,h)^0 = {phi = sum_(i,j) phi_(i j) hat(T)_(i j)^(p,0) in hat(W)_(p,h)^0 : phi_(1j) - phi_(0j) = 0 quad forall j },
    hat(V)_(p,h)^1 = {bold(v) = sum_(i j) v_(i j)^r hat(bold(T))_(i j)^(p,r) + sum_(i,j) v_(i j)^z hat(bold(T))_(i j)^(p,z) in hat(W)_(p,h)^1 : v_(0j)^r = v^z_(1j) - v_(0j)^z = 0 quad forall j },
    hat(V)_(p,h)^2 = {f = sum_(i,j) f_(i j) hat(T)_(i j)^(p,2) in hat(W)_(p,h)^2 : f_(0 j) = 0 quad forall j}
  ) \
  "toroidal" &
  cases(
    hat(V)_(theta,h)^1 = {psi = sum_(i,j) psi_(i j) hat(T)_(i j)^(theta,1) in hat(W)_(theta,h)^1 : psi_(0j) = psi_(1j) = 0 quad forall j},
    hat(V)_(theta,h)^2 = {bold(w) = sum_(i, j) w_(i j)^r hat(bold(T))_(i j)^(theta,r) + sum_(i,j) w_(i j)^z hat(bold(T))_(i j)^(theta,z) in hat(W)_(theta,h)^2 : w_(0j)^r = w_(1j)^r = w_(0j)^z = 0 quad forall j},
    hat(V)_(theta,h)^3 = {p = sum_(i,j) p_(i j) hat(T)_(i j)^(theta,3) in hat(W)_(theta,h)^3 : p_(0j) = 0 quad forall j}
  )
$
we have to enforce that $partial_r phi(0,z) = 0$ so
$
  partial_r phi = sum_(i,j) phi_(i j) B_j (z) partial_r B_i (r) = sum_(i, j) B_j (z) (phi_(i+1,j) - phi_(i j)) M_i (r)
$
using the fact that
$
  sum_(i=0)^(n_s - 1) phi_(i j) (M_(i-1) (r) - M_i (r)) =
  sum_(i=0)^(n_s-2) (phi_(i+1,j) - phi_(i j)) M_i (r)
$
then for
$
  partial_r phi (0, z) = sum_(i, j) B_j (z) (phi_(i+1,j) - phi_(i j)) M_i (0) = 0 \
  => phi_(1,j) - phi_(0,j) = 0
$
In the 2nd element of poloidal chain we have the transverse $E_r$ and axial $E_z$, so
$
  v^r (0,z) = 0, quad partial_r v^z (0,z) = 0 \
  => v_(0j)^r = v^z_(1j) - v_(0j)^z = 0 quad forall j
$

We can verify these constraints, first the poloidal sequence
$
  phi_(1 j) - phi_(0j) => partial_r phi (0,z) = 0 => (partial_r phi) (0,z) = 0, quad partial_r (partial_z phi) (0,z) = 0, quad f (0,z) = 0
$
toroidal
$
  psi_(0j) = psi_(1j) = 0 => psi (0, z) = partial_r psi(0,z) = 0 => partial_z partial_r psi(0,z) = 0 => partial_r partial_z psi(0,z) + partial_z partial_r psi(0,z) = 0
$
so all the constraints imply eachother.

== Projections

Now we define the projection operators:
$
  P_p^0 : hat(W)_(p,h)^0 -> hat(V)_(p,h)^0, quad cases(
    P_p^0 hat(T)_(0 j)^(p,0) := 0,
    P_p^0 hat(T)_(1 j)^(p,0) := hat(T)_(0 j)^(p,0) + hat(T)_(1 j)^(p,0),
    P_p^0 hat(T)_(i j)^(p,0) := hat(T)_(i j)^(p,0) quad forall (i >= 2),
  )
$
$
  P_p^1 : hat(W)_(p,h)^1 -> hat(V)_(p,h)^1, quad cases(
    P_p^1 hat(bold(T))_(0 j)^(p,r) := 0,
    P_p^1 hat(bold(T))_(i j)^(p,r) := hat(bold(T))_(i j)^(p,r) quad forall (i >= 1),
    P_p^1 hat(bold(T))_(0 j)^(p,z) := 0,
    P_p^1 hat(bold(T))_(1 j)^(p,z) := hat(bold(T))_(0 j)^(p,z) + hat(bold(T))_(1 j)^(p,z),
    P_p^1 hat(bold(T))_(i j)^(p,z) := hat(bold(T))_(i j)^(p,z) quad forall (i >= 2),
  )
$
$
  P_p^2 : hat(W)_(p,h)^2 -> hat(V)_(p,h)^2, quad cases(
    P_p^2 hat(T)_(0 j)^(p,2) := 0,
    P_p^2 hat(T)_(i j)^(p,2) := hat(T)_(i j)^(p,2),
  )
$
$
  P_theta^1 : hat(W)_(theta,h)^1 -> hat(V)_(theta,h)^1, quad cases(
    P_theta^1 hat(T)_(0 j)^(theta,1) := 0,
    P_theta^1 hat(T)_(1 j)^(theta,1) := - sum_(k=2)^(n_s-1) hat(T)_(k j)^(theta,1),
    P_theta^1 hat(T)_(i j)^(theta,1) := hat(T)_(i j)^(theta,1) quad forall (i >= 2),
  )
$
$
  P_theta^2 : hat(W)_(theta,h)^2 -> hat(V)_(theta,h)^2, quad cases(
    P_theta^2 hat(bold(T))_(0 j)^(theta,r) := 0,
    P_theta^2 hat(bold(T))_(1j)^(theta,r) := - sum_(k=2)^(n_s-1) hat(bold(T))_(k j)^(theta,r),
    P_theta^2 hat(bold(T))_(i j)^(theta,r) := hat(bold(T))_(i j)^(theta,r) quad forall (i >= 2),
    P_theta^2 hat(bold(T))_(0 j)^(theta,z) := 0,
    P_theta^2 hat(bold(T))_(i j)^(theta,z) := hat(bold(T))_(i j)^(theta,z) quad forall (i >= 1),
  )
$
$
  P_theta^3 : hat(W)_(theta,h)^3 -> hat(V)_(theta,h)^3, quad cases(
    P_theta^3 hat(T)_(0 j)^(theta,3) := 0,
    P_theta^3 hat(T)_(i j)^(theta,3) := hat(T)_(i j)^(theta,3) quad forall (i >= 1),
  )
$

We need
$
  (hat(nabla)_theta times) P_theta^1 = P_theta^2 (hat(nabla)_theta times)
$
note
$
  hat(nabla)_theta times hat(T)_(i j)^(theta,1) = vec(
    - partial_z hat(T)_(i j)^(theta,1),
    partial_r hat(T)_(i j)^(theta,1)
  ) = vec(
    - B_i (r) (M_(j-1) (z) - M_j (z)),
    (M_(i-1) (r) - M_i (r)) B_j (z)
  ) = - hat(bold(T))_(i,j-1)^(theta,r) + hat(bold(T))_(i,j)^(theta,r) + hat(bold(T))_(i-1,j)^(theta,z) - hat(bold(T))_(i,j)^(theta,z) \
$
so
$
  (hat(nabla)_theta times) P_theta^1 hat(T)_(1 j)^(theta,1) = 0 \
  hat(nabla)_theta times hat(T)_(1 j)^(theta,1) = - hat(bold(T))_(1,j-1)^(theta,r) + hat(bold(T))_(1,j)^(theta,r) + hat(bold(T))_(0,j)^(theta,z) - hat(bold(T))_(1,j)^(theta,z) \
  P_theta^2 (hat(nabla)_theta times hat(T)_(1 j)^(theta,1)) = - hat(bold(T))_(1,j)^(theta,z)
$
which doesn't commute. Here the projection operator can set that coefficient to 0 only by setting this to 0... the curl operator is obviously fixed so the 4 vectors we get from $hat(nabla)_theta times hat(T)_(1 j)^(theta,1)$ are fixed. The simplest fix for this is $P_theta^2 hat(bold(T))_(1,j)^(theta,z) := 0$. Then we check
$
  (hat(nabla)_theta times) P_theta^1 hat(T)_(0j)^(theta,1) = 0 \
  P_theta^2 (hat(nabla)_theta times hat(T)_(0j)^(theta,1)) =
  P_theta^2 (
    - hat(bold(T))_(0,j-1)^(theta,r)
    + hat(bold(T))_(0,j)^(theta,r)
    + hat(bold(T))_(-1,j)^(theta,z)
    - hat(bold(T))_(0,j)^(theta,z)
  ) = 0
$
$
  (hat(nabla)_theta times) P_theta^1 hat(T)_(2 j)^(theta,1) =
    - hat(bold(T))_(2,j-1)^(theta,r)
    + hat(bold(T))_(2,j)^(theta,r)
    + hat(bold(T))_(1,j)^(theta,z)
    - hat(bold(T))_(2,j)^(theta,z) \
  P_theta^2 (hat(nabla)_theta times hat(T)_(2 j)^(theta,1)) = 
    - hat(bold(T))_(2,j-1)^(theta,r)
    + hat(bold(T))_(2,j)^(theta,r)
    - hat(bold(T))_(2,j)^(theta,z)
$
so clearly if we increase the threshold for $i$ to cut off in $hat(bold(T))_(i,j)^(theta,z)$ we have to keep increasing it further, destroying the space. Thus we must solve this another way, for $hat(V)_(theta,h)^1$ anything past $i >= 2$ is unrestricted, so we can have $P_theta^1 hat(T)_(1 j)^(theta,1)$ equal bases from there. Say we just have
$
  P_theta^1 hat(T)_(1 j)^(theta,1) := hat(T)_(2 j)^(theta,1) \
  (hat(nabla)_theta times) P_theta^1 hat(T)_(1 j)^(theta,1) = 
    - hat(bold(T))_(2,j-1)^(theta,r)
    + hat(bold(T))_(2,j)^(theta,r)
    + hat(bold(T))_(1,j)^(theta,z)
    - hat(bold(T))_(2,j)^(theta,z) \
  P_theta^2 (hat(nabla)_theta times hat(T)_(1 j)^(theta,1)) =
  P_theta^2 (
    - hat(bold(T))_(1,j-1)^(theta,r)
    + hat(bold(T))_(1,j)^(theta,r)
    + hat(bold(T))_(0,j)^(theta,z)
    - hat(bold(T))_(1,j)^(theta,z)
  )
  = - hat(bold(T))_(1,j)^(theta,z)
$
try to fix this by messing with $P_theta^2$. We have no terms in common
$
  P_theta^2 hat(bold(T))_(1,j)^(theta,r) & := hat(bold(T))_(2,j)^(theta_r) \
  P_theta^2 hat(bold(T))_(0,j)^(theta,z) & := 2 hat(bold(T))_(1,j)^(theta,z) - hat(bold(T))_(2,j)^(theta,z)
$
$
  P_theta^1 hat(T)_(0j)^(theta,1) := hat(T)_(2j)^(theta,1) \
  (hat(nabla)_theta times) P_theta^1 hat(T)_(0j)^(theta,1) =
    - hat(bold(T))_(2,j-1)^(theta,r)
    + hat(bold(T))_(2,j)^(theta,r)
    + hat(bold(T))_(1,j)^(theta,z)
    - hat(bold(T))_(2,j)^(theta,z) \
  P_theta^2 (hat(nabla)_theta times hat(T)_(0j)^(theta,1)) =
  P_theta^2 (
    - hat(bold(T))_(0,j-1)^(theta,r)
    + hat(bold(T))_(0,j)^(theta,r)
    + hat(bold(T))_(-1,j)^(theta,z)
    - hat(bold(T))_(0,j)^(theta,z)
  ) =
    - hat(bold(T))_(1,j-1)^(theta,r)
    + hat(bold(T))_(1,j)^(theta,r)
    
$

it seems like we need some sort of telescoping. Let's again note that if we have
$
  P_theta^1 hat(T)_(1 j)^(theta,1) := hat(T)_(2 j)^(theta,1) \
  (hat(nabla)_theta times) P_theta^1 hat(T)_(1 j)^(theta,1) = 
    - hat(bold(T))_(2,j-1)^(theta,r)
    + hat(bold(T))_(2,j)^(theta,r)
    + hat(bold(T))_(1,j)^(theta,z)
    - hat(bold(T))_(2,j)^(theta,z) \
$
this has one term almost matching except for a sign, so let's just flip it:
$
  P_theta^1 hat(T)_(1 j)^(theta,1) := -hat(T)_(2 j)^(theta,1) \
  (hat(nabla)_theta times) P_theta^1 hat(T)_(1 j)^(theta,1) = 
      hat(bold(T))_(2,j-1)^(theta,r)
    - hat(bold(T))_(2,j)^(theta,r)
    - hat(bold(T))_(1,j)^(theta,z)
    + hat(bold(T))_(2,j)^(theta,z) \
$
here we can see that the last 2 terms are telescoping, so
$
  
  P_theta^1 hat(T)_(1 j)^(theta,1) := - sum_(k=2)^(n_s-1) hat(T)_(k j)^(theta,1) \
  (hat(nabla)_theta times) P_theta^1 hat(T)_(1 j)^(theta,1) =
      sum_(k=2)^(n_s-1) hat(bold(T))_(k,j-1)^(theta,r)
    - sum_(k=2)^(n_s-1) hat(bold(T))_(k,j)^(theta,r)
    - sum_(k=2)^(n_s-1) hat(bold(T))_(k-1,j)^(theta,z)
    + sum_(k=2)^(n_s-1) hat(bold(T))_(k,j)^(theta,z) \
    =
      sum_(k=2)^(n_s-1) hat(bold(T))_(k,j-1)^(theta,r)
    - sum_(k=2)^(n_s-1) hat(bold(T))_(k,j)^(theta,r)
    -hat(bold(T))_(1,j)^(theta,z)
$
Now
$
  P_theta^2 (hat(nabla)_theta times hat(T)_(1 j)^(theta,1)) = P_theta^2(
    - hat(bold(T))_(1,j-1)^(theta,r)
    + hat(bold(T))_(1,j)^(theta,r)
    + hat(bold(T))_(0,j)^(theta,z)
    - hat(bold(T))_(1,j)^(theta,z)
  )
$
we match terms to get
$
  P_theta^2 hat(bold(T))_(1j)^(theta,r) := - sum_(k=2)^(n_s-1) hat(bold(T))_(k j)^(theta,r) \
  P_theta^2(
    - hat(bold(T))_(1,j-1)^(theta,r)
    + hat(bold(T))_(1,j)^(theta,r)
    + hat(bold(T))_(0,j)^(theta,z)
    - hat(bold(T))_(1,j)^(theta,z)
  ) =
  sum_(k=2)^(n_s-1) hat(bold(T))_(k,j-1)^(theta,r) - sum_(k=2)^(n_s-1) hat(bold(T))_(k,j)^(theta,r) - hat(bold(T))_(1,j)^(theta,r) \
  (hat(nabla)_theta times) P_theta^1 hat(T)_(1j)^(theta,1) = P_theta^2 (hat(nabla)_theta times) hat(T)_(1j)^(theta,1)
$
then
$
  (hat(nabla)_theta times) P_theta^1 hat(T)_(0j)^(theta,1) = 0 \
  P_theta^2 (hat(nabla)_theta times) hat(T)_(1j)^(theta,1) = P_theta^2 (
    - hat(bold(T))_(0,j-1)^(theta,r)
    + hat(bold(T))_(0,j)^(theta,r)
    + hat(bold(T))_(-1,j)^(theta,z)
    - hat(bold(T))_(0,j)^(theta,z)
  ) = 0
$
and in general for $i >= 2$
$
  (hat(nabla)_theta times) P_theta^1 hat(T)_(i j)^(theta,1) = 
    - hat(bold(T))_(i,j-1)^(theta,r)
    + hat(bold(T))_(i,j)^(theta,r)
    + hat(bold(T))_(i-1,j)^(theta,z)
    - hat(bold(T))_(i,j)^(theta,z) \
  P_theta^2 (hat(nabla)_theta times) hat(T)_(i j)^(theta,1) =
    - hat(bold(T))_(i,j-1)^(theta,r)
    + hat(bold(T))_(i,j)^(theta,r)
    + hat(bold(T))_(i-1,j)^(theta,z)
    - hat(bold(T))_(i,j)^(theta,z)
$
then for the next set we must have
$
  (hat(nabla) dot) P_theta^2
  = P_theta^3 (hat(nabla) dot)
$

$
  hat(nabla) dot hat(bold(T))_(i j)^(theta,r) = hat(nabla) dot vec(B_i (r) M_j (z), 0) = (M_(i-1) (r) - M_i (r)) M_j (z) = hat(T)_(i-1,j)^(theta,3) - hat(T)_(i,j)^(theta,3)
$

$
  (hat(nabla) dot) P_theta^2 hat(bold(T))_(1j)^(theta,r) = - hat(nabla) dot (sum_(k=2)^(n_s-1) hat(bold(T))_(k j)^(theta,r)) = sum_(k=2)^(n_s-1) (hat(T)_(k,j)^(theta,3) - hat(T)_(k-1,j)^(theta,3)) = - hat(T)_(1,j)^(theta,3)
  \
  P_theta^3 (hat(nabla) dot) hat(bold(T))_(1j)^(theta,r) = P_theta^3 (hat(T)_(0,j)^(theta,3) - hat(T)_(1,j)^(theta,3)) = - hat(T)_(1,j)^(theta,3)
$
it can be verified that all the remaining cases commute! The original operator definitions have been updated.
