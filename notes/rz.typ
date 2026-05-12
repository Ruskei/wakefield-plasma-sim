#import("@preview/physica:0.9.8"): *

== De Rham Complexes

In axisymmetric cylindrical coordinates
$
  F(r, theta, z) = vec(r cos theta, r sin theta, z) \
  (D F (bold(xi)))_(i j) = pdv(x_i, xi_j) =
  mat(
    1, 0, 0;
    0, r, 0;
    0, 0, 1
  ), #h(1em) J_F (bold(xi)) = det(D F (bold(xi))) => r \
  N = (D F (bold(xi)))^(-T) =>
  mat(
    1, 0, 0;
    0, 1\/r, 0;
    0, 0, 1;
  ) \
  G_t = (D F)^T (D F) \
  G_n = N^T N \
  nabla f = partial_r f e_r + partial_z f e_z \
  nabla^2 f = partial_(r r) f + 1/r partial_r f + partial(z z) f \
  nabla dot bold(A) = 1/r partial_r (r A_r) + partial_z A_z \
  nabla times bold(A) = - partial_z A_0 bold(e)_r + (partial_z A_r - partial_r A_z) bold(e)_theta + 1/r partial_r (r A_theta) bold(e)_z \
$
now we need to find the discrete gradient, curl, and divergence matrices which satisfy:
$
  nabla_bold(xi) tilde(Lambda)^0 (bold(xi)) = tilde(bold(Lambda))^1 (bold(xi)) "G", #h(1em) nabla_bold(xi) times tilde(bold(Lambda))^1 (bold(xi)) = tilde(bold(Lambda))^2 (bold(xi)) "C", #h(1em) nabla_bold(xi) dot tilde(bold(Lambda))^2 (bold(xi)) = tilde(Lambda)^3 (bold(xi)) "D" \
  "G" in RR^(3 N_1 times N_0), #h(1em) "C" in RR^(3N_2 times 3N_1), #h(1em) "D" in RR^(N_3 times 3N_2) \
  "DC" = 0, "CG" = 0 \
$
with full axisymmetry we still have 3D space, since entities can move in $theta$ direction, so $"position" = (r, z)$ but $v = (r, theta, z)$, so we're in 2D3V. Anyway with the above definitions for $D F, J, N$, we have:
$
  tilde(bold(E)) = vec(E_r, r E_theta, E_z), #h(1em) tilde(bold(B)) = vec(r B_r, B_theta, r B_z)
$
$
  nabla_bold(xi) times tilde(bold(E)) = vec(
    - partial_z tilde(E)_theta,
    partial_z tilde(E)_r - partial_r tilde(E)_z,
    partial_r tilde(E)_theta
  )
$
full de Rham chain is, scalar potential, electric field, magnetic field, density spaces
$
  V^0 stretch(-->)^G V^1 stretch(-->)^C V^2 stretch(-->)^D V^3 \
  V^1 = V_p^1 plus.o V_theta^1 \
  V^2 = V_p^2 plus.o V_theta^2 \
$
where $V_p$ means poloidal
$
  V_p^1 = V_r^1 plus.o V_z^1 \
  V_p^2 = V_r^2 plus.o V_z^2 \
$
$
  V_p^1 : (tilde(E)_r, tilde(E)_z), V_theta^1 : tilde(E)_0 \
  V_p^2 : (tilde(B)_r, tilde(B)_z), V_theta^2 : tilde(B)_0 \
  C_t (tilde(E)_r, tilde(E)_z) = partial_z tilde(E)_r - partial_r tilde(E)_z \
  C_p (tilde(E)_theta) = (-partial_z tilde(E)_theta, partial r tilde(E)_theta)
$
where $C_t$ produces the toroidal magnetic field, and $C_p$ produces the poloidal field.
$
  V^0 stretch(-->)^(G_p) V_p^1 stretch(-->)^(C_t) V_theta^2
$
which means scalar $->$ poloidal E $->$ toroidal B, or
$
  phi -> (tilde(E)_r, tilde(E)_z) -> tilde(B)_theta
$

secondly:
$
  V_theta^1 stretch(-->)^(C_p) V_p^2 stretch(-->)^(D_p) V^3 \
  tilde(E)_theta -> (tilde(B)_r, tilde(B)_z) -> rho
$
since these two complexes are independent of eachother in the derivatives, we can choose different basis spaces for them. Spline spaces are:
$
  V^0 & = S_r^p times.o S_z^p \
  V_p^1 & = (S_r^(p-1) times.o S_z^p) dif r plus.o (S_r^p times.o S_z^(p-1)) dif z \
  V_t^2 & = S_r^(p-1) times.o S_z^(p-1) \
  V_t^1 & = S_r^p times.o S_z^p \
  V_p^2 & = (S_r^p times.o S_z^(p-1)) B_r plus.o (S_r^(p-1) times.o S_z^p) B_z \
  V^3 & = S_r^(p-1) times.o S_z^(p-1) \
$
now for Maxwell's equations, Faraday's law:
$
  nabla times bold(E) = - pdv(bold(B), t) \
$
$
  partial_t tilde(B)_r &= + partial_z tilde(E)_theta \
  partial_t tilde(B)_theta &= - partial_z tilde(E)_r + partial_r tilde(E)_z \
  partial_t tilde(B)_z &= -partial_r tilde(E)_theta \
$
$
  dot(b)_p = - C_p e_theta, #h(1em) dot(b)_theta = - C_t e_p
$
Gauss magnetic law:
$
  D_p b_p = 0
$
Ampere's law:
$
  pdv(bold(E), t) = nabla times bold(B)  - bold(J) \
  integral_Omega bold(phi) dot pdv(bold(E), t) dif bold(x) = integral_Omega nabla times bold(phi) dot bold(B) dif bold(x) - bold(J)(bold(phi))(t) #h(1em) bold(phi) in H("curl", Omega) \
  dv(, t) integral_tilde(Omega) N tilde(bold(phi)) dot N tilde(bold(E)) abs(J_F) dif bold(xi) = integral_tilde(Omega) (D F) / J_F nabla times tilde(bold(phi)) dot (D F) / J_F tilde(bold(B)) abs(J_F) dif bold(xi) - tilde(bold(J)) (N tilde(bold(phi))) \
$
Gauss electric law:
$
  - integral_tilde(Omega) N nabla tilde(psi) dot N tilde(bold(E)) abs(J_F) dif bold(xi) = tilde(rho) (tilde(psi))
$

we can separate these out by component for each complex.

== Splines

Now we need to find fitting basis functions that let us preserve de Rham structure and let us handle the $s=0$ singularity.

To extend the polar formulation we need to set $partial_theta = 0$ and eliminate the $theta$ component in practice. We then need to verify both our de Rham chains on these new spaces.

The extension to z-axis is likely trivial so we'll to the space reduction first. We should end up with a 1D spline space and a reduced polar singularity.

We have an open domain $Omega in RR^2$ which is a reduced logical cylinder $[0, L] times [0, W]$. We have a $C^1$ bijective mapping to the closure $overline(Omega)$ which is the image of this mapping
$
  F : hat(Omega) in.rev vec(s, z) |-> bold(x) in overline(Omega)
$

this is bijective rather than the paper's surjective mapping since there is only one point in each slice where $s=0$. We'll still write:
$
  overline(Omega)_0 := overline(Omega) \\ {bold(x)_0}, hat(Omega)_0 := hat(Omega) \\ {s=0} = F^(-1) (overline(Omega)_0)
$

now we can't convey the details of our physical space through a mapping, but looking at Definition 2.4, we can borrow ideas about the properties of first order polar singularities. The paper states that it's Jacobian matrix is of the form:
$
  J_F (s, theta)
  mat(
    C(theta) + cal(O)(s), s (C' (theta) + cal(O)(s));
    S(theta) + cal(O)(s), s (S' (theta) + cal(O)(s));
  ), #h(1em) det J_F (s, theta) = s (D(theta) + cal(O) (s))
$
well in our case we have, where $hat(s), hat(z)$ correspond to logical variables:
$
  J_F = jmat(s, z; hat(s), hat(z);big:#true) =
  mat(
    1, 0;
    0, 1;
  )
$

this means we can't fully pretend we're in 2D because we still have a $dif theta$ element, so from now on $Omega in RR^3$ and
$
  F : hat(Omega) in.rev vec(s, theta, z) |-> bold(x) in overline(Omega)
$
then our analytical mapping is:
$
  F : vec(s, theta, z) |-> vec(s cos theta, s sin theta, z)
$
if we didn't reduce our spline mapping would be:
$
  F : vec(s, theta, z) |-> vec(0, 0, z) + sum_(i=0)^(n_s-1) sum_(j=0)^(n_theta-1) bold(P)_(i j) B_i (s) circle(B)_j (theta)
$
but due to axisymmetry:
$
  F : vec(s, theta, z) |-> vec(0, 0, z) + sum_(i=0)^(n_s-1) bold(P)_i (theta) B_i (s) \
$
where
$
  bold(P)_i (theta) = vec(rho_i cos theta, rho_i sin theta, 0)
$
here $B_i$ is a B-spline of degree $p >= 1$ in the $s$ variable, which are defined by an open knot sequence
$
  s_0 = ... = s_p < s_(p+1) < ... < s_(n_s-1) < s_n_s = ... = s_(n_s+p) = L \
$
which just adds enough duplicates at start and end to ensure nice behavior. Also, Assumption 2.2 has a requirement on spline resolution in $theta$:
$
  n_theta = 4n'
$
for us this is just $n' = infinity$. For our analytical mapping the Jacobian is already known, written in the start in the cylindrical basis, but for our spline mapping we calculate:
$
  J_F = jmat(F_1, F_2, F_3; s, theta, z; big:#true) =
  mat(
    (cos theta) sum_(i=0)^(n_s-2) (rho_(i+1) - rho_i) M_(i) (s), -(rho_i sin theta) sum_(i=0)^(n_s - 1) B_i (s), 0;
    (sin theta) sum_(i=0)^(n_s-2) (rho_(i+1) - rho_i) M_(i) (s), (rho_i cos theta) sum_(i=0)^(n_s - 1) B_i (s), 0;
    0, 0, 1;
  )
$

Recall our De Rham chains:

$
  V^0 stretch(-->)^(G_p) V_p^1 stretch(-->)^(C_t) V_theta^2 & : phi -> (tilde(E)_r, tilde(E)_z) -> tilde(B)_theta \
  V_theta^1 stretch(-->)^(C_p) V_p^2 stretch(-->)^(D_p) V^3 & : tilde(E)_theta -> (tilde(B)_r, tilde(B)_z) -> p \
$
for our splines this means:
$
  hat(W)_h^0 := SS_(p, p) (hat(Omega)) stretch(->, size: #150%)^display(hat("grad")_p) hat(W)_(p, h)^1 :=
  vec(
    SS_(p - 1, p) (hat(Omega)),
    SS_(p, p - 1) (hat(Omega))
  )
  stretch(->, size: #150%)^display(hat("curl")_t) W_(theta, h)^2 := SS_(p-1, p-1) (hat(Omega)) \
  hat(W)_(theta, h)^1 := SS_(p, p) (hat(Omega)) stretch(->, size: #150%)^display(hat("curl")_p) hat(W)_(p, h)^2 :=
  vec(
    SS_(p, p - 1) (hat(Omega)),
    SS_(p - 1, p) (hat(Omega))
  ) stretch(->, size: #150%)^display(hat("div")_p) W_h^3 := SS_(p-1,p-1) (hat(Omega))
$

Here $SS_(p_1, p_2) (hat(Omega)) := SS_p_1 ([0, L)) times.o SS_(p_2) ([0, W])$. These look very similar, because the latter is just the rotated version of the former as pointed out at the start of the paper. For now we'll just consider the first chain. Our bases are:
$
  cases(
    "for" hat(W)_h^0: & hat(T)_(i j) (s, z) = B_i (s) B_j (z),
    "for" hat(W)_(p, h)^1 : & hat(bold(T))_(i j)^s (s, z) = vec(M_i (s) B_j (z), 0) "and"
    hat(bold(T))_(i j)^z (s, z) = vec(0, B_i (s) M_j (z)),
    "for" hat(W)_(theta, h)^2 : & hat(T)_(i j)^2 = M_i (s) M_j (z)
  )
$

On the physical reduced cylindrical domain $Omega$ we consider the pushforward operators:

$
  W_h^0 := cal(F)^0 hat(W)_h^0,
  quad
  W_(p, h)^1 := cal(F)^1 hat(W)_(p,h)^1,
  quad
  W_(theta,h)^2 := cal(F)^2 hat(W)_(theta,h)^2 \
  cases(
    cal(F)^0 : & hat(phi) |-> phi := hat(phi) compose F^(-1),
    cal(F)^1 : & hat(bold(v)) |-> bold(v) := (J_F^(-T) hat(bold(v))) compose F^(-1),
    cal(F)^2 : & hat(f) |-> f := (det J_F^(-1) hat(f)) compose F^(-1),
  )
$

$cal(F)^0$ makes sense because a simple scalar value at a coordinate doesn't care about distance from origin. $cal(F)^2$ relates in the first chain to toroidal vectors, which is scaled by $r$ when we go from logical to physical. Thus $cal(F)^2$ makes sense. For $cal(F)^1$ we're transforming whole vectors, which are fully 3D, which means $J_F$ is definitely 3D, which checks out. Our spline spaces are only 2D, but we remember that the first sequence only stores $tilde(E)_r, tilde(E)_z$ fields so this makes sense. But this means the way we've currently defined these pushforward operators is acting on the full 3D scheme and not our split scheme, so the above is wrong.

Then for the first chain
$
  cal(F)^0 : hat(phi) |-> phi := hat(phi) compose F^(-1)
$
this is clearly correct due to the above reasoning. Then for the gradient of this, $hat(W)_(p,h)^1$, contains just poloidal components, none of which change with cylindrical coordinates. This part preserves line integrals, which if we evaluate without a $theta$ component in non-reduced cylindrical coordinates, since it'll be in the form:
$
  integral bold(v) dif s, integral bold(v) dif v
$

thus
$
  cal(F)^1 : & hat(bold(v))_p |-> bold(v) := hat(bold(v))_p compose F^(-1),
$
now for the toroidal flux term, $tilde(B)_theta$, we want to preserve face integrals, so for this everything stays simple in that
$
  cal(F)^2 : & hat(f) |-> f := (1/r hat(f)) compose F^(-1),
$

$
  cases(
    cal(F)^0 : hat(f) |-> f := hat(phi) compose F_(r z)^(-1),
    cal(F)^1 : hat(bold(v))_p |-> bold(v)_p := hat(bold(v))_p compose F_(r z)^(-1),
    cal(F)^2 : hat(f) |-> f := (hat(f)\/r) compose F_(r z)^(-1)
  )
$

so now we can define the physical spline basis
$
  cases(
    W_h^0 = "Span" ({T_(i j) := cal(F)^0 hat(T)_(i j) : 0 <= i <= n_s, 0 <= j <= n_z}),
    W_(p, h)^1 = "Span" ({bold(T)_(i j)^s := cal(F)^1 hat(bold(T))_(i j)^s : 0 <= i < n_s - 1, 0 <= j < n_z} union {bold(T)_(i j)^s := cal(F)^1 hat(bold(T))_(i j)^z : 0 <= i < n_s, 0 <= j < n_z - 1}),
    W_(theta, h)^2 = "Span" ({T_(i j)^2 := cal(F)^2 hat(T)_(i j)^2 : 0 <= i < n_s - 1, 0 <= j < n_z - 1})
  )
$

we have inverse operators called pullbacks
$
  cases(
    cal(B)^0 : phi |-> hat(phi) := phi compose F_(r z),
    cal(B)^1 : bold(v)_p |-> hat(bold(v))_p := bold(v) compose F_(r z),
    cal(B)^2 : f |-> hat(f) := r (f compose F_(r z))
  )
$
where again scalar fields the poloidal suvectors remain untouched but the flux is given a factor of $r$. If we try to us ethis formulation for the full field, it clearly shows itself to only work for the punctured domain and explodes if we try to include all of $Omega$. We fix this with discrete projection operators that map arbitrary tensor-product splines onto the conforming subspaces. These projections will be local, so only coefficients close to the pole will be modified, and commuting in the de Rham sequence. We will also have local conforming projections of the form
$
  P_Z^cal(l) : W_h^cal(l) -> Z_h^cal(l)
$

with the paper's context, we write the relations
$
  cases(
    V_h^0 = {phi = sum_(i, j) phi_(i j) T_(i, j)^0 in W_h^0 : phi_(0j) = gamma_j quad forall j},
    V_(p,h)^1 = {bold(v) = sum_(i, j) v_(i j)^s bold(T)_(i, j)^s + sum_(i,j) v_(i, j)^z bold(T)_(i,j)^z in W_(p,h)^1},
    V_theta^2 = {f = sum_(i,j) f_(i j) T_(i,j)^2 in W_(theta,h)^2 : f_(0j) = 0 quad forall j}
  )
$
we have the restriction that each point on the polar singularity is the same from every angle. this guarantees also that points between nodes share the same value so everything is consistent. since all our vectors are poloidal in this sequence then there's no issues with singularities. then for the final $theta$ flux through the center we make it 0.
