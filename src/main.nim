import strutils
import random

type
  f64 = float64
  Vec2[T] = object
    r, z: T
  Vec2D = Vec2[f64]
  Vec3[T] = object
    r, θ, z: T
  Vec3D = Vec3[f64]
  Matrix[NR, NZ: static int; T] = object
    data: array[NR * NZ, T]
  MatrixRef[NR, NZ: static int; T] = ref Matrix[NR, NZ, T]
  SimulationSettings = object
    NR: int # number of nodes in r
    NZ: int # number of nodes in z
    P: int # b-spline order
    H: f64 # node spacing
  Simulation[NR, NZ: static int] = object
    # free field coefficients
    E_r: Matrix[NR - 2, NZ, f64]
    rE_θ: Matrix[NR - 2, NZ, f64]
    E_z: Matrix[NR - 1, NZ - 1, f64]

    rB_r: Matrix[NR - 2, NZ - 1, f64]
    B_θ: Matrix[NR - 2, NZ - 1, f64]
    rB_z: Matrix[NR - 2, NZ, f64]
  ParticleSpecies = object
    q, m: f64
  MacroParticle = object
    species: ParticleSpecies
    n: f64
    p: Vec2D

const ε = 1e-11

proc idx[NR, NZ: static int; T](m: MatrixRef[NR, NZ, T]; r_idx, z_idx: int): int =
  NZ * r_idx + z_idx

proc `[]`[NR, NZ: static int; T](m: MatrixRef[NR, NZ, T]; r_idx, z_idx: int): var T =
  m.data[m.idx(r_idx, z_idx)]

proc `[]=`[NR, NZ: static int; T](m: MatrixRef[NR, NZ, T]; r_idx, z_idx: int; v: T) =
  m[m.idx(r_idx, z_idx)] = v

proc knot[N, P: static int](i: int; h: f64): f64 =
  if (i < P + 1): result = 0
  elif (i >= N): result = h * (N - 1)
  else: result = h * (i - P).float64

proc find_span[N, P: static int](h, x: f64): int =
  ## Open uniform knot vector, N basis functions, degree P.
  ## Returns span k such that xi[k] <= x < xi[k+1].
  ## Valid returned range: P .. N-1.
  when N <= P:
    {.error: "Need N > P".}

  let xMax = h * f64(N - P)

  if x <= 0.0:
    return P
  if x >= xMax:
    return N - 1

  result = P + int(x / h) # truncates toward zero; x > 0 here

  if result < P:
    result = P
  elif result > N - 1:
    result = N - 1

proc b_spline_local[N, P: static int](h, x: f64): tuple[first: int, vals: array[P + 1, f64]] =
  var vals: array[P + 1, f64]
  let xMax = h * f64(N - P)

  if x < 0.0 or x > xMax:
    return (0, vals) # all zero

  let span = find_span[N, P](h, x)
  var left, right: array[P + 1, f64]

  vals[0] = 1.0
  for j in 1 .. P:
    left[j] = x - knot[N, P](span + 1 - j, h)
    right[j] = knot[N, P](span + j, h) - x

    var saved = 0.0
    for r in 0 ..< j:
      let denom = right[r + 1] + left[j - r]
      let temp = if denom != 0.0: vals[r] / denom else: 0.0

      vals[r] = saved + right[r + 1] * temp
      saved = left[j - r] * temp

    vals[j] = saved

  return (span - P, vals)

proc ρ_from_particles[S: static SimulationSettings](
  particles: openArray[MacroParticle]
): MatrixRef[S.NR, S.NZ, f64] =
  var ρ = MatrixRef[S.NR, S.NZ, f64]()
  for particle in particles:
    let (r_first, r_vals) = b_spline_local[S.NR, S.P](S.H, particle.p.r)
    let (z_first, z_vals) = b_spline_local[S.NZ, S.P](S.H, particle.p.z)

    for rk in 0 .. S.P:
      let ri = r_first + rk
      if (ri < 0 or ri >= S.NR): continue
      if (r_vals[rk] < ε): continue
      for zk in 0 .. S.P:
        let zi = z_first + zk
        if (zi < 0 or zi >= S.NZ): continue
        if (z_vals[zk] < ε): continue
        ρ[ri, zi] += r_vals[rk] * z_vals[zk]
  result = ρ

proc write_npy_f64_2d(path: string; data: openArray[f64]; nx, ny: int) =
  assert data.len == nx * ny
  let f = open(path, fmWrite)
  defer: f.close()

  let magic = "\x93NUMPY"
  f.write(magic)
  f.write(char(1))
  f.write(char(0))

  var header = "{'descr': '<f8', 'fortran_order': False, 'shape': (" &
               $ny & ", " & $nx & "), }"
  let preamble_len = 10
  let pad_len = (16 - ((preamble_len + header.len + 1) mod 16)) mod 16
  header &= repeat(' ', pad_len) & "\n"
  let hlen = uint16(header.len)
  f.write(char(hlen and 0xff))
  f.write(char((hlen shr 8) and 0xff))

  f.write (header)

  for v in data:
    var x = v
    discard f.writeBuffer(addr x, sizeof(f64))



proc main() =
  randomize()

  const S = SimulationSettings(
    NR: 518,
    NZ: 518,
    P: 3,
    H: 1
  )
  
  let e_species = ParticleSpecies(q: 1, m: 1)
  var particles: seq[MacroParticle] = @[]

  for _ in 0 ..< 1_000_000:
    particles.add MacroParticle(
      species: e_species,
      n: 1,
      p: Vec2D(r: rand(0.0 .. 512.0), z: rand(0.0 .. 512.0))
    )

  let ρ = ρ_from_particles[S](particles)
  var s = 0.0
  for r in 0 ..< S.NR:
    for z in 0 ..< S.NZ:
      s += ρ[r, z]
  echo s

  write_npy_f64_2d("rho.npy", ρ.data, S.NR, S.NZ)

main()
