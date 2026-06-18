import std/random
import std/math

import short_names
import legendre_quadrature
import math
import spline
import python

type
  MacroParticle = object
    constituents: int
    q: f64
    position: Vec2[f64]
  SimulationSettings = object
    ns, nz: int
    base_spline_order: int
    spacing: f64
  Simulation[settings: static SimulationSettings] = ref object
    ρ: DenseMatrix[settings.ns - 1, settings.nz, f64]
  SpaceKind = enum
    kd_poloidal_0, kd_poloidal_1, kd_poloidal_2,
    kd_toroidal_1, kd_toroidal_2, kd_toroidal_3,
  ReducedMassMatrix[
    space: static SpaceKind,
    settings: static SimulationSettings
  ] = object
    # b_a = a factor for b component
    when space == kd_poloidal_0:
      s_factors: BandMatrix[settings.ns - 1, settings.base_spline_order + 1, f64]
      z_factors: BandMatrix[settings.nz, settings.base_spline_order + 1, f64]
    elif space == kd_poloidal_1:
      s_s_factors: BandMatrix[settings.ns - 2, settings.base_spline_order, f64]
      s_z_factors: BandMatrix[settings.nz, settings.base_spline_order + 1, f64]
      z_s_factors: BandMatrix[settings.ns - 1, settings.base_spline_order + 1, f64]
      z_z_factors: BandMatrix[settings.nz - 1, settings.base_spline_order, f64]
    elif space == kd_poloidal_2:
      s_factors: BandMatrix[settings.ns - 2, settings.base_spline_order, f64]
      z_factors: BandMatrix[settings.nz - 1, settings.base_spline_order, f64]
    elif space == kd_toroidal_1:
      s_factors: BandMatrix[settings.ns - 2, settings.base_spline_order + 1, f64]
      z_factors: BandMatrix[settings.nz, settings.base_spline_order + 1, f64]
    elif space == kd_toroidal_2:
      s_s_factors: BandMatrix[settings.ns - 2, settings.base_spline_order + 1, f64]
      s_z_factors: BandMatrix[settings.nz - 1, settings.base_spline_order, f64]
      z_s_factors: BandMatrix[settings.ns - 2, settings.base_spline_order, f64]
      z_z_factors: BandMatrix[settings.nz, settings.base_spline_order + 1, f64]
    else:
      s_factors: BandMatrix[settings.ns - 2, settings.base_spline_order, f64]
      z_factors: BandMatrix[settings.nz - 1, settings.base_spline_order, f64]
  ReducedMassMatrixRef[
    space: static SpaceKind,
    settings: static SimulationSettings
  ] = ref ReducedMassMatrix[space, settings]
  ReducedMassMatrixEntry[space: static SpaceKind] = object
    when space == kd_poloidal_1 or space == kd_toroidal_2:
      s, z: f64
    else:
      v: f64

proc reduced_mass_matrix_weight(space: SpaceKind; point: f64): f64 =
  case space:
  of kd_poloidal_0: point
  of kd_poloidal_1: point
  of kd_poloidal_2: point
  of kd_toroidal_1: 1'f64 / point
  of kd_toroidal_2: 1'f64 / point
  of kd_toroidal_3: 1'f64 / point

proc z_m_mass_matrix_factors[
  settings: static SimulationSettings
](): BandMatrix[settings.nz - 1, settings.base_spline_order, f64] =
  const z_spline = Spline(h: settings.spacing, n: settings.nz, p: settings.base_spline_order)
  for k in z_spline.p ..< settings.nz:
    let offset = knot[z_spline](k + 1)
    let knot_offset = k - z_spline.p

    for i in 0 ..< z_spline.p:
      let point = (z_spline.h * quadrature_root(z_spline.p, i)) / 2 + offset + z_spline.h / 2
      let weight = z_spline.h / 2 * quadrature_weight(z_spline.p, i)

      let values = evaluate_m_spline[z_spline](point)
      for a in 0 ..< z_spline.p:
        for b in a ..< z_spline.p:
          result[knot_offset + a, knot_offset + b] +=
            weight * values[a] * values[b]

proc z_b_mass_matrix_factors[
  settings: static SimulationSettings
](): BandMatrix[settings.nz, settings.base_spline_order + 1, f64] =
  const z_spline = Spline(h: settings.spacing, n: settings.nz, p: settings.base_spline_order)
  for k in z_spline.p ..< settings.nz:
    let offset = knot[z_spline](k + 1)
    let knot_offset = k - z_spline.p

    for i in 0 .. z_spline.p:
      let point = (z_spline.h * quadrature_root(z_spline.p + 1, i)) / 2 + offset + z_spline.h / 2
      let weight = z_spline.h / 2 * quadrature_weight(z_spline.p + 1, i)

      let values = evaluate_b_spline[z_spline](point)
      for a in 0 .. z_spline.p:
        for b in a .. z_spline.p:
          result[knot_offset + a, knot_offset + b] +=
            weight * values[a] * values[b]

proc s_rbeq_mass_matrix_factors[
  space: static SpaceKind;
  settings: static SimulationSettings
](): BandMatrix[settings.ns - 1, settings.base_spline_order + 1, f64] =
  const s_spline = Spline(h: settings.spacing, n: settings.ns, p: settings.base_spline_order)
  for k in s_spline.p ..< settings.ns:
    let offset = knot[s_spline](k + 1)
    let knot_offset = k - s_spline.p

    for i in 0 .. s_spline.p:
      let point = (s_spline.h * quadrature_root(s_spline.p + 1, i)) / 2 + offset + s_spline.h / 2
      let weight = s_spline.h / 2 * quadrature_weight(s_spline.p + 1, i) *
        reduced_mass_matrix_weight(space, point)

      let values = evaluate_b_spline[s_spline](point)
      if knot_offset == 0:
        var reduced_values: array[s_spline.p, f64]
        reduced_values[0] = values[0] + values[1]
        for i in 1 ..< reduced_values.len:
          reduced_values[i] = values[i + 1]

        for a in 0 ..< s_spline.p:
          for b in a ..< s_spline.p:
            result[a, b] += weight * reduced_values[a] *
              reduced_values[b]
      else:
        for a in 0 .. s_spline.p:
          for b in a .. s_spline.p:
            result[knot_offset + a - 1, knot_offset + b - 1] +=
              weight * values[a] * values[b]

proc s_rm0_mass_matrix_factors[
  space: static SpaceKind;
  settings: static SimulationSettings
](): BandMatrix[settings.ns - 2, settings.base_spline_order, f64] =
  const s_spline = Spline(h: settings.spacing, n: settings.ns, p: settings.base_spline_order)
  for k in s_spline.p ..< settings.ns:
    let offset = knot[s_spline](k + 1)
    let knot_offset = k - s_spline.p
    let start_idx = if knot_offset == 0: 1 else: 0

    for i in 0..<s_spline.p:
      let point = (s_spline.h * quadrature_root(s_spline.p, i)) / 2 + offset + s_spline.h / 2
      let weight = s_spline.h / 2 * quadrature_weight(s_spline.p, i) *
        reduced_mass_matrix_weight(space, point)

      let values = evaluate_m_spline[s_spline](point)
      for a in start_idx ..< s_spline.p:
        for b in a ..< s_spline.p:
          result[knot_offset + a - 1, knot_offset + b - 1] +=
            weight * values[a] * values[b]

proc s_rb00_mass_matrix_factors[
  space: static SpaceKind;
  settings: static SimulationSettings
](): BandMatrix[settings.ns - 2, settings.base_spline_order + 1, f64] =
  const s_spline = Spline(h: settings.spacing, n: settings.ns, p: settings.base_spline_order)
  for k in s_spline.p ..< settings.ns:
    let offset = knot[s_spline](k + 1)
    let knot_offset = k - s_spline.p
    let start_idx = if knot_offset < 2: 2 - knot_offset else: 0

    for i in 0 .. s_spline.p:
      let point = (s_spline.h * quadrature_root(s_spline.p + 1, i)) / 2 + offset + s_spline.h / 2
      let weight = s_spline.h / 2 * quadrature_weight(s_spline.p + 1, i) *
        reduced_mass_matrix_weight(space, point)

      let values = evaluate_b_spline[s_spline](point)
      for a in start_idx .. s_spline.p:
        for b in a .. s_spline.p:
          result[knot_offset + a - 2, knot_offset + b - 2] +=
            weight * values[a] * values[b]

proc build_reduced_mass_matrix(
  space: static SpaceKind;
  settings: static SimulationSettings
): ReducedMassMatrixRef[space, settings] =
  new result
  when space == kd_poloidal_0:
    result.s_factors = s_rbeq_mass_matrix_factors[space, settings]()
    result.z_factors = z_b_mass_matrix_factors[settings]()
  elif space == kd_poloidal_1:
    result.s_s_factors = s_rm0_mass_matrix_factors[space, settings]()
    result.s_z_factors = z_b_mass_matrix_factors[settings]()
    result.z_s_factors = s_rbeq_mass_matrix_factors[space, settings]()
    result.z_z_factors = z_m_mass_matrix_factors[settings]()
  elif space == kd_poloidal_2:
    result.s_factors = s_rm0_mass_matrix_factors[space, settings]()
    result.z_factors = z_m_mass_matrix_factors[settings]()
  elif space == kd_toroidal_1:
    result.s_factors = s_rb00_mass_matrix_factors[space, settings]()
    result.z_factors = z_b_mass_matrix_factors[settings]()
  elif space == kd_toroidal_2:
    result.s_s_factors = s_rb00_mass_matrix_factors[space, settings]()
    result.s_z_factors = z_m_mass_matrix_factors[settings]()
    result.z_s_factors = s_rm0_mass_matrix_factors[space, settings]()
    result.z_z_factors = z_b_mass_matrix_factors[settings]()
  else: # kd_toroidal_3
    result.s_factors = s_rm0_mass_matrix_factors[space, settings]()
    result.z_factors = z_m_mass_matrix_factors[settings]()

proc `[]`[
  space: static SpaceKind;
  settings: static SimulationSettings
](matrix: ReducedMassMatrixRef[space, settings]; μ, ν: (int, int)): ReducedMassMatrixEntry[space] =
  when space == kd_poloidal_1 or space == kd_toroidal_2:
    result.s = 2 * PI *
      matrix.s_s_factors[μ[0], ν[0]] * matrix.s_z_factors[μ[1], ν[1]]
    result.z = 2 * PI *
      matrix.z_s_factors[μ[0], ν[0]] * matrix.z_z_factors[μ[1], ν[1]]
  else:
    result.v = 2 * PI *
      matrix.s_factors[μ[0], ν[0]] * matrix.z_factors[μ[1], ν[1]]

# needs to be refactored/checked to properly handle the fact that it's a reduced representation
proc deposit_ρ[settings: static SimulationSettings](
  ρ: var DenseMatrix[settings.ns - 1, settings.nz, f64];
  particles: seq[MacroParticle]
) =
  const s_spline = Spline(h: settings.spacing, n: settings.ns - 1, p: settings.base_spline_order)
  const z_spline = Spline(h: settings.spacing, n: settings.nz, p: settings.base_spline_order)

  for particle in particles:
    let s = particle.position.s
    let z = particle.position.z

    let s_offset = knot_offset[s_spline](s)
    let z_offset = knot_offset[z_spline](z)

    let s_spline_values = evaluate_b_spline[s_spline](s)
    let z_spline_values = evaluate_b_spline[z_spline](z)

    for s_idx in 0 ..< s_spline_values.len:
      for z_idx in 0 ..< z_spline_values.len:
        ρ[s_offset + s_idx, z_offset + z_idx] +=
          particle.constituents.f64 * particle.q *
          s_spline_values[s_idx] * z_spline_values[z_idx]

proc ρ_test() =
  randomize()

  var particles: seq[MacroParticle] = @[]
  for _ in 1 .. 1_000:
    particles.add MacroParticle(
      constituents: 1,
      q: 1,
      position: Vec2[f64](
        s: rand(10'f64 .. 100'f64),
        z: rand(10'f64 .. 100'f64),
      )
    )

  const settings = SimulationSettings(
    ns: 200, nz: 200,
    base_spline_order: 5,
    spacing: 1
  )

  var simulation = Simulation[settings]()
  deposit_ρ[settings](simulation.ρ, particles)

  var sum = 0'f64
  for entry in simulation.ρ.data: sum += entry
  echo sum

  "out.npy".write_npy_dense_matrix_f64(simulation.ρ)

proc mass_matrix_test() =
  const space = kd_toroidal_1
  const settings = SimulationSettings(
    ns: 7, nz: 7,
    base_spline_order: 3,
    spacing: 1
  )
  let matrix = build_reduced_mass_matrix(space, settings)
  echo matrix[(2, 4), (3, 5)]

mass_matrix_test()
