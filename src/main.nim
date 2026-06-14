import std/random

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
  MultiIndex = Vec2[int]
  SpaceKind = enum
    kd_poloidal_0, kd_poloidal_1, kd_poloidal_2,
    kd_toroidal_1, kd_toroidal_2, kd_toroidal_3,
  ReducedMassMatrix[
    space: static SpaceKind,
    settings: static SimulationSettings
  ] = object
    # b_a = a factor for b component
    when space == kd_poloidal_1:
      s_s_factors: BandMatrix[settings.ns - 2, settings.base_spline_order, f64]
      s_z_factors: BandMatrix[settings.nz, settings.base_spline_order, f64]
    else:
      discard

proc build_reduced_mass_matrix(
  space: static SpaceKind,
  settings: static SimulationSettings
): ReducedMassMatrix[space, settings] =
  when space == kd_poloidal_1:
    const s_spline = Spline(h: settings.spacing, n: settings.ns, p: settings.base_spline_order)
    # loop over knot-spans with nonzero length
    for k in s_spline.p ..< settings.ns:
      let offset = knot[s_spline](k + 1)
      let knot_offset = k - s_spline.p
      let start_idx = if knot_offset == 0: 1 else: 0
      for i in 0..<s_spline.p:
        let point = (s_spline.h * quadrature_root(s_spline.p, i)) / 2 + offset + s_spline.h / 2
        let weight = s_spline.h / 2 * quadrature_weight(s_spline.p, i)

        let values = evaluate_m_spline[s_spline](point)
        for a in start_idx ..< s_spline.p:
          for b in a ..< s_spline.p:
            result.s_s_factors[knot_offset + a - 1, knot_offset + b - 1] +=
              weight * values[a] * values[b] * point

    const z_spline = Spline(h: settings.spacing, n: settings.nz, p: settings.base_spline_order)
    for k in z_spline.p ..< settings.nz:
      let offset = knot[z_spline](k + 1)
      let knot_offset = k - z_spline.p
      for i in 0..z_spline.p:
        let point = (z_spline.h * quadrature_root(z_spline.p + 1, i)) / 2 + offset + z_spline.h / 2
        let weight = z_spline.h / 2 * quadrature_weight(z_spline.p + 1, i)

        let values = evaluate_b_spline[z_spline](point)
        for a in 0 .. z_spline.p:
          for b in a .. z_spline.p:
            result.s_z_factors[knot_offset + a, knot_offset + b] +=
              weight * values[a] * values[b]
  else:
    {.fatal: "TODO".}

# to get mass matrix elements at runtime we need to store
# the subparts at runtime. 
# proc reduced_mass_matrix_poloidal_1_s_entry[
#   settings: static SimulationSettings
# ](μ, ν: MultiIndex): f64 =
#
#

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
  const space = kd_poloidal_1
  const settings = SimulationSettings(
    ns: 512, nz: 512,
    base_spline_order: 3,
    spacing: 1
  )
  const matrix = build_reduced_mass_matrix(space, settings)
  "rmm.npy".write_npy_band_matrix_f64(matrix.s_s_factors)

mass_matrix_test()
