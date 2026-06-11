# now we need to represent simulation
# for this we need to represent the simulation fields
# simulation is defined by ns, nz, which control the number of nodes in each direction
# since ρ is a zero-form, we just use full dimensions for it
# however we want to only store reduced ρ, which has different dimensions
# ρ is a functional of a basis function that returns how much of the charge is contributed to by that function

import std/random

import short_names
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
    ρ: Matrix[settings.ns - 1, settings.nz, f64]

proc deposit_ρ[settings: static SimulationSettings](
  ρ: var Matrix[settings.ns - 1, settings.nz, f64];
  particles: seq[MacroParticle]
) =
  const s_spline = Spline(h: settings.spacing, n: settings.ns - 1, p: settings.base_spline_order)
  const z_spline = Spline(h: settings.spacing, n: settings.nz, p: settings.base_spline_order)

  for particle in particles:
    let s = particle.position.s
    let z = particle.position.z

    let s_offset = knot_offset[s_spline](s)
    let z_offset = knot_offset[z_spline](z)

    let s_spline_values = evaluate[s_spline](s)
    let z_spline_values = evaluate[z_spline](z)

    for s_idx in 0 ..< s_spline_values.len:
      for z_idx in 0 ..< z_spline_values.len:
        ρ[s_offset + s_idx, z_offset + z_idx] +=
          particle.constituents.f64 * particle.q *
          s_spline_values[s_idx] * z_spline_values[z_idx]

proc main() =
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

  "out.npy".write_npy_matrix_f64(simulation.ρ)

main()
