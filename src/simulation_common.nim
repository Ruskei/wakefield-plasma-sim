import math
import short_names

type
  SpaceKind* = enum
    kd_poloidal_0, kd_poloidal_1, kd_poloidal_2,
    kd_toroidal_1, kd_toroidal_2, kd_toroidal_3,
  SimulationSettings* = object
    ns*, nz*: int
    base_spline_order*: int
    spacing*: f64
  FieldDimensions*[space: static SpaceKind] = object
    when space == kd_poloidal_1 or space == kd_toroidal_2:
      s_s*, s_z*, z_s*, z_z*: int
    else:
      s*, z*: int

proc field_dimensions*[
  space: static SpaceKind,
  settings: static SimulationSettings
](): FieldDimensions[space] =
  when space == kd_poloidal_0:
    result.s = settings.ns - 1
    result.z = settings.nz
  elif space == kd_poloidal_1:
    result.s_s = settings.ns - 2
    result.s_z = settings.nz
    result.z_s = settings.ns - 1
    result.z_z = settings.nz - 1
  elif space == kd_poloidal_2:
    result.s = settings.ns - 2
    result.z = settings.nz - 1
  elif space == kd_toroidal_1:
    result.s = settings.ns - 2
    result.z = settings.nz
  elif space == kd_toroidal_2:
    result.s_s = settings.ns - 2
    result.s_z = settings.nz - 1
    result.z_s = settings.ns - 2
    result.z_z = settings.nz
  else:
    result.s = settings.ns - 2
    result.z = settings.nz - 1

type
  Field*[
    space: static SpaceKind,
    settings: static SimulationSettings
  ] = object
    when space == kd_poloidal_1 or space == kd_toroidal_2:
      s*: DenseMatrix[
        field_dimensions[space, settings].s_s,
        field_dimensions[space, settings].s_z,
        f64]
      z*: DenseMatrix[
        field_dimensions[space, settings].z_s,
        field_dimensions[space, settings].z_z,
        f64]
    else:
      v*: DenseMatrix[
        field_dimensions[space, settings].s,
        field_dimensions[space, settings].z,
        f64]
