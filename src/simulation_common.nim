import short_names

type
  SpaceKind* = enum
    kd_poloidal_0, kd_poloidal_1, kd_poloidal_2,
    kd_toroidal_1, kd_toroidal_2, kd_toroidal_3,
  SimulationSettings* = object
    ns*, nz*: int
    base_spline_order*: int
    spacing*: f64
