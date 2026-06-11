import std/math

import short_names

type
  Spline* = object
    ## these variables are as defined in [notes_sangalli.pdf]
    h*: f64
    n*, p*: int

template low(spline: static Spline): f64 = 0
template high(spline: static Spline): f64 = (spline.n - spline.p).f64 * spline.h

proc knot[spline: static Spline](i: int): f64 =
  ((i - spline.p - 1).f64 * spline.h).clamp(low(spline), high(spline))

proc knot_offset*[spline: static Spline](x: f64): int =
  floor(x / spline.h).int + 1

proc evaluate*[spline: static Spline](x: f64): array[spline.p + 1, f64] =
  if x notin low(spline) .. high(spline):
    return

  let knot_offset = knot_offset[spline](x)
  result[spline.p] = 1
  for i in 1 .. spline.p:
    let knot_ip1 = knot[spline](knot_offset + spline.p + 1)
    let knot_i1 = knot[spline](knot_offset + spline.p - i + 1)

    result[spline.p - i] = result[spline.p - i + 1] *
      (knot_ip1 - x) /
      (knot_ip1 - knot_i1)
    
    for j in (spline.p - i + 1) .. (spline.p - 1):
      let knot_j = knot[spline](knot_offset + j)
      let knot_jp = knot[spline](knot_offset + j + i)

      let knot_jp1 = knot[spline](knot_offset + j + i + 1)
      let knot_j1 = knot[spline](knot_offset + j + 1)

      result[j] = result[j] *
        (x - knot_j) /
        (knot_jp - knot_j) +
        result[j + 1] *
        (knot_jp1 - x) /
        (knot_jp1 - knot_j1)

    let knot_k = knot[spline](knot_offset + spline.p)
    let knot_kp = knot[spline](knot_offset + spline.p + i)

    result[spline.p] *= (x - knot_k) /
      (knot_kp - knot_k)
