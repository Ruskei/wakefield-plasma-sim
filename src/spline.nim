import std/strutils
import std/macros
import std/math

import short_names

type
  Spline* = object
    ## these variables are as defined in [notes_sangalli.pdf]
    ## more accurately, this represents a spline space, as both
    ## B-splines and M-splines are evaluated with this object
    h*: f64
    n*, p*: int

template low(spline: Spline): f64 = 0
template high(spline: Spline): f64 = (spline.n - spline.p).f64 * spline.h

# 1-indexed
proc knot*[spline: static Spline](i: int): f64 =
  max(min(((i - spline.p - 1).f64 * spline.h), high(spline)), low(spline))

proc knot_offset*[spline: static Spline](x: f64): int =
  floor(x / spline.h).int + 1

# not handling p=0 case rn
# not clamping input
# x is expected to be of type f64
# result is expected to be of type array[p + 1, f64]
macro generate_b_spline_evaluation(spline: static Spline; p: static int; x: typed; res: typed): untyped =
  do_assert spline.p >= 0
  do_assert spline.n > spline.p
  do_assert p in 0 .. spline.p

  result = new_stmt_list()

  let low_lit = new_lit(low(spline))
  let high_lit = new_lit(high(spline))

  let uniform_statement_list = new_stmt_list()

  # generates code for the simple case in the middle of the spline space
  block:
    proc multiply_linear(poly: seq[f64], constant, u_factor: f64): seq[f64] =
      result = new_seq[f64](poly.len + 1)
      for i, coefficient in poly:
        result[i] += constant * coefficient
        result[i + 1] += u_factor * coefficient

    proc add_poly(a, b: seq[f64]): seq[f64] =
      result = new_seq[f64](max(a.len, b.len))
      for i, coefficient in a:
        result[i] += coefficient
      for i, coefficient in b:
        result[i] += coefficient

    let u_node = ident("u")
    let scaled_x_node = ident("scaled_x")
    let inv_h_value = new_lit(1'f64 / spline.h)

    let scaled_x_initialization = quote do:
      let `scaled_x_node` = `x` * `inv_h_value`
    uniform_statement_list.add scaled_x_initialization

    let u_initialization = quote do:
      let `u_node` = `scaled_x_node` - floor(`scaled_x_node`)

    uniform_statement_list.add u_initialization

    var polynomials = new_seq[seq[f64]](p + 2)
    polynomials[p + 1] = @[1'f64]

    for l in 1 .. p:
      block:
        let source = polynomials[p + 2 - l]
        polynomials[p + 1 - l] = multiply_linear(source, 1'f64, -1'f64)

      for j in (p - l + 2) .. p:
        let left = multiply_linear(polynomials[j], (p + 1 - j).f64, 1'f64)
        let right = multiply_linear(polynomials[j + 1], (j + l - p).f64, -1'f64)
        polynomials[j] = add_poly(left, right)

      block:
        let source = polynomials[p + 1]
        polynomials[p + 1] = multiply_linear(source, 0'f64, 1'f64)

    let final_factor = new_lit(1 / fac(p))
    
    for i in 1 .. (p + 1):
      let index = new_lit(i - 1)
      var coefficients = polynomials[i]
      for coefficient in coefficients.mitems:
        coefficient *= final_factor.float_val

      var polynomial_expression = new_lit(coefficients[^1])
      if coefficients.len > 1:
        for coefficient_index in countdown(coefficients.high - 1, 0):
          let coefficient = new_lit(coefficients[coefficient_index])
          polynomial_expression = quote do:
            `polynomial_expression` * `u_node` + `coefficient`

      let stmt = quote do:
        `res`[`index`] = `polynomial_expression`
      uniform_statement_list.add stmt

  let boundary_statement_list = new_stmt_list()

  # generate intersecting boundary case
  block:
    let k_identifier = ident("k")
    let h_value = new_lit(spline.h)
    let spline_p_value = new_lit(spline.p)
    let k_initialization = quote do:
      let `k_identifier` = floor(`x` / `h_value`).int + `spline_p_value` + 1
    boundary_statement_list.add k_initialization

    for i in 2 .. (2 * p + 1):
      let knot_identifier = ident("knot_" & $i)
      let offset_lit = new_lit(p + spline.p + 2 - i)
      let knot_initialization = quote do:
        let `knot_identifier` = min(max((`k_identifier` - (`offset_lit`)).f64 * `h_value`, `low_lit`), `high_lit`)
      boundary_statement_list.add knot_initialization

    for i in 2 .. (p + 1):
      let coefficient_identifier = ident("left_coefficient_" & $i)
      let coefficient_value = ident("knot_" & $i)
      let coefficient_initialization = quote do:
        let `coefficient_identifier` = `x` - `coefficient_value`
      boundary_statement_list.add coefficient_initialization

    for i in (p + 2) .. (2 * p + 1):
      let coefficient_identifier = ident("right_coefficient_" & $i)
      let coefficient_value = ident("knot_" & $i)
      let coefficient_initialization = quote do:
        let `coefficient_identifier` = `coefficient_value` - `x`
      boundary_statement_list.add coefficient_initialization

    let constant_identifier = ident("b_$#-0" % $(p + 1))
    let constant_initialization = quote do:
      const `constant_identifier` = 1'f64
    boundary_statement_list.add constant_initialization

    for l in 1 .. p:
      for i in (p - l + 2) .. (p + 1):
        let denom_identifier = ident("d_$#-$#" % [$i, $l])
        let right_identifier = ident("right_coefficient_" & $(i + l))
        let left_identifier = ident("left_coefficient_" & $i)
        let denom_initialization = quote do:
          let `denom_identifier` = `right_identifier` + `left_identifier`
        boundary_statement_list.add denom_initialization

      block:
        let b_identifier = ident("b_$#-$#" % [$(p + 1 - l), $l])
        let numer_identifier = ident("right_coefficient_" & $(p + 2))
        let denom_identifier = ident("d_$#-$#" % [$(p - l + 2), $l])
        let factor_identifier = ident("b_$#-$#" % [$(p + 2 - l), $(l - 1)])
        let b_statement = quote do:
          let `b_identifier` = `numer_identifier` / `denom_identifier` * `factor_identifier`
        boundary_statement_list.add b_statement

      for j in (p - l + 2) .. p:
        let b_identifier = ident("b_$#-$#" % [$j, $l])
        let left_numer_identifier = ident("left_coefficient_" & $j)
        let left_denom_identifier = ident("d_$#-$#" % [$j, $l])
        let left_factor_identifier = ident("b_$#-$#" % [$j, $(l - 1)])
        let right_numer_identifier = ident("right_coefficient_" & $(j + l + 1))
        let right_denom_identifier = ident("d_$#-$#" % [$(j + 1), $l])
        let right_factor_identifier = ident("b_$#-$#" % [$(j + 1), $(l - 1)])
        let b_statement = quote do:
          let `b_identifier` = `left_numer_identifier` / `left_denom_identifier` * `left_factor_identifier` +
            `right_numer_identifier` / `right_denom_identifier` * `right_factor_identifier`
        boundary_statement_list.add b_statement

      block:
        let b_identifier = ident("b_$#-$#" % [$(p + 1), $l])
        let numer_identifier = ident("left_coefficient_" & $(p + 1))
        let denom_identifier = ident("d_$#-$#" % [$(p + 1), $l])
        let factor_identifier = ident("b_$#-$#" % [$(p + 1), $(l - 1)])
        let b_statement = quote do:
          let `b_identifier` = `numer_identifier` / `denom_identifier` * `factor_identifier`
        boundary_statement_list.add b_statement

    for i in 1 .. (p + 1):
      let index = new_lit(i - 1)
      let b = ident("b_$#-$#" % [$i, $p])
      let stmt = quote do:
        `res`[`index`] = `b`
      boundary_statement_list.add stmt

  # boundaries for where spline evaluation is uniform
  let low_uniform = new_lit((p.f64 - 1'f64) * spline.h)
  let high_uniform = new_lit((spline.n - spline.p - p + 1).f64 * spline.h)

  do_assert low_uniform.float_val <= high_uniform.float_val

  let if_statement = quote do:
    if likely(`x` == min(max(`x`, `low_uniform`), `high_uniform`)):
      `uniform_statement_list`
    elif `x` == min(max(`x` , `low_lit`), `high_lit`):
      `boundary_statement_list`
  result.add if_statement

proc evaluate_b_spline*[spline: static Spline](x: f64): array[spline.p + 1, f64] =
  generate_b_spline_evaluation(spline, spline.p, x, result)

proc evaluate_m_spline*[spline: static Spline](x: f64): array[spline.p + 0, f64] =
  generate_b_spline_evaluation(spline, spline.p - 1, x, result)
