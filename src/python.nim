import strutils

import short_names
import math

proc write_npy_2d_f64(path: string; data: openArray[f64]; rows, cols: int) =
  doAssert rows * cols == data.len
  let f = open(path, fmWrite)
  defer: f.close()

  const magic = "\x93NUMPY"

  f.write magic
  f.write '\x01'
  f.write '\x00'

  var header = "{'descr':'<f8', 'fortran_order': False, 'shape': (" &
    $rows & ", " & $cols & "), }"

  f.write char(header.len.uint16 and 0xFF)
  f.write char((header.len.uint16 shr 8) and 0xFF)

  let header_len = header.len

  header &= " ".repeat(64 - (magic.len + 4 + header_len + 1) mod 64)
  header &= '\n'

  f.write header

  discard f.writeBuffer(addr data[0], data.len * 8)

proc write_npy_matrix_f64*[rows, cols: static int](path: string, m: Matrix[rows, cols, f64]) =
  write_npy_2d_f64(path, m.data, rows, cols)
