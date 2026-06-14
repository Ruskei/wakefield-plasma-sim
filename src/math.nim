type
  Vec2*[T] = object
    s*, z*: T
  DenseMatrix*[rows, cols: static int; T] = object
    data*: array[rows * cols, T]
  BandMatrix*[rows, bandwidth: static int; T] = object
    data*: array[rows * bandwidth, T]

proc `[]`*[rows, cols: static int; T](m: var DenseMatrix[rows, cols, T]; r, c: int): var T =
  m.data[r * cols + c]

proc `[]`*[rows, bandwidth: static int; T](m: var BandMatrix[rows, bandwidth, T], r, c: int): var T =
  let a = min(r, c)
  let b = max(r, c)
  m.data[a * bandwidth + b - a]
