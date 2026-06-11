type
  Vec2*[T] = object
    s*, z*: T
  Matrix*[rows, cols: static int; T] = object
    data*: array[rows * cols, T]

proc `[]`*[rows, cols: static int; T](m: var Matrix[rows, cols, T]; r, c: int): var T =
  m.data[r * cols + c]
