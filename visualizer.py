import numpy as np
import matplotlib.pyplot as plt

data = np.load("rmm2.npy")

plt.imshow(
    data,
    origin="lower",
    aspect="auto",
    cmap="plasma",
    interpolation="nearest"
)
plt.show()
