import numpy as np
import matplotlib.pyplot as plt

data = np.load("out.npy")

plt.imshow(
    data,
    origin="lower",
    aspect="auto",
    cmap="plasma"
)
plt.show()
