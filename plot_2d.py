import numpy as np
import matplotlib.pyplot as plt

rho = np.load("rho.npy")

plt.figure()
plt.imshow(
    rho,
    origin="lower",
    aspect="auto",
    interpolation="nearest",
)
plt.colorbar(label="rho")
plt.xlabel("x index")
plt.ylabel("y indey")
plt.title("rho")
plt.show()
