import numpy as np
import matplotlib.pyplot as plt

def U(r, ε, σ):
	return 4*ε*((σ/r)**12 - (σ/r)**6)

r_ray = np.linspace(0.9, 3, 100)
U_ray = U(r_ray, 1, 1)

plt.figure(figsize=(8, 6))
plt.plot(r_ray, U_ray)
plt.title("Lennard-Jones potential")
plt.xlabel("r [m]")
plt.ylabel("U [J]")
plt.savefig("figs/oppg1ai.png", dpi=400)
