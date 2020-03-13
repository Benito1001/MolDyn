import matplotlib.pyplot as plt
from oppg2ai import simulate
from oppg2ai_Atom import Atom
from vector import Vector3

def main(distance, name):
	dt = 0.01
	length = 5

	atoms = [
		Atom(Vector3(0, 0, 0)),
		Atom(Vector3(distance, 0, 0))
	]

	t_list, r_list = simulate(atoms, dt, length, "chromer")
	plt.figure(figsize=(8, 6))
	plt.plot(t_list, r_list, color="#11966e", label="distance")
	plt.legend()

	min_r, max_r = min(r_list), max(r_list)
	padding = (max_r-min_r)/5
	plt.axis([0, length, min_r-padding, max_r+padding])

	plt.title(f"Two atom simulation, distance={distance}$\\sigma$")
	plt.xlabel("t$'$")
	plt.ylabel("r$'$")
	plt.savefig(f"figs/{name}.png", dpi=400)
main(1.5, "oppg2bii")
