import matplotlib.pyplot as plt
import itertools
from functools import reduce
from oppg2ai_Atom import Atom
from vector import Vector3
from oppg2ai import deep_copy, step

def U(r_sqrd):
	return 4*(r_sqrd**-6 - r_sqrd**-3)

def get_energy(atoms):
	# Calculate potential energy:
	potential_energy = 0
	for atom1, atom2 in itertools.combinations(atoms, 2):
		between_vec = atom1.pos - atom2.pos
		potential_energy += U(between_vec.get_length_sqrd())

	# Calculate kinetic energy:
	kinetic_energy = reduce(lambda accumelator, atom: accumelator + 0.5*atom.vel.get_length_sqrd(), atoms, 0)

	return potential_energy, kinetic_energy, potential_energy+kinetic_energy

def simulate(atoms, dt, t_max, update_func):
	t_list = [0]
	kin_list, pot_list, tot_list = ([] for i in range(3))
	start_atoms = deep_copy(atoms)

	while t_list[-1] < t_max:
		pot, kin, tot = get_energy(atoms)

		kin_list.append(kin)
		pot_list.append(pot)
		tot_list.append(tot)

		step(atoms, dt, update_func)
		t_list.append(t_list[-1] + dt)

	atoms = start_atoms
	return t_list[:-1], kin_list, pot_list, tot_list


def main(distance, ax):
	dt = 0.01
	length = 5

	atoms = [
		Atom(Vector3(0, 0, 0)),
		Atom(Vector3(distance, 0, 0))
	]

	t_list, kin_list, pot_list, tot_list = simulate(atoms, dt, length, "chromer")
	data_list = [
		("kinetic", "r", kin_list),
		("potential", "g", pot_list),
		("total", "k", tot_list),
	]
	for name, color, values in data_list:
		ax.plot(t_list, values, color=color, label=name)
	ax.legend()

	min_value = min(pot_list)
	max_value = max(kin_list)
	buffer = (max_value - min_value)/5
	ax.axis([0, length, min_value-buffer, max_value+buffer])
	ax.set_title(f"Two atom simulation, distance={distance}$\\sigma$")
	ax.set_xlabel("t$'$")
	ax.set_ylabel("r$'$")

fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
main(1.5, axs[0])
main(0.95, axs[1])
axs[0].set_xlabel("")
plt.savefig("figs/oppg2ci.png", dpi=400)
