import random
import math
from vector import Vector3
from atom import Atom
from oppg3eii import box_positions
from oppg4ai import simulate

def create_atoms(atom_count, d, temperature):
	n = int((atom_count**(1/3))/(4**(1/3)))
	L = d*n

	atoms = []
	normal = lambda: random.gauss(0, math.sqrt(temperature/119.7))
	positions = box_positions(n, d)
	velocities = [Vector3(normal(), normal(), normal()) for i in range(atom_count)]
	for position, velocity in zip(positions, velocities):
		atoms.append(Atom(position, velocity))

	return L, atoms

def main(filename):
	dt = 0.01
	length = 5

	L, atoms = create_atoms(108, 1.7, 300)

	t_list, pot_list, kin_list, tot_list, tmp_list = simulate(atoms, dt, length, filename, "verlet", L)
	with open("data/"+filename+".energy", "w") as outfile:
		for t, pot, kin, tot, temp in zip(t_list, pot_list, kin_list, tot_list, tmp_list):
			outfile.write(f"{t} {pot} {kin} {tot} {temp}\n")

if __name__ == "__main__":
	main("data4aii")
