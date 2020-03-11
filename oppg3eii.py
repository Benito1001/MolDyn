import itertools
import time
import random
from functools import reduce
from vector import Vector3
from atom import Atom

def deep_copy(list):
	return [item.copy() for item in list]

def U(r_sqrd):
	return 4*(r_sqrd**-6 - r_sqrd**-3) + 2912/531441

def get_force(between_vec, r_sqrd):
	direction_vec = between_vec/r_sqrd
	force = 24*(2*r_sqrd**-6 - r_sqrd**-3)

	return direction_vec*force

def periodic_boundry(between_vec, L):
	dx = between_vec.x
	dx = dx - round(dx/L)*L
	dy = between_vec.y
	dy = dy - round(dy/L)*L
	dz = between_vec.z
	dz = dz - round(dz/L)*L

	direction_vec = Vector3(dx, dy, dz)
	r_sqrd = dx**2 + dy**2 + dz**2

	return direction_vec, r_sqrd

def get_energy(atoms, L):
	# Calculate potential energy:
	potential_energy = 0
	for atom1, atom2 in itertools.combinations(atoms, 2):
		between_vec = atom1.pos - atom2.pos
		direction_vec, r_sqrd = periodic_boundry(between_vec, L)
		potential_energy += U(r_sqrd)

	# Calculate kinetic energy:
	kinetic_energy = reduce(lambda accumelator, atom: accumelator + 0.5*atom.vel.get_length_sqrd(), atoms, 0)

	return potential_energy, kinetic_energy, potential_energy+kinetic_energy


def step(atoms, dt, update_func, L, datafile):
	# Add the force acting on the particles efficiently using pairs
	for atom1, atom2 in itertools.combinations(atoms, 2):
		between_vec = atom1.pos - atom2.pos
		direction_vec, r_sqrd = periodic_boundry(between_vec, L)

		if r_sqrd < 3*3:
			force = get_force(direction_vec, r_sqrd)
			atom1.force += force
			atom2.force -= force

	datafile.write(f"{len(atoms)}\ntype x y z\n")
	for atom in atoms:
		# Save current atom positions to file
		atom.save_state(datafile)
		# Update atom positions using given method
		atom.update(dt, L, update_func)

def simulate(atoms, dt, t_max, filename, update_func, L):
	# Declare variables and save atoms
	t_list = [0]
	pot_list = []
	kin_list = []
	tot_list = []
	start_atoms = deep_copy(atoms)

	datafile = open("data/"+filename+".xyz", "w")

	start_time = time.time()
	while t_list[-1] < t_max:
		# Fancy progress indicator
		if int(t_list[-1]/dt) % t_max == 0:
			print(f"\r{t_list[-1]*100/t_max:.0f} %", end="")

		# Get and store energy
		pot, kin, tot = get_energy(atoms, L)
		pot_list.append(pot)
		kin_list.append(kin)
		tot_list.append(tot)

		step(atoms, dt, update_func, L, datafile)
		t_list.append(t_list[-1] + dt)

	# Print total elapsed time
	print(f"\n{update_func:>7}: {time.time()- start_time:.3g} s")

	# Save final atom positions
	for atom in atoms:
		atom.save_state(datafile)

	datafile.close()

	atoms = start_atoms
	return t_list[:-1], pot_list, kin_list, tot_list


def box_positions(n, d):
	positions = []
	for i in range(0, n):
		for j in range(0, n):
			for k in range(0, n):
				positions.append(Vector3(i, j, k)*d)
				positions.append(Vector3(i, 0.5 + j, 0.5 + k)*d)
				positions.append(Vector3(0.5 + i, j, 0.5 + k)*d)
				positions.append(Vector3(0.5 + i, 0.5 + j, k)*d)
	return positions

def create_atoms(atom_count, d):
	n = int((atom_count**(1/3))/(4**(1/3)))
	L = d*n

	random.seed(1)
	atoms = []
	for position in box_positions(n, d):
		atoms.append(Atom(position, (Vector3(random.random(), random.random(), random.random())*2)-1))

	return L, atoms


def main(filename):
	dt = 0.01
	length = 5

	L, atoms = create_atoms(108, 1.7)

	t_list, pot_list, kin_list, tot_list = simulate(atoms, dt, length, filename, "verlet", L)
	with open("data/"+filename+".energy", "w") as outfile:
		for t, pot, kin, tot in zip(t_list, pot_list, kin_list, tot_list):
			outfile.write(f"{t} {pot} {kin} {tot}\n")

if __name__ == "__main__":
	main("data3eii")
