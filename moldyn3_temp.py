import itertools
import time
import random
import math
from functools import reduce
from vector import Vector3
from moldyn import Atom

def deep_copy(list):
	return [item.copy() for item in list]

def periodic_boundry(between_vec, L):
	sign = lambda x: 1 if x > 0 else -1 if x < 0 else 0
	dx = between_vec.x
	dx = dx - sign(dx)*round(abs(dx)/L)*L
	dy = between_vec.y
	dy = dy - sign(dy)*round(abs(dy)/L)*L
	dz = between_vec.z
	dz = dz - sign(dz)*round(abs(dz)/L)*L

	direction_vec = Vector3(dx, dy, dz)
	r_sqrd = dx**2 + dy**2 + dz**2

	return direction_vec, r_sqrd

def U(r_sqrd):
	return 4*(r_sqrd**-6 - r_sqrd**-3) + 2912/531441

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

def get_force(between_vec, r_sqrd):
	direction_vec = between_vec/r_sqrd
	force = 24*(2*r_sqrd**-6 - r_sqrd**-3)

	return direction_vec*force

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
	# Update the atoms position using given method
	for atom in atoms:
		atom.update(dt, update_func)
		atom.save_state(datafile)

def simulate(atoms, dt, t_max, update_func, filename, L, completion, total_runs):
	t_list = [0]
	pot_list = []
	kin_list = []
	tot_list = []
	tmp_list = []
	start_atoms = deep_copy(atoms)

	datafile = open("data/"+filename+".xyz", "w")

	while t_list[-1] < t_max:
		if int(t_list[-1]/dt) % t_max == 0:
			print(f"\r{(100*completion + t_list[-1]*100/t_max)/total_runs:.0f} %", end="")

		pot, kin, tot = get_energy(atoms, L)
		temperature = (2/(3*len(atoms)))*kin

		pot_list.append(pot)
		kin_list.append(kin)
		tot_list.append(tot)
		tmp_list.append(temperature)
		step(atoms, dt, update_func, L, datafile)
		t_list.append(t_list[-1] + dt)

	datafile.close()

	atoms = start_atoms
	return t_list[:-1], pot_list, kin_list, tot_list, tmp_list

def main(name, simulate_count):
	dt = 0.01
	length = 5

	sum_lists = [[0]*(int(length/dt)+1) for i in range(4)]

	start_time = time.time()
	for i in range(simulate_count):
		L, atoms = create_atoms(108, 1.7, 180)

		t_list, pot_list, kin_list, tot_list, tmp_list = simulate(atoms, dt, length, "verlet", name, L, i, simulate_count)
		for i, (pot, kin, tot, temp) in enumerate(zip(pot_list, kin_list, tot_list, tmp_list)):
			sum_lists[0][i] += pot/simulate_count
			sum_lists[1][i] += kin/simulate_count
			sum_lists[2][i] += tot/simulate_count
			sum_lists[3][i] += temp/simulate_count
	print(f"\ntime: {time.time() - start_time:.3g}")

	with open("data/"+name+".energy", "w") as outfile:
		for t, pot, kin, tot, temp in zip(t_list, *sum_lists):
			outfile.write(f"{t} {pot} {kin} {tot} {temp}\n")

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


if __name__ == "__main__":
	main("data4a3", 10)
