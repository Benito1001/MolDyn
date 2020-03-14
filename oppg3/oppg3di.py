import itertools
from functools import reduce
import time
from vector import Vector3
from oppg2di_Atom import Atom

def U(r_sqrd):
	return 4*(r_sqrd**-6 - r_sqrd**-3) + 2912/531441

def get_force(between_vec, r_sqrd):
	direction_vec = between_vec/r_sqrd
	force = 24*(2*r_sqrd**-6 - r_sqrd**-3)

	return direction_vec*force

def get_energy(atoms):
	# Calculate potential energy:
	potential_energy = 0
	for atom1, atom2 in itertools.combinations(atoms, 2):
		between_vec = atom1.pos - atom2.pos
		r_sqrd = between_vec.get_length_sqrd()
		if r_sqrd < 3*3:
			potential_energy += U(r_sqrd)

	# Calculate kinetic energy:
	kinetic_energy = reduce(lambda accumelator, atom: accumelator + 0.5*atom.vel.get_length_sqrd(), atoms, 0)

	return potential_energy, kinetic_energy, potential_energy+kinetic_energy


def step(atoms, dt, update_func, datafile):
	# Add the force acting on the particles efficiently using pairs
	for atom1, atom2 in itertools.combinations(atoms, 2):
		between_vec = atom1.pos - atom2.pos
		r_sqrd = between_vec.get_length_sqrd()

		if r_sqrd < 3*3:
			force = get_force(between_vec, r_sqrd)
			atom1.force += force
			atom2.force -= force

	datafile.write(f"{len(atoms)}\ntype x y z\n")
	for atom in atoms:
		# Save current atom positions to file
		atom.save_state(datafile)
		# Update atom positions using given method
		atom.update(dt, update_func)

def simulate(atoms, dt, t_max, update_func, filename):
	# Declare variables
	t_list = [0]
	pot_list = []
	kin_list = []
	tot_list = []

	datafile = open("data/"+filename+".xyz", "w")

	start_time = time.time()
	while t_list[-1] < t_max:
		# Fancy progress indicator
		if int(t_list[-1]/dt) % t_max == 0:
			print(f"\r{t_list[-1]*100/t_max:.0f} %", end="")

		# Get and store energy
		pot, kin, tot = get_energy(atoms)
		pot_list.append(pot)
		kin_list.append(kin)
		tot_list.append(tot)

		step(atoms, dt, update_func, datafile)
		t_list.append(t_list[-1] + dt)

	# Print total elapsed time
	print(f"\n{update_func:>7}: {time.time()- start_time:.3g} s")

	datafile.close()

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

def main(filename):
	dt = 0.01
	length = 5

	atoms = []
	for position in box_positions(4, 1.7):
		atoms.append(Atom(position))

	t_list, pot_list, kin_list, tot_list = simulate(atoms, dt, length, "verlet", filename)
	with open("data/"+filename+".energy", "w") as outfile:
		for t, pot, kin, tot in zip(t_list, pot_list, kin_list, tot_list):
			outfile.write(f"{t} {pot} {kin} {tot}\n")
main("data3di")
