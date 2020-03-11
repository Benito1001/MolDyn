import itertools
import time
import random
import math
from functools import reduce
from vector import Vector3

class Atom:
	def __init__(self, pos, vel=None):
		self.force = Vector3()
		self.acc = Vector3()
		if vel is None:
			self.vel = Vector3()
		else:
			self.vel = vel
		self.vel0 = vel
		self.pos = pos

	def update(self, dt, func="chromer"):
		if func == "chromer":
			self.update_chromer(dt)
		elif func == "euler":
			self.update_euler(dt)
		elif func == "verlet":
			self.update_verlet(dt)
		else:
			print(">:[")

	def update_euler(self, dt):
		"""
		r[i+1] = r[i] + v[i]*dt
		v[i+1] = v[i] + a[i]*dt
		"""
		acc = self.force
		self.pos += self.vel*dt
		self.vel += acc*dt
		self.force.set(0, 0, 0)

	def update_chromer(self, dt):
		"""
		v[i+1] = v[i] + a[i]*dt
		r[i+1] = r[i] + v[i+1]*dt
		"""
		acc = self.force
		self.vel += acc*dt
		self.pos += self.vel*dt
		self.force.set(0, 0, 0)

	def update_verlet(self, dt):
		"""
		v[i] = v[i-1] + 0.5*(a[i-1] + a[i])*dt
		r[i+1] = v[i]*dt + 0.5*a[i]*dt^2
		"""
		acc_prev = self.acc
		acc = self.force
		self.vel += 0.5*(acc_prev + acc)*dt
		self.pos += self.vel*dt + 0.5*acc*dt**2
		self.acc = acc
		self.force = Vector3()

	def save_state(self, file):
		file.write(f"Ar {self.pos.x:f} {self.pos.y:f} {self.pos.z:f}\n")

	def length_to(self, atom):
		return (self.pos - atom.pos).length

	def copy(self):
		return Atom(self.pos.copy())

	def __repr__(self):
		return str(self.pos)

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
	vac_list = []
	start_atoms = deep_copy(atoms)

	datafile = open("data/"+filename+".xyz", "w")

	while t_list[-1] < t_max:
		if int(t_list[-1]/dt) % t_max == 0:
			print(f"\r{(100*completion + t_list[-1]*100/t_max)/total_runs:.0f} %", end="")

		pot, kin, tot = get_energy(atoms, L)
		temperature = (2/(3*len(atoms)))*kin

		vel_autocor = 0
		for atom in atoms:
			vel_ac_emum = atom.vel.dot(atom.vel0)
			vel_ac_denom = atom.vel0.get_length_sqrd()
			vel_autocor += vel_ac_emum/vel_ac_denom
		vel_autocor = vel_autocor/len(atoms)

		pot_list.append(pot)
		kin_list.append(kin)
		tot_list.append(tot)
		tmp_list.append(temperature)
		vac_list.append(vel_autocor)
		step(atoms, dt, update_func, L, datafile)
		t_list.append(t_list[-1] + dt)

	datafile.close()

	atoms = start_atoms
	return t_list[:-1], pot_list, kin_list, tot_list, tmp_list, vac_list

def main(name, simulate_count):
	dt = 0.01
	length = 5

	sum_lists = [[0]*(int(length/dt)+1) for i in range(5)]

	start_time = time.time()
	for i in range(simulate_count):
		L, atoms = create_atoms(108, 1.7, 180)

		t_list, pot_list, kin_list, tot_list, tmp_list, vac_list = simulate(atoms, dt, length, "verlet", name, L, i, simulate_count)
		for i, (pot, kin, tot, temp, vac) in enumerate(zip(pot_list, kin_list, tot_list, tmp_list, vac_list)):
			sum_lists[0][i] += pot/simulate_count
			sum_lists[1][i] += kin/simulate_count
			sum_lists[2][i] += tot/simulate_count
			sum_lists[3][i] += temp/simulate_count
			sum_lists[4][i] += vac/simulate_count
	print(f"\ntime: {time.time() - start_time:.3g}")

	with open("data/"+name+".energy", "w") as outfile:
		for t, pot, kin, tot, temp, vac in zip(t_list, *sum_lists):
			outfile.write(f"{t} {pot} {kin} {tot} {temp} {vac}\n")

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
	main("data4b1", 10)
