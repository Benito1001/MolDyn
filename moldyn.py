import itertools
import time
# import matplotlib.pyplot as plt
from vector import Vector3

class Atom:
	def __init__(self, pos, vel=None):
		self.force = Vector3()
		self.acc = Vector3()
		if vel is None:
			self.vel = Vector3()
		else:
			self.vel = vel
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

atoms = [
	Atom(Vector3(0, 0, 0)),
	Atom(Vector3(1.5, 0, 0))
]

def deep_copy(list):
	return [item.copy() for item in list]

def get_force(atom1, atom2):
	between_vec = atom1.pos - atom2.pos
	r_sqrd = between_vec.get_length_sqrd()

	direction_vec = between_vec/r_sqrd
	force = 24*(2*r_sqrd**-6 - r_sqrd**-3)

	return direction_vec*force

def step(dt, update_func, datafile=None):
	global atoms

	# Add the force acting on the particles efficiently using pairs
	for atom1, atom2 in itertools.combinations(atoms, 2):
		force = get_force(atom1, atom2)
		atom1.force += force
		atom2.force -= force

	if datafile is not None:
		datafile.write(f"{len(atoms)}\ntype x y z\n")
	# Update the atoms position using given method
	for atom in atoms:
		atom.update(dt, update_func)
		if datafile is not None:
			atom.save_state(datafile)

def simulate(dt, t_max, update_func):
	global atoms

	t_list = [0]
	r_list = []
	start_atoms = deep_copy(atoms)

	datafile = open("data.xyz", "w")

	start_time = time.time()
	while t_list[-1] < t_max:
		r_list.append(atoms[0].length_to(atoms[1]))
		step(dt, update_func, datafile)
		t_list.append(t_list[-1] + dt)
	print(f"{update_func:>7}: {time.time()- start_time:.3g} s")

	datafile.close()

	atoms = start_atoms
	print(min(r_list))
	return t_list[:-1], r_list


def main():
	dt = 0.001
	length = 5
	plt.plot(*simulate(dt, length, "chromer"), color="k", label="chromer")
	plt.legend()
	plt.axis([0, length, 0, 3])
	plt.show()

if __name__ == "__main__":
	main()
