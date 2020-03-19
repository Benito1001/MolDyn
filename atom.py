from vector import Vector3

class Atom:
	def __init__(self, pos, vel=None):
		self.force = Vector3()
		self.acc = Vector3()
		if vel is None:
			self.vel = Vector3()
		else:
			self.vel = vel
		self.vel0 = self.vel.copy()
		self.pos = pos
		self.pos0 = self.pos.copy()

		self.teleports = Vector3()

	def update(self, dt, L, func):
		# Update position and velocity
		self.update_verlet(dt)

		# Teleport if outside of box
		if self.pos.x > L:
			self.pos.x = self.pos.x - L
			self.teleports.x += 1
		elif self.pos.x < 0:
			self.pos.x = self.pos.x + L
			self.teleports.x -= 1

		if self.pos.y > L:
			self.pos.y = self.pos.y - L
			self.teleports.y += 1
		elif self.pos.y < 0:
			self.pos.y = self.pos.y + L
			self.teleports.y -= 1

		if self.pos.z > L:
			self.pos.z = self.pos.z - L
			self.teleports.z += 1
		elif self.pos.z < 0:
			self.pos.z = self.pos.z + L
			self.teleports.z -= 1

	def dist_traveled(self, L):
		return self.pos + self.teleports*L

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
