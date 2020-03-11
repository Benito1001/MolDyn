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

	def update(self, dt, L, func="chromer"):
		# Update position and velocity
		if func == "chromer":
			self.update_chromer(dt)
		elif func == "euler":
			self.update_euler(dt)
		elif func == "verlet":
			self.update_verlet(dt)
		else:
			print(">:[")

		# Teleport if outside of box
		if self.pos.x > L:
			self.pos.x = self.pos.x - L
		elif self.pos.x < 0:
			self.pos.x = L - self.pos.x

		if self.pos.y > L:
			self.pos.y = self.pos.y - L
		elif self.pos.y < 0:
			self.pos.y = L - self.pos.y

		if self.pos.z > L:
			self.pos.z = self.pos.z - L
		elif self.pos.z < 0:
			self.pos.z = L - self.pos.z

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
