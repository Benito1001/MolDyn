from vector import Vector3

class Atom:
	def __init__(self, pos, vel=None):
		self.force = Vector3()
		self.acc = Vector3()
		if vel is None:
			self.vel = Vector3()
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
		acc = self.force.copy()
		self.vel += 0.5*(acc_prev + acc)*dt
		self.pos += self.vel*dt + 0.5*acc*dt**2
		self.acc = acc
		self.force.set(0, 0, 0)

	def length_to(self, atom):
		return (self.pos - atom.pos).length

	def copy(self):
		return Atom(self.pos.copy())

	def __repr__(self):
		return str(self.pos)
