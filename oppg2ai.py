import itertools

def deep_copy(list):
	return [item.copy() for item in list]

def get_force(atom1, atom2):
	between_vec = atom1.pos - atom2.pos
	r_sqrd = between_vec.get_length_sqrd()

	direction_vec = between_vec/r_sqrd
	force = 24*(2*r_sqrd**-6 - r_sqrd**-3)

	return direction_vec*force

def step(atoms, dt, update_func):
	# Add the force acting on the particles efficiently using pairs
	for atom1, atom2 in itertools.combinations(atoms, 2):
		force = get_force(atom1, atom2)
		atom1.force += force
		atom2.force -= force

	# Update the atoms position using given method
	for atom in atoms:
		atom.update(dt, update_func)

def simulate(atoms, dt, t_max, update_func):
	t_list = [0]
	r_list = []
	start_atoms = deep_copy(atoms)

	while t_list[-1] < t_max:
		r_list.append(atoms[0].length_to(atoms[1]))
		step(atoms, dt, update_func)
		t_list.append(t_list[-1] + dt)

	atoms = start_atoms
	return t_list[:-1], r_list
