import time
from oppg3eii import step, get_energy, deep_copy
from oppg4aii import create_atoms
from oppg4biii import ez_simulate

def simulate(atoms, dt, t_max, filename, update_func, L, completion, total_runs):
	# Declare variables and save atoms
	t_list = [0]
	pot_list = []
	kin_list = []
	tot_list = []
	tmp_list = []
	vac_list = []
	msd_list = []
	start_atoms = deep_copy(atoms)

	datafile = open("data/"+filename+".xyz", "w")

	while t_list[-1] < t_max:
		# Fancy progress indicator, now even more complicated
		if int(t_list[-1]/dt) % t_max == 0:
			print(f"\r{' '*7}{completion+1}/{total_runs} : {(100*completion + t_list[-1]*100/t_max)/total_runs:3.0f} %", end="")

		# Get and store energy
		pot, kin, tot = get_energy(atoms, L)
		pot_list.append(pot)
		kin_list.append(kin)
		tot_list.append(tot)

		# Calculate and store temperature
		temperature = (2/(3*len(atoms)))*kin
		tmp_list.append(temperature)

		# Calculate and store velocity autocorrelation
		vac = 0
		for atom in atoms:
			vac_emum = atom.vel.dot(atom.vel0)
			vac_denom = atom.vel0.get_length_sqrd()
			vac += vac_emum/vac_denom
		vac = vac/len(atoms)
		vac_list.append(vac)

		# Calculate and store mean squared displacement
		msd = sum((atom.dist_traveled(L) - atom.pos0).get_length_sqrd() for atom in atoms)/len(atoms)
		msd_list.append(msd)

		step(atoms, dt, update_func, L, datafile)
		t_list.append(t_list[-1] + dt)

	datafile.close()

	atoms = start_atoms
	return t_list[:-1], pot_list, kin_list, tot_list, tmp_list, vac_list, msd_list

def create_equalibrium_atoms(atom_count, d, temperature, dt, t_max, update_func):
	L, start_atoms = create_atoms(atom_count, d, temperature)
	warm_atoms = ez_simulate(start_atoms, dt, t_max, update_func, L)
	for atom in warm_atoms:
		atom.vel0 = atom.vel.copy()
		atom.pos0 = atom.pos.copy()

	return L, warm_atoms

def main(filename, simulate_count):
	dt = 0.01
	length = 5

	sum_lists = [[0]*(int(length/dt)+1) for i in range(6)]

	start_time = time.time()
	for i in range(simulate_count):
		L, atoms = create_equalibrium_atoms(864, 1.7, 180, dt, length, "verlet")

		list_list = simulate(atoms, dt, length, filename, "verlet", L, i, simulate_count)
		t_list, pot_list, kin_list, tot_list, tmp_list, vac_list, msd_list = list_list
		for i, (pot, kin, tot, temp, vac, msd) in enumerate(zip(pot_list, kin_list, tot_list, tmp_list, vac_list, msd_list)):
			sum_lists[0][i] += pot/simulate_count
			sum_lists[1][i] += kin/simulate_count
			sum_lists[2][i] += tot/simulate_count
			sum_lists[3][i] += temp/simulate_count
			sum_lists[4][i] += vac/simulate_count
			sum_lists[5][i] += msd/simulate_count
	print(f"\ntime: {time.time() - start_time:.3g}")

	with open("data/"+filename+".energy", "w") as outfile:
		for t, pot, kin, tot, temp, vac, msd in zip(t_list, *sum_lists):
			outfile.write(f"{t} {pot} {kin} {tot} {temp} {vac} {msd}\n")

if __name__ == '__main__':
	main("null", 1)
	pass
