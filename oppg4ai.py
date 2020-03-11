import time
from oppg3eii import deep_copy, get_energy, step, create_atoms

def simulate(atoms, dt, t_max, filename, update_func, L):
	# Declare variables and save atoms
	t_list = [0]
	pot_list = []
	kin_list = []
	tot_list = []
	tmp_list = []
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

		# Calculate and store temperature
		temperature = (2/(3*len(atoms)))*kin
		tmp_list.append(temperature)

		step(atoms, dt, update_func, L, datafile)
		t_list.append(t_list[-1] + dt)
	print(f"\n{update_func:>7}: {time.time()- start_time:.3g} s")

	datafile.close()

	atoms = start_atoms
	return t_list[:-1], pot_list, kin_list, tot_list, tmp_list

def main(filename):
	dt = 0.01
	length = 5

	L, atoms = create_atoms(108, 1.7)

	t_list, pot_list, kin_list, tot_list, tmp_list = simulate(atoms, dt, length, filename, "verlet", L)
	with open("data/"+filename+".energy", "w") as outfile:
		for t, pot, kin, tot, temp in zip(t_list, pot_list, kin_list, tot_list, tmp_list):
			outfile.write(f"{t} {pot} {kin} {tot} {temp}\n")

if __name__ == "__main__":
	main("data4ai")
