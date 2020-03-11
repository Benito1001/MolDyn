import numpy as np
import matplotlib.pyplot as plt

def plotinator(datafile, draw_list=["energy", "temp", "vac"], size=(10, 8)):
	t_list, pot_list, kin_list, tot_list, tmp_list, vac_list, msd_list = [[] for i in range(7)]
	with open("data/"+datafile, "r") as infile:
		for line in infile:
			t, pot, kin, tot, *temp_vac_msd = (float(v) for v in line.split())
			t_list.append(t)
			pot_list.append(pot)
			kin_list.append(kin)
			tot_list.append(tot)
			if len(temp_vac_msd) > 2:
				msd_list.append(temp_vac_msd[2])
			if len(temp_vac_msd) > 1:
				vac_list.append(temp_vac_msd[1])
			if len(temp_vac_msd) > 0:
				tmp_list.append(temp_vac_msd[0])

	data_list = {
		"energy": [[
			("kinetic", "r", kin_list),
			("potential", "g", pot_list),
			("total", "k", tot_list),
		], "Energy plot", "Time", "Energy"],
		"temp": [[
			("temperature", "r", tmp_list),
		], "Temperature plot", "Time", "Temperature [K]"],
		"vac": [[
			("velocity autocorrelation", "g", vac_list),
		], "Velocity autocorrelation plot", "Time", "Velocity autocorrelation"],
		"msd": [[
			("mean squared displacemen", "b", msd_list),
		], "Mean squared displacemen plot", "Time", "Mean squared displacemen"],
	}
	fig, axs = plt.subplots(len(draw_list), 1, figsize=(size[0], size[1]), sharex=True, gridspec_kw={'hspace': 0.2})
	if not hasattr(axs, "__getitem__"):
		axs = [axs]
	for i, type in enumerate(draw_list):
		data, title, xlabel, ylabel = data_list[type]
		for name, color, values in data:
			axs[i].plot(t_list, values, color=color, label=name)
		axs[i].set_title(title)
		axs[i].set_ylabel(ylabel)
		if i == len(draw_list)-1:
			axs[i].set_xlabel(xlabel)
	for ax in axs:
		ax.legend()

	plt.show()

def plot_vac(energy_file):
	t_list, pot_list, kin_list, tot_list, tmp_list, vac_list = [[] for i in range(6)]
	with open("data/"+energy_file, "r") as infile:
		for line in infile:
			t, pot, kin, tot, temp, vac = (float(v) for v in line.split())
			t_list.append(t)
			pot_list.append(pot)
			kin_list.append(kin)
			tot_list.append(tot)
			tmp_list.append(temp*119.7)
			vac_list.append(vac)

	plt.plot(t_list, vac_list)
	plt.title("Velocity autocorrelation plot")
	plt.show()

def plot_energy_temp_vac(energy_file):
	t_list, pot_list, kin_list, tot_list, tmp_list, vac_list = [[] for i in range(6)]
	with open("data/"+energy_file, "r") as infile:
		for line in infile:
			t, pot, kin, tot, temp, vac = (float(v) for v in line.split())
			t_list.append(t)
			pot_list.append(pot)
			kin_list.append(kin)
			tot_list.append(tot)
			tmp_list.append(temp*119.7)
			vac_list.append(vac)

	data_list = [
		("kinetic", "r", kin_list, 0),
		("potential", "g", pot_list, 0),
		("total", "k", tot_list, 0),
		("temperature", "r", tmp_list, 1),
		("velocity autocorrelation", "r", vac_list, 2)
	]
	fig, axs = plt.subplots(3, 1, figsize=(8, 12), sharex=True, gridspec_kw={'hspace': 0.11})
	for name, color, values, index in data_list:
		axs[index].plot(t_list, values, color=color, label=name)
	for ax in axs:
		ax.legend()
	axs[0].set_title("Energy plot")
	axs[0].set_ylabel("Energy")
	axs[1].set_title("Temperature plot")
	axs[1].set_ylabel("Temperature [K]")
	axs[2].set_title("Velocity autocorrelation plot")
	axs[2].set_ylabel("Velocity autocorrelation")
	axs[2].set_xlabel("Time")
	plt.show()

def plot_energy_temp(energy_file):
	t_list, pot_list, kin_list, tot_list, tmp_list = [[] for i in range(5)]
	with open("data/"+energy_file, "r") as infile:
		for line in infile:
			t, pot, kin, tot, temp = (float(v) for v in line.split())
			t_list.append(t)
			pot_list.append(pot)
			kin_list.append(kin)
			tot_list.append(tot)
			tmp_list.append(temp*119.7)

	print(f"Equilibrium ~{np.average(tmp_list[int(len(tmp_list)/2):]):.1f} Kelvin")
	data_list = [
		("kinetic", "r", kin_list, 0),
		("potential", "g", pot_list, 0),
		("total", "k", tot_list, 0),
		("temperature", "r", tmp_list, 1)
	]
	fig, axs = plt.subplots(2, 1, figsize=(8, 12), sharex=True, gridspec_kw={'hspace': 0.1})
	for name, color, values, index in data_list:
		axs[index].plot(t_list, values, color=color, label=name)
	for ax in axs:
		ax.legend()
	axs[0].set_title("Energy plot")
	axs[0].set_ylabel("Energy")
	axs[1].set_title("Temperature plot")
	axs[1].set_ylabel("Temperature [K]")
	axs[1].set_xlabel("Time")
	plt.show()

def plot_energy(energy_file):
	t_list, pot_list, kin_list, tot_list = [[] for i in range(4)]
	with open("data/"+energy_file, "r") as infile:
		for line in infile:
			t, pot, kin, tot = (float(v) for v in line.split())
			t_list.append(t)
			pot_list.append(pot)
			kin_list.append(kin)
			tot_list.append(tot)

	data_list = [
		("kinetic", "r", kin_list),
		("potential", "g", pot_list),
		("total", "k", tot_list),
	]
	plt.figure(figsize=(8, 6))
	for name, color, values in data_list:
		plt.plot(t_list, values, color=color, label=name)
	plt.legend()
	plt.xlabel("Time")
	plt.ylabel("Energy")
	plt.title("Energy plot")
	plt.show()

# plot_energy("data3eii.energy")
plotinator("data4ci.energy", ["msd"], (8, 6))
