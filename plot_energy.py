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
				tmp_list.append(temp_vac_msd[0]*119.7)

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

def plotinator_comparinator(datafiles, draw_list=["energy", "temp", "vac"], size=(10, 8)):
	fig, axs = plt.subplots(len(draw_list), 1, figsize=(size[0], size[1]), sharex=True, gridspec_kw={'hspace': 0.2})

	for datafile in datafiles:
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
					tmp_list.append(temp_vac_msd[0]*119.7)

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

		if not hasattr(axs, "__getitem__"):
			axs = [axs]
		for i, type in enumerate(draw_list):
			data, title, xlabel, ylabel = data_list[type]
			for name, color, values in data:
				axs[i].plot(t_list, values, label=name + " ("+datafile+")")
			axs[i].set_title(title)
			axs[i].set_ylabel(ylabel)
			if i == len(draw_list)-1:
				axs[i].set_xlabel(xlabel)
		for ax in axs:
			ax.legend()

	plt.show()

# plot_energy("data3eii.energy")
plotinator_comparinator(["null.energy", "data4ci.energy", "data4ci_cpp.energy"], ["energy", "temp", "vac", "msd"], (16, 32))
# plotinator("null.energy", ["energy", "temp", "vac", "msd"], (8, 12))
