#include "../moldyn_functions.hpp"

double get_equalibrium(int simulate_count, float temperature) {
	double dt = 0.01;
	int length = 5;

	double equalibrium_temp = 0;
	mutex mutex;

	#pragma omp parallel for
	for (int i = 0; i < simulate_count; i++) {
		auto [L, atoms, atom_combinations] = create_equalibrium_atoms(108, 1.7, temperature, dt, length);
		auto [t_list, pot_list, kin_list, tot_list, tmp_list, vac_list, msd_list] = simulate(atoms, atom_combinations, dt, length, "null", L, i, simulate_count);

		mutex.lock();

		double avg = 0;
		for (size_t i = 0; i < tmp_list.size(); i++) {
			avg += tmp_list[i];
		}
		equalibrium_temp += (avg/tmp_list.size())/simulate_count;

		mutex.unlock();
	}

	return equalibrium_temp;
}

int main(int argc, char const *argv[]) {
	int Tmin = 80;
	int Tmax = 280;
	vector<double> start_temps;
	vector<double> equalibrium_temps;

	double start_time = get_time();
	for (double T = Tmin; T <= Tmax; T += 5) {
		cout << "\n" << T << ":" << "\n";
		start_temps.push_back(T);
		equalibrium_temps.push_back(get_equalibrium(512, T)*119.7);
	}
	printf("\ntime: %.3g s\n", get_time() - start_time);

	FILE* datafile = fopen("data/null", "w");
	for (size_t i = 0; i < equalibrium_temps.size(); i++) {
		fprintf(datafile, "%g -> %f\n", start_temps[i], equalibrium_temps[i]);
	}
	return 0;
}
