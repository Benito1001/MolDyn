#include "../moldyn_functions.hpp"

vector<string> split(string s, string delimiter) {
	vector<string> elements;
	size_t pos = 0;
	string token;
	while ((pos = s.find(delimiter)) != string::npos) {
		token = s.substr(0, pos);
		elements.push_back(token);
		s.erase(0, pos + delimiter.length());
	}
	elements.push_back(s);

	return elements;
}

void run(string filename, int simulate_count) {
	double dt = 0.01;
	int length = 5;

	vector<vector<double>> sum_lists(6, vector<double>((int) (length/dt)));
	mutex mutex;
	vector<double> t_list;

	double start_time = get_time();
	#pragma omp parallel for
	for (int i = 0; i < simulate_count; i++) {
		auto [L, atoms, atom_combinations] = create_equalibrium_atoms(864, 1.7, 180, dt, length);

		auto [_t_list, pot_list, kin_list, tot_list, tmp_list, vac_list, msd_list] = simulate(atoms, atom_combinations, dt, length, filename, L, i, simulate_count);
		mutex.lock();
		t_list = move(_t_list);
		for (int j = 0; j < t_list.size(); j++) {
			sum_lists[0][j] += pot_list[j]/((double) simulate_count);
			sum_lists[1][j] += kin_list[j]/((double) simulate_count);
			sum_lists[2][j] += tot_list[j]/((double) simulate_count);
			sum_lists[3][j] += tmp_list[j]/((double) simulate_count);
			sum_lists[4][j] += vac_list[j]/((double) simulate_count);
			sum_lists[5][j] += msd_list[j]/((double) simulate_count);
		}
		mutex.unlock();
	}
	printf("\ntime: %.3g s\n", get_time() - start_time);

	FILE* datafile = fopen(("data/"+filename+".energy").c_str(), "w");

	for (int i = 0; i < t_list.size(); i++) {
		fprintf(datafile, "%f %f %f %f %f %f %f\n", t_list[i], sum_lists[0][i], sum_lists[1][i], sum_lists[2][i], sum_lists[3][i], sum_lists[4][i], sum_lists[5][i]);
	}
	fclose(datafile);
}

int main(int argc, char const *argv[]) {
	run("null", 1);

	return 0;
}

// SANIC: -O3 -march=native -fopenmp -ffast-math

// Preformance improvements: (864 atoms, i7-7700K)
// python          :~1120 s
// pypy            : 112 s
// default + -03   : 29 s
// + -march=native : 28.2 s
// + -ffast-math   : 4.41 s
// + -fopenmp      : 0.88 (simulate_count = cpu threads)
