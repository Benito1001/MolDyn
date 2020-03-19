#include <unordered_map>
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

double run(double temperature, int simulate_count) {
	double dt = 0.01;
	int length = 5;

	vector<double> avg_msd_list((int) (length/dt), 0);
	vector<double> t_list;
	mutex mutex;

	#pragma omp parallel for
	for (int i = 0; i < simulate_count; i++) {
		auto [L, atoms, atom_combinations] = create_equalibrium_atoms(864, 1.7, temperature, dt, length);

		auto [_t_list, pot_list, kin_list, tot_list, tmp_list, vac_list, msd_list] = simulate(atoms, atom_combinations, dt, length, "null", L, i, simulate_count);
		mutex.lock();
		t_list = move(_t_list);
		for (size_t i = 0; i < avg_msd_list.size(); i++) {
			avg_msd_list[i] += msd_list[i]/simulate_count;
		}
		mutex.unlock();
	}

	// Find derivative
	vector<double> diff_list(avg_msd_list.size()-1);
	for (size_t i = 0; i < diff_list.size(); i++) {
		diff_list[i] = (avg_msd_list[i+1] - avg_msd_list[i])/(t_list[i+1] - t_list[i]);
	}

	// Calculate average
	double D = 0;
	for (size_t i = 0; i < diff_list.size(); i++) {
		D += diff_list[i]/diff_list.size();
	}

	return D/6.0;
}

int main(int argc, char const *argv[]) {
	int Tmin = 85;
	int Tmax = 220;
	vector<int> T_list;
	for (int T = Tmin; T <= Tmax; T += 5) {
		T_list.push_back(T);
	}

	unordered_map<int, string> T2eqT;
	string line;
	ifstream temp_map_file("data/temp2equalibrium.dat");
	while (getline(temp_map_file, line)) {
		vector<string> split_line = split(line, " -> ");
		T2eqT[stoi(split_line[0])] = split_line[1];
	}
	temp_map_file.close();

	double start_time = get_time();
	vector<double> D_list;
	for (int T : T_list) {
		cout << "\n" << T << ":\n";
		D_list.push_back(run(T, 288));
	}
	printf("\ntime: %.3g s\n", get_time() - start_time);

	FILE* datafile = fopen("data/oppg4ciii.dat", "w");
	for (int i = 0; i < D_list.size(); i++) {
		fprintf(datafile, "%s %f\n", T2eqT[T_list[i]].c_str(), D_list[i]);
	}
	fclose(datafile);

	return 0;
}
