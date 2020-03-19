#include <algorithm>
#include "../moldyn_functions.hpp"

vector<double> calc_g(vector<double> &r_list, double dr, vector<Atom> &atoms, double L) {
	size_t atom_count = atoms.size();
	vector<vector<double>> atom_length_lists(atom_count, vector<double>(atom_count - 1));

	for (size_t i = 0; i < atom_count; i++) {
		vector<double> length_list;
		length_list.reserve(atom_count);

		for (size_t j = 0; j < atom_count; j++) {
			if (i != j) {
				length_list.push_back(sqrt(get_length_sqrd(atoms[i].pos-atoms[j].pos)));
			}
		}

		sort(length_list.begin(), length_list.end());
		atom_length_lists[i] = move(length_list);
	}

	// cout << "\n";
	// for (vector<double> length_list : atom_length_lists) {
	// 	for (double length : length_list) {
	// 		cout << length << "  ";
	// 	}
	// 	cout << "\n";
	// }

	vector<double> g_list;
	g_list.reserve(r_list.size());

	for (double r : r_list) {
		double avg_particles = 0;

		for (size_t i = 0; i < atom_length_lists.size(); i++) {
			// Get particles between r and r+dr
			int neighbour_count = 0;
			for (double length : atom_length_lists[i]) {
				if (length > r) {
					if (length > r + dr) {
						atom_length_lists[i].erase(atom_length_lists[i].begin(), atom_length_lists[i].begin()+neighbour_count);
						break;
					}
					neighbour_count += 1;
				}
			}
			avg_particles += neighbour_count/((double) atom_count);
		}

		double numerator = pow(L, 3) * avg_particles;
		double denominator = atom_count*4*M_PI*pow(r, 2)*dr;
		g_list.push_back(numerator/denominator);
	}

	return g_list;
}


void run(string filename, int simulate_count, int timestep_count) {
	double dt = 0.01;
	int length = 50;

	int atom_count = 1372;
	double d = 1.7;

	double r_max = ceil(d*cbrt(atom_count/4.0))*1.5;
	double dr = 0.01;
	vector<double> r_list;
	for (double r = dr; r <= r_max+dr; r += dr) {
		r_list.push_back(r);
	}

	vector<double> avg_avg_g_list(r_list.size());
	mutex mutex;

	double start_time = get_time();

	#pragma omp parallel for
	for (int i = 0; i < simulate_count; i++) {
		FILE* datafile = fopen("data/null.xyz", "w");
		auto [L, atoms, atom_combinations] = create_equalibrium_atoms(atom_count, d, 180, dt, length);
		simulate(atoms, atom_combinations, dt, length, "null", L, i, simulate_count);

		vector<double> avg_g_list(r_list.size());

		for (int j = 0; j < timestep_count; j++) {
			vector<double> g_list = calc_g(r_list, dr, atoms, L);
			printf("\r                  %i/%i", j+1, timestep_count);
			fflush(stdout);
			for (size_t k = 0; k < r_list.size(); k++) {
				avg_g_list[k] += g_list[k]/timestep_count;
			}
			step(atoms, atom_combinations, dt, L, datafile);
		}
		mutex.lock();
		for (size_t j = 0; j < r_list.size(); j++) {
			avg_avg_g_list[j] += avg_g_list[j]/simulate_count;
		}
		mutex.unlock();
	}

	FILE* datafile = fopen(("data/"+filename+".g").c_str(), "w");

	for (size_t i = 0; i < r_list.size(); i++) {
		fprintf(datafile, "%f %f\n", r_list[i], avg_avg_g_list[i]);
	}

	fclose(datafile);

	printf("\ntime: %.3g s\n", get_time() - start_time);
}

int main(int argc, char const *argv[]) {
	run("null", 96, 100);
	return 0;
}
