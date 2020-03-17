#include <algorithm>
#include "../moldyn_functions.hpp"

tuple<vector<double>, vector<double>> calc_g(vector<Atom> &atoms, double dr, double r_max, double L) {
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

	vector<double> r_vec;
	r_vec.reserve((int) r_max/dr);
	vector<double> g_vec;
	g_vec.reserve((int) r_max/dr);


	for (double r = dr; r <= r_max+dr; r += dr) {
		double avg_particles = 0;

		for (size_t i = 0; i < atom_length_lists.size(); i++) {
			// Get particles between r and r+dr
			int neighbour_count = 0;
			for (double length : atom_length_lists[i]) {
				if (length > r) {
					if (length > r + dr) {
						break;
					}
					neighbour_count += 1;
				}
			}
			avg_particles += neighbour_count/((double) atom_count);
		}

		double numerator = pow(L, 3) * avg_particles;
		double denominator = atom_count*4*M_PI*pow(r, 2)*dr;
		g_vec.push_back(numerator/denominator);
		r_vec.push_back(r);
	}

	return tuple<vector<double>, vector<double>>(r_vec, g_vec);
}


void run(string filename, int simulate_count, int timestep_count) {
	double dt = 0.01;
	int length = 50;

	int atom_count = 864;
	double d = 1.7;

	double r_max = ceil(d*cbrt(atom_count/4.0))*2;
	double dr = 0.01;

	vector<vector<double>> sum_lists(2, vector<double>((int) r_max/dr));
	mutex mutex;


	double start_time = get_time();

	#pragma omp parallel for
	for (int i = 0; i < simulate_count; i++) {
		FILE* datafile = fopen("data/null.xyz", "w");
		auto [L, atoms, atom_combinations] = create_equalibrium_atoms(atom_count, d, 180, dt, length);
		ez_simulate(atoms, atom_combinations, dt, length, L);

		mutex.lock();
		for (int i = 0; i < timestep_count; i++) {
			auto [r_vec, g_vec] = calc_g(atoms, dr, r_max, L);
			for (size_t j = 0; j < r_vec.size(); j++) {
				sum_lists[0][j] += r_vec[j]/(simulate_count*timestep_count);
				sum_lists[1][j] += g_vec[j]/(simulate_count*timestep_count);
			}
			step(atoms, atom_combinations, dt, L, datafile);
		}
		mutex.unlock();
	}

	FILE* datafile = fopen(("data/"+filename+".g").c_str(), "w");

	for (int i = 0; i < sum_lists[0].size(); i++) {
		fprintf(datafile, "%f %f\n", sum_lists[0][i], sum_lists[1][i]);
	}

	fclose(datafile);

	printf("\ntime: %.3g s\n", get_time() - start_time);
}

int main(int argc, char const *argv[]) {
	run("data4di", 16, 200);
	return 0;
}
