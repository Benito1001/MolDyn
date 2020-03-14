#include "../moldyn_functions.hpp"

vector<double> calc_g(double dr, double r_max, double L) {
	for (Atom &atom1 : atoms) {
		vector<double> neighbours;
		neighbours.reserve(atoms.size());
		for (Atom &atom2 : atoms) {
			if (atom1 != atom2) {
				neighbours.push_back(sqrt(get_length_sqrd(atom1.pos-atom2.pos)));
			}
		}
		sort(neighbours);
	}

	vector<double> r_ray;
	r_vec.reserve((int) r_max/dr);
	vector<double> g_vec;
	g_vec.reserve((int) r_max/dr);
	for (double r = 0; r < r_max; r += dr) {
		double avg_neighbours = 0;
		for (Atom &atom : atoms) {
			int n = 0;
			for (double length : atom.neighbours) {
				if (length > r) {
					if (length > r + dr) {
						break;
					}
					n += 1;
				}
			}
			avg_neighbours += n/((double) atoms.size());
		}
		double numerator = pow(L, 3) * avg_neighbours;
		double denominator = atoms.size()*4*pi*pow(r, 2)*dr;
		g_vec.push_back(numerator/denominator);
		r_vec.push_back(r);
	}

	return r_vec, g_vec;
}


void run(string filename, int simulate_count) {
	double dt = 0.01;
	int length = 5;

	int atom_count = 864;
	double d = 1.7;

	double start_time = get_time();

	auto [L, atoms, atom_combinations] = create_equalibrium_atoms(atom_count, d, 180, dt, length);
	auto [_t_list, pot_list, kin_list, tot_list, tmp_list, vac_list, msd_list] = simulate(atoms, atom_combinations, dt, length, filename, L, i, simulate_count);

	r_vec, g_vec = calc_g(L);

	FILE* datafile = fopen(("data/"+filename+".g").c_str(), "w");

	for (int i = 0; i < r_vec.size(); i++) {
		fprintf(datafile, "%f %f\n", r_vec[i], g_vec[i]);
	}

	fclose(datafile);

	printf("\ntime: %.3g s\n", get_time() - start_time);
}

int main(int argc, char const *argv[]) {
	run("data4di", 100)
	return 0;
}
