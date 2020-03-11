#include <iostream>
#include <chrono>
#include <string>
#include <random>
#include <cmath>
#include <tuple>
#include <list>
#include <vector>
#include <valarray>
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>

using namespace std;
using vec3 = glm::vec<3, double>;

class Atom;
void step(vector<Atom> &atoms, vector<vector<Atom*>> &atom_combinations, double dt, double L, FILE *datafile);


class Atom {
	public:
		vec3 force;
		vec3 acc;
		vec3 vel;
		vec3 vel0;
		vec3 pos;
		vec3 pos0;
		vec3 teleports;
		Atom(vec3 pos, vec3 vel=vec3 {0, 0, 0}) {
			this->pos = pos;
			this->pos0 = pos;
			this->vel = vel;
			this->vel0 = vel;
			this->force = (vec3) {0, 0, 0};
			this->acc = (vec3) {0, 0, 0};
			this->teleports = (vec3) {0, 0, 0};
		}

		void update(double dt, double L) {
			// Update position and velocity
			update_verlet(dt);

			// Teleport if outside of box
			if (pos.x > L) {
				pos.x = pos.x - L;
				teleports.x += 1;
			} else if (pos.x < 0) {
				pos.x = L - pos.x;
				teleports.x -= 1;
			}

			if (pos.y > L) {
				pos.y = pos.y - L;
				teleports.y += 1;
			} else if (pos.y < 0) {
				pos.y = L - pos.y;
				teleports.y -= 1;
			}

			if (pos.z > L) {
				pos.z = pos.z - L;
				teleports.z += 1;
			} else if (pos.z < 0) {
				pos.z = L - pos.z;
				teleports.z -= 1;
			}
		}

		void update_verlet(double dt) {
			vec3 acc_prev = acc;
			acc = force;
			vel += 0.5*(acc_prev + acc)*dt;
			pos += vel*dt + ((double) 0.5)*acc*pow(dt, 2);
			force = (vec3) {0, 0, 0};
		}

		vec3 dist_traveled(double L) {
			return pos + teleports*L;
		}

		void save_state(FILE *file) {
			fprintf(file, "Ar %f %f %f\n", pos.x, pos.y, pos.z);
		}
};

double get_time() {
	return (double) chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()).count()/1000.0;
}

void printray(vec3 ray) {
	cout << ray.x << " " << ray.y << " " << ray.z << "\n";
}


double get_length_sqrd(vec3 ray) {
	return pow(ray.x, 2) + pow(ray.y, 2) + pow(ray.z, 2);
}

void get_atom_combinations(vector<vector<Atom*>> &combinations, vector<Atom> &atoms) {
	string bitmask(2, 1); // 2 leading 1's
	bitmask.resize(atoms.size(), 0); // N-2 trailing 0's

	int i = 0;
	do {
		for (int j = 0; j < atoms.size(); j++) {
				if (bitmask[j]) combinations[i].push_back(&atoms[j]);
		}
		i++;
	} while (prev_permutation(bitmask.begin(), bitmask.end()));
}


double U(double r_sqrd) {
	return 4*(pow(r_sqrd, -6) - pow(r_sqrd, -3)) + 2912/531441;
}

vec3 get_force(vec3 &between_vec, double r_sqrd) {
	vec3 direction_vec = between_vec/r_sqrd;
	double force = 24*(2*pow(r_sqrd, -6) - pow(r_sqrd, -3));

	return direction_vec*force;
}


tuple<vec3, double> periodic_boundry(vec3 &between_ray, double L) {
	double dx = between_ray[0];
	dx = dx - round(dx/L)*L;
	double dy = between_ray[1];
	dy = dy - round(dy/L)*L;
	double dz = between_ray[2];
	dz = dz - round(dz/L)*L;

	vec3 direction_ray = {dx, dy, dz};
	double r_sqrd = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);

	return tuple<vec3, double>(direction_ray, r_sqrd);
}

tuple<double, double, double> get_energy(vector<Atom> &atoms, vector<vector<Atom*>> &atom_combinations, double L) {
	// Calculate potential energy:
	double potential_energy = 0;
	for (vector<Atom*> &atom_combination : atom_combinations) {
		Atom *atom1 = atom_combination[0];
		Atom *atom2 = atom_combination[1];

		vec3 between_vec = atom1->pos - atom2->pos;
		auto [direction_vec, r_sqrd] = periodic_boundry(between_vec, L);
		potential_energy += U(r_sqrd);
	}

	// Calculate kinetic energy:
	double kinetic_energy = 0;
	for (Atom &atom : atoms) {
		kinetic_energy += 0.5*get_length_sqrd(atom.vel);
	}

	return tuple<double, double, double>(potential_energy, kinetic_energy, potential_energy+kinetic_energy);
}


void ez_simulate(vector<Atom> &atoms, vector<vector<Atom*>> &atom_combinations, double dt, double t_max, double L) {
	vector<double> t_list((int) t_max/dt, 0);

	FILE* datafile = fopen("data/null", "w");

	for (size_t i = 0; i < t_list.size(); i++) {
		t_list[i] = i*dt;
		// Fancy progress indicator
		if ((int) (t_list[i]*6/t_max) % 1 == 0) {
			cout << "\r";
			for (int j = 0; j < t_list[i]*6/t_max; j++) {
				cout << ".";
			}
			cout << flush;
		}
		step(atoms, atom_combinations, dt, L, datafile);
	}

	fclose(datafile);
}

tuple<vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>> simulate(
	vector<Atom> &atoms, vector<vector<Atom*>> &atom_combinations, double dt, double t_max, string filename, double L, int completion, int total_runs
) {
	// Declare variables
	int n = t_max/dt;
	vector<double> t_list(n, 0);
	vector<double> pot_list(n, 0);
	vector<double> kin_list(n, 0);
	vector<double> tot_list(n, 0);
	vector<double> tmp_list(n, 0);
	vector<double> vac_list(n, 0);
	vector<double> msd_list(n, 0);

	FILE *datafile = fopen(("data/"+filename+".xyz").c_str(), "w");

	for (size_t i = 0; i < t_list.size(); i++) {
		t_list[i] = i*dt;
		// Fancy progress indicator, now even more complicated
		if ((int) (t_list[i]*100/t_max) % 1 == 0) {
			printf("\r...... %i/%i : %3.0f %%", completion+1, total_runs, (100*completion + t_list[i]*100/t_max)/total_runs);
			fflush(stdout);
		}

		// Get and store energy
		auto [pot, kin, tot] = get_energy(atoms, atom_combinations, L);
		pot_list[i] = pot;
		kin_list[i] = kin;
		tot_list[i] = tot;

		// Calculate and store temperature
		double temperature = (2/(3*atoms.size()))*kin;
		tmp_list[i] = temperature;

		// Calculate and store velocity autocorrelation
		double vac = 0;
		for (Atom &atom : atoms) {
			double vac_emum = glm::dot(atom.vel, atom.vel0);
			double vac_denom = get_length_sqrd(atom.vel0);
			vac += vac_emum/vac_denom;
		}
		vac = vac/atoms.size();
		vac_list[i] = vac;

		// Calculate and store mean squared displacement
		double msd = 0;
		for (Atom &atom : atoms) {
			msd += get_length_sqrd(atom.dist_traveled(L) - atom.pos0);
		}
		msd = msd/atoms.size();
		msd_list[i] = msd;

		step(atoms, atom_combinations, dt, L, datafile);
	}
	cout << "\n";

	fclose(datafile);
	return make_tuple(t_list, pot_list, kin_list, tot_list, tmp_list, vac_list, msd_list);
}

void step(vector<Atom> &atoms, vector<vector<Atom*>> &atom_combinations, double dt, double L, FILE *datafile) {
	// Add the force acting on the particles efficiently using pairs
	for (vector<Atom*> &atom_combination : atom_combinations) {
		Atom *atom1 = atom_combination[0];
		Atom *atom2 = atom_combination[1];

		vec3 between_ray = atom1->pos - atom2->pos;
		auto [direction_ray, r_sqrd] = periodic_boundry(between_ray, L);

		if (r_sqrd < 3*3) {
			vec3 force = get_force(direction_ray, r_sqrd);
			atom1->force += force;
			atom2->force -= force;
		}
	}

	fprintf(datafile, "%zi\ntype x y z\n", atoms.size());
	for (Atom &atom : atoms) {
		// Save current atom positions to file
		atom.save_state(datafile);
		// Update atom positions using given method
		atom.update(dt, L);
	}
}


vector<vec3> box_positions(int n, double d) {
	vector<vec3> positions;
	positions.reserve(4*pow(n, 3));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				positions.push_back(vec3(i, j, k)*d);
				positions.push_back(vec3(i, 0.5 + j, 0.5 + k)*d);
				positions.push_back(vec3(0.5 + i, j, 0.5 + k)*d);
				positions.push_back(vec3(0.5 + i, 0.5 + j, k)*d);
			}
		}
	}
	return positions;
}

tuple<double, vector<Atom>, vector<vector<Atom*>>> create_atoms(int atom_count, double d, double temperature) {
	int n = cbrt(atom_count)/cbrt(4);
	double L = d*n;

	vector<Atom> atoms;
	atoms.reserve(atom_count);

	default_random_engine generator;
	normal_distribution<double> normal_temperature(0, sqrt(temperature/119.7));

	vector<vec3> positions = box_positions(n, d);
	for (int i = 0; i < atom_count; i++) {
		vec3 velocitiy = vec3(normal_temperature(generator), normal_temperature(generator), normal_temperature(generator));
		atoms.push_back(Atom(positions[i], velocitiy));
	}

	int atom_combination_count = (atom_count*(atom_count-1))/2;
	vector<vector<Atom*>> atom_combinations(atom_combination_count);
	get_atom_combinations(atom_combinations, atoms);

	return tuple<double, vector<Atom>, vector<vector<Atom*>>>(L, move(atoms), move(atom_combinations));
}

tuple<double, vector<Atom>, vector<vector<Atom*>>> create_equalibrium_atoms(int atom_count, double d, double temperature, double dt, double t_max) {
	auto [L, warm_atoms, atom_combinations] = create_atoms(atom_count, d, temperature);
	ez_simulate(warm_atoms, atom_combinations, dt, t_max, L);
	for (Atom &atom : warm_atoms) {
		atom.vel0 = atom.vel;
		atom.pos0 = atom.pos;
	}

	return tuple<double, vector<Atom>, vector<vector<Atom*>>>(L, move(warm_atoms), move(atom_combinations));
}

void run(string filename, int simulate_count) {
	double dt = 0.01;
	int length = 5;

	vector<vector<double>> sum_lists(6, vector<double>((int) (length/dt)+1));
	vector<double> t_list;

	double start_time = get_time();
	for (int i = 0; i < simulate_count; i++) {
		auto [L, atoms, atom_combinations] = create_equalibrium_atoms(108, 1.7, 180, dt, length);

		auto [t_list, pot_list, kin_list, tot_list, tmp_list, vac_list, msd_list] = simulate(atoms, atom_combinations, dt, length, filename, L, i, simulate_count);
		for (int j; j < t_list.size(); j++) {
			sum_lists[0][j] += pot_list[j]/simulate_count;
			sum_lists[1][j] += kin_list[j]/simulate_count;
			sum_lists[2][j] += tot_list[j]/simulate_count;
			sum_lists[3][j] += tmp_list[j]/simulate_count;
			sum_lists[4][j] += vac_list[j]/simulate_count;
			sum_lists[5][j] += msd_list[j]/simulate_count;
		}
	}
	printf("time: %.3g s\n", get_time() - start_time);

	FILE* datafile = fopen(("data/"+filename+".energy").c_str(), "w");

	for (int i = 0; i < t_list.size(); i++) {
		fprintf(datafile, "%f %f %f %f %f %f\n", sum_lists[0][i], sum_lists[1][i], sum_lists[2][i], sum_lists[3][i], sum_lists[4][i], sum_lists[5][i]);
	}
	fclose(datafile);
}

int main(int argc, char const *argv[]) {
	run("null", 1);

	return 0;
}
