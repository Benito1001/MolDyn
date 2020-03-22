#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <random>
#include <utility>
#include <cmath>
#include <tuple>
#include <list>
#include <vector>
#include <valarray>
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>
#include <omp.h>
#include <mutex>

using namespace std;
using vec3 = glm::vec<3, double>;

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
				pos.x = pos.x + L;
				teleports.x -= 1;
			}

			if (pos.y > L) {
				pos.y = pos.y - L;
				teleports.y += 1;
			} else if (pos.y < 0) {
				pos.y = pos.y + L;
				teleports.y -= 1;
			}

			if (pos.z > L) {
				pos.z = pos.z - L;
				teleports.z += 1;
			} else if (pos.z < 0) {
				pos.z = pos.z + L;
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

double get_time();
void printray(vec3 ray);
double get_length_sqrd(vec3 ray);

void get_atom_combinations(vector<pair<Atom*, Atom*>> &combinations, vector<Atom> &atoms);
double U(double r_sqrd);
vec3 get_force(vec3 &between_vec, double r_sqrd);
tuple<vec3, double> periodic_boundry(vec3 &between_ray, double L);
tuple<double, double, double> get_energy(vector<Atom> &atoms, vector<pair<Atom*, Atom*>> &atom_combinations, double L);

void ez_simulate(vector<Atom> &atoms, vector<pair<Atom*, Atom*>> &atom_combinations, double dt, double t_max, double L);
tuple<vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>> simulate(
	vector<Atom> &atoms, vector<pair<Atom*, Atom*>> &atom_combinations, double dt, double t_max, string filename, double L, int completion, int total_runs
);
void step(vector<Atom> &atoms, vector<pair<Atom*, Atom*>> &atom_combinations, double dt, double L, FILE *datafile);
vector<vec3> box_positions(int n, double d);
tuple<double, vector<Atom>, vector<pair<Atom*, Atom*>>> create_atoms(int atom_count, double d, double temperature);
tuple<double, vector<Atom>, vector<pair<Atom*, Atom*>>> create_equalibrium_atoms(int atom_count, double d, double temperature, double dt, double t_max);
