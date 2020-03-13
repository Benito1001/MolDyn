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

double get_length_sqrd();
tuple<vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>> simulate();
void step(vector<Atom> &atoms, vector<vector<Atom*>> &atom_combinations, double dt, double L, FILE *datafile);
tuple<double, vector<Atom>, vector<vector<Atom*>>> create_equalibrium_atoms();
