all:


oppg4ci: oppg4/oppg4ci.cpp moldyn_functions.cpp
	g++ $^ -o oppg4/$@ -std=c++17 -O3 -march=native -fopenmp -ffast-math

oppg4ciii_get_equalibrium: oppg4/oppg4ciii_get_equalibrium.cpp moldyn_functions.cpp
	g++ $^ -o oppg4/$@ -std=c++17 -O3 -march=native -fopenmp -ffast-math

oppg4ciii: oppg4/oppg4ciii.cpp moldyn_functions.cpp
	g++ $^ -o oppg4/$@ -std=c++17 -O3 -march=native -fopenmp -ffast-math
