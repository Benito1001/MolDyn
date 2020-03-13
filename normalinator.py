import math
import random

with open("data/normal_speeds.dat", "w") as outfile:
	normal = lambda: random.gauss(0, math.sqrt(180/119.7))
	for i in range(108):
		outfile.write(f"{normal()}, {normal()}, {normal()}\n")
