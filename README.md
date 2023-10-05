# 3DToystarSPH
Code to simulate the collapse of a polytropic gas cloud into a more dense object resembling a star. For gravitational modelling it uses a simple linear acceleration -lambda * (x, y, z) approximation, no particle-particle gravity is used. A dampening force -nu * (vx, vy, vz) is used to model energy dissipation.

## Code Manual
The shell script runStarSim.sh should handle everything by itself. Alternativly one can set up the folder structure seen in the shell script and then compile and execute the main.cpp script. It uses the random_distribution.dat file as input data and outputs a txt file with the results. By running plot.py you can generate the plots needed for the animation with ffmpeg. By changing the initial data of particles in the input file different particle distributions can be generated. The physical parameters specified in main.cpp can be changed to adjust the behaviour of the simulated system.
