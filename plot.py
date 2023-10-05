
# # File 2 while manually executing
# ## Plots the solution
# This script is used to plot the results from the cpp file, which has to be executed first. It generates the png files used by ffmpeg to create the movies.


import numpy as np
import matplotlib.pyplot as plt

# We begin by defining a function to load the output data of the c++ script.
def load(sourcefile):
    x, y, z, rho = np.loadtxt(sourcefile, delimiter=' ', unpack=True)
    return x, y, z, rho


# Now we use this to load and plot the distribution and the density curve for all particles. These plots are now saved in the Plots/star/ and Plots/rho/ dir
sources = []
for i in range(1, 201):
    sources.append(f'results/out{i}.txt')

i=1000
for source in sources:
    x, y, z, rho = load(source)
    plt.scatter(x, y, label=source, c=rho)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.colorbar(label=r'$\rho$')
    plt.xlim(-0.8, 0.8)
    plt.ylim(-0.8, 0.8)
    plt.savefig('Plots/star/'+str(i)+'.png')
    plt.close()
    dist = (x**2+y**2+z**2)**0.5
    plt.scatter(dist, rho, label='Simulation')
    plt.plot(np.linspace(0, 0.8, 1000), 2.012/(0.2*2.0)*(0.75**2-np.linspace(0, 0.8, 1000)**2), label='analytic')
    plt.ylim(0, 3.5)
    plt.xlabel('r')
    plt.ylabel(r'$\rho$')
    plt.legend()
    plt.savefig('Plots/rho/' + str(i) + '.png')
    plt.close()
    i += 1


