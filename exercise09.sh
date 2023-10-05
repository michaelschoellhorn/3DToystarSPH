# a complete shell script to simulate the formation of a toy star via SPH and generate movies from it
# Should handle everthing by itself
# !compiler or python3 command might need to be changed!

# Can also be done manually by 
#1. setting up folder structure(Plots/rho/ Plots/star/ and results/)
#2. compiling main.cpp and running it
#3. running plot.py
#4. running the two ffmpeg instructions in the block below

mkdir results
mkdir Plots
cd Plots/
mkdir rho
mkdir star
cd ..
g++ main.cpp
./a.out
python3 plot.py
ffmpeg -framerate 30 -pattern_type glob -i "Plots/star/*.png" -c:v libx264 -pix_fmt yuv420p star.mp4
ffmpeg -framerate 30 -pattern_type glob -i "Plots/rho/*.png" -c:v libx264 -pix_fmt yuv420p rho.mp4
rm -r results
rm -r Plots
rm a.out