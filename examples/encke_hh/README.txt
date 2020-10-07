To compile this example, type the following:

make clean
make

To run the encke code with an initial stepsize of 40 days and adaptive stepsize
with epsilon=1e-7 step size, type:

./rebound 40 1e7 1e-7


To run the encke code with a fixed 40-day step size, type:

./rebound 40 1e7 -1


Currently the initial conditions are set to those corresponding to Fig. 1.  
