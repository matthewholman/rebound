To compile this example, type the following:

make clean
make

To run EnckeHH with a fixed 50-day step size, for 5*10^7 days, type:

./rebound 50 5e7 -1

This outputs the cpu time and the energy error, and can be used to plot a point in Figure 1 of Hernandez & Holman (2020).  The initial conditions correspond to the outer Solar System plus the Sun.  These initial conditions are the same as those in the example problem, outer_solar_system, except for the Pluto contribution.

EnckeHH can also be run with adaptive steps, but this is not usually recommended.  To use EnckeHH with an initial stepsize of 50 days, adaptive step parameter epsilon=1e-7, with the same initial conditions for 5*10^7 days, type:

./rebound 50 5e7 1e-7 
