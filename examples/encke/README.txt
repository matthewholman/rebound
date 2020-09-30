To use nbody.c please place all the files in the same directory and type ’make’ in that directory.  

To run an N-body simulation please type ‘./nbody a b’.  a is the time step in years and b is the runtime in years.

Currently the initial conditions are set to those corresponding to Fig. 1.  For example, to generate a point that would go in Fig. 1, you could type ‘time ./nbody .002 10000’ NOTE: the runtime in the paper should be 10,000, not 100,000 years.  You should get something like this:

time ./nbody .002 10000
dE/E=1.64996e-06

real	0m8.655s
user	0m8.638s
sys	0m0.013s

To change the initial conditions, you can edit subroutine ‘infig1.’  Alternatively, I have prepared subroutines with initial conditions that correspond to other figures in the paper: ‘infig2’,’infig3’,’infig5’,and ’infig6.’  By editing nbody.c you can use them.  Note there are sometimes special instructions for using the initial conditions, described in the Paper text.  