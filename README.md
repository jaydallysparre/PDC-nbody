# PDC-nbody
The goal of this assignment was to take a naive N-body solver written in C with MPI and optimize its memory usage.

# Usage
A makefile is provided to build each iteration of the program. `make all` will build all three; `nbody1a`, `nbody1b`, and `nbodyref`.

To run, in terminal use

`mpirun -n <processes> ./<program> <number of particles> <number of timesteps> <size of timestep> <output frequency> <g|i>`

where `g` generates the initial conditions, and `i` allows them to be read from stdin.
