EXECS=nbody1a nbody1b nbodyref
MPICC?=mpicc

all: ${EXECS}

mpi_part1a: part1a.c
	${MPICC} -o nbody1a part1a.c -lm

mpi_part1b: part1b.c
	${MPICC} -o nbody1b part1b.c -lm

reference: mpi_nbody_basic.c
	${MPICC} -o nbodyref mpi_nbody_basic.c -lm

clean:
	rm -f ${EXECS}
