EXECS=nbody1a nbody1b nbodyref
MPICC?=mpicc

# going to redirect for tests
DIR ?= .

all: ${EXECS}

nbody1a: part1a.c
	${MPICC} -o ${DIR}/nbody1a part1a.c -lm

nbody1b: part1b.c
	${MPICC} -o ${DIR}/nbody1b part1b.c -lm

nbodyref: mpi_nbody_basic.c
	${MPICC} -o ${DIR}/nbodyref mpi_nbody_basic.c -lm

clean:
	rm -f ${EXECS}
