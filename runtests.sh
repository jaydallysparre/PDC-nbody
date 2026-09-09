#!/bin/bash

mkdir .temp
make all DIR=.temp

# test for a number of different process counts
procN=(1 2 3 4 6 12)

eps="0.000000001" # acceptable error

for test in testcases/*.test; do
    name=$(basename "$test" .test)

    for n in "${procN[@]}"; do
        echo "-----------------------------------------------------------------"
        echo "TEST INPUT: $test ON $n PROCESSES"
        echo "-----------------------------------------------------------------"
        mpirun -n "$n" .temp/nbodyref $(cat "$test") < testcases/12_particle.ic > ".temp/${name}_${n}_ref"
        mpirun -n "$n" .temp/nbody1a $(cat "$test") < testcases/12_particle.ic > ".temp/${name}_${n}_1a"

        python3 epsdiff.py ".temp/${name}_${n}_ref" ".temp/${name}_${n}_1a" "$eps"
        statusA=$? # exit status for ref vs 1a

        if [ "$statusA" -eq 0 ]; then
            echo "A PASS"
        else
            echo "A FAILED"
        fi

        # placed after because its easier to catch mpi errors
        mpirun -n "$n" .temp/nbody1b $(cat "$test") < testcases/12_particle.ic > ".temp/${name}_${n}_1b"

        python3 epsdiff.py ".temp/${name}_${n}_ref" ".temp/${name}_${n}_1b" "$eps"
        statusB=$? # exit status for ref vs 1b

        if [ "$statusB" -eq 0 ]; then
            echo "B PASS"
        else
            echo "B FAILED"
        fi
    done
done

rm -rf .temp

