#!/bin/bash
mpirun -np 1 ./project.x 1000 500 500 10 1 2 1000 1000 simulation
mpirun -np 2 ./project.x 1000 500 500 10 1 2 1000 1000 simulation
mpirun -np 4 ./project.x 1000 500 500 10 1 2 1000 1000 simulation
mpirun -np 8 ./project.x 1000 500 500 10 1 2 1000 1000 simulation
mpirun -np 16 ./project.x 1000 500 500 10 1 2 1000 1000 simulation
mpirun -np 32 ./project.x 1000 500 500 10 1 2 1000 1000 simulation
mpirun -np 1 ./project.x 2000 1000 1000 10 1 2 1000 1000 simulation
mpirun -np 2 ./project.x 2000 1000 1000 10 1 2 1000 1000 simulation
mpirun -np 4 ./project.x 2000 1000 1000 10 1 2 1000 1000 simulation
mpirun -np 8 ./project.x 2000 1000 1000 10 1 2 1000 1000 simulation
mpirun -np 16 ./project.x 2000 1000 1000 10 1 2 1000 1000 simulation
mpirun -np 32 ./project.x 2000 1000 1000 10 1 2 1000 1000 simulation
mpirun -np 1 ./project.x 4000 2000 2000 10 1 2 1000 1000 simulation
mpirun -np 2 ./project.x 4000 2000 2000 10 1 2 1000 1000 simulation
mpirun -np 4 ./project.x 4000 2000 2000 10 1 2 1000 1000 simulation
mpirun -np 8 ./project.x 4000 2000 2000 10 1 2 1000 1000 simulation
mpirun -np 16 ./project.x 4000 2000 2000 10 1 2 1000 1000 simulation
mpirun -np 32 ./project.x 4000 2000 2000 10 1 2 1000 1000 simulation
mpirun -np 1 ./project.x 8000 4000 4000 10 1 2 1000 1000 simulation
mpirun -np 2 ./project.x 8000 4000 4000 10 1 2 1000 1000 simulation
mpirun -np 4 ./project.x 8000 4000 4000 10 1 2 1000 1000 simulation
mpirun -np 8 ./project.x 8000 4000 4000 10 1 2 1000 1000 simulation
mpirun -np 16 ./project.x 8000 4000 4000 10 1 2 1000 1000 simulation
mpirun -np 32 ./project.x 8000 4000 4000 10 1 2 1000 1000 simulation
mpirun -np 1 ./project.x 16000 8000 8000 10 1 2 1000 1000 simulation
mpirun -np 2 ./project.x 16000 8000 8000 10 1 2 1000 1000 simulation
mpirun -np 4 ./project.x 16000 8000 8000 10 1 2 1000 1000 simulation
mpirun -np 8 ./project.x 16000 8000 8000 10 1 2 1000 1000 simulation
mpirun -np 16 ./project.x 16000 8000 8000 10 1 2 1000 1000 simulation
mpirun -np 32 ./project.x 16000 8000 8000 10 1 2 1000 1000 simulation