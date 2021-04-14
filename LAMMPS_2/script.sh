#!/bin/sh
#!/bin/bash
#SBATCH -p long
#SBATCH -J 0.6N0.5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=14-00
#SBATCH --output=energy.out
#SBATCH --cpus-per-task=16
export OMP_NUM_THREADS=1

mpirun -np 16 ./lmp -in in.complexation 


