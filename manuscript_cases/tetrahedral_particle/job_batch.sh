#!/bin/bash -l
#SBATCH --job-name=tetra=0.0016
#SBATCH --output=tetra=0.0016.log
#Number of processors
#SBATCH --time=3-00:00:0
#SBATCH --nodes=1
#Number of cores
#SBATCH --ntasks-per-node=1
##SBATCH --exclusive
#SBATCH --mail-user=didarula@buffalo.edu
#SBATCH --mem=12G
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc
./MemDynamics
