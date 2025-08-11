#!/bin/bash -l
#SBATCH --job-name=nv=252
#SBATCH --output=nv=252.log
#Number of processors
#SBATCH --time=10-00:00:0
#SBATCH --nodes=1
#Number of cores
#SBATCH --ntasks-per-node=1
##SBATCH --exclusive
#SBATCH --mail-user=didarula@buffalo.edu
#SBATCH --mem=12G
#SBATCH --partition=xinyong
#SBATCH --qos=xinyong
#SBATCH --cluster=faculty
./MemDynamics
