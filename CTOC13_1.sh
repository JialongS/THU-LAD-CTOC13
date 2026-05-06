#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1  
#SBATCH -n 1  
#SBATCH -c 64 
#SBATCH -J SongJL
srun inplane