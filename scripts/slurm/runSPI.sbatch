#!/bin/sh

#SBATCH --job-name=runSPI
#SBATCH --output=runSPI_%j.out
#SBATCH --error=runSPI_%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2G  # NOTE DO NOT USE THE --mem=OPTION 
#SBATCH array XXXXX

# When running a large number of tasks simultaneously, it may be
# necessary to increase the user process limit.
# ulimit -u 10000

module load julia-1.7.2

outprefix=

srun julia \
    --project=@. \
    --threads=20 \
    ./runners/runSPI.jl \
        --inputfile ...\
        --outputdir ...\
        --nboot 100