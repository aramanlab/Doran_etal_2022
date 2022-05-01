#!/bin/sh

#SBATCH --job-name=testarrayjob
#SBATCH --output=/scratch/midway2/bend/projects/Doran_etal_2022/_research/logs/testarrayjob_%j_%t.out
#SBATCH --error=/scratch/midway2/bend/projects/Doran_etal_2022/_research/logs/testarrayjob_%j_%t.err
#SBATCH --partition=broadwl
#SBATCH --array=1-288%20
# #SBATCH --ntasks=20
#SBATCH --cpus-per-task=10
#SBATCH --time=0-00:02:00
#SBATCH --mem-per-cpu=2G  # NOTE DO NOT USE THE --mem=OPTION 

# When running a large number of tasks simultaneously, it may be
# necessary to increase the user process limit.
ulimit -u 10000

module load julia/1.7.2

projdir="/scratch/midway2/bend/projects/Doran_etal_2022"
sourcefiles=(ls $projdir/data/sims/MSAs/*)
inputfile="${sourcefiles[$SLURM_ARRAY_TASK_ID]}"
name=$(basename $inputfile .phy)
outputdir="${projdir}/_research/runSPI/${name}"

# julia --threads=5 scripts/slurm/runners/runSPI_test.jl \
#     -i $inputfile \
#     -o $outputdir \
#     --nboot 100 > $outputdir/runSPI.out

echo "datadir: " $datadir
echo "inputfile: " $inputfile
echo "name: " $name
echo "outputdir: " $outputdir