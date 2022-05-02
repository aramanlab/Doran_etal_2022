#!/bin/sh

#SBATCH --job-name=testarrayjob
#SBATCH --output=/scratch/midway2/bend/projects/Doran_etal_2022/_research/logs/testarrayjob_%j_%t.out
#SBATCH --error=/scratch/midway2/bend/projects/Doran_etal_2022/_research/logs/testarrayjob_%j_%t.err
#SBATCH --partition=broadwl
#SBATCH --array=1-5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:00:10
# #SBATCH --mem-per-cpu=1G  # NOTE DO NOT USE THE --mem=OPTION 

# When running a large number of tasks simultaneously, it may be
# necessary to increase the user process limit.
# ulimit -u 10000

module load julia/1.7.2

projdir="/scratch/midway2/bend/projects/Doran_etal_2022"
sourcefiles=(ls $projdir/data/sims/MSAs/*)
inputfile="${sourcefiles[$SLURM_ARRAY_TASK_ID]}"
name=$(basename $inputfile .phy)
outputdir="${projdir}/_research/runSPI/${name}"

echo "datadir: " $datadir
echo "inputfile: " $inputfile
echo "name: " $name
echo "outputdir: " $outputdir
