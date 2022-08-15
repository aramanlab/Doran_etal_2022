#!/bin/sh
#SBATCH --job-name=runSPI
#SBATCH --output=/scratch/midway2/bend/projects/Doran_etal_2022/_research/logs/runSPI/runSPI_%A_%a.out
#SBATCH --error=/scratch/midway2/bend/projects/Doran_etal_2022/_research/logs/runSPI/runSPI_%A_%a.err
#SBATCH --partition=broadwl
#SBATCH --array=1-288%20
#SBATCH --cpus-per-task=10
#SBATCH --time=0-00:05:00
#SBATCH --mem-per-cpu=2G  # NOTE DO NOT USE THE --mem=OPTION 
#SBATCH --mail-type=START,END
#SBATCH --mail-user=bend@uchicago.edu

# When running a large number of tasks simultaneously, it may be
# necessary to increase the user process limit.
# ulimit -u 10000

module load julia/1.7.2

projdir="/scratch/midway2/bend/projects/Doran_etal_2022"
sourcefiles=(ls $projdir/data/sims/MSAs/*)
inputfile="${sourcefiles[$SLURM_ARRAY_TASK_ID]}"
name=$(basename $inputfile .phy)
outputdir="${projdir}/_research/runSPI/${name}"
mkdir -p $outputdir

echo "projdir: " $projdir
echo "inputfile: " $inputfile
echo "name: " $name
echo "outputdir: " $outputdir

start=`date +%s`
srun julia --threads=$SLURM_CPUS_PER_TASK scripts/slurm/runners/runSPI.jl \
    -i $inputfile \
    -o $outputdir \
    --nboot 100 > $outputdir/runSPI.out
end=`date +%s`
timeelapsed=$((end-start))
echo "time_elapsed: " $timeelapsed
echo "time_elapsed: " $timeelapsed >> $outputdir/runSPI.out
