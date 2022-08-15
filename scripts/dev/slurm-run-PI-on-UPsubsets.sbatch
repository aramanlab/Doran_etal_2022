#!/bin/sh
#SBATCH --output=/scratch/midway2/bend/projects/Doran_etal_2022/_research/logs/emapperUP8881-%j.out
#SBATCH --error=/scratch/midway2/bend/projects/Doran_etal_2022/_research/logs/emapperUP8881-%j.err
#SBATCH --partition=broadwl
#SBATCH --time=4:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28

module load parallel
module load julia
ulimit -u 10000

projdir="/home/bend/scratchspace/projects/Doran_etal_2022/"
indir="$projdir/data/exp_pro/UPsubsetMSAs"
arg_list=$(ls $indir) 
outdir="$projdir/_research/UPsubsetMSAs"
runnerscripts=('runSPI', 'runFastTree', 'runFastME', 'runPhyML', 'runRAxML', 'runMrBayes')
mkdir -p $outdir

srun="srun --exclusive -N1 -n1 --cpus-per-task $SLURM_CPUS_PER_TASK"
parallel="parallel --dryrun --delay 0.2 -j $SLURM_NNODES --joblog $outdir/SPItask.log --tmpdir=$TMPDIR --resume" # --retry-failed --resume-failed
$parallel "$srun julia --threads=$SLURM_CPUS_PER_TASK $projdir/scripts/slurm/runners/{2}.jl \
    -i {1} \
    -o $outdir \
    -m JC69 \
    --nboot 100 > $outdir/{2}.out \
    " ::: $arg_list ::: $runnerscripts
