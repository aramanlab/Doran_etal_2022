#!/bin/sh
#SBATCH --output=/scratch/midway2/bend/projects/Doran_etal_2022/_research/logs/emapperUP8881-%j.out
#SBATCH --error=/scratch/midway2/bend/projects/Doran_etal_2022/_research/logs/emapperUP8881-%j.err
#SBATCH --partition=broadwl
#SBATCH --time=36:00:00
#SBATCH --nodes=100
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28

module load parallel
module load python
ulimit -u 10000

activate /project2/araman/eggnog-mapper/emapperenv
export PATH=/project2/araman/eggnog-mapper:/project2/araman/eggnog-mapper/eggnogmapper/bin:"$PATH"
export EGGNOG_DATA_DIR=/project2/araman/eggnog-mapper/data

arg_list=$(cat /project2/araman/UP202202/UP8881proteomefiles.txt) 
outdir="/home/bend/scratchspace/projects/Doran_etal_2022/_research/emapperUP8881"
mkdir -p $outdir

srun="srun --exclusive -N1 -n1 --cpus-per-task $SLURM_CPUS_PER_TASK"
parallel="parallel --delay 0.2 -j $SLURM_NNODES --joblog $outdir/task.log --tmpdir=$TMPDIR --resume" # --retry-failed --resume-failed

$parallel "$srun /project2/araman/eggnog-mapper/emapperenv/bin/python3 /project2/araman/eggnog-mapper/emapper.py \
    --cpu 0 -m diamond --evalue 0.001 --seed_ortholog_score 60.0 --query_cover 20.0 --subject_cover 0.0 --tax_scope 2 --tax_scope_mode inner_narrowest \
    --data_dir $EGGNOG_DATA_DIR --dmnd_db $EGGNOG_DATA_DIR/bacteria.dmnd \
    --temp_dir $TMPDIR --override --dmnd_ignore_warnings \
    -i {1} -o {2} --output_dir $outdir" ::: $arg_list :::+ $(basename --multiple $arg_list)
