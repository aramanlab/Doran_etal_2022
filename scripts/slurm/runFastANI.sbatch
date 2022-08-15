#!/bin/sh
#SBATCH --output=/scratch/midway2/bend/projects/Doran_etal_2022/_research/FastANI/fastANI-%j.out
#SBATCH --error=/scratch/midway2/bend/projects/Doran_etal_2022/_research/FastANI/fastANI-%j.err
#SBATCH --partition=broadwl
# #SBATCH --time=02:00:00
#SBATCH --nodes=11
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --mail-type=START,END
#SBATCH --mail-user=bend@uchicago.edu

# Load the default version of GNU parallel.
module load parallel

# setup inputs and outputs
outdir="/home/bend/scratchspace/projects/Doran_etal_2022/_research/FastANI"
fastafileslist="/project2/araman/UP202202/UP8881genomefiles.txt"
shardlist=$(ls /project2/araman/UP202202/genomeinputshards/*)
shardnamelist=$(basename --multiple $shardlist)


## setup env to use julia artifact binary
alias fastANI=/home/bend/.julia/artifacts/3804aeeedd511a74e3216b7551fda62a70cfc8a2/bin/fastANI
export PATH="/home/bend/.julia/artifacts/3908767681eedc038a8c56d4c4ba2e1ca633e690/bin:/home/bend/.julia/artifacts/3804aeeedd511a74e3216b7551fda62a70cfc8a2/bin:/software/slurm-current-el7-x86_64/bin:/software/parallel-latest-el7-x86_64/bin:/software/julia-1.7.2-el7-x86_64/bin:/software/curl-7.50-el7-x86_64/bin:/software/slurm-current-el7-x86_64/bin:/bin:/home/bend/.vscode-server/bin/c3511e6c69bb39013c4a4b7b9566ec1ca73fc4d5/bin/remote-cli:/bin:/software/python-anaconda-2020.02-el7-x86_64/bin:/software/python-anaconda-2020.02-el7-x86_64/condabin:/bin:/software/git-2.10-el7-x86_64/bin:/software/subversion-1.9.4-el7-x86_64/bin:/software/bin:/srv/adm/bin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/lpp/mmfs/bin:/home/bend/bin:/home/bend/bin"
export LD_LIBRARY_PATH="/home/bend/.julia/artifacts/3908767681eedc038a8c56d4c4ba2e1ca633e690/lib:/software/julia-1.7.2-el7-x86_64/lib/julia:/software/julia-1.7.2-el7-x86_64/bin/../lib/julia:/software/julia-1.7.2-el7-x86_64/bin/../lib/julia:/software/julia-1.7.2-el7-x86_64/bin/../lib:/software/julia-1.7.2-el7-x86_64/lib/julia:/software/julia-1.7.2-el7-x86_64/lib:/software/curl-7.50-el7-x86_64/lib:/software/slurm-current-el7-x86_64/lib64:/software/slurm-current-el7-x86_64/lib:/software/git-2.10-el7-x86_64/lib64:/software/subversion-1.9.4-el7-x86_64/lib:/software/serf-1.3.9-el7-x86_64/lib"
export fastANI

# setup slurm commands
srun="srun --exclusive -N$SLURM_NNODES -n$SLURM_NTASKS --cpus-per-task $SLURM_CPUS_PER_TASK"
parallel="parallel --env PATH --env LD_LIBRARY_PATH --env fastANI --delay 0.2 -j $SLURM_NNODES --joblog $outdir/task.log --resume --resume-failed"

# Run the parallel command.
$parallel "$srun fastANI \
    --threads $SLURM_CPUS_PER_TASK \
    --ql $fastafileslist --rl {1} \
    --output $outdir/{2}shard.tsv" ::: $shardlist :::+ $shardnamelist
