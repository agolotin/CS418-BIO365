#!/bin/bash

#SBATCH --time=167:58:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=7168 # memory per CPU core
#SBATCH -J "mapper_optimal"   # job name
#SBATCH --mail-user=artem.golotin@gmail.com   # email address
#SBATCH --mail-type=END

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
chunk=
kmer_size=

while test $# -gt 0
do
	case $1 in
		-c)
			chunk=$2
			shift
			;;
		-k)
			kmer_size=$2
			shift
			;;
		*) 
			echo >&2 "Invalid argument: $1"
			;;
	esac
	shift
done

./dna_mapper_optimal.py shotgun_small/chr19.small.error.fastq shotgun_small/chunk$chunk\.chr19.small.error.reads.fastq $kmer_size $chunk
#./dna_mapper_threadpool.py shotgun_small/chr19.small.error.fastq shotgun_small/chr19.small.error.reads.fastq 10 errors
#./dna_mapper.py shotgun_small/chr19.small.error.fastq shotgun_small/chr19.small.error.reads.fastq 10 errors #12 16 21
#./dna_mapper_threadded.py shotgun_large/chr19.chr.fasta shotgun_large/chr19.reads.fastq 16 32
