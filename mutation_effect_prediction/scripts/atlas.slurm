#!/bin/bash
## Job Name
#SBATCH --job-name=atlas
## Allocation Definition
## The account and partition options should be the same except in a few cases (e.g. ckpt queue and genpool queue).
#SBATCH --account=spe
#SBATCH --partition=gpu-a40
#SBATCH --gres=gpu:1
## Resources
## Total number of Nodes
#SBATCH --nodes=1
## Number of cores per node
#SBATCH --ntasks-per-node=1
## Walltime. Do not specify a walltime substantially more than your job needs.
#SBATCH --time=23:00:00
## Memory per node. It is important to specify the memory since the default memory is very small.
## For mox, --mem may be more than 100G depending on the memory of your nodes.
## For ikt, --mem may be 58G or more depending on the memory of your nodes.
## See above section on "Specifying memory" for choices for --mem.
#SBATCH --mem=64G
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/spe/gvisan01/ProteinMPNN-copy/mutation_effect_prediction/scripts
##turn on e-mail notification
## export all your environment variables to the batch job session
#SBATCH --export=all
#SBATCH -e atlas.err
#SBATCH -o atlas.log

bash atlas.sh

