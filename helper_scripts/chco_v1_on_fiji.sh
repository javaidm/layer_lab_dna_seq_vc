#!/bin/sh 
#SBATCH --job-name=nf_chco_v1
#SBATCH --output=nf_chco_v1.out
#SBATCH --error=nf_chco_v1.err
#SBATCH --time=120:00:00 
#SBATCH -p highmem 
#SBATCH --cpus=1 
#SBATCH --mem=16gb
nextflow -log logs/log run ../main.nf -profile fiji
