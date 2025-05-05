#!/bin/bash
#SBATCH --job-name=extract_sv
#SBATCH --output=/scratch/tmh6573/CSE556-2/logs/extract_sv_%j.out
#SBATCH --error=/scratch/tmh6573/CSE556-2/logs/extract_sv_%j.err
#SBATCH --time=20:00:00
#SBATCH --mem=100GB

input_folder="/scratch/tmh6573/CSE556-2/VCF"
output_folder="/scratch/tmh6573/CSE556-2/analysis"

#--- NO CHANGE AFTER THIS LINE ---

python /scratch/tmh6573/CSE556-2/scripts/extract_svcall.py --base_dir "$input_folder" --output_dir "$output_folder"
