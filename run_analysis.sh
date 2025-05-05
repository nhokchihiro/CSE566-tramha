#!/bin/bash
#SBATCH --job-name=analysis_sv
#SBATCH --output=/scratch/tmh6573/CSE556-2/logs/analysis_sv_%j.out
#SBATCH --error=/scratch/tmh6573/CSE556-2/logs/analysis_sv_%j.err
#SBATCH --time=20:00:00
#SBATCH --mem=100GB

all_file="/scratch/tmh6573/CSE556-2/analysis/all_svcall.tsv"
simulated_file="/scratch/tmh6573/CSE556-2/analysis/simulated_svcall.tsv"
out_dir="/scratch/tmh6573/CSE556-2/analysis/tolerance_2000"
tolerance_base=2000

#--- NO CHANGE AFTER THIS LINE ---

python /scratch/tmh6573/CSE556-2/scripts/analysis.py \
  --all_svcall "$all_file" \
  --simulated_svcall "$simulated_file" \
  --output_dir "$out_dir" \
  --tolerance "$tolerance_base"

