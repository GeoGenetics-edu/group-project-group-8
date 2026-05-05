#!/bin/bash
#SBATCH --job-name=extract_hq_mags
#SBATCH --output=/maps/projects/course_1/scratch/group8/group-project-group-8/week18-preprocessing/logs/extract_hq_mags_%x_%j.out   # stdout log
#SBATCH --error=/maps/projects/course_1/scratch/<group_#>/logs/extract_hq_mags_%x_%j.err    # stderr log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00

set -euo pipefail
