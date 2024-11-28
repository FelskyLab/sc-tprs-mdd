#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=megan.ding@mail.utoronto.ca
#SBATCH --job-name=modify_immune_models
#SBATCH --output=output/slurm-impute_abcd-%A.%a.out
#SBATCH --error=error/slurm-impute_abcd-%A.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=8192MB
#SBATCH --array=1-22
#SBATCH --time=0-400:00:00


#Load modules
#module load PYTHON/3.8.5-Anaconda3-2021.03
#conda activate imlabtools
module purge
source /external/rprshnas01/netdata_kcni/dflab/.bashrc.d/miniforge.bashrc

#paths to env variables
export MODELS=~/sc-tprs-mdd/immune_cell_models
export PATH=~/sc-tprs-mdd
sed $'s/\r$//' ./rsid_from_pos_models.R > ./rsid_from_pos_models.unix.R

for filename in $MODELS/*.db; do
  #Run rsid_from_pos_models.R
  $PATH/rsid_from_pos_models.unix.R $filename
  
done
