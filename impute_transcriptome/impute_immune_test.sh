#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=megan.ding@mail.utoronto.ca
#SBATCH --job-name=Impute_abcd_expression_CMC
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
mamba activate metaxcan

#paths to env variables
export METAXCAN=~/sc-tprs-mdd/MetaXcan/software
export MODELS=/genome/scratch/nikolova/share/fsantos/PrediXcan/MetaXcan/software/data/weights/Downloads/CMC_DLPFC-PrediXcan_Models/
export DATA=~/sc-tprs-mdd/abcd_data
export RESULTS=~/sc-tprs-mdd/abcd_imputed

#Run MetaXcan
python3 $METAXCAN/Predict.py \
--model_db_path $MODELS/DLPFC_newMetax.db \
--vcf_genotypes $DATA/abcd_genotypes.vcf \
--vcf_mode genotyped \
--prediction_output $RESULTS/with_CMC/abcd_pred_${SLURM_ARRAY_TASK_ID}.txt  \
--prediction_summary_output $RESULTS/with_CMC/summary_abcd_pred_${SLURM_ARRAY_TASK_ID}.txt \
--verbosity 9 \
--throw
