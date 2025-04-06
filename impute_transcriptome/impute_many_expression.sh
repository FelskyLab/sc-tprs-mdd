#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=megan.ding@mail.utoronto.ca
#SBATCH --job-name=Impute_abcd_expression_immune
#SBATCH --output=/external/rprshnas01/external_data/abcd/workspace/mding/output/slurm-impute_abcd-%A.%a.out
#SBATCH --error=/external/rprshnas01/external_data/abcd/workspace/mding/error/slurm-impute_abcd-%A.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=8192MB
#SBATCH --array=1
#SBATCH --time=0-400:00:00


module purge
source /external/rprshnas01/netdata_kcni/dflab/.bashrc.d/miniforge.bashrc
mamba activate metaxcan

#paths to env variables
export METAXCAN=~/sc-tprs-mdd/MetaXcan/software
export MODELS=~/abcd/workspace/mding/immune_cell_models
export DATA=~/abcd/workspace/mding/abcd_data
export RESULTS=~/abcd/workspace/mding/abcd_imputed

for filename in $MODELS/*.db; do
  fn=`basename $filename .db`
  mkdir -p $RESULTS/with_sc_blood38/$fn
  echo $filename
  #Run MetaXcan
  python3 $METAXCAN/Predict.py \
  --model_db_path $filename \
  --vcf_genotypes $DATA/abcd_genotypes38_varid_tab.vcf \
  --vcf_mode genotyped \
  --prediction_output $RESULTS/with_sc_blood38/$fn/abcd_pred_${SLURM_ARRAY_TASK_ID}.txt  \
  --prediction_summary_output $RESULTS/with_sc_blood38/$fn/summary_abcd_pred_${SLURM_ARRAY_TASK_ID}.txt \
  --verbosity 9 \
  --throw
done