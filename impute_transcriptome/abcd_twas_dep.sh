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
#SBATCH --array=0-21
#SBATCH --time=0-400:00:00


module purge
source /external/rprshnas01/netdata_kcni/dflab/.bashrc.d/miniforge.bashrc
mamba activate metaxcan

#paths to env variables
export METAXCAN=~/sc-tprs-mdd/MetaXcan/software
export MODELS=~/abcd/workspace/mding/immune_cell_models
export DATA=~/abcd/workspace/mding/abcd_data
export RESULTS=~/abcd/workspace/mding/abcd_imputed

arr=($MODELS/*.db)

fn=`basename ${arr[${SLURM_ARRAY_TASK_ID}]} .db`
# mkdir -p $RESULTS/with_sc_blood38/$fn
echo ${arr[${SLURM_ARRAY_TASK_ID}]}
#Run MetaXcan
python3 $METAXCAN/PrediXcanAssociation.py \
--expression_file $RESULTS/with_sc_blood38/$fn/abcd_pred_adtrain.txt \
--input_phenos_file $DATA/abcd_ad_train.txt \
--input_phenos_column anxdep_t \
--covariates_file $DATA/abcd_ad_train.txt \
--covariates genetic_pc_1	genetic_pc_2	genetic_pc_3	genetic_pc_4	genetic_pc_5	genetic_pc_6	\
genetic_pc_7	genetic_pc_8	genetic_pc_9	genetic_pc_10 age male site22 site02 site03 site04 site05 \
site06 site07 site08 site09 site10 site11 site12 site13 site14 site15 site16 site17 site18 site19 site20 site21 \
--output $RESULTS/with_sc_blood38/$fn/abcd_anxdep_assoc_covars.txt \
--verbosity 9 \
--throw

python3 $METAXCAN/PrediXcanAssociation.py \
--expression_file $RESULTS/with_sc_blood38/$fn/abcd_pred_wdtrain.txt \
--input_phenos_file $DATA/abcd_wd_train.txt \
--input_phenos_column withdep_t \
--covariates_file $DATA/abcd_wd_train.txt \
--covariates genetic_pc_1	genetic_pc_2	genetic_pc_3	genetic_pc_4	genetic_pc_5	genetic_pc_6	\
genetic_pc_7	genetic_pc_8	genetic_pc_9	genetic_pc_10 age male site22 site02 site03 site04 site05 \
site06 site07 site08 site09 site10 site11 site12 site13 site14 site15 site16 site17 site18 site19 site20 site21 \
--output $RESULTS/with_sc_blood38/$fn/abcd_withdep_assoc_covars.txt \
--verbosity 9 \
--throw




# for filename in $MODELS/*.db; do
#   rm /external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/with_sc_blood38/`basename $filename .db`/abcd_anxdep_assoc_1.txt
# done