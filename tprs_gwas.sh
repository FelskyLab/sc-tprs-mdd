#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=megan.ding@mail.utoronto.ca
#SBATCH --job-name=tprs_gwas
#SBATCH --output=/external/rprshnas01/external_data/abcd/workspace/mding/output/slurm-impute_abcd-%A.%a.out
#SBATCH --error=/external/rprshnas01/external_data/abcd/workspace/mding/error/slurm-impute_abcd-%A.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=8192MB
#SBATCH --array=1
#SBATCH --time=0-400:00:00


module purge

#paths to env variables
export MODELS=~/abcd/workspace/mding/immune_cell_models
i=1
for filename in $MODELS/*.db; do
  fn=`basename $filename .db`
  i=`expr $i + 1`
  /usr/bin/time /external/rprshnas01/kcni/mding/plink2 \
  --pfile /external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/abcd_geno_binaries/abcd_geno_fileset \
  --pheno-col-nums $i \
  --glm allow-no-covars \
  --pheno /external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/tprs_pheno.txt iid-only --no-psam-pheno \
  --pfilter 0.05 \
  --out /external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/gwas/no_covar \
  --no-input-missing-phenotype \
  --extract $MODELS/snps/${fn}_uw_from_gene.txt

  /usr/bin/time /external/rprshnas01/kcni/mding/plink2 \
  --pfile /external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/abcd_geno_binaries/abcd_geno_fileset \
  --pheno-col-nums $i \
  --glm hide-covar --covar /external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/tprs_pheno_genpc.txt iid-only \
  --pheno /external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/tprs_pheno.txt iid-only --no-psam-pheno \
  --pfilter 0.05 \
  --out /external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/gwas/w_covar \
  --no-input-missing-phenotype \
  --extract $MODELS/snps/${fn}_uw_from_gene.txt
done
# /external/rprshnas01/kcni/mding/plink2 \
# --vcf /external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/abcd_genotypes38_varid_tab.vcf \
# --out /external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/abcd_geno_fileset \
# --make-pgen \
# --psam /external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/abcd_geno_fileset-temporary.psam