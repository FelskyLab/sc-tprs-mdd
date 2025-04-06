#!/bin/bash
#SBATCH --mail-type=all
#SBATCH --mail-user=megan.ding@mail.utoronto.ca
#SBATCH --job-name=modify_abcd
#SBATCH --output=/external/rprshnas01/kcni/mding/abcd/workspace/mding/output/mod_abcd-%A.%a.out
#SBATCH --error=/external/rprshnas01/kcni/mding/abcd/workspace/mding/error/mod_abcd-%A.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8192MB
#SBATCH --array=1
#SBATCH --time=0-400:00:00


module purge


#paths to env variables
export DATA=/external/rprshnas01/kcni/mding/abcd/workspace/mding

#format: chr21_17290620_C_T_b38

#see a single row: sed '32q;d' /external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/abcd_genotypes38_varid_tab2.vcf

# get a single column: awk '$3 == $3 {printf "%s\n",$3}' abcd_genotypes.vcf > abcd_rsids.txt
#zcat snp151Seq.txt.gz | grep -Fwf abcd_data/abcd_rsids.txt > abcd_data/abcd_pos.txt
#----------------------------------------------

/usr/bin/wget --directory-prefix=$DATA/abcd_data ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/snp151.txt.gz
/usr/bin/awk 'BEGIN { OFS="\t" }{print $5,$4"\n"}' <(/usr/bin/zcat $DATA/abcd_data/snp151.txt.gz) > $DATA/abcd_data/rsid_pos38.txt

/external/rprshnas01/kcni/mding/plink_files/plink \
--bfile $DATA/abcd_data/genotype_QCed/ABCD_release_3.0_QCed \
--update-map $DATA/abcd_data/rsid_pos38.txt \
--make-bed \
--out $DATA/abcd_data/genotype_QCed/ABCD_release_3.0_QCed38

/external/rprshnas01/kcni/mding/plink_files/plink \
--bfile $DATA/abcd_data/genotype_QCed/ABCD_release_3.0_QCed38 \
--recode vcf --out $DATA/abcd_data/abcd_genotypes38

/usr/bin/awk 'NR>31 { $3 = sprintf("chr%s_%s_%s_%s_b38",$1,$2,$4,$5) } 1' $DATA/abcd_data/abcd_genotypes38.vcf > $DATA/abcd_data/abcd_genotypes38_varid.vcf
/usr/bin/sed 's/[[:space:]]\+/\t/g' $DATA/abcd_data/abcd_genotypes38_varid.vcf > $DATA/abcd_data/abcd_genotypes38_varid_tab.vcf
#------------------------------------------------