#!/usr/bin/env bash

MY_PWD=$1
cd $MY_PWD"/data"

mkdir raw_data
cd raw_data


curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/058/SRR10530358/SRR10530358_1.fastq.gz -o SRR10530358_GSM4192305_Wild_type_Rep_2_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/058/SRR10530358/SRR10530358_2.fastq.gz -o SRR10530358_GSM4192305_Wild_type_Rep_2_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/051/SRR10530351/SRR10530351_1.fastq.gz -o SRR10530351_GSM4192298_WT_Ino80_Rapamycin_90_min_Rep_1_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/051/SRR10530351/SRR10530351_2.fastq.gz -o SRR10530351_GSM4192298_WT_Ino80_Rapamycin_90_min_Rep_1_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/052/SRR10530352/SRR10530352_1.fastq.gz -o SRR10530352_GSM4192299_WT_Ino80_Rapamycin_90_min_Rep_2_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/052/SRR10530352/SRR10530352_2.fastq.gz -o SRR10530352_GSM4192299_WT_Ino80_Rapamycin_90_min_Rep_2_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/057/SRR10530357/SRR10530357_1.fastq.gz -o SRR10530357_GSM4192304_Wild_type_Rep_1_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/057/SRR10530357/SRR10530357_2.fastq.gz -o SRR10530357_GSM4192304_Wild_type_Rep_1_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/073/SRR10530373/SRR10530373_1.fastq.gz -o SRR10530373_GSM4192320_TKO_Rep_2_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/073/SRR10530373/SRR10530373_2.fastq.gz -o SRR10530373_GSM4192320_TKO_Rep_2_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/046/SRR10530346/SRR10530346_1.fastq.gz -o SRR10530346_GSM4192293_TKO_Ino80_Rapamycin_90_min_Rep_1_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/046/SRR10530346/SRR10530346_2.fastq.gz -o SRR10530346_GSM4192293_TKO_Ino80_Rapamycin_90_min_Rep_1_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/048/SRR10530348/SRR10530348_1.fastq.gz -o SRR10530348_GSM4192295_TKO_Ino80_Rapamycin_90_min_Rep_2_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/048/SRR10530348/SRR10530348_2.fastq.gz -o SRR10530348_GSM4192295_TKO_Ino80_Rapamycin_90_min_Rep_2_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/072/SRR10530372/SRR10530372_1.fastq.gz -o SRR10530372_GSM4192319_TKO_Rep_1_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/072/SRR10530372/SRR10530372_2.fastq.gz -o SRR10530372_GSM4192319_TKO_Rep_1_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
