#!/usr/bin/env bash

MY_PWD=$1
cd $MY_PWD"/data"

mkdir raw_data
cd raw_data


curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/058/SRR10530358/SRR10530358_1.fastq.gz -o SRR10530358_GSM4192305_Wild_type_Rep_2_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/058/SRR10530358/SRR10530358_2.fastq.gz -o SRR10530358_GSM4192305_Wild_type_Rep_2_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/040/SRR10530340/SRR10530340_1.fastq.gz -o SRR10530340_GSM4192287_TKO_Rpb1_Rapamycin_120_min_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/040/SRR10530340/SRR10530340_2.fastq.gz -o SRR10530340_GSM4192287_TKO_Rpb1_Rapamycin_120_min_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/042/SRR10530342/SRR10530342_1.fastq.gz -o SRR10530342_GSM4192289_TKO_Rpb1_Ino80_Rapamycin_120_min_Rep_1_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/042/SRR10530342/SRR10530342_2.fastq.gz -o SRR10530342_GSM4192289_TKO_Rpb1_Ino80_Rapamycin_120_min_Rep_1_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/039/SRR10530339/SRR10530339_1.fastq.gz -o SRR10530339_GSM4192286_TKO_Rpb1_Vehicle_120_min_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/039/SRR10530339/SRR10530339_2.fastq.gz -o SRR10530339_GSM4192286_TKO_Rpb1_Vehicle_120_min_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/044/SRR10530344/SRR10530344_1.fastq.gz -o SRR10530344_GSM4192291_TKO_Rpb1_Ino80_Rapamycin_120_min_Rep_2_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/044/SRR10530344/SRR10530344_2.fastq.gz -o SRR10530344_GSM4192291_TKO_Rpb1_Ino80_Rapamycin_120_min_Rep_2_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/041/SRR10530341/SRR10530341_1.fastq.gz -o SRR10530341_GSM4192288_TKO_Rpb1_Ino80_Vehicle_120_min_Rep_1_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/041/SRR10530341/SRR10530341_2.fastq.gz -o SRR10530341_GSM4192288_TKO_Rpb1_Ino80_Vehicle_120_min_Rep_1_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/043/SRR10530343/SRR10530343_1.fastq.gz -o SRR10530343_GSM4192290_TKO_Rpb1_Ino80_Vehicle_120_min_Rep_2_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/043/SRR10530343/SRR10530343_2.fastq.gz -o SRR10530343_GSM4192290_TKO_Rpb1_Ino80_Vehicle_120_min_Rep_2_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/057/SRR10530357/SRR10530357_1.fastq.gz -o SRR10530357_GSM4192304_Wild_type_Rep_1_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/057/SRR10530357/SRR10530357_2.fastq.gz -o SRR10530357_GSM4192304_Wild_type_Rep_1_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/072/SRR10530372/SRR10530372_1.fastq.gz -o SRR10530372_GSM4192319_TKO_Rep_1_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/072/SRR10530372/SRR10530372_2.fastq.gz -o SRR10530372_GSM4192319_TKO_Rep_1_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/073/SRR10530373/SRR10530373_1.fastq.gz -o SRR10530373_GSM4192320_TKO_Rep_2_Saccharomyces_cerevisiae_MNase-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/073/SRR10530373/SRR10530373_2.fastq.gz -o SRR10530373_GSM4192320_TKO_Rep_2_Saccharomyces_cerevisiae_MNase-Seq_2.fastq.gz
