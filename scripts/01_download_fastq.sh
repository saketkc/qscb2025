#!/bin/bash
set -eox pipefail
wget --content-disposition -c "http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR173/020/SRR17375020/SRR17375020_1.fastq.gz" -o GSM5764349_1.fastq.gz
wget --content-disposition -c "http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR173/020/SRR17375020/SRR17375020_2.fastq.gz" -o GSM5764349_2.fastq.gz


#wget --content-disposition -c "http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR173/050/SRR17375050/SRR17375050_1.fastq.gz" -o GSM5764393_1.fastq.gz
#wget --content-disposition -c "http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR173/050/SRR17375050/SRR17375050_2.fastq.gz" -o GSM5764393_2.fastq.gz
