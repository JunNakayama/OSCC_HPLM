# QCcheck

fastqc -t 8 ${sample}_S1_L001_R1_001.fastq.gz -o QC

# kallisto
kallisto index -i hs_kallisto_index Homo_sapiens.GRCh38.dna.alt.fa
kallisto quant -t 10 -i hs_kallisto_index -o ${sample} ./FASTQ/${sample}_S1_L001_R1_001.fastq.gz ./FASTQ/${sample}_S1_L001_R1_001.fastq.gz


