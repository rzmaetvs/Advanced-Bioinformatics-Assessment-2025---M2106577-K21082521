#!/bin/bash
#SECTION2.1
# --------------- set up environment -----------

#create and activate conda environment
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
PREFIX=/home/ubuntu/miniconda
export PATH="$HOME/miniconda/bin:$PATH"

#create conda environment for pipeline
conda create -y -n ngs_pipeline python=3.10

#install tools
conda install -y -c bioconda fastqc trimmomatic bwa samtools picard freebayes bedtools snpeff
conda install annovar

#Create data directory
mkdir -p ~/ngs_project/data && cd ~/ngs_project/data

#Download the FastQ files and annotation BED file
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

#rename files and unzip
mv NGS0001.R1.fastq.qz NGS0001.R1.fastq.gz
mv NGS0001.R2.fastq.qz NGS0001.R2.fastq.gz
gunzip NGS0001.R1.fastq.gz
gunzip NGS0001.R2.fastq.gz


#SECTION 2.2
#Step 1. Quality check raw reads
mkdir -p qc_reports
fastqc NGS0001.R1.fastq NGS0001.R2.fastq -o qc_reports

#Step 2: Trimming using Trimmomatic
trimmomatic PE \
  NGS0001.R1.fastq NGS0001.R2.fastq \
  NGS0001.R1.trimmed.fastq NGS0001.R1.unpaired.fastq \
  NGS0001.R2.trimmed.fastq NGS0001.R2.unpaired.fastq \
  ILLUMINACLIP:/home/ubuntu/anaconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
#Step 3: Quality check trimmed reads
mkdir -p qc_trimmed_reports
fastqc NGS0001.R1.trimmed.fastq NGS0001.R2.trimmed.fastq -o qc_trimmed_reports

#SECTION 2.3
#------------------ ALIGNMENT ------------------

#Step 1: Download and index hg19 reference
cd ~/ngs_project/data
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
bwa index hg19.fa
samtools faidx hg19.fa

#Create sequence dictionary for Picard
picard CreateSequenceDictionary R=hg19.fa O=hg19.dict

#Step 2: BWA MEM alignment (include read group)
bwa mem -M -R "@RG\tID:sample1\tSM:NGS0001\tPL:ILLUMINA" hg19.fa \
  NGS0001.R1.trimmed.fastq NGS0001.R2.trimmed.fastq > aln.sam

#Step 3: Convert SAM to BAM and sort
samtools view -bS aln.sam | samtools sort -o aln.sorted.bam

#Step 4: Mark duplicates
picard MarkDuplicates I=aln.sorted.bam O=aln.dedup.bam M=dup_metrics.txt

#Step 5: Filter by quality
samtools view -q 20 -b aln.dedup.bam > aln.filtered.bam

#Step 6: Index the final BAM
samtools index aln.filtered.bam

#Step 7: Alignment statistics
samtools flagstat aln.filtered.bam > flagstat.txt
samtools idxstats aln.filtered.bam > idxstats.txt
samtools depth aln.filtered.bam > depth.txt
picard CollectInsertSizeMetrics I=aln.filtered.bam O=insert_size_metrics.txt H=insert_size_histogram.pdf M=0.5

#SECTION 2.4
#------------------ VARIANT CALLING ------------------
#Step 1: Call variants using FreeBayes restricted to BED file
freebayes -f hg19.fa -t annotation.bed aln.filtered.bam > raw.vcf

#Step 2: Quality filtering
vcffilter -f "QUAL > 20 & DP > 10" raw.vcf > filtered.vcf

#SECTION 2.5
# ------------------ VARIANT ANNOTATION ------------------
# Step 1: ANNOVAR annotation
cd annovar
chmod +x annotate_variation.pl
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/

ls /home/ubuntu/annovar.latest/annovar/
export PATH=$PATH:/home/ubuntu/annovar.latest/annovar
source ~/.bashrc  
chmod +x convert2annovar.pl
cd ~/ngs_project/data
convert2annovar.pl -format vcf4 filtered.vcf > filtered.avinput

cd /home/ubuntu/annovar.latest/annovar
perl annotate_variation.pl -downdb -buildver hg19 refGene humandb/
annotate_variation.pl --buildver hg19 --downdb seq humandb/hg19_seq
retrieve_seq_from_fasta.pl humandb/hg19_refGene.txt -seqdir humandb/hg19_seq -format refGene -outfile humandb/hg19_refGeneMrna.fa
chmod +x /home/ubuntu/annovar.latest/annovar/retrieve_seq_from_fasta.pl
retrieve_seq_from_fasta.pl humandb/hg19_refGene.txt -seqdir humandb/hg19_seq -format refGene -outfile humandb/hg19_refGeneMrna.fa


mkdir ~/ngs_project/annovar
cd ~/ngs_project/annovar
perl /home/ubuntu/annovar.latest/annovar/annotate_variation.pl -downdb -buildver hg19 refGene /home/ubuntu/ngs_project/annovar/humandb/
export PATH=$PATH:/home/ubuntu/annovar.latest/annovar

#dbNSFP (functional predictions for missense variants)
#ClinVar (clinical significance)

cd ~/ngs_project/data
chmod +x /home/ubuntu/annovar.latest/annovar/*.pl
perl /home/ubuntu/annovar.latest/annovar/table_annovar.pl \
filtered.avinput \
/home/ubuntu/annovar.latest/annovar/humandb/ \
-buildver hg19 \
-out annovar_output \
-remove \
-protocol refGene,dbnsfp31a_interpro,clinvar_20180603 \
-operation g,f,f \
-nastring .

# Step 2: snpEFF annotation
conda install -c bioconda snpeff
#System Java 11 needed to run snpEff 
snpEff download -v hg19
snpEff hg19 filtered.vcf > snpeff_annotated.vcf

# Step 3: Variant prioritisation
grep -v "dbSNP" annovar_output.hg19_multianno.txt | grep "exonic" > prioritized_variants.txt


