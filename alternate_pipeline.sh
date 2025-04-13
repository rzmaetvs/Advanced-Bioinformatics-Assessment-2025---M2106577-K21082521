SECTION 2.6
# ------------------ ALTERNATIVE TOOL: HISAT2 ------------------

conda install -y -c bioconda hisat2
#This command installs the HISAT2 tool using Conda from the Bioconda channel, which is a popular and reliable source for bioinformatics software. Conda ensures that dependencies are handled and provides an easy way to install HISAT2 without worrying about compatibility or version issues.
#Installing HISAT2 via Conda is convenient for reproducibility, as Conda can handle specific versions of software and dependencies, ensuring that your pipeline will work as expected on different systems or environments. If you are using the same version of HISAT2 across different machines, the results should be consistent.

# Build HISAT2 index
hisat2-build hg19.fa hg19_index
#The hisat2-build command is used to generate an index of the reference genome (in this case, hg19.fa), which is necessary for efficient read alignment. HISAT2 uses an index to quickly align the reads to the genome, and building the index is a one-time process (unless the reference genome changes).
#Building an index is crucial for the alignment process. If the index is not properly generated or if you use a different reference genome or annotation file, the alignments could be incorrect. Additionally, if the reference genome is mismatched (e.g., hg19 vs hg38), the results will be impacted due to differences in sequence and annotations.


# Align reads
hisat2 -x hg19_index -1 NGS0001.R1.trimmed.fastq -2 NGS0001.R2.trimmed.fastq -S hisat2.sam
#The -x hg19_index option specifies the reference genome index file that HISAT2 will use for alignment (which was built earlier).
#The -1 NGS0001.R1.trimmed.fastq and -2 NGS0001.R2.trimmed.fastq options specify the input files for the paired-end reads.
#The -S hisat2.sam option tells HISAT2 to output the results in the SAM file format, which is a standard for storing alignment data.
#The paired-end reads (-1 for the forward read, -2 for the reverse read) are aligned against the reference genome using HISAT2’s algorithm. HISAT2 is known for being very fast and efficient with aligning large datasets, especially when there are spliced alignments (e.g., RNA-Seq data). If the reference genome is poorly indexed or mismatched, this could lead to lower-quality alignments or misaligned reads. The choice of SAM format allows for easier post-processing, but using other formats like BAM directly can avoid the need for a conversion step.

# Convert to BAM and continue rest of pipeline
samtools view -bS hisat2.sam | samtools sort -o hisat2.sorted.bam
#The samtools view -bS hisat2.sam command converts the SAM file into the compressed BAM format. The -b option tells samtools to output a BAM file, and -S specifies that the input is in SAM format.
#The samtools sort -o hisat2.sorted.bam command sorts the resulting BAM file by position in the genome, which is necessary for downstream analysis (e.g., variant calling, visualization).
#Converting SAM to BAM is crucial because BAM files are compressed and indexed, making them faster to work with in subsequent steps (e.g., variant calling).
#Sorting the BAM file ensures that the reads are in the correct order (by genomic coordinates), which is essential for various downstream analyses like variant calling and viewing alignments in genome browsers.
#If you skip sorting, tools like bcftools or GATK (for variant calling) may not function properly, as they rely on sorted BAM files for efficient processing.

#Pipeline with integrated alternate tool
#!/bin/bash
#SECTION 2.1
# --------------- set up environment -----------

#create and activate conda environment
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
PREFIX=/home/ubuntu/miniconda
export PATH="$HOME/miniconda/bin:$PATH"

#create conda environment for pipeline
conda create -y -n ngs_pipeline python=3.10

#install tools
conda install -y -c bioconda fastqc trimmomatic bwa samtools picard freebayes bedtools snpeff hisat2
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

# Build HISAT2 index
hisat2-build hg19.fa hg19_index
#The hisat2-build command generates an index for HISAT2 to align reads against.

#Step 2: HISAT2 alignment (using the trimmed reads)
hisat2 -x hg19_index -1 NGS0001.R1.trimmed.fastq -2 NGS0001.R2.trimmed.fastq -S hisat2.sam
#This aligns the trimmed paired-end reads using the HISAT2 index generated for hg19.

#Step 3: Convert SAM to BAM and sort
samtools view -bS hisat2.sam | samtools sort -o hisat2.sorted.bam
#This converts the SAM file generated by HISAT2 to BAM format and sorts the BAM file by genomic coordinates.

#Step 4: Mark duplicates
picard MarkDuplicates I=hisat2.sorted.bam O=hisat2.dedup.bam M=dup_metrics.txt

#Step 5: Filter by quality
samtools view -q 20 -b hisat2.dedup.bam > hisat2.filtered.bam

#Step 6: Index the final BAM
samtools index hisat2.filtered.bam

#Step 7: Alignment statistics
samtools flagstat hisat2.filtered.bam > flagstat.txt
samtools idxstats hisat2.filtered.bam > idxstats.txt
samtools depth hisat2.filtered.bam > depth.txt
picard CollectInsertSizeMetrics I=hisat2.filtered.bam O=insert_size_metrics.txt H=insert_size_histogram.pdf M=0.5

#SECTION 2.4
#------------------ VARIANT CALLING ------------------
#Step 1: Call variants using FreeBayes restricted to BED file
freebayes -f hg19.fa -t annotation.bed hisat2.filtered.bam > raw.vcf

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


