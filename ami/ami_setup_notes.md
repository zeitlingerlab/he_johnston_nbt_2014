# Amazon EC2 AMI setup notes

We start with the Ubuntu 14.04 EBS-backed AMI for the us-east-1 region.

## Boot AMI

AMI was started as type **c3.4xlarge** with a 150GB EBS boot volume. This instance type has 16 CPU cores and ~30GB of memory.

## Set up data directory

```bash
sudo mkdir /data
sudo chown -R ubuntu:ubuntu /data
```

## Install R 3.1.1 and Bioconductor 2.15

Create /etc/apt/sources.list.d/cran.list and add the line:

```
deb http://cran.rstudio.com/bin/linux/ubuntu trusty/
```

Execute the following to add the appropriate signing key:

```bash
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
```

Update packages and install R:

```bash
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install r-base r-base-dev
sudo apt-get install libxml2-dev libcurl3-dev # for Rsamtools
sudo apt-get install htop samtools pigz fastx-toolkit sra-toolkit parallel git bowtie bedtools python-pip python-numpy haskell-platform procmail
sudo cabal update
sudo cabal install pandoc --global
```

Create ~/.Rprofile:

```S
options(repos=c(CRAN="http://cran.rstudio.com"))

q <- function(save = "no", status = 0, runLast = TRUE) {
  .Internal(quit(save, status, runLast))
}

update_bioc = function() {
  library(BiocInstaller)
  update.packages(repos=biocinstallRepos(), ask=FALSE)
}

update_all = function() {
  message("Updating R packages...")
  update.packages(ask=F)
  message("Updating Bioconductor packages...")
  update_bioc()
}
````

Install R and Bioconductor packages by sourcing the following script from a root R session:

```S
update.packages(ask=F)
install.packages(c("gtools", "reshape", "ggplot2", "xtable", 
                   "data.table", "matrixStats", "stringr", "optparse", "xts",
                   "dplyr", "knitr", "knitcitations"))

source("http://bioconductor.org/biocLite.R")
biocLite()

bioc_packages <- c("GenomicRanges", 
                   "BSgenome.Dmelanogaster.UCSC.dm3", 
                   "BSgenome.Hsapiens.UCSC.hg19",
                   "TxDb.Hsapiens.UCSC.hg19.knownGene",
                   "TxDb.Dmelanogaster.UCSC.dm3.ensGene",
                   "rtracklayer", 
                   "ShortRead",
                   "Rsamtools",
                   "chipseq",
                   "VariantAnnotation",
                   "Gviz",
                   "seqLogo")

biocLite(bioc_packages)
````

## Download analysis code

```bash
cd /home/ubuntu
git clone https://github.com/zeitlingerlab/TODO analysis_code
```

## Install cutadapt

```bash
mkdir /data/software
cd /data/software
git clone https://github.com/marcelm/cutadapt.git
cd cutadapt
git checkout v1.3
sudo apt-get install cython
sudo python setup.py install
```

## Install SRA tools

```bash
mkdir /data/software/sra
cd /data/software/sra
wget 'http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.4.2-1/sratoolkit.2.4.2-ubuntu64.tar.gz'
tar xf sratoolkit.2.4.2-ubuntu64.tar.gz
/data/software/sra/sratoolkit.2.4.2-ubuntu64/bin/vdb-config -i
```

## Build UCSC dm3 reference genome for Bowtie

```bash
mkdir /data/genomes
cd /data/genomes
mkdir dm3
cd dm3
wget 'http://hgdownload.cse.ucsc.edu/goldenPath/dm3/bigZips/chromFa.tar.gz'
tar xf chromFa.tar.gz
rm chromFa.tar.gz
cat *.fa > dm3.fasta
rm *.fa
bowtie-build -o 1 dm3.fasta dm3
```

## Build UCSC hg19 reference genome for Bowtie

```bash
cd /data/genomes
mkdir hg19
cd hg19
wget 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz'
tar xf chromFa.tar.gz
rm chromFa.tar.gz
cat *.fa > hg19.fasta
rm *.fa
bowtie-build -o 1 hg19.fasta hg19
pigz -vp 8 hg19.fasta
```

## Download ChIP-exo TBP K562 data (Pugh lab)

```bash
mkdir /data/sra
cd /data/sra
wget 'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR770/SRR770743/SRR770743.sra'
wget 'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR770/SRR770744/SRR770744.sra'
parallel -uj 1 /data/software/sra/sratoolkit.2.4.2-ubuntu64/bin/fastq-dump --gzip /data/sra/{} ::: *.sra
rm *.sra
ln -s SRR770743.fastq.gz hsap_k562_tbp_chipexo_01.fastq.gz
ln -s SRR770744.fastq.gz hsap_k562_tbp_chipexo_02.fastq.gz
```

## Preprocess ChIP-nexus FASTQ reads

```bash
mkdir /data/preprocessed_fastq
cd /data/preprocessed_fastq
parallel -uj 6 Rscript /data/analysis_code/pipeline/preprocess_fastq.r -p 3 -f {} -k 22 -b CTGA -r 5 -o \`basename {} .fastq.gz\`_processed.fastq.gz ::: /data/raw_fastq/chipnexus/*.fastq.gz
```

## Oregon-R reference genome modification

Align Dorsal and Twist embryo ChIP-seq samples with a maximum of 3 mismatches to identify as many SNPs as possible:

```bash
mkdir /data/genomes/dm3_oregonr
cd $_
zcat /data/raw_fastq/chipseq/dmel_embryo* | pigz -p 16 > embryo_combined_chipseq.fastq.gz
/data/analysis_code/pipeline/align_chipseq_3mm.sh embryo_combined_chipseq.fastq.gz /data/genomes/dm3/dm3
/data/analysis_code/pipeline/sort_bam embryo_combined_chipseq.bam
```

Call variants using **samtools**:

```bash
samtools mpileup -uD -f /data/genomes/dm3/dm3.fasta embryo_combined_chipseq.bam | bcftools view -vcg - | gzip > variants.vcf.gz
```

Apply variants to reference:

```bash
Rscript /data/analysis_code/pipeline/apply_variants_to_reference.r -f variants.vcf.gz
```

Build **bowtie** index:

```bash
bowtie-build -o 1 dm3_oregonr.fasta dm3_oregonr
```

## Align ChIP-nexus reads

```bash
cd /data/bam

/data/analysis_code/pipeline/align_chipnexus.sh /data/preprocessed_fastq/dmel_embryo_dl_chipnexus_01_processed.fastq.gz /data/genomes/dm3_oregonr/dm3_oregonr
/data/analysis_code/pipeline/align_chipnexus.sh /data/preprocessed_fastq/dmel_embryo_dl_chipnexus_02_processed.fastq.gz /data/genomes/dm3_oregonr/dm3_oregonr
/data/analysis_code/pipeline/align_chipnexus.sh /data/preprocessed_fastq/dmel_embryo_twi_chipnexus_01_processed.fastq.gz /data/genomes/dm3_oregonr/dm3_oregonr
/data/analysis_code/pipeline/align_chipnexus.sh /data/preprocessed_fastq/dmel_embryo_twi_chipnexus_02_processed.fastq.gz /data/genomes/dm3_oregonr/dm3_oregonr
/data/analysis_code/pipeline/align_chipnexus.sh /data/preprocessed_fastq/dmel_s2_max_chipnexus_01_processed.fastq.gz /data/genomes/dm3/dm3
/data/analysis_code/pipeline/align_chipnexus.sh /data/preprocessed_fastq/dmel_s2_max_chipnexus_02_processed.fastq.gz /data/genomes/dm3/dm3
/data/analysis_code/pipeline/align_chipnexus.sh /data/preprocessed_fastq/dmel_s2_myc_chipnexus_01_processed.fastq.gz /data/genomes/dm3/dm3
/data/analysis_code/pipeline/align_chipnexus.sh /data/preprocessed_fastq/dmel_s2_myc_chipnexus_02_processed.fastq.gz /data/genomes/dm3/dm3

/data/analysis_code/pipeline/align_chipnexus.sh /data/preprocessed_fastq/hsap_k562_tbp_chipnexus_01_processed.fastq.gz /data/genomes/hg19/hg19
/data/analysis_code/pipeline/align_chipnexus.sh /data/preprocessed_fastq/hsap_k562_tbp_chipnexus_02_processed.fastq.gz /data/genomes/hg19/hg19
parallel -uj 4 /data/analysis_code/pipeline/sort_bam {} ::: *chipnexus*.bam
```

## Align ChIP-exo reads

```bash
cd /data/bam
/data/analysis_code/pipeline/align_chipexo.sh /data/sra/hsap_k562_tbp_chipexo_01.fastq.gz /data/genomes/hg19/hg19
/data/analysis_code/pipeline/align_chipexo.sh /data/sra/hsap_k562_tbp_chipexo_02.fastq.gz /data/genomes/hg19/hg19
parallel -uj 1 /data/analysis_code/pipeline/sort_bam {} ::: hsap*chipexo*.bam
```

## Process aligned ChIP-nexus reads

```bash
mkdir /data/rdata
cd /data/rdata
parallel -uj 2 Rscript /data/analysis_code/pipeline/process_bam.r -f {} -n \`basename {} .bam\` ::: ../bam/*chipnexus*.bam
```

## Process aligned ChIP-exo reads

```bash
cd /data/rdata
parallel -uj 2 Rscript /data/analysis_code/pipeline/bam_to_granges.r -f {} ::: ../bam/*chipexo*.bam
```

## Generate ChIP-nexus and ChIP-exo BigWigs

```bash
cd /data/rdata
parallel -uj 6 Rscript /data/analysis_code/pipeline/split_granges_by_strand.r -r {} ::: *.rds
mv *.bw /data/bigwigs
```

## Align ChIP-seq samples

Align all ChIP-seq samples:

```bash
cd /data/bam
/data/analysis_code/pipeline/align_chipseq.sh /data/raw_fastq/chipseq/dmel_embryo_dl_chipseq_01.fastq.gz /data/genomes/dm3_oregonr/dm3_oregonr
/data/analysis_code/pipeline/align_chipseq.sh /data/raw_fastq/chipseq/dmel_embryo_twi_chipseq_01.fastq.gz /data/genomes/dm3_oregonr/dm3_oregonr

/data/analysis_code/pipeline/align_chipseq.sh /data/raw_fastq/chipseq/dmel_s2_max_chipseq_01.fastq.gz /data/genomes/dm3/dm3

/data/analysis_code/pipeline/align_chipseq.sh /data/raw_fastq/chipseq/hsap_k562_tbp_chipseq_01.fastq.gz /data/genomes/hg19/hg19

parallel -uj 2 /data/analysis_code/pipeline/sort_bam {} ::: *chipseq*.bam
````

## R object generation for ChIP-seq samples

Generate R objects for all aligned samples:

```bash
cd /data/bam
Rscript /data/analysis_code/pipeline/bamtoolsr.r -f dmel_embryo_dl_chipseq_01.bam -e 136 -n dmel_embryo_dl_chipseq_01
Rscript /data/analysis_code/pipeline/bamtoolsr.r -f dmel_embryo_twi_chipseq_01.bam -e 124 -n dmel_embryo_twi_chipseq_01
Rscript /data/analysis_code/pipeline/bamtoolsr.r -f dmel_s2_max_chipseq_01.bam    -e 83  -n dmel_s2_max_chipseq_01
Rscript /data/analysis_code/pipeline/bamtoolsr.r -f hsap_k562_tbp_chipseq_01.bam  -e 74  -n hsap_k562_tbp_chipseq_01
mv *.bw /data/bigwigs
mv *.rds /data/rdata
```

## DNA shape

Align all 1,024 pentamers to reference:

```bash
mkdir /data/dna_shape
cd /data/dna_shape
cp /data/analysis_code/internal_data/pentamers.fasta .
cp /data/analysis_code/internal_data/[HMPR]*.csv .
/data/analysis_code/pipeline/align_pentamers.sh
```

Build DNA shape BigWigs

```bash
Rscript /data/analysis_code/pipeline/build_shape_table.r 
Rscript /data/analysis_code/pipeline/build_pentamers_granges.r
Rscript /data/analysis_code/pipeline/build_dna_shape_bigwigs.r
```

## Additional R processing

```bash
cd /data/analysis_code
Rscript preprocess_samples.r
```

## Peak calling

Install MACS2

```bash
sudo pip install MACS2==2.0.10.20131216
```

Generate BED files

```bash
mkdir /data/macs
cd /data/macs
parallel -uj 8 Rscript /data/analysis_code/pipeline/granges_to_bed.r -r -f {} -b \`basename {} granges.rds\`bed -s 1 -i ::: /data/rdata/*.granges.rds
```

Run MACS 

```bash
parallel -uj 4 macs2 callpeak -f BED -t {} -g hs -n \`basename {} .bed.bgz\` --call-summits --keep-dup=all ::: hsap*chipnexus*.bed.bgz
parallel -uj 1 macs2 callpeak -f BED -t {} -g hs -n \`basename {} .bed.bgz\` --call-summits ::: hsap*chipseq*.bed.bgz

parallel -uj 4 macs2 callpeak -f BED -t {} -g dm -n \`basename {} .bed.bgz\` --call-summits --keep-dup=all ::: dmel*chipnexus*.bed.bgz
parallel -uj 4 macs2 callpeak -f BED -t {} -g dm -n \`basename {} .bed.bgz\` --call-summits ::: dmel*chipseq*.bed.bgz
gzip *_summits.bed
```


