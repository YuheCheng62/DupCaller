# DupCaller (IN DEVELOPMENT. PLEASE DON'T USE FOR PRODUCTION PURPOSE)
DupCaller is a somatic variant caller for DNA duplex-sequencing.

## Installation and dependencies
The scripts requires **Python 3** to run abd has been tested with Python3.7+. After cloning the repository to your local directory, go to the repository and run the following command to install the scripts and dependencies:
```
python setup.py install
```
The utility scripts should be automatically add to your path after installation.

The pipeline also requires **samtools** and **bwa** to be installed.

## Usage
### Trimming of fastq files
The pipeline of DupCaller starts from the raw paired-end fastq files from the sequencing platform. Adaptor trimmed fastqs are NOT preferred and has not been tested. Use the follwoing command to trim the fastqs and save barcoding informations:
```
DupCallerTrim -i {sample}_1.fastq -i2 {sample}_2.fastq -o {sample}_trimmed -p XXXNNNN
```
The command will output the trimmed fastq files as *{sample}_trimmed_1.fastq* and *{sample}_trimmed_2.fastq*. Indicate the pattern of your barcoding with the option -p, where X represents a barcode base and N represents a skipped base. In the example, the first three bases will be the barcode and the next four bases will be skipped.
### Aligning with bwa
Align the trimmed fastq files with bwa and sorted by coordinate:
```
bwa -p threads hg38.fa {sample}_trimmed_1.fastq {sample}_trimmed_2.fastq | samtools sort -@ {threads} > {sample}.bam
```
### Somatic Mutation Calling from aligned sample bam and matched normal
After alignment is finished, run the following command to call mutations:
```
DupCallerCall -b {sample}.bam -g af-only-gnomad.hg38.vcf.gz -f hg38.fa -o PD43272 -p 28 -n {normal}.bam -m {noise}.bed -rl 151
```
The command will output a vcf file and a table of the number of duplex groups for each duplex read counts from 2 to 30.


