
# <u>RNA-Seq Analysis Phase I: Quality Control
## <i>This Notebook Performs QC on Paired-End RNA-Seq Data</i>
#### Last Revision: July  2017
#### Author: Charles David and Dan Jones
#### This analysis done by Dan Jones and Karmun Chooi

This version of the workflow has been simplified for inclusion into a manuscript. It contains all code required to replicate the experiment, but does not contain information that relates to configuration specific to Plant and Food Research's system. A number of non-critical sanity checks have also been removed for clarity. However, the full workflow is also available in the relevant [github repository](https://github.com/PlantandFoodResearch/bioinf_Vitis_Nicotiana_RNAseq/)

## <u>Step I: Establish Data Management Structure on PowerPlant</u>


```bash
# Define the user as a variable
USER="hradxj"
PROJECTNAME="Vitis_Nicotiana_experiment"
# Define the project directory and temp subdirectory as a variable
PROJECT="/workspace/$USER/$PROJECTNAME"

# Define the location of various subdirectories within PROJECT

RAW=$PROJECT/000.raw
FASTQC_RAW=$PROJECT/001.fastqc_raw
SORTMERNA=$PROJECT/002.SMRNA
FASTQC_SORTMERNA=$PROJECT/003.fastqc_smrna
TRIMMOMATIC=$PROJECT/004.trimmomatic
FASTQC_TRIMMOMATIC=$PROJECT/005.fastqc_trim
TEMP="$PROJECT/TEMP"
```

### Create analysis directories

At this point, we have not actually created any directories... only defined what the directory is going to be called __when__ we create it. We still need to actually create the directories.

We use the Unix shell command `mkdir` to create the directories. The switch `-p` suppresses error messages if the directory already exists.


```bash
# Create the project directory

mkdir -p $PROJECT

# Create project subdirectories

mkdir -p $RAW
mkdir -p $FASTQC_RAW
mkdir -p $SORTMERNA
mkdir -p $FASTQC_SORTMERNA
mkdir -p $TRIMMOMATIC
mkdir -p $FASTQC_TRIMMOMATIC
mkdir -p $TEMP

```


```bash
# Create symlinks of all input fastq files and put them in $RAW

ln -s /input/genomic/viral/metagenomic/170621_150PE_HS4K2A/Almeida/*.fastq.gz $RAW
```


```bash
# Change filenames to comply with naming requirements for downstream steps
FILENAMES=$(ls $RAW)
#echo $FILENAMES
for FILE in $FILENAMES
do
NEWFILENAME=$(echo $FILE | sed 's/_001//g')
mv $RAW/$FILE $RAW/$NEWFILENAME
done
```


```bash
ls -s $RAW

```

    total 18
    2 RACP005_11_S11_L002_R1.fastq.gz  2 RACP005_1_S8_L002_R1.fastq.gz
    2 RACP005_11_S11_L002_R2.fastq.gz  2 RACP005_1_S8_L002_R2.fastq.gz
    2 RACP005_12_S12_L002_R1.fastq.gz  2 RACP005_5_S9_L002_R1.fastq.gz
    2 RACP005_12_S12_L002_R2.fastq.gz  2 RACP005_5_S9_L002_R2.fastq.gz
    2 RACP005_13_S13_L002_R1.fastq.gz  2 RACP005_8_S10_L002_R1.fastq.gz
    2 RACP005_13_S13_L002_R2.fastq.gz  2 RACP005_8_S10_L002_R2.fastq.gz


## <u>Step II Part 1: FastQC RAW Data</u>
- The input for this step is the raw data from the provider in FASTQ format
- The output from this step are the HTML FASTQC Reports


```bash
# Define the location for the QC reports:
OUT="${PROJECT}/001.fastqc_raw"
LOG="${OUT}/logs"

mkdir -p $LOG

# Define the list of files to process:
FILES=`ls ${RAW}/*.gz`

# Load the FastQC module:
module load FastQC

for file in $FILES
    do
        COMMAND="fastqc --nogroup -q -t 2 -o ${OUT} ${file}"
        bsub -o ${LOG}/FQC.out -e ${LOG}/FQC.err -J FASTQC -n 2 $COMMAND
    done

```


```bash
# Create multiQC report of FastQC results
module load MultiQC;
multiqc $OUT -o $OUT
```

The output MultiQC report is available at: https://github.com/PlantandFoodResearch/bioinf_Vitis_Nicotiana_RNAseq/blob/master/multiqc/raw_multiqc_report.html


Download the file, save as `.html`, and open in any web browser.

## <u>Step II Part 2: SortMeRNA</u>
* In this step we will remove any rRNA contamination by comparing our reads to 6 databases of known rRNA's
* We will capture the rRNA reads in case further investigation is needed
* We will output the filtered reads to use for our workflow

### Merge (interleave) paired fastq files
* The input to this step is the raw data in FASTQ format

__NOTE THAT THIS WILL WORK ONLY IF PAIRED FILES END IN "_R1.fastq.gz" and "_R2.fastq.gz"__

* SortMeRNA requires that paired-end files are merged (intereaved) prior to execution
* There is a Bash shell script that does this called `merge-paired-reads.sh` in the `/scripts` subdirectory
* NOTE! This script can NOT process zipped files!!!
  * Therefore, all files must be de-compressed prior to use...
  * This is best done with process substitution, using `<(zcat ...)`
* The output from this step are interleaved fastq files


```bash
# Define the location of the various QC programs we will be using,
# and the location of SortMeRNA rRNA databases

SMRNA="/workspace/cflcyd/software/sortmerna-2.1b"
SCRIPTS="${SMRNA}/scripts"
DB="${SMRNA}/rRNA_databases"
INDEX="${SMRNA}/index"
SORTMERNADB="${DB}/silva-bac-16s-id90.fasta,\
${INDEX}/silva-bac-16s-db:\
${DB}/silva-bac-23s-id98.fasta,\
${INDEX}/silva-bac-23s-db:\
${DB}/silva-arc-16s-id95.fasta,\
${INDEX}/silva-arc-16s-db:\
${DB}/silva-arc-23s-id98.fasta,\
${INDEX}/silva-arc-23s-db:\
${DB}/silva-euk-18s-id95.fasta,\
${INDEX}/silva-euk-18s-db:\
${DB}/silva-euk-28s-id98.fasta,\
${INDEX}/silva-euk-28s:\
${DB}/rfam-5s-database-id98.fasta,\
${INDEX}/rfam-5s-db:\
${DB}/rfam-5.8s-database-id98.fasta,\
${INDEX}/rfam-5.8s-db"
```


```bash
# Define the location for the merged files:
OUT="${PROJECT}/002.SMRNA"
LOG="${OUT}/logs"

# Define a set of unique names of the paired files,
# but excluding the _R1.fastq.gz and _R2.fastq.gz.
# This means that the variable "FILES" will consist of a unique name
# for each PAIR of paired fastq files. We can then append the 
# _R1.fastq.gz and _R2.fastq.gz suffix within the loop, to ensure
# that each iteration of the loop is working on two correctly paired files.

FILES=`basename -a ${RAW}/*.gz | sed 's/_R[1,2].fastq.gz//g'|sort -u `

for file in $FILES
     do

        file1=${file}_R1.fastq.gz
        file2=${file}_R2.fastq.gz
        COMMAND="${SCRIPTS}/Merge.sh \
                <(zcat $RAW/${file1}) \
                <(zcat $RAW/${file2}) \
                ${OUT}/${file}_MERGED.fastq"
        #echo "$COMMAND"
        bsub -o ${LOG}/MERGE.out -e ${LOG}/MERGE.err -J MERGE bash -c "${COMMAND}"
     done

### Note that the bash -c is needed to open a proper bash shell
### (instead of a bourne shell) for the processes substitution to work with OpenLava ###
```

### Run the Main SortMeRNA Program
* This is what actually does the sorting
* The input to this step are the merged fastq files
* The output are the rRNA matches and the filtered raw reads in interleaved fastq format


```bash
# # Define the location for the input and output files:
IN="${PROJECT}/002.SMRNA"
OUT="${PROJECT}/002.SMRNA"
FILTERED="${OUT}/filtered/merged"
rRNA="${OUT}/rRNA"
LOG="${OUT}/logs"

mkdir -p $rRNA
mkdir -p $LOG
mkdir -p $FILTERED

FILES=`ls ${IN}/*_MERGED.fastq`

for file in $FILES
    do
        NAME=`basename $file`
        COMMAND="${SMRNA}/sortmerna --ref ${SORTMERNADB} --reads ${file} \
                --paired_in -a 4 -m 3911 -v --log --fastx \
                --aligned ${rRNA}/${NAME}_rRNA \
                --other ${FILTERED}/${NAME}_sortmerna"
        bsub -o ${LOG}/${NAME}.out -e ${LOG}/${NAME}.err -J SMRNA -n 4 $COMMAND
     done
```
 Results:
    Total reads = 26,184,138
    Total reads passing E-value threshold = 225,953 (0.86%)
    Total reads failing E-value threshold = 25,958,185 (99.14%)
    Minimum read length = 100
    Maximum read length = 100
    Mean read length = 100
 By database:
    /workspace/cflcyd/software/sortmerna-2.1b/rRNA_databases/silva-bac-16s-id90.fasta		    0.20%
    /workspace/cflcyd/software/sortmerna-2.1b/rRNA_databases/silva-bac-23s-id98.fasta		    0.07%
    /workspace/cflcyd/software/sortmerna-2.1b/rRNA_databases/silva-arc-16s-id95.fasta		    0.00%
    /workspace/cflcyd/software/sortmerna-2.1b/rRNA_databases/silva-arc-23s-id98.fasta		    0.00%
    /workspace/cflcyd/software/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95.fasta		    0.34%
    /workspace/cflcyd/software/sortmerna-2.1b/rRNA_databases/silva-euk-28s-id98.fasta		    0.25%
    /workspace/cflcyd/software/sortmerna-2.1b/rRNA_databases/rfam-5s-database-id98.fasta		0.00%
    /workspace/cflcyd/software/sortmerna-2.1b/rRNA_databases/rfam-5.8s-database-id98.fasta		0.00%

### Un-merge (de-interleave) the filtered fastq files
* The other programs in our workflow use standard non-interleaved files, so we unmerge them.
* The input is the merged fastq files
* The output are the unmerged fastq files


```bash
IN="${PROJECT}/002.SMRNA/filtered/merged"
OUT="${PROJECT}/002.SMRNA/filtered"

FILENAMES=`ls ${IN}/*sortmerna*`

for file in $FILENAMES
     do
        # echo $file
        NAME=`basename $file`
        PREFIX=`echo $NAME | awk -F'[. ]' '{print $1}'`
        #echo $PREFIX
        file1=${OUT}/${PREFIX}_R1.fastq
        file2=${OUT}/${PREFIX}_R2.fastq
        COMMAND="${SCRIPTS}/Unmerge.sh $file $file1 $file2"
       # echo $COMMAND
        bsub -J UNMERGE -n 3 ${COMMAND}
     done

```


```bash
# Remove merged files
rm -rf ${PROJECT}/002.SMRNA/filtered/merged/*
```

## <u>Step II Part 3: FastQC SortMeRNA Filtered Output</u>
* We now verify that we did not break anything and re-check the quality of our reads after sorting
* The input for this step is the filtered data from SortMeRNA in FASTQ format
* The output from this step are the HTML FastQC Reports


```bash
IN="${PROJECT}/002.SMRNA/filtered"
OUT="${PROJECT}/003.fastqc_smrna"
LOG="${OUT}/logs"

# Get the files to check:
FILES=`ls ${IN}/*.fastq`
#echo $FILES

# Load the FastQC module:
module load FastQC

for file in $FILES
    do
        COMMAND="fastqc --nogroup -q -t 2 -o ${OUT} ${file}"
        bsub -o ${LOG}/FQC.out -e ${LOG}/FQC.err -J FASTQC -n 2 $COMMAND
    done

```


```bash
# Create multiQC report of FastQC results
module load MultiQC;
multiqc $OUT -o $OUT
```

The output MultiQC report is available at: https://github.com/PlantandFoodResearch/bioinf_Vitis_Nicotiana_RNAseq/blob/master/multiqc/rRNA_removed_multiqc_report.html


Download the file, save as `.html`, and open in any web browser.

## <u>Step II Part 4: TRIMMOMATIC</u>
* Now that the reads are filtered, we will remove adapters, over-represented sequences, and poor-quality bases from the reads
* The command specifies that bases with quality scores less than 30 will be clipped
* Also, after clipping, the min length for a read will be 50 bp
* The `Illumina.fa` file contains the TruSeq adapter sequences and homo-polymer sequences to clip
  * This file needs to be edited to contain the appropriate sequences.
* The input for this step are the SortMeRNA filtered reads
* The output are the trimmed reads


```bash
# Run the Trimmomatic program on the filtered data to remove Illumina adapters, homo-polymers, and low quality reads:
  # Note that to do this, it is necessary to edit the file containing the adapter sequences
  # to include all sequences that you wish to remove:
  # This file is called Illumina.fa and is in the 000.raw directory.

IN="${PROJECT}/002.SMRNA/filtered"
OUT="${PROJECT}/004.trimmomatic"
UNPAIRED="${OUT}/unpaired"
LOG="${OUT}/logs"

mkdir -p $IN
mkdir -p $OUT
mkdir -p $UNPAIRED
mkdir -p $LOG

# Set the path to the adapter file:
CLIP="${PROJECT}/Illumina.fa"

# Get the files to trim:
# Use echo statements to be sure that the results from awk are what you really want...
FILES=`basename -a ${IN}/*.fastq | sed 's/_R[1,2].fastq//g'|sort -u `

module load Trimmomatic

for FILE in $FILES
     do
        In_File1=${IN}/${FILE}_R1.fastq
        In_File2=${IN}/${FILE}_R2.fastq
        Out_PAIRED_1=${OUT}/${FILE}_trimmomatic_R1.fastq
        Out_UNPAIRED_1=${UNPAIRED}/${FILE}_trimmomatic_unpaired_1.fastq
        Out_PAIRED_2=${OUT}/${FILE}_trimmomatic_R2.fastq
        Out_UNPAIRED_2=${UNPAIRED}/${FILE}_trimmomatic_unpaired_2.fastq
        COMMAND="java -jar -Xms8G -Xmx8G \
                 ${TRIMMOMATIC} PE -threads 3 \
                 ${In_File1} ${In_File2} \
                 ${Out_PAIRED_1} ${Out_UNPAIRED_1} ${Out_PAIRED_2} ${Out_UNPAIRED_2} \
                 ILLUMINACLIP:${CLIP}:2:30:10 SLIDINGWINDOW:5:20 MINLEN:50"
        bsub -o ${LOG}/${PREFIX}.out -e ${LOG}/${PREFIX}.err -J TRIM -n 3 $COMMAND
     done

# It is critical to set the -X settings for Java for the program to run correctly
# Here, the VM is instantiated with 8GB of heap space, with a max of 8GB...

```

### Results Summary:
* ILLUMINACLIP: Using 1 prefix pairs, 8 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
* Quality encoding detected as phred33
* Input Read Pairs: 12,928,498 
* Both Surviving: 11,702,894 (90.52%) 
* Forward Only Surviving: 570,865 (4.42%) 
* Reverse Only Surviving: 181,229 (1.40%) 
* Dropped: 473,510 (3.66%)
* So we have over 90% of reads passing our criteria!


## <u>Step II Part 5: FASTQC of TRIMMED READS</u>
* We now verify that we did not break anything and re-check the quality of our reads after trimming
* The input for this step are the filtered trimmed reads in FASTQ format
* The output from this step are the HTML Reports


```bash
IN="${PROJECT}/004.trimmomatic"
OUT="${PROJECT}/005.fastqc_trim"
LOG="${OUT}/logs"

# Get the files to check:
FILES=`ls ${IN}/*trimmomatic*`

# Load the FastQC module:
module load FastQC

for file in $FILES
    do
        COMMAND="fastqc --nogroup -q -t 2 -o ${OUT} ${file}"
        bsub -o ${LOG}/FQC.out -e ${LOG}/FQC.err -J FASTQC -n 2 $COMMAND
    done
```


```bash
# # Create multiQC report of FastQC results
module load MultiQC;
multiqc $OUT -o $OUT;
module unload MultiQC;
```

The output MultiQC report is available at: https://github.com/PlantandFoodResearch/bioinf_Vitis_Nicotiana_RNAseq/blob/master/multiqc/rRNA_removed_trimmed_multiqc_report.html

Download the file, save as `.html`, and open in any web browser.


Minor adaptor contamination in R1 reads. This was successfully removed, which you can see when comparing the pre- and post- trimmomatic MultiQC reports.


```bash
# Render notebook to html and markdown
module load pfr-python3
jupyter nbconvert --to markdown /workspace/$USER/bioinf_Vitis_Nicotiana_RNAseq/01_RNAseq_Quality_Control_Vitis_Nicotiana-Manuscript.ipynb
jupyter nbconvert --to html /workspace/$USER/bioinf_Vitis_Nicotiana_RNAseq/01_RNAseq_Quality_Control_Vitis_Nicotiana-Manuscript.ipynb
module unload pfr-python3
```

    [NbConvertApp] Converting notebook /workspace/hradxj/bioinf_Vitis_Nicotiana_RNAseq/01_RNAseq_Quality_Control_Vitis_Nicotiana-Manuscript.ipynb to markdown
    [NbConvertApp] Writing 14295 bytes to /workspace/hradxj/bioinf_Vitis_Nicotiana_RNAseq/01_RNAseq_Quality_Control_Vitis_Nicotiana-Manuscript.md
    [NbConvertApp] Converting notebook /workspace/hradxj/bioinf_Vitis_Nicotiana_RNAseq/01_RNAseq_Quality_Control_Vitis_Nicotiana-Manuscript.ipynb to html
    [NbConvertApp] Writing 293356 bytes to /workspace/hradxj/bioinf_Vitis_Nicotiana_RNAseq/01_RNAseq_Quality_Control_Vitis_Nicotiana-Manuscript.html

