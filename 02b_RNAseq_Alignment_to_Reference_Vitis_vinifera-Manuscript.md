
# <u>RNA-Seq Analysis Phase IIa: Alignment to Reference Genome</u>
## This Notebook illustrates how to align paired-end RNA-Seq reads <br>that have already been processed through the QC pipeline.
#### Last Revision: July 2017
#### Author: Charles David
#### This analysis done by Dan Jones and Karmun Chooi

## <u>Step I: Establish Data Management Structure on PowerPlant (Continuing from the QC part)</u>

Already completed in the QC notebook

## <u>Step II, Part 1: Get the Genome and Annotation Files to be Used in the Alignment Process</u>
##### Note that best results are obtained if the reference is good quality and closely related, with annotations

### Define Project Variables:
* Note that we are using the latest version of STAR: 2.5.2b
* We are also using the latest version of Picard Tools: 2.9.4


```bash
# Define the user as a variable
USER="hradxj"
PROJECTNAME="Vitis_Nicotiana_experiment"
# Define the project directory and temp subdirectory as a variable

PROJECT="/workspace/$USER/$PROJECTNAME"
TEMP="$PROJECT/TEMP"

mkdir -p ${PROJECT}/009.STAR
mkdir -p ${PROJECT}/009.STAR/logs
mkdir -p ${PROJECT}/009.STAR/annotation
mkdir -p ${PROJECT}/009.STAR/genome
mkdir -p ${PROJECT}/009.STAR/index
mkdir -p ${PROJECT}/009.STAR/index/logs
mkdir -p ${PROJECT}/009.STAR/Single_Pass_Results
mkdir -p ${PROJECT}/009.STAR/Two_Pass_Results

ln -s /workspace/ComparativeDataSources/Vitis/vinifera/Genoscope_12X/Genes/Vitis_vinifera_annotation.gff3 ${PROJECT}/009.STAR/annotation
ln -s /workspace/ComparativeDataSources/Vitis/vinifera/Genoscope_12X/Genome/reference.fasta ${PROJECT}/009.STAR/genome

ANNOT=${PROJECT}/009.STAR/annotation/Vitis_vinifera_annotation.gff3
ANNOTGTF=${PROJECT}/009.STAR/annotation/Vitis_vinifera_annotation.gtf
GENOME=${PROJECT}/009.STAR/genome/reference.fasta
```

## <u>Step II, Part 2: Index the Genome Using STAR</u>
* The inputs to this step are the genome as a multi FASTA file and the annotations as a GFF or GTF file
* The outputs include the genome index files used in the alignment steps

Convert gff3 to gtf for use with STAR


```bash
/software/bioinformatics/cufflinks-2.2.1/gffread $ANNOT -g $GENOME -T -o $ANNOTGTF
```


```bash
COMMAND="module load STAR; \
STAR \
--runMode genomeGenerate \
--limitGenomeGenerateRAM 240000000000 \
--runThreadN 32 \
--genomeFastaFiles $GENOME \
--sjdbGTFfile ${ANNOTGTF} \
--sjdbGTFtagExonParentTranscript transcript_id \
--sjdbGTFtagExonParentGene gene_id \
--genomeDir ${PROJECT}/009.STAR/index; \
module unload STAR"
echo $COMMAND;
bsub \
-J STAR_Dan \
-o ${PROJECT}/009.STAR/index/logs/%J_STAR_index.out \
-e ${PROJECT}/009.STAR/index/logs/%J_STAR_index.err \
-n 32 \
-q lowpriority \
$COMMAND
```

This uses genome version (1.68.5) and model version (15...34)


```bash
#Find input trimmed rRNA removed reads
ls ${PROJECT}/004.trimmomatic
```

    logs
    RACP005_11_S11_L002_MERGED_trimmomatic_R1.fastq
    RACP005_11_S11_L002_MERGED_trimmomatic_R2.fastq
    RACP005_12_S12_L002_MERGED_trimmomatic_R1.fastq
    RACP005_12_S12_L002_MERGED_trimmomatic_R2.fastq
    RACP005_13_S13_L002_MERGED_trimmomatic_R1.fastq
    RACP005_13_S13_L002_MERGED_trimmomatic_R2.fastq
    RACP005_1_S8_L002_MERGED_trimmomatic_R1.fastq
    RACP005_1_S8_L002_MERGED_trimmomatic_R2.fastq
    RACP005_5_S9_L002_MERGED_trimmomatic_R1.fastq
    RACP005_5_S9_L002_MERGED_trimmomatic_R2.fastq
    RACP005_8_S10_L002_MERGED_trimmomatic_R1.fastq
    RACP005_8_S10_L002_MERGED_trimmomatic_R2.fastq
    unpaired



```bash
TRIMMED=${PROJECT}/004.trimmomatic
INDEX=${PROJECT}/009.STAR/index
OUT_STAR=${PROJECT}/009.STAR/Single_Pass_Results
mkdir -p ${OUT_STAR}
PREFIXLIST=`basename -a ${TRIMMED}/*.fastq | sed 's/_R[1,2].fastq//g'|sort -u `
echo $PREFIXLIST

for PREFIX in ${PREFIXLIST}
do
echo $PREFIX
R1=${TRIMMED}/${PREFIX}_R1.fastq
R2=${TRIMMED}/${PREFIX}_R2.fastq


COMMAND="module load STAR; \
            STAR \
            --runThreadN 8 \
            --genomeDir ${INDEX} \
            --readFilesIn ${R1} ${R2} \
            --sjdbGTFfile ${ANNOTGTF} \
            --sjdbGTFtagExonParentGene gene_id \
            --sjdbGTFtagExonParentTranscript transcript_id \
            --outFileNamePrefix ${OUT_STAR}/${PREFIX}.bam \
            --quantMode TranscriptomeSAM GeneCounts \
            --outStd BAM_SortedByCoordinate"
bsub \
-J STAR \
-o ${LOG}/%J_STAR_index.out \
-e ${LOG}/%J_STAR_index.err \
-n 8 \
-q lowpriority \
$COMMAND
done
```


```bash
# # Run STAR in two-pass mode
# # Create list of files containing splice junctions
SJLIST=$(ls ${PROJECT}/009.STAR/Single_Pass_Results/*SJ.out.tab)
```


```bash
TRIMMED=${PROJECT}/004.trimmomatic
INDEX=${PROJECT}/009.STAR/index
OUT_STAR=${PROJECT}/009.STAR/Two_Pass_Results
mkdir -p ${OUT_STAR}
PREFIXLIST=`basename -a ${TRIMMED}/*.fastq | sed 's/_R[1,2].fastq//g'|sort -u `
echo $PREFIXLIST

for PREFIX in ${PREFIXLIST}
do
echo $PREFIX
R1=${TRIMMED}/${PREFIX}_R1.fastq
R2=${TRIMMED}/${PREFIX}_R2.fastq


COMMAND="module load STAR; \
            STAR \
            --runThreadN 8 \
            --genomeDir ${INDEX} \
            --readFilesIn ${R1} ${R2} \
            --sjdbGTFfile ${ANNOTGTF} \
            --sjdbGTFtagExonParentGene gene_id \
            --sjdbGTFtagExonParentTranscript transcript_id \
            --outFileNamePrefix ${OUT_STAR}/${PREFIX}.bam \
            --sjdbFileChrStartEnd ${SJLIST} \
            --quantMode GeneCounts \
            --outStd BAM_SortedByCoordinate"
bsub \
-J STAR \
-o ${LOG}/%J_STAR_index.out \
-e ${LOG}/%J_STAR_index.err \
-n 32 \
-q lowpriority \
$COMMAND
done
```


```bash
# Find location of read count files
ls ${PROJECT}/009.STAR/Two_Pass_Results | grep ReadsPerGene
```

    RACP005_1_S8_L002_MERGED_trimmomatic.bamReadsPerGene.out.tab
    RACP005_5_S9_L002_MERGED_trimmomatic.bamReadsPerGene.out.tab
    RACP005_8_S10_L002_MERGED_trimmomatic.bamReadsPerGene.out.tab



```bash
# First create a new file with one column of all gene names. Do this for two files, one for sense
# and one for antisense
mkdir -p ${PROJECT}/011.edgeR_Vv;

cat ${PROJECT}/009.STAR/Two_Pass_Results/RACP005_1_S8_L002_MERGED_trimmomatic.bamReadsPerGene.out.tab |awk '{print $1}'\
> ${PROJECT}/011.edgeR_Vv/GRLaV3_Vv_EdgeR-genenames.tab;
```


```bash
# Now add the appropriate column
# Create a list of read count file names
READCOUNTFILELIST=$(ls ${PROJECT}/009.STAR/Two_Pass_Results | grep ReadsPerGene)


for READCOUNTFILE in $READCOUNTFILELIST
do
awk '{print $3}' < ${PROJECT}/009.STAR/Two_Pass_Results/${READCOUNTFILE} > ${PROJECT}/011.edgeR_Vv/${READCOUNTFILE}.col3;
done

paste ${PROJECT}/011.edgeR_Vv/GRLaV3_Vv_EdgeR-genenames.tab \
${PROJECT}/011.edgeR_Vv/RACP005_1_S8_L002_MERGED_trimmomatic.bamReadsPerGene.out.tab.col3 \
${PROJECT}/011.edgeR_Vv/RACP005_5_S9_L002_MERGED_trimmomatic.bamReadsPerGene.out.tab.col3 \
${PROJECT}/011.edgeR_Vv/RACP005_8_S10_L002_MERGED_trimmomatic.bamReadsPerGene.out.tab.col3 \
> ${PROJECT}/011.edgeR_Vv/GRLaV3_Vv_EdgeR_with_unmapped.tab

```

For this analysis, there are three samples:

RACP005_1_S8_L0022                 Grape-Healthy
RACP005_5_S9_L002                  Grape-Infected-1
RACP005_8_S10_L002                 Grape-Infected-2

Read counts have been generated by STAR. 
It is neccessary to create a single tab-delimited text file for import into R.

In this case, STAR produces a table of read counts that has the following characteristics:

1) The first four lines of the table contain information about the number of unmapped and multimapped reads.An example is shown below:

These four lines are relevant but not used in the DE analysis, so should be removed.

2) The columns are ordered as follows

Gene name | Unstranded read counts | Sense strand read counts | Antisense strand read counts
--- | --- | --- | ---

We are interested in columns 1 and 3 (and column 4 if we choose to investigate antisense transcripts).

Therefore we need to create a file that contains the gene names in column 1, and the Sense strand read counts from 3 samples in columns 2,3,4.



```bash
head  ${PROJECT}/011.edgeR_Vv/GRLaV3_Vv_EdgeR_with_unmapped.tab
```

    N_unmapped	5053374	11312084	3028748
    N_multimapping	1917215	2934146	514757
    N_noFeature	4601554	10738708	1439012
    N_ambiguous	78822	288253	24490
    GSVIVG01012261001	14	37	1
    GSVIVG01012259001	6	18	1
    GSVIVG01012257001	311	1151	151
    GSVIVG01012255001	755	1819	234
    GSVIVG01012253001	1	0	0
    GSVIVG01012250001	0	0	0



```bash
wc -l  ${PROJECT}/011.edgeR_Vv/GRLaV3_Vv_EdgeR_with_unmapped.tab
```

    19763 /workspace/hradxj/karmun_awesome_experiment/011.edgeR_Vv/GRLaV3_Vv_EdgeR_with_unmapped.tab



```bash
tail -19759 ${PROJECT}/011.edgeR_Vv/GRLaV3_Vv_EdgeR_with_unmapped.tab > ${PROJECT}/011.edgeR_Vv/GRLaV3_Vv_EdgeR.tab
```


```bash
sed -i '1s/^/Gene\tGrape-Healthy\tGrape-Infected-1\tGrape-Infected-2\n/' ${PROJECT}/011.edgeR_Vv/GRLaV3_Vv_EdgeR.tab;
```


```bash
head ${PROJECT}/011.edgeR_Vv/GRLaV3_Vv_EdgeR.tab;
```

    Gene	Grape-Healthy	Grape-Infected-1	Grape-Infected-2
    GSVIVG01012261001	14	37	1
    GSVIVG01012259001	6	18	1
    GSVIVG01012257001	311	1151	151
    GSVIVG01012255001	755	1819	234
    GSVIVG01012253001	1	0	0
    GSVIVG01012250001	0	0	0
    GSVIVG01012249001	0	0	0
    GSVIVG01012247001	32	32	6
    GSVIVG01012246001	0	0	0


We have derived a set of raw read counts. These are the inputs to the downstream RMarkdown notebooks.


```bash
# Render notebook to html and markdown
module load pfr-python3
jupyter nbconvert --to markdown /workspace/$USER/bioinf_Vitis_Nicotiana_RNAseq/02b_RNAseq_Alignment_to_Reference_Vitis_vinifera-Manuscript.ipynb
jupyter nbconvert --to html /workspace/$USER/bioinf_Vitis_Nicotiana_RNAseq/02b_RNAseq_Alignment_to_Reference_Vitis_vinifera-Manuscript.ipynb
module unload pfr-python3
```

    [NbConvertApp] Converting notebook /workspace/hradxj/bioinf_Vitis_Nicotiana_RNAseq/02b_RNAseq_Alignment_to_Reference_Vitis_vinifera-Manuscript.ipynb to markdown
    [NbConvertApp] Writing 9332 bytes to /workspace/hradxj/bioinf_Vitis_Nicotiana_RNAseq/02b_RNAseq_Alignment_to_Reference_Vitis_vinifera-Manuscript.md
    [NbConvertApp] Converting notebook /workspace/hradxj/bioinf_Vitis_Nicotiana_RNAseq/02b_RNAseq_Alignment_to_Reference_Vitis_vinifera-Manuscript.ipynb to html
    [NbConvertApp] Writing 276881 bytes to /workspace/hradxj/bioinf_Vitis_Nicotiana_RNAseq/02b_RNAseq_Alignment_to_Reference_Vitis_vinifera-Manuscript.html

