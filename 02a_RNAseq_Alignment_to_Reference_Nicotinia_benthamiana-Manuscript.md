
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
RAW=$PROJECT/000.raw
PROJECT="/workspace/$USER/$PROJECTNAME"
TEMP="$PROJECT/TEMP"
LOG=

TEMP="${PROJECT}/TEMP"
ANNOT="${PROJECT}/007.STAR/annotation/Niben101_annotation.gene_models.gff"
ZIPPEDGENOME="/workspace/ComparativeDataSources/Nicotiana/benthamiana/Genome/Niben.v1.0.1/assemblies/Niben.genome.v1.0.1.scaffolds.nrcontigs.fasta.gz"
INDEX="${PROJECT}/007.STAR/index"

PICARD="/workspace/cflcyd/software/picard/picard.jar"
LOG=$PROJECT/log
mkdir -p $LOG
```

#### Make appropriate directories and symlinks to files


```bash
mkdir ${PROJECT}/TEMP
mkdir ${PROJECT}/007.STAR
mkdir ${PROJECT}/007.STAR/logs
mkdir ${PROJECT}/007.STAR/annotation
mkdir ${PROJECT}/007.STAR/genome
mkdir ${PROJECT}/007.STAR/index
mkdir ${PROJECT}/007.STAR/index/logs
mkdir ${PROJECT}/007.STAR/Single_Pass_Results
mkdir ${PROJECT}/007.STAR/Two_Pass_Results
mkdir ${PROJECT}/008.MBA
mkdir ${PROJECT}/008.MBA/logs
mkdir ${PROJECT}/008.MBA/Single_Pass_Results
mkdir ${PROJECT}/008.MBA/Two_Pass_Results

ln -s  /workspace/ComparativeDataSources/Nicotiana/benthamiana/Genome/Niben.v1.0.1/annotation/Niben101/Niben101_annotation.gene_models.gff ${ANNOT}
```


```bash
# Obtain Grapevine leafroll-associated virus 3 and annotations
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/851/885/GCF_000851885.1_ViralProj14906/GCF_000851885.1_ViralProj14906_genomic.fna.gz -O $RAW/GCF_000851885.1_ViralProj14906_genomic.fna.gz;
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/851/885/GCF_000851885.1_ViralProj14906/GCF_000851885.1_ViralProj14906_genomic.gff.gz -O $RAW/GCF_000851885.1_ViralProj14906_genomic.gff.gz;
gunzip $RAW/*.gz;

```

## <u>Step II, Part 2: Index the Genome Using STAR</u>
* The inputs to this step are the genome as a multi FASTA file and the annotations as a GFF or GTF file
* The outputs include the genome index files used in the alignment steps


```bash
# The genome is a gzipped file and STAR requires unzipped files. Unzip

zcat ${ZIPPEDGENOME} > ${PROJECT}/007.STAR/genome/Niben.genome.v1.0.1.scaffolds.nrcontigs.fasta
```


```bash
GENOME=${PROJECT}/007.STAR/genome/Niben.genome.v1.0.1.scaffolds.nrcontigs.fasta;
echo $GENOME
echo $ANNOT
```

    /workspace/hradxj/karmun_awesome_experiment/007.STAR/genome/Niben.genome.v1.0.1.scaffolds.nrcontigs.fasta
    /workspace/hradxj/karmun_awesome_experiment/007.STAR/annotation/Niben101_annotation.gene_models.gff



```bash
# We need to examine the annotation file to set 3 parameters: the chromosome prefix and the text used
# to identify the parent-child relationship between exons and genes.
head -4 $ANNOT
# So the chromosome prefix is Niben101Ctg

```

    Niben101Ctg00001	.	contig	1	500	.	.	.	ID=Niben101Ctg00001
    Niben101Ctg00002	.	contig	1	500	.	.	.	ID=Niben101Ctg00002
    Niben101Ctg00003	.	contig	1	500	.	.	.	ID=Niben101Ctg00003
    Niben101Ctg00004	.	contig	1	500	.	.	.	ID=Niben101Ctg00004



```bash
# Let's check that the genome contains the right records
grep -A9 "Niben101Ctg00054" < $GENOME
# YES, this FASTA does contain a record exactly matching the annotation.
```

    >Niben101Ctg00054 cov=9.0
    AAAACCCAATTATTCTCTGAATCAATTCTCCTTCTTCCCTCTACCTCTCCTTTTCACTAA
    AAACCTAAACATTTTTCAATATCTCTCTATTAACCCATTTATACATAAATCTACAACGCA
    GTTCAGTTTGTTAAAGTTATTGCACTGTCTAAAAAAAAGAGCATCAATGGCTGAGCCAAC
    AACACCAAACTCAGAATCAGAATCTAGTAGTTATAACTCTTGTTCTCTTTCTTCTACAAT
    TTCATCTTCTTCTGTTCTTATAAAGAATATCAACTCGAAAAACCGACTCAAGAGATGCCG
    TGAAGTAGCAGAAGAAAATGATGGTCAGACAAATATTAGTAAGTGCAGTAATAGCAGTAA
    GAAAAATGGTGGCAACAAAACTAGCGACGGTTCAAAACACCCATCGTACGTTGGTGTACG
    AAAGAGGGCATGGGGAAAATGGGTGTCCGAAATTCGTGAACCGAAGAAGAAATCAAGAAT
    CTGGTTAGGTACTTTCGCCAC


From the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf):
>2.2.3 Annotations in GFF format.
In addition to the aforementioned options, for GFF3 formatted annotations you need to use
--sjdbGTFtagExonParentTranscript Parent. In general, for --sjdbGTFfile files STAR only
processes lines which have --sjdbGTFfeatureExon (=exon by default) in the 3rd field (column).
The exons are assigned to the transcripts using parent-child relationship defined by the
--sjdbGTFtagExonParentTranscript (=transcript id by default) GTF/GFF attribute.


```bash
# First let's look for the first line with "gene" in it
cat $ANNOT | grep gene | head -1
# We can see that the gene has an "ID=". This text identifies the gene. 
# The subcomponents of the gene will refer to this as a "parent"
```

    Niben101Ctg00054	maker	gene	167	487	.	+	.	ID=Niben101Ctg00054g00001;Alias=snap_masked-Niben101Ctg00054-processed-gene-0.0
    grep: write error
    cat: write error: Broken pipe



```bash
# Let's examine all lines that contain this "ID=" tag.
grep "ID=Niben101Ctg00054g00001" < $ANNOT
# So the tag that identifies what transcript an exon comes from is "ID"
# and the tag that identifies what gene an mRNA or CDS comes from is also "ID" 

```

    Niben101Ctg00054	maker	gene	167	487	.	+	.	ID=Niben101Ctg00054g00001;Alias=snap_masked-Niben101Ctg00054-processed-gene-0.0
    Niben101Ctg00054	maker	mRNA	167	487	.	+	.	ID=Niben101Ctg00054g00001.1;Parent=Niben101Ctg00054g00001;Alias=snap_masked-Niben101Ctg00054-processed-gene-0.0-mRNA-1;Note="Ethylene-responsive transcription factor 7";Ontology_term=GO:0003677,GO:0006355,GO:0003700
    Niben101Ctg00054	maker	exon	167	487	.	+	.	ID=Niben101Ctg00054g00001.1:exon:001;Parent=Niben101Ctg00054g00001.1;Alias=snap_masked-Niben101Ctg00054-processed-gene-0.0-mRNA-1:exon:4812
    Niben101Ctg00054	maker	CDS	167	487	.	+	0	ID=Niben101Ctg00054g00001.1:cds:001;Parent=Niben101Ctg00054g00001.1;Alias=snap_masked-Niben101Ctg00054-processed-gene-0.0-mRNA-1:cds


Convert gff3 to gtf for use with STAR


```bash
/software/bioinformatics/cufflinks-2.2.1/gffread $ANNOT -T -o ${PROJECT}/007.STAR/annotation/Niben101_annotation.gene_models.gtf

```


```bash
# Convert the gtf to use gene_ID as a tag with a find-and-replace.
sed 's/geneID/gene_ID/g' < $ANNOTGTF > ${PROJECT}/007.STAR/annotation//Niben101_annotation.gene_models_fixedg05002.gtf;

```


```bash
# Reset the $ANNOTGTF variable
ANNOTGTF=${PROJECT}/007.STAR/annotation/Niben101_annotation.gene_models_fixedg05002.gtf
```


```bash
mkdir ${PROJECT}/007.STAR;
mkdir ${PROJECT}/007.STAR/logs;
mkdir ${PROJECT}/007.STAR/annotation;
mkdir ${PROJECT}/007.STAR/genome;
mkdir ${PROJECT}/007.STAR/index;
mkdir ${PROJECT}/007.STAR/index/logs;
mkdir ${PROJECT}/007.STAR/Single_Pass_Results;
mkdir ${PROJECT}/007.STAR/Two_Pass_Results;

GENOME=${PROJECT}/007.STAR/genome/Niben.genome.v1.0.1.scaffolds.nrcontigs.fasta;


COMMAND="module load STAR; \
STAR \
--runMode genomeGenerate \
--limitGenomeGenerateRAM 240000000000 \
--runThreadN 32 \
--genomeFastaFiles $GENOME \
--genomeDir ${PROJECT}/007.STAR/index \
--sjdbGTFfile ${ANNOTGTF} \
--sjdbGTFtagExonParentGene gene_id \
--sjdbGTFtagExonParentTranscript transcript_id; \
module unload STAR"
echo $COMMAND;
bsub \
-J STAR_Dan \
-o ${PROJECT}/007.STAR/index/logs/%J_STAR_index.out \
-e ${PROJECT}/007.STAR/index/logs/%J_STAR_index.err \
-n 32 \
$COMMAND
```

The inputs for the next step are the final outputs of the QC notebook: trimmed, rRNA removed FASTQ files.


```bash
TRIMMED=${PROJECT}/004.trimmomatic
OUT_STAR="${PROJECT}/007.STAR/Single_Pass_Results"
INDEX="${PROJECT}/007.STAR/index"
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
echo $COMMAND
            
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
SJLIST=$(ls ${PROJECT}/007.STAR/Single_Pass_Results/*SJ.out.tab)
```


```bash
TRIMMED=${PROJECT}/004.trimmomatic
INDEX=${PROJECT}/007.STAR/index
OUT_STAR=${PROJECT}/007.STAR/Two_Pass_Results
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

For this analysis, there are three samples:

RACP005_11_S11_L002                 Benth-Healthy-1
RACP005_12_S12_L002                 Benth-Healthy-2
RACP005_13_S13_L002                 Benth-Infected

Read counts have been generated by STAR. 
It is neccessary to create a single tab-delimited text file for import into R.

In this case, STAR produces a table of read counts that has the following characteristics:

1) The first four lines of the table contain information about the number of unmapped and multimapped reads.An example is shown below:

```
head RACP005_13_S13_L002_MERGED_trimmomatic.bamReadsPerGene.out.tab
N_unmapped      24827661        24827661        24827661
N_multimapping  6684546 6684546 6684546
N_noFeature     7429471 7663155 34125351
N_ambiguous     717277  256966  6802
```
These four lines are relevant but not used in the DE analysis, so should be removed.

2) The columns are ordered as follows

Gene name | Unstranded read counts | Sense strand read counts | Antisense strand read counts
--- | --- | --- | ---

We are interested in columns 1 and 3 (and column 4 if we choose to investigate antisense transcripts).

Therefore we need to create a file that contains the gene names in column 1, and the Sense strand read counts from 3 samples in columns 2,3,4.



```bash
# Find location of read count files
ls ${PROJECT}/007.STAR/Two_Pass_Results | grep ReadsPerGene
```

    RACP005_11_S11_L002_MERGED_trimmomatic.bamReadsPerGene.out.tab
    RACP005_12_S12_L002_MERGED_trimmomatic.bamReadsPerGene.out.tab
    RACP005_13_S13_L002_MERGED_trimmomatic.bamReadsPerGene.out.tab



```bash
# First create a new file with one column of all gene names.

mkdir -p ${PROJECT}/010.edgeR_Nb;

cat ${PROJECT}/007.STAR/Two_Pass_Results/RACP005_11_S11_L002_MERGED_trimmomatic.bamReadsPerGene.out.tab |awk '{print $1}'\
> ${PROJECT}/010.edgeR_Nb/GRLaV3_Nb_EdgeR-genenames.tab;

```


```bash
# Now add the appropriate column
# Create a list of read count file names
READCOUNTFILELIST=$(ls ${PROJECT}/007.STAR/Two_Pass_Results | grep ReadsPerGene)


for READCOUNTFILE in $READCOUNTFILELIST
do
awk '{print $3}' < ${PROJECT}/007.STAR/Two_Pass_Results/${READCOUNTFILE} > ${PROJECT}/010.edgeR_Nb/${READCOUNTFILE}.col3;
done

paste ${PROJECT}/010.edgeR_Nb/GRLaV3_Nb_EdgeR-genenames.tab \
${PROJECT}/010.edgeR_Nb/RACP005_11_S11_L002_MERGED_trimmomatic.bamReadsPerGene.out.tab.col3 \
${PROJECT}/010.edgeR_Nb/RACP005_12_S12_L002_MERGED_trimmomatic.bamReadsPerGene.out.tab.col3 \
${PROJECT}/010.edgeR_Nb/RACP005_13_S13_L002_MERGED_trimmomatic.bamReadsPerGene.out.tab.col3 \
> ${PROJECT}/010.edgeR_Nb/GRLaV3_Nb_EdgeR_with_unmapped.tab

```


```bash
head ${PROJECT}/010.edgeR_Nb/GRLaV3_Nb_EdgeR_with_unmapped.tab

```

    N_unmapped	11575247	1193884	24827661
    N_multimapping	6457575	575648	6684546
    N_noFeature	6314912	701699	7663155
    N_ambiguous	219454	21609	256966
    Niben101Ctg00054g00001	55	23	79
    Niben101Ctg00074g00004	0	0	0
    Niben101Ctg00116g00002	65	17	133
    Niben101Ctg00129g00001	0	0	0
    Niben101Ctg00141g00002	1	0	1
    Niben101Ctg00174g00001	0	0	0



```bash
tail -59814  ${PROJECT}/010.edgeR_Nb/GRLaV3_Nb_EdgeR_with_unmapped.tab >  ${PROJECT}/010.edgeR_Nb/GRLaV3_Nb_EdgeR.tab

```


```bash
# Add header line

sed -i '1s/^/Gene\tBenth-Healthy-1\tBenth-Healthy-2\tBenth-Infected\n/' ${PROJECT}/010.edgeR_Nb/GRLaV3_Nb_EdgeR.tab;

```


```bash
# Check final formatting of file
head ${PROJECT}/010.edgeR_Nb/GRLaV3_Nb_EdgeR.tab
# We are now ready to analyse this in EdgeR
```

    Gene	Benth-Healthy-1	Benth-Healthy-2	Benth-Infected
    Niben101Ctg00054g00001	55	23	79
    Niben101Ctg00074g00004	0	0	0
    Niben101Ctg00116g00002	65	17	133
    Niben101Ctg00129g00001	0	0	0
    Niben101Ctg00141g00002	1	0	1
    Niben101Ctg00174g00001	0	0	0
    Niben101Ctg00192g00001	0	0	0
    Niben101Ctg00219g00002	3	1	4
    Niben101Ctg00228g00001	43	5	53


We have derived a set of raw read counts. These are the inputs to the downstream RMarkdown notebooks.


```bash
# Render notebook to html and markdown
module load pfr-python3
jupyter nbconvert --to markdown /workspace/$USER/bioinf_Vitis_Nicotiana_RNAseq/02a_RNAseq_Alignment_to_Reference_Nicotinia_benthamiana-Manuscript.ipynb
jupyter nbconvert --to html /workspace/$USER/bioinf_Vitis_Nicotiana_RNAseq/02a_RNAseq_Alignment_to_Reference_Nicotinia_benthamiana-Manuscript.ipynb
module unload pfr-python3
```

    [NbConvertApp] Converting notebook /workspace/hradxj/bioinf_Vitis_Nicotiana_RNAseq/02a_RNAseq_Alignment_to_Reference_Nicotinia_benthamiana-Manuscript.ipynb to markdown
    [NbConvertApp] Writing 14187 bytes to /workspace/hradxj/bioinf_Vitis_Nicotiana_RNAseq/02a_RNAseq_Alignment_to_Reference_Nicotinia_benthamiana-Manuscript.md
    [NbConvertApp] Converting notebook /workspace/hradxj/bioinf_Vitis_Nicotiana_RNAseq/02a_RNAseq_Alignment_to_Reference_Nicotinia_benthamiana-Manuscript.ipynb to html
    [NbConvertApp] Writing 288296 bytes to /workspace/hradxj/bioinf_Vitis_Nicotiana_RNAseq/02a_RNAseq_Alignment_to_Reference_Nicotinia_benthamiana-Manuscript.html

