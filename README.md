# UMI-DSBseq
"""1. align to reference sequence
2. make unmapped bam
3. merge mapped and umapped bam
4. annotate reads with UMIs
5. merge UMI annotated BAMs for this run and the previous run
6. group reads by UMI
7. generate consensus reads 
8. make Fastqs 
9. join Fastqs """

module load picard/2.8.3
module load jdk/8.111 
module load samtools
module load bwa/0.7.13
module load fgbio/1.1.0

# 1. demultiplex raw run  
## bcl2fastq: 
### --barcode-mismatches 0
### --use-bases-mask Y151,I8Y9,I8,Y151
### --ignore-missing-bcl
### --no-lane-splitting
### --mask-short-adapter-reads 0

bcl2fastq --input-dir /BaseCalls/ \
    --runfolder-dir run_folder/ \
    --output-dir /home/labs/alevy/Collaboration/Novaseq_May16_2021/demultiplex_data/ \
    --sample-sheet /path_to_sample/Sample_Sheet.csv \
    --barcode-mismatches 0 \
    --use-bases-mask Y151,I8Y9,I8,Y151 \
    --create-fastq-for-index-reads \
    --ignore-missing-bcl \
    --no-lane-splitting \
    --mask-short-adapter-reads 0

# 2. Mapped Bam: Align demultiplexed Fastq paired end files to to reference made from Intact amplicon sequence for each target seperately 
## bwa mem: bwa/0.7.13
bwa mem -t 8 $reference $fastq_R1 $fastq_R3 > $mapped_bam

# 3. Umapped Bam: create unmapped bam from fastq files 
java -Xmx4g -jar /apps/RH7U2/general/picard/2.8.3/picard.jar FastqToSam F1=$fastq_R1 F2=$fastq_R3 \
    O=$unmapped_bam \
    SM=${sample}
# 4. Merged Bam: Merge mapped and unmapped Bam files 
java -Xmx4g -jar /apps/RH7U2/general/picard/2.8.3/picard.jar MergeBamAlignment ALIGNED=$mapped_bam UNMAPPED=$unmapped_bam OUTPUT=$output REFERENCE_SEQUENCE=$reference MAX_GAPS=-1 SORT_ORDER=coordinate CREATE_INDEX=true CLIP_ADAPTERS=true
# 5. Annotate Bam: Annotate Bam with UMI's from UMI read (R3) File
java -Xmx10g -jar $FGBIO_JAR_FILE AnnotateBamWithUmis -i $merged_bam -f $umi_file -o $merged_bam_UMI
# 6. Group reads by UMI: Agroup reads by UMI without filtering
java -Xmx20g -jar $FGBIO_JAR_FILE GroupReadsByUmi \
    --input=$merged_bam_UMI \
    --output=$merged_bam_UMI_group \
    --family-size-histogram=${sample} \
    --strategy=adjacency --edits=1 --min-map-q=0 \
    --assign-tag=MI
# 7. Make Concensus sequences based on min base quality=20 and minumum reads=2
java -Xmx20g -jar $FGBIO_JAR_FILE CallMolecularConsensusReads \
    --input=$merged_bam_UMI_group \
    --output=$merged_bam_UMI_group_2readconsensus \
    --min-reads=2 \
    --rejects=$merged_bam_UMI_group_2readconsensus_rejects \
    --min-input-base-quality=20 \
    --read-group-id=${sample}
# 8. Convert consensus Bam to Fastq
java -Xmx4g -jar /apps/RH7U2/general/picard/2.8.3/picard.jar SamToFastq \
    INPUT=$merged_bam_UMI_group_2readconsensus \
    FASTQ=$merged_bam_UMI_group_2readconsensus_fastq1 \
    SECOND_END_FASTQ=$merged_bam_UMI_group_2readconsensus_fastq2
# 9. join fastqs for further analysis
/apps/RH7U2/gnu/ea-utils/1.1.2/bin/fastq-join -o fastq_join/${sample}.%.fastq $merged_bam_UMI_group_2readconsensus_fastq1 $merged_bam_UMI_group_2readconsensus_fastq2

## UMI-DSBseq Analysis WorkFlow
### 1.User Input
####  target information:  
###### name of target (target_name) - for naming output
###### directory containing joined-fastqs and Sample_data.xlxs (my_dir), 
###### 5bp end of any known PCR contamination (primer)

#### Sample_data.xlxs file: must be placed in fastq-join folder
##### Required columns, rows for each sample in time-course:
###### 'library name': name of fastq-join file
###### 'Sample Name': name of sample
###### 'Guide sequence': sgRNA 20bp sequence (no PAM)
###### 'Amplicon sequence': expected sequence of wt intact molecule (primer--> restriction site:
###### 'forward_primer': sequene of target-specific amplification primer
###### 'restricted_end': expected 5 bp  end of amplicon at end-repaired restr. site 
### 2. Create 100bp reference sequence
##### reads in Sample_data.xlxs file and extracts reference sequence from WT_ref (50bp to each side of DSB using Guide sequence and Amplicon sequence)
### 3. Convert fastq-join to table of id, read sequence, wt sequence 
###### Parse fastq-joined files using SeqIO, make simple table for further processessing
### 4. Differentiate state of read: Intact vs. DSB
##### align the last 12bp at each end of the 100bp WT reference window (ali1 and ali) to each read
###### intact = both ali1 and ali2 >10/12 match to pass
###### DSB = only ali1 >10/12 match to pass
###### filter out = only ali2 >10/12 match to pass
### 5. Filter out Intact reads resulting from contamination of PCR products - ensure that the ends correspond to the restriction site
##### uses the restriced_end expected at each target - can filter using user input of contaminated end (known 5bp sequence from amplicon) entered in step 1
##### grouping reads by sequence within the 100 bp window of interest
###### adding a column 'count'- serves to reduce number of sequences needing alignment in the next steps
### 6. Determining read Type: Calling Indels and Putative DSBs vs WT reads, determines type of indels, and name of indel
#### Type: WT, Ins, Del, Precise DSB,  Extended DSB, Resected DSB
##### Uses Mut_col function to name the indels and DSBs according to their missing or added bases, places those in columns 'Mut name' and 'DSB name' respectively, identifies any microhomology around the deletions, groups and counts the reads
##### changing window sizes requires adjusting the expected window length has to be changed 
### 7. Group reads and count 
##### retains columns sequence, length, sate, type, mut name, DSB
### 8. Identify read groups with microhomology associated deletions of 2bps or more
#### cut-off can be manually adjusted
### 9. Output all_grouped_reads.csv of fully characterized grouped, counted reads for each sample 
#### saved in directory /save_dir/target_name_all_grouped_reads_csv/'
### 10. Produce individual summary tables of counts for each of  for States, Types, Indels, and DSBs 
#### states: Intact, Putative DSB, NHEJ
#### types: Insertion, Deletion, Precise DSB, Extended DSB, Resected DSB
#### Makes Tables of counts and percentages for all samples:
###### [target_name]_States_df.csv: Indels, Intact (WT), DSB 
###### [target_name]_Types_MH_df.csv: input to kinetics pipeline 
###### [target_name]_indel.csv: error-prone repair footprints
###### [target_name]_DSB.csv:specific DSBs by position and other characteristics 
