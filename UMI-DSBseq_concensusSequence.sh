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