# Generating Consensus Sequences
requires:
bcl2fastq, picard/2.8.3, jdk/8.111, bwa/0.7.13, fgbio/1.1.0

Design of the UMI-DSBseq adaptors was adapted from on the P7 tail of the xGen UDIUMI
adaptors from IDT (https://eu.idtdna.com/pages/products/next-generationsequencing/
adapters/xgen-udi-umi-adapters), containing an 8bp index for barcoding
from the IDT8_i7 index list, and a 9bp unique molecular identifier, for analysis with the
recommended pipeline for single molecule sequencing.

UMI-DSBseq libraries are sequenced using 150bp paired-end.
The run settings: 17 bp index 1 and 8 bp index 2, with 149-151 bp for each read 1 and read 2. Depending on the kit. Example here seen with 151 bp reads.

We recommend a minimum of between 2-5 million reads per sample: 

Demultiplexing is done using bcl2fasq, splitting the index 1 read into the UMI file and index file. Data processing pipeline was adapted from the IDT pipeline (https://www.youtube.com/watch?v=68sca_jsqg8&ab_channel=IntegratedDNATechnologies).
for building consensus sequences using FGBIO (https://github.com/fulcrumgenomics/fgbio). 

Fastqs are aligned to a reference of the target sequence using bwa-mem
Unmapped BAM files are generated using picard FastqToSam. 
The unmapped and mapped BAM files are merged using picard MergeBamAlignment (MAX_GAPS=-1, CLIP_ADAPTORS=true)
and annotated with UMIs using FGBIO.
Reads are grouped by UMI using FGBIOs GroupReadsByUmi (strategy=adjacency, edits=1, min-map-q=0,assign-tag=MI). Finally,Consensus sequences are generated using FGBIO CallMolecularConsensusReads with min-reads=2 and minimum input base quality set to 20. Final BAM files with consensusreads are converted to Fastqs using SamToFastq from picard, and joined using ea-utils fastq-join.

### 1. demultiplex raw run  
    bcl2fastq --input-dir /BaseCalls/ \
        --runfolder-dir run_folder/ \
        --output-dir /demultiplex_data/ \
        --sample-sheet /path_to_sample/Sample_Sheet.csv \
        --barcode-mismatches 0 \
        --use-bases-mask Y151,I8Y9,I8,Y151 \
        --create-fastq-for-index-reads \
        --ignore-missing-bcl \
        --no-lane-splitting \
        --mask-short-adapter-reads 0
   
### 2. Mapped Bam: Align demultiplexed Fastq paired end files to to reference made from Intact amplicon sequence for each target seperately with bwa mem: version used- bwa/0.7.13
    bwa mem -t 8 reference fastq_R1 fastq_R3 > mapped_bam

### 3. Umapped Bam: create unmapped bam from fastq files 
    java -Xmx4g -jar /picard/2.8.3/picard.jar FastqToSam F1=fastq_R1 F2=fastq_R3 O=unmapped_bam SM=sample_name
### 4. Merged Bam: Merge mapped and unmapped Bam files 
    java -Xmx4g -jar /picard/2.8.3/picard.jar MergeBamAlignment ALIGNED=mapped_bam UNMAPPED=unmapped_bam OUTPUT=output REFERENCE_SEQUENCE=reference MAX_GAPS=-1 SORT_ORDER=coordinate CREATE_INDEX=true CLIP_ADAPTERS=true
### 5. Annotate Bam: Annotate Bam with UMI's from UMI read (R3) File
    java -Xmx10g -jar FGBIO_JAR_FILE AnnotateBamWithUmis -i merged_bam -f umi_file -o merged_bam_UMI
### 6. Group reads by UMI: Agroup reads by UMI without filtering
    java -Xmx20g -jar FGBIO_JAR_FILE GroupReadsByUmi \
        --input=merged_bam_UMI \
        --output=merged_bam_UMI_group \
        --family-size-histogram=sample_name_hist \
        --strategy=adjacency --edits=1 --min-map-q=0 \
        --assign-tag=MI
### 7. Make Concensus sequences based on min base quality=20 and minumum reads=2
    java -Xmx20g -jar FGBIO_JAR_FILE CallMolecularConsensusReads \
        --input=merged_bam_UMI_group \
        --output=merged_bam_UMI_group_2readconsensus \
        --min-reads=2 \
        --rejects=merged_bam_UMI_group_2readconsensus_rejects \
        --min-input-base-quality=20 \
        --read-group-id=[ENTER sample_name]
### 8. Convert consensus Bam to Fastq
    java -Xmx4g -jar /apps/RH7U2/general/picard/2.8.3/picard.jar SamToFastq \
        INPUT=merged_bam_UMI_group_2readconsensus \
        FASTQ=merged_bam_UMI_group_2readconsensus_fastq1 \
        SECOND_END_FASTQ=merged_bam_UMI_group_2readconsensus_fastq2
### 9. join fastqs for further analysis
    /ea-utils/1.1.2/bin/fastq-join -o fastq_join/${sample}.%.fastq merged_bam_UMI_group_2readconsensus_fastq1 merged_bam_UMI_group_2readconsensus_fastq2
    
# UMI-DSBseq Consensus Sequence Analysis WorkFlow
Provided in the Jupyter notebook: UMI-DSBseq_analysis
## 1.User Input
    name of target (target_name) - for naming output
    directory containing joined-fastqs and Sample_data.xlxs (my_dir), 
    Sample_data.xlxs file: must be placed in fastq-join folder
        Required columns, rows for each sample in time-course:
             'library name': name of fastq-join file
             'Sample Name': name of sample
             'Guide sequence': sgRNA 20bp sequence (no PAM)
             'Amplicon sequence': expected sequence of wt intact molecule (primer--> restriction site:
             'forward_primer': sequene of target-specific amplification primer
             'restricted_end': expected 5 bp end of amplicon at end-repaired restriction site 
#### EXAMPLE
     target_name = 'Psy1'
     restricted_end = 'GTCCG'#enter 5bp end of amplicon at the restriction site
     imprecise_end='' #enter sequence for direct filtering
     my_dir = '/fastq_join/'
## 2. Create 100bp reference sequence
    reads in Sample_data.xlxs file
    extracts reference sequence from WT_ref (50bp to each side of DSB using Guide sequence and Amplicon sequence)
## 3. Convert fastq-join to table of id, read sequence, wt sequence 
    Parse fastq-joined files using SeqIO, make simple table for further processessing
## 4. Define state of consensus sequences as either Intact or unrepaired DSB
     align the last 12bp at each end of the 100bp WT reference window (ali1 and ali) to each read
         intact = both ali1 and ali2 >10/12 match to pass
         DSB = only ali1 >10/12 match to pass
         filter out = only ali2 >10/12 match to pass
## 5. Filter Intact reads
    uses the restriced_end expected at each target 
    grouping reads by sequence within the 100 bp window of interest
    adding a column 'count'- serves to reduce number of sequences needing alignment in the next steps
## 6. Characterizing consensus sequences
    Uses local alignment to the WT reference (given in the sample sheet) to characterize the consensus sequences
    Defines each Intact sequence as either WT or Indel (insertion/deletion)
        indels are further divided into insertions (Ins) and deletions (Del)
    Defines unrepaired DSBs
        Precise DSB: exactly at the expected position with a length of 50 bp
        guide-side DSB: position is shifted towards the sgRNA sequence
        PAM-side DSB: position is shifted towards the PAM sequence
        Extended DSB: longer than expected DSB of >50 bp with extension that does not match the reference sequence
    Calling indels and DSB names using Mut_col function according to their missing or added bases 
        Places those in columns 'Mut name' and 'DSB name' respectively
    Identifies any microhomology around the deletions, 
    groups and counts the sequences by all parameters.
        *changing window sizes requires adjusting the expected window length

#### Mut_col functions 
    1. Mut_col_reverse_guide_only: in the case that the primer is oriented correctly, but the sgRNA is in reverse
    2. Mut_col_reverse: in the case that both the primer and sgRNA are in reverse
    3. Mut_col_Forward: in the case that the primer and sgRNA are oriented correctly

## 7. Group reads and count 
     retains columns sequence, length, state, type, mut name, DSB
## 8. Identify microhomology associated deletions
    cut-off can be manually adjusted, default is 2 bp minimum microhomology
## 9. Output all_grouped_reads.csv of fully characterized grouped, counted reads for each sample 
    saved in directory /save_dir/target_name_all_grouped_reads_csv/'
## 10. Produce individual summary tables of counts for each of  for States, Types, Indels, and DSBs 
    states: Intact, Putative DSB, NHEJ
    types: Insertion, Deletion, Precise DSB, Extended DSB, Resected DSB
        Makes Tables of counts and percentages for all samples:
      [target_name]_States_df.csv: Indels, Intact (WT), DSB 
      [target_name]_Types_MH_df.csv: input to kinetics pipeline 
      [target_name]_indel.csv: error-prone repair footprints
      [target_name]_DSB.csv:specific DSBs by position and other characteristics 
