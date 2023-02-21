# UMI-DSBseq Consensus Sequence Analysis WorkFlow
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
             'restricted_end': expected 5 bp end of amplicon at end-repaired restrictioin site 
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
        ### Functions
#### Mut_col functions 
    1. Mut_col_reverse_guide_only: in the case that the primer is oriented correctly, but the sgRNA is in reverse
    2. Mut_col_reverse: in the case that both the primer and sgRNA are in reverse
    3. Mut_col_Forward: in the case that the primer and sgRNA are oriented correctly

## 7. Group reads and count 
     retains columns sequence, length, sate, type, mut name, DSB
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
