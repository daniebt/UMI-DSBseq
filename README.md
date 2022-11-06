# UMI-DSBseq
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
