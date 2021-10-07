# iDSBindel
DNA double strand breaks are one of the most deleterious DNA lesions. The repair by nonhomologous end joining is error-prone, which would indced small insertions or  remove damaged or mismatched nucleotides.

To detect the products of DNA double strand break repair, we developed a novel techonology -- Break-Ins. In accompany with this technology, we also developed iDSBindel, which would sensentively detect short insertion, short deletion, short deletion and short insertion, and indel-free products at DNA double strans break site. iDSBindel can also calculate the frequency of each indel repair events.


# Description
iDSBindel: The software can detect small insertions, small deletions, no changes at DNA double strand break sites. Using iDSBindel, we can also measure the efficiency and accuracy of non-homologous end joining.

# Feature 
	 -- Ignore the sequence errors that happened at two side of MATA regions
	 -- Measure the read counts of each short insertions events, short deletion events and no change events
	 -- Define the best sequence to represent each event.

# Usage
```
Usage: sh DSBsindel.sh -a SampleID -i ReadsWithLargeInsertion -f HighQuality Forward Read -r HighQuality Reverse Read -b WorkingDirectory -o OutputFolder -p SoftwareDirectory [Options]



Request Parameters:
	-a Sample Id (Example: yYY398-B_S10)
	-i Reads With Large Insertion events (The first column is read ID)
	-f High quality forward reads
	-r High quality reverse reads
	-b The working directory, where the raw read stored and we perform the analyses of the large insertion
	-o Name of Output folder
	-p Software installed Path, Attation: This required to install BWA in the same folder (Default:)


Optional Parameters:

	iDSBindel: Overall requirements:
	-n Number of threads (Default: 15)
	-gs Genome sequence, we only used chrIII fasta file (Must corrected with chromosome ID)

	iDSBindel: Trimm the index:
	-si The size of left custom index (Default 3)
	-sr The size of right custom index (Default 3)

	iDSBindel: Define the MATA information:
	-il minimum length of large insertion (Default 10bp)
	-ms Total size of whole MATA region (Default 84)
	-mc Mapped chromosme of MATA reference position (Default chrIII)
	-ms Mapped start site of MATA reference position (Default 294300)
	-me Mapped end site of MATA reference position (Default 294500)
	-ml Size of left MATA region, here we allowed 4 nucleotide shift considering the microhomology or sequence errors (Default 45)
	-mr Size of left MATA region, here we allowed 4 nucleotide shift considering the microhomology or sequence errors (Default 51)

	iDSBindel: Define the Primer information:

	-ps The size of left primer, extended 5bp if there are deletion on the primer (default 25)
	-pf The size of right primer, extended 5bp if there are deletion on the primer (default 22)
	iDSBindel: Define the junction size to determine the unique of events:
	-u The collect upstream and downstream for the unique of deletion or insertions from the raw reads (default 5bp)
	-uc The read counts cut-off to determine confident indel(default 5)


	-h help
	
```


Alternatively, you could run the pipeline step by step:
	1. Extract the reads that did not contain the large insertions (Optional)
	2. Map the reads to the genome
	3. Detect and Measure the short indel events

# Output
Output fold contains several files to explain the results

1. Fasta files: These files includes the representive sequence for short insertion (outpout_insertion.fasta) and deletion events (outpout_deletion.fasta). These files also contain the information of sequenceID, mapping CIGAR, mapping junction index, read counts, indel size. 
2. Final Statistic: The read number of short insertion, short deleton, or no changes.
3. Indel stat file: The mapping junction index, read counts, indel size and detailed mapping information of each indel events
		
For more detail information, please feel free to contact: xin.wang@childrens.harvard.edu


