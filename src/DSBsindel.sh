#!/bin/bash
##################################################################
# Author: Xin Wang                                                   
# Email: xin.wang@childrens.harvard.edu                            
# Copyright (c) 2021 Dr. Kaifu lab                                   
# PI: Kaifu Chen                                                   
# Description: 
#  DSBins contains four modules: 
#    iDSBquality, iDSBdetection, iDSBdeduplication and iDSBannotation. 
################################################################


### Usage function tells users how to run the software
helpFunction()
{
	echo "*********************************** how to use iDSBindel ***********************************"
	echo "Usage: sh $0 -a SampleID -i ReadsWithLargeInsertion -f HighQuality Forward Read -r HighQuality Reverse Read -b WorkingDirectory -o OutputFolder -p SoftwareDirectory [Options]"
	echo ""
	echo "iDSBindel: The software can detect small insertions/ small deletions/ no changes/ at DNA double strand break sites"
	echo -e "\t Features:"
	echo -e "\t -- Ignore the sequence errors that happened at two side of MATA regions"
	echo -e "\t -- Measure the read counts of each short insertions events, short deletion events and no change events"
	echo -e "\t -- Define the best sequence to represent each event."
	echo ""
	echo ""
	echo -e "Request Parameters:"
	echo -e "\t-a Sample Id (Example: yYY398-B_S10)"
	echo -e "\t-i Reads With Large Insertion events (The first column is read ID)"
	echo -e "\t-f High quality forward reads"
	echo -e "\t-r High quality reverse reads"
	echo -e "\t-b The working directory, where the raw read stored and we perform the analyses of the large insertion"
	echo -e "\t-o Name of Output folder"
	echo -e "\t-p Software installed Path, Attation: This required to install BWA in the same folder (Default:"")" 
	echo ""
	echo ""
	echo -e "Optional Parameters:"
	echo "" 
	echo -e "Optional Parameters -- overall requirements:"
	echo -e "\t-n Number of threads (Default: 15)"
	echo -e "\t-gs Genome sequence, we only used chrIII fasta file (Must corrected with chromosome ID) "
	echo "" 
	echo -e "Optional Parameters of Mapping:"
	echo -e "\tOutsoftware will use bwa mem ... "
	
	echo -e "iDSBindel: Trimm the index:"
	echo -e "\t-si The size of left custom index (Default 3)"  
	echo -e "\t-sr The size of right custom index (Default 3)"
	echo "" 

	echo -e "iDSBindel: Define the MATA information:"
	echo -e "\t-il minimum length of large insertion (Default 10bp)"  
	echo -e "\t-ms Total size of whole MATA region (Default 84)"
	echo -e "\t-mc Mapped chromosme of MATA reference position (Default chrIII)"
	echo -e "\t-ms Mapped start site of MATA reference position (Default 294300)"
	echo -e "\t-me Mapped end site of MATA reference position (Default 294500)"
	echo -e "\t-ml Size of left MATA region, here we allowed 4 nucleotide shift considering the microhomology or sequence errors (Default 45)"
	echo -e "\t-mr Size of left MATA region, here we allowed 4 nucleotide shift considering the microhomology or sequence errors (Default 51)"
	echo "" 
	echo -e "iDSBindel: Define the Primer information:"
	echo "" 
	echo -e "\t-ps The size of left primer, extended 5bp if there are deletion on the primer (default 25)"
	echo -e "\t-pf The size of right primer, extended 5bp if there are deletion on the primer (default 22)"
	
	echo -e "iDSBindel: Define the junction size to determine the unique of events:"
	echo -e "\t-u The collect upstream and downstream for the unique of deletion or insertions from the raw reads (default 5bp)"
	echo -e "\t-uc The read counts cut-off to determine confident indel(default 5)"
	echo "" 
	echo ""
	
	echo "Alternatively, you could run the pipeline step by step: "
	echo -e "\t\t Extract the reads that did not contain the large insertions (Optional)"
	echo -e "\t\t Map the reads to the genome "
	echo -e "\t\t Detect and Measure the short indel events "
	echo "" 
	echo ""
	echo -e "\t-h help"
	
	echo "For more detail information, please feel free to contact: xin.wang@childrens.harvard.edu"
	echo "**************************************************"
   exit 1 # Exit script after printing help
}


# Get absolute path for scripts and check if required scripts exist
in=$PWD


### Set default options
# some general paramaters including mata locus thread, and software path
# default threads
nproc=15
# default software path 
softwarepath=''

# default the customed index size
Findexsize=3
Rindexsize=3

# default whole MAT length
Matasize=84
# default Mata chromosome 
Matachr="ChrIII"
# default Mata start site
Matastart=294300
# default Mata end site
Mataend=294500

# default minimum insertion length
Ilength=10

# default parameter for first cuting-edge deduplication
### the size of primer (We allowed 5bp frameshift)
LeftPrimerSize=25
RightPrimerSize=22


### the size of left and right MAT size (we also allowed 5bp frameshift)
LeftMataSize=45
RightMataSize=51

### The read count to define the unique events
ReadCount=5

### The collect upstream and downstream for the unique of deletion or insertions from the raw reads (5bp)
Cutstream=5


while getopts "a:i:f:r:o:p:gs:n:si:sr:il:ms:mc:mb:me:ml:mr:ps:pf:u:uc" opt
do
   case "$opt" in
      a ) SampleID="$OPTARG" ;;
	  b ) IN="$OPTARG" ;;
      i ) INSERT="$OPTARG" ;;
	  f ) FORWARD="$OPTARG" ;;
	  r ) REVERSE="$OPTARG" ;;
      o ) OUTPUT="$OPTARG" ;;
	  p ) softpath="$OPTARG" ;;
  	 
	  
      gs) genomeseq="$OPTARG" ;;
      n ) nproc="$OPTARG" ;;
	  si) Findexsize="$OPTARG" ;;
	  sr) Rindexsize="$OPTARG" ;;
	  
	  il ) Ilength="$OPTARG" ;;
      ms ) Matasize="$OPTARG" ;;
	  mc ) Matachr="$OPTARG" ;;
      mb ) Matastart="$OPTARG" ;;
	  me ) Mataend="$OPTARG" ;;
	  ml ) LeftMataSize="$OPTARG" ;;
  	  mr ) RightMataSize="$OPTARG" ;;
	    

      ps ) LeftPrimerSize="$OPTARG" ;;
	  pf ) RightMataSize="$OPTARG" ;;
      u ) Cutstream="$OPTARG" ;;
	  
	 
	  uc ) ReadCount="$OPTARG" ;;

      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


# defualt script path

srcDir=${softpath}/iDSBindel/src
# default genome sequence
Chr3=${softpath}/iDSBindel/Database/chr3.fasta


# Print helpFunction in case parameters are empty
if [ -z "${SampleID}" ] 
then
   echo "*** error: input Sample ID must be provided ***";
   helpFunction
fi

if [ -z "${INSERT}" ] 
then
   echo "*** error: input insertion file must be provided ***";
   helpFunction
fi


if [ -z "${FORWARD}" ] 
then
   echo "*** error: input forward read must be provided ***";
   helpFunction
fi

if [ -z "${REVERSE}" ] 
then
   echo "*** error: input reverse read must be provided ***";
   helpFunction
fi

if [ -z "${OUTPUT}" ] 
then
   echo "*** error: Output folder must be defined ***";
   helpFunction
fi

if [ -z "${softpath}" ] 
then
   echo "*** error: input software path must be provided ***";
   helpFunction
fi


# Begin script in case all parameters are correct
echo "${SampleID}"
echo "${in}"
echo "${INSERT}"
echo "${FORWAR}"
echo "${REVERSE}"
echo "${softpath}"

# 

echo "All the paramter are sucessfully provided, now let's detect the large insertion event"


#### Here is the last version of large insertion events


echo "The job array of indel detection is started ..."
echo ""
date



#softpath=/lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/Updated_V0903/UsedMutants/TestSoftWare/Software

### Set up path file:

echo "Change to the WorkingPath/SampleID, where you store your raw fastq reads"

cd ${in}
cd  ${SampleID}

echo -e "Create a Small indel folder that store all the results of indel detection"
mkdir SmallIndel
cd SmallIndel

### Extract the non large inserted reads, Supposing we have already detected whether each reads contained large inserted sequence.

#cp ${in}/${SampleID}/Detection/${SampleID}_detected_final/${SampleID}_detected_combined.highqual.txt ./

echo -e "Extract the highquality reads that did not have large inserted sequence, we also trimmed two sides of index sequences"

perl ${softpath}/iDSBindel/src/Extract_insert_bothfasta_fastq_forindels_v2.pl  -i ${in}/${INSERT} -f ${in}/${FORWARD} -r ${in}/${REVERSE} -o ${SampleID}_noninsertion -m ${Findexsize} -n ${Rindexsize}

## bwa version 0.7.17-r1188
#bwa index ${Chr3}

echo  -e "Mapping the reads aganst chrIII using BWA mem..."

bwa mem ${Chr3} -t ${nproc} ${SampleID}_noninsertion.R1.fastq ${SampleID}_noninsertion.R2.fastq >${SampleID}_noninsertion.sam

echo -e  "Detect the short deletion, short insertion events..."

perl ${softpath}/iDSBindel/src/Indel_detection_fromsam_uniq_4.pl -i ${SampleID}_noninsertion.sam -o ${OUTPUT} -q ${SampleID} -m ${Matasize} -n ${Ilength} -c ${ReadCount} -a ${Matachr} -b ${Matastart} -e ${Mataend} -f ${LeftPrimerSize} -r ${RightPrimerSize} -l ${LeftMataSize} -d ${RightMataSize} -u ${Cutstream}

echo "Congratulation! Short indel detection job array is finished !"



