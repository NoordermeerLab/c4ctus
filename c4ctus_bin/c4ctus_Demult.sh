#!/bin/sh


#   PIPELINE TO DEMULTIPLEX FastQ READS AS DONE IN HTSSTATION
#   author=Benoit Moindrot
#   credit goes to bbcf/bbcflib
#
#  V1   -  28/01/2019
#  V1.1 -  31/01/2019
#    * Boolean to trim or not the reads
#    * Run Exonerate in parallel: bug corrected in the variable i -> iii
#    * Find Exonerate output: bug corrected in the regex
#  V1.2 - 13/02/2019
#    * change variable i to variable x
#
#  V2 - 20-22/02/2019
#    * change Clean-Up and file naming (some files had the .fastq extension and where not fastq files)
#    * change of the argument parser and usage
#    * Add demultiplexing statistics and Report.pdf  (need R, with ggplot2 and grid libraries)
#    * compress output (fastq files, Exonerate output)
#  V2.1 - 25/02/2019
#    * Python function called here: their name has been changed


usage ()
{
	echo 'c4ctus_Demult.sh -f <FastQ> -p <PrimerFile> -g <GroupName>'
	echo '                 [--Ncpu <Nbr cpu>  --trim <True|False>]'
	echo '                 [options]'
	echo ''
	echo 'MANDATORY ARGUMENTS'
	echo '   -f/--fastq       fastq file'
	echo '                    Note: Need to be decompressed'
	echo ''
	echo '   -p/--primer      Primer File'
	echo '                    Note: see https://github.com/bbcf/bbcflib/blob/'
	echo '                              master/doc/tutorial_demultiplexing.rst'
	echo ''
	echo '   -g/--group       Prefix for file Naming after demultiplexing'
	echo ''
	echo 'DEMUTLIPLEXING PARAMETERS [OPTIONAL / DEFAULT VALUE ARE SHOWN]'
	echo '   -c/--Ncpu     Maximum Number of core to use simultaneously'
	echo '                     [DEFAULT: --Ncpu 15]'
	echo ''
	echo '   -n            Search the primer from base i'
	echo '                     [DEFAULT: -n 1]'
	echo ''
	echo '   -x            Search the primer in the next n bps of the reads [i to i+n]'
	echo '                     [DEFAULT: -x 22]'
	echo ''
	echo '   -s            Minimum score for Exonerate'
	echo '                 MinScore can be evaluated as minScore = (5 x L) - (9 x M)'
	echo '                 Where L is the length of the prefix of the primer sequence given'
	echo '                 in the primer file and M is the number of mismatches allowed'
	echo '                     [DEFAULT: -s 77]'
	echo ''
	echo '   -l            Length of the reads to align'
	echo '                     [DEFAULT: -l 30]'
	echo ''
	echo '   -t/--trim     Trim reads? Can be True or False'
	echo '                     [DEFAULT: --trim True]'
	echo ''
	exit
}



#################################
#--    Parse Arguments        --#
#################################

while test -n "$1"; do
    case "$1" in
      -f|--fastq)
          FastQReads=$2
          shift 2
          ;;  
      -p|--primer)
          PrimerFile=$2
          shift 2
          ;;  
      -g|--group)
          Group=$2
          shift 2
          ;;
      -c|--Ncpu)
          Ncpu=$2
          shift 2
          ;;
      -n)
          n=$2
          shift 2
          ;;
      -x)
          x=$2
          shift 2
          ;;
      -s)
          s=$2
          shift 2
          ;;
      -l)
          l=$2
          shift 2
          ;;
      -t|--trim)
          Trim=$2
          [[ "$Trim" = "True" || "$Trim" = "False" ]] || usage
          shift 2
          ;;
      *)
          echo "Unknown argument: $1"
          usage
          exit 1
          ;;  
    esac
done 

# for mandatory arguments
if [ -z "$FastQReads" ]; then echo "FastQReads is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$PrimerFile" ]; then echo "PrimerFile is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$Group" ]; then echo "Group is not set => abort"; echo "####"; usage; exit 1;fi


# for optional aguments
if [ -z "$Ncpu" ]; then echo "'  Max number of core to use' not set => Using default value (15)"; Ncpu=15;fi
if [ -z "$n" ]; then echo "'  Search the primer from base i' not set => Using default value (1)"; n=1;fi
if [ -z "$x" ]; then echo "'  Search the primer in the next n bps' not set => Using default value (22)"; x=22;fi
if [ -z "$s" ]; then echo "'  Minimum score for Exonerate' not set => Using default value (77)"; s=77;fi
if [ -z "$l" ]; then echo "'  Length of the reads to align' not set => Using default value (30)"; l=30;fi
if [ -z "$Trim" ]; then echo "'  Trim read?' not set => Using default value (True)"; Trim="True";fi
echo "======="

# Check for files existence. Die if files are not found
if [ ! -e "$FastQReads" ]; then echo "FastQ Reads: ${FastQReads} => FILE NOT FOUND"; exit 1;fi
if [ ! -e "$PrimerFile" ]; then echo "Primer File: ${PrimerFile} => FILE NOT FOUND"; exit 1;fi


echo "Demultiplexing Parameters"
echo "FastQReads=${FastQReads}"
echo "PrimerFile=${PrimerFile}"
echo "Group=${Group}"
echo "Ncpu=$Ncpu   n=$n   x=$x    s=$s    l=$l   Trim=$Trim"
echo "-----";echo " "; echo " "











##################################
#-- Step 1: Convert FastQ To Fasta, using a dedicated script. 
##################################
# FastQ Reads
FastQBaseName=$(basename "$FastQReads")
FastQShortName="${FastQBaseName%.*}"

# The Fasta header contains the information relative to the read name, and the sequence quality
echo "############ Processing $FastQReads ###############"
echo " "
echo "----- Converting to fasta : $(date +'%T') -----"
MyFastqToFasta.py \
		--input $FastQReads \
		--output ${FastQShortName}.fasta \
		--start $n \
		--length $x \




##################################
#-- Step 2: Assign to a given Primer, using the primer File 
#           Rely on exonerate software
##################################
echo " "
echo "----- Alignment to $PrimerFile : $(date +'%T') -----"

#How many job to run. By default here, the jobs are s
Nreads2="$(( $(wc -l ${FastQShortName}.fasta | awk '{print $1}') / 2))" 
Njob=$( perl -e 'use POSIX; print ceil($ARGV[0] / 5e6)' "$Nreads2" )
echo "$Nreads2 reads => splitting into $Njob jobs"


#Run up to Ncpu jobs in parallel, if reads number > 5e6
(
job=1
while [ "$job" -le "$Njob" ]; do
   ((iii=iii%Ncpu)); ((iii++==0)) && wait
   JobLeadingZero=$( printf "%03d" $job)
   echo "     - running Exonerate job $job / $Njob"
   exonerate --showalignment no --model affine:local \
         -o -4 -e -12 \
         -s 10 \
         --querychunkid $job  --querychunktotal $Njob \
         --query ${FastQShortName}.fasta \
         --target $PrimerFile > ${FastQShortName}_RezExonerate.$JobLeadingZero.txt &
   job=$(($job+1))
done
wait
)



#combine Exonerate output into one file and remove the intermediate files
find . -maxdepth 1 -type f -regextype sed -regex ".*${FastQShortName}_RezExonerate.[0-9]\{1,3\}.txt*" -exec cat {} \; > ${FastQShortName}_RezExonerate_all.txt
find . -maxdepth 1 -type f -regextype sed -regex ".*${FastQShortName}_RezExonerate.[0-9]\{1,3\}.txt*" -exec rm -f {} \;
echo "     -> All Exonerate jobs done"
echo " "



##################################
#-- Step 3: Parse Exonerate Output
##################################
echo "----- Parse Exonerate output : $(date +'%T') -----"
ParseExonerate.py \
		--file ${FastQShortName}_RezExonerate_all.txt \
		--group $Group \
		--minScore $s \
		--length $l \
		--FromBase $n \
		--trim ${Trim} > ${FastQShortName}_DemultiplexingFileList.txt

echo " "




##################################
#-- Step 4: Cleanup + statistics on demultiplexing
##################################
echo "----- Clean-Up, Statistics and compress : $(date +'%T') -----"

mkdir -p Demultiplexing/debug

# Where read count are going to be stored
TEMP_Count=$(mktemp ReadCount.XXXXXXXX)

# List of files that will be gzip
ListFileToCompress=""


# Reading the FileList file: wc and move files
while IFS=$'\t' read -r key File
do
	Viewpoint="$( echo $key | sed -E 's/key=(.*)/\1/' )"
	IndFastQ="$( echo $File | sed -E 's/file=(.*)/\1/' )"
	
	#echo "   > doing VP=${Viewpoint}  file=${IndFastQ}"
	
	#Calculate Read count
	if [[ "$Viewpoint" = "ambiguous_fastq" ]]			# fastQ: Ambiguous reads in a fastQ format
	then
		ReadCount="$(( $(wc -l ${IndFastQ} | awk '{print $1}') / 4))"
		mv $IndFastQ Demultiplexing/${IndFastQ}.fastq
		ListFileToCompress="$ListFileToCompress Demultiplexing/${IndFastQ}.fastq"
		
	elif [[ "$Viewpoint" = "unaligned" ]]				# Vulgar: Unaligned reads - alignement score < Minimum score for Exonerate (-s)
	then
		ReadCount="$(( $(wc -l ${IndFastQ} | awk '{print $1}') / 1))"
		mv $IndFastQ Demultiplexing/${IndFastQ}.Exonerate.txt
		ListFileToCompress="$ListFileToCompress Demultiplexing/${IndFastQ}.Exonerate.txt"
		
	elif [[ "$Viewpoint" = "discarded" ]]				# fastQ: Score OK, but length < l/2
	then
		ReadCount="$(( $(wc -l ${IndFastQ} | awk '{print $1}') / 4))"
		mv $IndFastQ Demultiplexing/${IndFastQ}.fastq
		ListFileToCompress="$ListFileToCompress Demultiplexing/${IndFastQ}.fastq"
		
	elif [[ "$Viewpoint" = "ambiguous" ]]				# ???: Ambiguous reads in an unknown format
	then
		mv $IndFastQ Demultiplexing/${IndFastQ}
		ListFileToCompress="$ListFileToCompress Demultiplexing/${IndFastQ}"
		
	else												# fastQ: Reads assigned to a viewpoint
		ReadCount="$(( $(wc -l ${IndFastQ} | awk '{print $1}') / 4))"
		mv $IndFastQ Demultiplexing/${IndFastQ}.fastq
		ListFileToCompress="$ListFileToCompress Demultiplexing/${IndFastQ}.fastq"
		
	fi
	
	# Print in file, if $Viewpoint is not ambiguous (number already counted in the fastq format)
	if [[ "$Viewpoint" != "ambiguous" ]]
	then
		echo -e "${IndFastQ}\t${Viewpoint}\t${ReadCount}" >> ${TEMP_Count}
	fi
	
done < ${FastQShortName}_DemultiplexingFileList.txt


#Plot Demultplexing Statistics
PlotDemultStat.r ${TEMP_Count} ${Nreads2} ${Group}
mv ${Group}_DemutliplexingReport.pdf Demultiplexing/${Group}_DemutliplexingReport.pdf
mv ${TEMP_Count} Demultiplexing/${TEMP_Count}
echo "    > Statistics (Read Cound) done"

#gzip files (multiprocess)
echo $ListFileToCompress | xargs -P 10 -n 1 sh -c 'gzip $1 ; sleep 2' {}
echo "    > Files Compression done"


#Move other files
mv ${FastQShortName}_RezExonerate_all.txt Demultiplexing/debug/${FastQShortName}_RezExonerate_all.txt
mv ${FastQShortName}.fasta Demultiplexing/debug/${FastQShortName}.fasta
mv ${FastQShortName}_DemultiplexingFileList.txt Demultiplexing/${FastQShortName}_DemultiplexingFileList.txt


echo "############ All Done ###############"
echo "############ $(date +'%T') ###############"