#!/bin/sh


# Script to streamline the Mapping of 4C-seq reads (batch mode)
#
#   author=Benoit Moindrot
#
# 21/02/2019 => Now adapted to use process fastq.gz files
# 25/02/2019 => Change usage (now works for dm6 & mm10)
# 04/03/2019 : Change usage and argument parsing




usage ()
{
	echo 'c4ctus_MultiMapping.sh -p <PrimerFile> -g <GroupName> --genome <Assembly>'
	echo '                       [--Ncpu <NumTHREADS>]'
	echo ''
	echo '             This script will process the demultiplexed files named as follow:'
	echo '             Demultiplexing/<GroupName>_<PrimerName>.fastq.gz'
	echo ''
	echo '    -p/--primer         Primer File (same as used for demultiplexing'
	echo '                        Note: see https://github.com/bbcf/bbcflib/blob/'
	echo '                              master/doc/tutorial_demultiplexing.rst'
	echo ''
	echo '    -g/--group          Prefix USED during demultiplexing step for naming'
	echo '                        files (demultiplexed Files should be now named as'
	echo '                        follow: Demultiplexing/<GroupName>_<Primer>.fastq.gz)'
	echo ''
	echo '    -a/--genome         Genome Assembly'
	echo '                        ex: dm6, mm10, hg19...'
	echo '                        Note: Only works for dm6 & mm10 until now'
	echo ''
	echo '   [-c/--Ncpu    ]      Optional: Max Number of THREADS to run in parallel]'
	echo '                           [DEFAULT: --Ncpu 8]'
	echo '   [-m/--maxhits ]      Optional: Maximum amounts of hits tolerated for'
	echo '                        the reads => criteria to filter the bam file'
	echo '                           [default = 5]'
	exit
}




while test -n "$1"; do
    case "$1" in
      -p|--primer)
          PrimerFile=$2
          shift 2
          ;;  
      -g|--group)
          GroupName=$2
          shift 2
          ;;
      -a|--genome)
          Assembly=$2
          [[ "$Assembly" = "dm6" || "$Assembly" = "mm10" ]] || usage
          shift 2
          ;;
      -c|--Ncpu)
          Ncpu=$2
          shift 2
          ;;
      -m|--maxhits)
          maxhits=$2
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
if [ -z "$PrimerFile" ]; then echo "PrimerFile is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$GroupName" ]; then echo "GroupName is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$Assembly" ]; then echo "Genome Assembly is not set => abort"; echo "####"; usage; exit 1;fi

# for optional aguments
if [ -z "$Ncpu" ]; then echo "'  Max number of core to use' not set => Using default value (8)"; Ncpu=8;fi
if [ -z "$maxhits" ]; then echo "'  maxhits number' not set => Using default value (5)"; maxhits=5;fi

# Check for files existence. Die if files are not found
if [ ! -e "$PrimerFile" ]; then echo "Primer File: ${PrimerFile} => FILE NOT FOUND"; exit 1;fi





echo "----- Starting : $(date +'%T') -----"
echo "LOG FILES WILL BE FOUND IN THE MAPPING FOLDER"
echo ''


##################################
#-- Parse Primer File
#   to get a list of viewpoints
##################################
#ListAllViewpoint=""
ArrayAllViewpoint=()

while read -r line; do
	if [[ $line =~ ">" ]]
	then
		CurrentViewpoint=$(echo $line | perl -ne 'm/\>(.*?)\|/i; print $1')
		#ListAllViewpoint="$ListAllViewpoint $CurrentViewpoint"
		ArrayAllViewpoint+=($CurrentViewpoint)
	fi
done < "$PrimerFile"

echo "List of All Viewpoints to process"
#echo $ListAllViewpoint
echo "${ArrayAllViewpoint[*]}"
echo ""


##################################
#-- Run Mapping in parallel
##################################
mkdir -p Mapping

#echo $ListAllViewpoint | xargs -P$Ncpu -n1 bash -c 'echo "$1 started"; ./TestSleep.sh' {}
#echo $ListAllViewpoint | xargs -P $Ncpu -n 1 sh -c 'echo Processing Demultiplexing/'"$GroupName"'_$1.fastq.gz ...' {}
#echo $ListAllViewpoint | xargs -P $Ncpu -n 1 sh -c 'echo Processing Demultiplexing/'"$GroupName"'_$1.fastq.gz ... ; c4ctus_Mapping.sh Demultiplexing/'"$GroupName"'_$1.fastq.gz '"$Assembly"'&> Mapping/'"$GroupName"'_Mapping_$1.log ; sleep 2' {}


(
job=1
while [ "$job" -le "${#ArrayAllViewpoint[@]}" ]; do
   ((iii=iii%Ncpu)); ((iii++==0)) && wait
   echo "Job ${job} / ${#ArrayAllViewpoint[@]} RUNNING : "
   Viewpoint=${ArrayAllViewpoint[$job-1]}
   echo "     - running Mapping module for Demultiplexing/${GroupName}_${Viewpoint}.fastq.gz"
   c4ctus_Mapping.sh -f Demultiplexing/${GroupName}_${Viewpoint}.fastq.gz --genome ${Assembly} --maxhits ${maxhits} &> Mapping/${GroupName}_Mapping_${Viewpoint}.log &
   #./TestSleep.sh &>> _TUTU_${Viewpoint}_${Assembly}_4Cseq.log &
   job=$(($job+1))
done
wait
)



echo "----- Finished : $(date +'%T') -----"

