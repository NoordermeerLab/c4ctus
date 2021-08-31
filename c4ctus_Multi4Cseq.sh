#!/bin/bash


# Script to streamline the 4C-seq module (batch mode)
#
#   author=Benoit Moindrot
#
# 31/08/2021 => Change the way arguments are set-up



usage ()
{
	echo 'c4ctus_Multi4Cseq <4CLib> <PrimerFile> <GroupName> <Assembly> <NumTHREADS>'
	echo ''
	echo '    -l/--lib4C          4C-Library in a bed format'
	echo '                        Note: This library need to be processed by the CreateLib'
	echo '                             module. It can be repeat masked or not'
	echo ''
	echo '    -p/--primer         Primer File (same as used for demultiplexing'
	echo '                        Note: see https://github.com/bbcf/bbcflib/blob/'
	echo '                              master/doc/tutorial_demultiplexing.rst'
	echo ''
	echo '    -g/--group          Prefix USED during demultiplexing step for naming'
	echo '                        files'
	echo ''
	echo '    -a/--genome         Genome Assembly'
	echo '                        ex: dm6, mm10, hg19...'
	echo '                        Note: Only works for dm6 until now'
	echo ''
	echo '    [-w/--window ]      Window size for the rolling mean (smoothing of 4C data)'
	echo '                           [DEFAULT: --window 11]'
	echo ''
	echo '    [-c/--Ncpu ]        Optional: Max Number of THREADS to run in parallel'
	echo '                           [DEFAULT: --Ncpu 8]'
	exit
}




##################################
#-- Die if wrong number of arguments
#   or if BASH_VERSION < 4
##################################
if [ "${BASH_VERSINFO}" -lt 4 ]; then echo "Sorry, you need at least bash-4.0 to run this script." >&2; exit 1; fi




##################################
#-- Step 0: Parse Arguments
##################################


while test -n "$1"; do
    case "$1" in
      -l|--lib4C)
          Lib4C=$2
          shift 2
          ;;
      -p|--primer)
          PrimerFile=$2
          shift 2
          ;;
      -a|--genome)
          Assembly=$2
          shift 2
          ;;
      -g|--group)
          GroupName=$2
          shift 2
          ;;
      -w|--window)
          RollMeanSize=$2
          shift 2
          ;;
      -c|--Ncpu)
          Ncpu=$2
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
if [ -z "$Lib4C" ]; then echo "4C-Library is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$PrimerFile" ]; then echo "PrimerFile is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$GroupName" ]; then echo "GroupName is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$Assembly" ]; then echo "Genome Assembly is not set => abort"; echo "####"; usage; exit 1;fi

# for optional aguments
if [ -z "$Ncpu" ]; then echo "'  Max number of core to use' not set => Using default value (8)"; Ncpu=8;fi
if [ -z "$RollMeanSize" ]; then echo "'  window size' not set => Using default value (11)"; RollMeanSize=11;fi


# Check for file existence. Die if not found
if [ ! -e "$Lib4C" ]; then echo "4C-Library: ${Lib4C} => FILE NOT FOUND"; exit 1;fi
if [ ! -e "$PrimerFile" ]; then echo "Primer File: ${PrimerFile} => FILE NOT FOUND"; exit 1;fi



##################################
#-- Starting
##################################
echo "----- Starting : $(date +'%T') -----"
echo "## LOG FILES WILL BE FOUND IN THE 4Cseq FOLDER"
echo "## 4CLib used = ${Lib4C}"
echo "## Running up to ${Ncpu} jobs in parallel"
echo '##########';echo '';echo ''





##################################
#-- Parse Primer File
#   to get a list of viewpoints and Exclude region in an associative array
##################################
declare -A ListVPandExclude


while read -r line; do
	if [[ $line =~ ">" ]]
	then
		CurrentViewpoint=$(echo $line | perl -ne 'm/\>(.*?)\|/i; print $1')
		CurrentExclude=$(echo $line | perl -ne 'm/.*\|Exclude=(.*?)[\||\s]/i; print $1')
		ListVPandExclude[$CurrentViewpoint]=$CurrentExclude
	fi
done < "$PrimerFile"

#List the content of the associative array
echo "The viewpoint that will be processed: "
for elem in ${!ListVPandExclude[*]} ; do
        echo "  - VP \"${elem}\", RegToExclude : "${ListVPandExclude[${elem}]} ;
done
echo "";echo ""







mkdir -p 4Cseq

#Run up to Ncpu jobs in parallel
Njob=${#ListVPandExclude[@]}
KEYS=(${!ListVPandExclude[@]})

(
job=1
while [ "$job" -le "$Njob" ]; do
   ((iii=iii%Ncpu)); ((iii++==0)) && wait
   Viewpoint=${KEYS[$job-1]}
   ExcludeRegion=${ListVPandExclude[${KEYS[$job-1]}]}
   echo "Job ${job} / ${Njob} RUNNING : "
   echo "     - running 4Cseq module for Viewpoint=${Viewpoint} with Exclude=${ExcludeRegion}"
   Density_file_Nickname="Mapping/${GroupName}_${Viewpoint}_${Assembly}_Density"
   echo "       using:"
   echo "         Density_fwd=${Density_file_Nickname}_fwd.txt"
   echo "         Density_rev=${Density_file_Nickname}_rev.txt"
   echo "";echo ""
   c4ctus_4Cseq.sh --lib4C ${Lib4C} --genome ${Assembly} --densF ${Density_file_Nickname}_fwd.txt --densR ${Density_file_Nickname}_rev.txt --Excl ${ExcludeRegion} --window ${RollMeanSize} &>> 4Cseq/${GroupName}_${Viewpoint}_${Assembly}_4Cseq.log &
   #./TestSleep.sh &>> 4Cseq/_TUTU_${Viewpoint}_${Assembly}_4Cseq.log &
   sleep 1
   job=$(($job+1))
done
wait
)


echo "----- Finished : $(date +'%T') -----"