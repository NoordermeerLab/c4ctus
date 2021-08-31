#!/bin/sh


# Script to re-scale 4C tracks.
#    Based on a Left and Right Coordinates corresponding to 5 TAD + 5 restriction sites
#    Signal within this interval is 
#         - scaled to 1e6
#         - smoothed using a rolling mean
#
# Authors: Benoit Moindrot
#
# List of changes
#	v1.1	27/08/2020 => Typo corrected in usage
# v1.2  03/09/2020 => check for file existence




usage ()
{
	echo 'c4ctus_Normalize.sh -f <SQL> --Excl <chrA:pos1-pos2>'
	echo '                    --LeftTAD <Coordinate>'
	echo '                    --rightTAD <Coordinate>'
	echo '                    [-o <output file>]'
	echo ''
	echo '     -f/--sqlFile        4C-seq SQL file'
	echo ''
	echo '     -e/--Excl           Genomic coordinate to Exclude'
	echo '                         Ex: chr3R:16730708-16737272'
	echo ''
	echo '      -w|--window        Window size for the rolling mean (smoothing of 4C data)'
	echo '                         Ex: 11'
	echo ''
	echo 'Definition of Normalization boundaries. Usually consists in '
	echo '* 2 TAD + 5 restriction sites upstream'
	echo '* 2 TAD + 5 restriction sites downstream'
	echo '     -l/--LeftTAD        Left Coordinate for the normalization window'
	echo '     -r/--RightTAD       Right Coordinate for the normalization window'
	echo '                         (sum of scores will be set to 1e6 in this window)'
	echo ''
	echo '     [-o/--output]       [Optional: output file prefix]'
	echo '                         if provided, file will be named' 
	echo '                           <output>_<CHR>_smoothed_<window>FragPerWin_Norm.bedgraph'
	echo '                         if not provided, file will be named after the SQL file'
	exit
}


# -------------------------------
#   Parse Arguments
# -------------------------------
while test -n "$1"; do
    case "$1" in
      -f|--sqlFile)
          SQLFile=$2
          shift 2
          ;;  
      -e|--Excl)
          RegToExlude=$2
          shift 2
          ;;  
      -w|--window)
          RollMeanSize=$2
          shift 2
          ;;
      -l|--LeftTAD)
          ExcludeLeft=$2
          shift 2
          ;;
      -r|--RightTAD)
          ExcludeRight=$2
          shift 2
          ;;
      -o|--output)
          OutFile=$2
          shift 2
          ;;
      *)
          echo "Unknown argument: $1"
          usage
          exit 1
          ;;  
    esac
done



# Check if mandatory arguments are set - otherwise die
if [ -z "$SQLFile" ]; then echo "SQLFile is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$RegToExlude" ]; then echo "RegToExlude is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$RollMeanSize" ]; then echo "RollMeanSize is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$ExcludeLeft" ]; then echo "Left TAD border is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$ExcludeRight" ]; then echo "Right TAD border is not set => abort"; echo "####"; usage; exit 1;fi


# Optional arguments
if [ -z "$OutFile" ]
then 
	SQLBaseName=$(basename "$SQLFile")
	OutFile=$(echo $SQLBaseName | sed 's/_all.sql//g')
	OutFile=$(echo $OutFile | sed 's/segToFrag_//g')
fi




# -------------------------------
#   Test For file existence
# -------------------------------
if [ ! -e "$SQLFile" ]; then echo "SQLfile not found"; exit; fi



Norm4C.r ${SQLFile} \
         ${RegToExlude} \
         ${RollMeanSize} \
         ${ExcludeLeft} \
         ${ExcludeRight} \
         ${OutFile}