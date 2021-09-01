#!/bin/sh


#   PIPELINE TO COMPUTE 4C-SCORES AS DONE IN HTSSTATION
#   author=Benoit Moindrot
#   credit goes to bbcf/bbcflib
#
#   V1
#
#   V2 - 22/02/2019
#      * change of the argument parser and usage
#      * check for file existence before computing
#   V2.1 - 28/02/2019
#      * Assembly name corrected in the SQL file
#   V2.2 - 16/08/2020
#	   * remove the presorted assumption for bedtools intersect


usage ()
{
	echo 'c4ctus_4Cseq.sh --lib4C <4CLib> --genome <Assembly>'
	echo '                -f <MapDensityFwd> -r <MapDensityRev>'
	echo '                -e <RegToExclude> -w <RollWindowSize>'
	echo '                [-n <NickName>]'
	echo ''
	echo '    -c/--lib4C           4C-Library in a bed format'
	echo '                       Note: This library need to be processed by the CreateLib'
	echo '                             module. It can be repeat masked or not'
	echo ''
	echo '    -g/--genome        Genome Assembly (dm6, mm9, mm10, hg19, hg38)'
	echo ''
	echo '    -f/--densF         Mapping density FORWARD (BED format)'
	echo '                       Note: Obtained by the Mapping module'
	echo ''
	echo '    -r/--densR         Mapping density REVERSE (BED format)'
	echo '                       Note: Obtained by the Mapping module'
	echo ''
	echo '    -e/--Excl          Genomic coordinate to Exclude'
	echo '                       Ex: chr3R:16730708-16737272'
	echo ''
	echo '    -w|--window        Window size for the rolling mean (smoothing of 4C data)'
	echo '                       Ex: 11'
	echo ''
	echo '    -n|--NickName      Optional <Sample>_<Primer>_<Assembly> NickName for'
	echo '                       file naming   Ex: S2_rep1_Ubx_dm6'
	exit
	
}



##################################
#-- Step 0: Parse Arguments
##################################


while test -n "$1"; do
    case "$1" in
      -c|--lib4C)
          Lib4C=$2
          shift 2
          ;;  
      -g|--genome)
          Assembly=$2
          shift 2
          ;;  
      -f|--densF)
          MapDensityFwd=$2
          shift 2
          ;;  
      -r|--densR)
          MapDensityRev=$2
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
      -n|--NickName)
          NickName=$2
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
if [ -z "$Lib4C" ]; then echo "Lib4C is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$Assembly" ]; then echo "Genome assembly is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$MapDensityFwd" ]; then echo "MapDensityFwd is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$MapDensityRev" ]; then echo "MapDensityRev is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$RegToExlude" ]; then echo "RegToExlude is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$RollMeanSize" ]; then echo "RollMeanSize is not set => abort"; echo "####"; usage; exit 1;fi


# Optional arguments
if [ -z "$NickName" ]; then NickName=0;fi


# Guessing file Nickname or Parsing the optional NickName argument
DensityFileName=$(basename "$MapDensityFwd")
if [[ $NickName = "0" ]]
then
	re="(^[[:graph:]]+)_Density_(fwd|rev|FWD|REV)"
	if [[ $DensityFileName =~ $re ]]
	then
		NickName=${BASH_REMATCH[1]}
	else
		echo "Cannot recognized the sample/Primer name"
		echo "Provide a <Sample>_<Primer>_<Assembly> NickName as optional argument"
		echo ''
		usage
		exit
	fi
fi


# Check for files existence. Die if files are not found
if [ ! -e "$Lib4C" ]; then echo "Lib4C: ${Lib4C} => FILE NOT FOUND"; exit 1;fi
if [ ! -e "$MapDensityFwd" ]; then echo "Density Forward: ${MapDensityFwd} => FILE NOT FOUND"; exit 1;fi
if [ ! -e "$MapDensityRev" ]; then echo "Density Reverse: ${MapDensityRev} => FILE NOT FOUND"; exit 1;fi



#Rep Count [as in HTSstation]  ### !!!! EDIT THE FOLLOWING TWO LINES TO UPDATE THE PATH TO THE LogCount.cnt file !!! ###
count=$(cat /store/EQUIPES/CHRODY/Benoit/Lab/c4ctus/c4ctus_bin/LogCount.cnt)
echo $(($count+1)) > /store/EQUIPES/CHRODY/Benoit/Lab/c4ctus/c4ctus_bin/LogCount.cnt


echo "############ 4C-seq module ###############"
echo "## Files will be named with the following prefix: $NickName"
echo "## Genome Assembly = $Assembly"
echo "## The two density files are: $MapDensityFwd  and  $MapDensityRev"
echo "## The 4C-Lib used if $Lib4C"
echo "## The region to exclude is: $RegToExlude"
echo "## rep count = ${count}"
echo "## Window size = ${RollMeanSize} for smoothing the RAW 4C data"
echo ""








##################################
#-- Step 1: Generate Intersect between 4Clib and mapping density
#           Score kept = number of base of overlap
##################################
echo "----- Generating intersect between 4Clib and mapping density : $(date +'%T') -----"

bedtools intersect -wao \
		-a $Lib4C \
		-b $MapDensityFwd $MapDensityRev | \
		awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$9*$10}' | \
		sort -k1,1 -k2,2n > ${NickName}_4CLibIntersect.txt


ParseIntersect_4CLib_Density.pl ${NickName}_4CLibIntersect.txt \
		> ${NickName}_4CLibIntersect_Parsed.txt
echo "> done"; echo ''




##################################
#-- Step 2: Combine scores of the two ends of a given restriction fragment
#           Parse and Sort
#           Output = Bed format of all restriction fragments
#                  AND
#                    Bed format of restriction fragments with score > 0
##################################
echo "----- Combining scores from both restriction fragment ends : $(date +'%T') -----"

segToFrag.awk -v reg2Excl=${RegToExlude} ${NickName}_4CLibIntersect_Parsed.txt > ${NickName}_awkOutput.txt
ParseAWK_4Cdata_V1.py ${NickName}_awkOutput.txt > ${NickName}_awkOutput_Parsed.txt
bedtools sort -i ${NickName}_awkOutput_Parsed.txt > segToFrag_${NickName}_rep${count}_all.bedraph

#Create bed file with restriction fragment having score > 0
awk '$4 > 0 {print}' segToFrag_${NickName}_rep${count}_all.bedraph > segToFrag_${NickName}_rep${count}_NonZeros.bedraph
echo "> done"; echo ''




##################################
#-- Step 3: Generate SQL file
#           NB: The alignement is a bit messy, but sqlite3 does not tolerate tabs...
##################################
echo "----- Generating SQL file with results : $(date +'%T') -----"


#List of unique chromosomes
ListCHR=($(cut -f 1 segToFrag_${NickName}_rep${count}_all.bedraph | sort | uniq))

#Create Temp directory
TempDir=$(openssl rand -hex 3)
TempDir="TempSplit${TempDir}"
mkdir -p $TempDir

#split BED files per chromosome  (Works, but is rather slow...)
for chr in "${ListCHR[@]}"
do
  grep -w $chr segToFrag_${NickName}_rep${count}_all.bedraph | awk '{print $2,$3,$4}' OFS='\t' \
  		> $TempDir/$chr.output.bed
done


#Generate empty SQL file
if [ -e "segToFrag_${NickName}_rep${count}_all.sql" ]; then
  echo "SQL exists - the corresponding file has been deleted to be replaced"
  rm "segToFrag_${NickName}_rep${count}_all.sql"
fi 
sqlite3 segToFrag_${NickName}_rep${count}_all.sql ""


#Fill-in SQL file with 4C-seq data (one table per chromosome)
for i in "${ListCHR[@]}"
do
sqlite3 segToFrag_${NickName}_rep${count}_all.sql <<!
CREATE TABLE $i (start INTEGER PRIMARY KEY, end INTEGER, score REAL);
.separator "\t"
.import $TempDir/$i.output.bed $i
!
rm $TempDir/$i.output.bed
done


#Delete Temp directory
rm -r $TempDir


#Add chrNames table to the SQL file   !!! line 254: Update path to a directory with ${Assembly}.chrom.sizes files
tmpfile=$(mktemp /tmp/abc-script.XXXXXX)
awk -v var=${Assembly} 'BEGIN {OFS="\t"} {print $0, var}' /store/EQUIPES/CHRODY/COMMON/_GENOME_SIZES/UCSC_${Assembly}.chrom.sizes > $tmpfile
sqlite3 segToFrag_${NickName}_rep${count}_all.sql <<!
CREATE TABLE chrNames (name TEXT, length INTEGER, assembly TEXT);
.separator "\t"
.import $tmpfile chrNames
!
rm $tmpfile




####################################################################
#-- STEP 4: Smoothing of RAW data (Running Mean)
#           Only performed on the chromosome containing the viewpoint
####################################################################
mkdir -p 4Cseq/debug

Smooth4C.r segToFrag_${NickName}_rep${count}_all.sql \
           ${RegToExlude} \
           ${RollMeanSize} \
           4Cseq/${NickName}_rep${count}


####################################################################
#-- STEP 4: 
#   Clean-up
####################################################################


mv ${NickName}_4CLibIntersect.txt 4Cseq/debug/${NickName}_4CLibIntersect.txt
mv ${NickName}_4CLibIntersect_Parsed.txt 4Cseq/debug/${NickName}_4CLibIntersect_Parsed.txt
mv ${NickName}_awkOutput.txt 4Cseq/debug/${NickName}_awkOutput.txt
mv ${NickName}_awkOutput_Parsed.txt 4Cseq/debug/${NickName}_awkOutput_Parsed.txt

mv segToFrag_${NickName}_rep${count}_all.bedraph 4Cseq/segToFrag_${NickName}_rep${count}_all.bedraph
mv segToFrag_${NickName}_rep${count}_all.sql 4Cseq/segToFrag_${NickName}_rep${count}_all.sql
mv segToFrag_${NickName}_rep${count}_NonZeros.bedraph 4Cseq/segToFrag_${NickName}_rep${count}_NonZeros.bedraph


echo "############ All Done ###############"
echo "############ $(date +'%T') ###############"
