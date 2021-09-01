#!/bin/sh


# Script to Compute domainogram analysis from 4C-seq data
#
#   author=Benoit Moindrot
#   adapted from Marion Leuleu





usage ()
{
	echo 'c4ctus_Domainogram.sh -i <SQL_file>'
	echo '                 [options]'
	echo ''
	echo 'MANDATORY ARGUMENTS'
	echo '   -i/--sqlfile     SQL file'
	echo ''
	echo '   -o/--outPrefix   Prefix used for file naming'
	echo ''
	echo ''
	echo 'DOMAINOGRAM PARAMETERS [OPTIONAL / DEFAULT VALUE ARE SHOWN]'
	echo '   -c/--Ncpu       Maximum Number of core to use simultaneously'
	echo '                     [DEFAULT: --Ncpu 8]'
	echo ''
	echo '   -d/--outDir     Output Directory, where files will be saved'
	echo '                     [DEFAULT: -n Domainogram]'
	echo ''
	exit
}



#################################
#--    Parse Arguments        --#
#################################

while test -n "$1"; do
    case "$1" in
      -i|--sqlfile)
          inSQLfile=$2
          shift 2
          ;;  
      -o|--outPrefix)
          outNamePrefix=$2
          shift 2
          ;;
      -c|--Ncpu)
          Ncpu=$2
          shift 2
          ;;
      -d|--outDir)
          outDir=$2
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
if [ -z "$inSQLfile" ]; then echo "inSQLfile is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$outNamePrefix" ]; then echo "outNamePrefix is not set => abort"; echo "####"; usage; exit 1;fi

# for optional aguments
if [ -z "$Ncpu" ]; then echo "'  Max number of core to use' not set => Using default value ()"; Ncpu=8;fi
if [ -z "$outDir" ]; then echo "'  Output directory' not set => Using default value (Domainogram)"; outDir="Domainogram";fi

# Check for files existence. Die if files are not found
if [ ! -e "$inSQLfile" ]; then echo "SQL file: ${inSQLfile} => FILE NOT FOUND"; exit 1;fi




echo "Starting: `date +%Y-%m-%d:%H:%M:%S`"
echo "Computing domainograms for ${inSQLfile} (${Ncpu} chromosomes simultaneously)"
echo "   Files will have the following prefix: ${outNamePrefix}_"
echo "   They will be found in the directory :${outDir}"
echo ""




#################################################################################
#--    Convert SQL to bedgraph file                                           --#
#--    [without chr*_random (unlocalized) and chrUn_* (unplaced clone contigs)] #
#################################################################################

#Get List of chromosomes [without chr*_random (unlocalized) and chrUn_* (unplaced clone contigs)]
#Selected chromosomes have the following name:  "chr_" or "chr__" where _ can be any character
ChrList=$(sqlite3 -noheader ${inSQLfile} "select name from chrNames WHERE name LIKE 'chr_' OR name LIKE 'chr__' ORDER BY name")
echo "Computing domainograms for the following chromosomes:"
echo "$ChrList"
echo "---------------------";echo ""


#Covert SQL to bedgraph for every chromosomes in ${ChrList}
TEMP_RAW_ALL=$(mktemp ${inSQLfile}.Raw4CAll.XXXXXXXX)

for chrom in ${ChrList[*]}
do
TEMP_RAW=$(mktemp Raw4CSingle${chrom}.XXXXXXXX)
#echo "     doing chr=${chrom} in ${TEMP_RAW}"
sqlite3 ${inSQLfile} <<!
.headers off
.mode tabs
.output "${TEMP_RAW}"
select "$chrom",* from $chrom order by start;
!
cat ${TEMP_RAW} >> ${TEMP_RAW_ALL}
rm ${TEMP_RAW}
done




###########################################################
#--    Compute domainograms                             --#
#--    up to ${Ncpu} chromosomes simultaneously         --#
###########################################################

arrCHR=($ChrList)
Njob=${#arrCHR[@]}

mkdir -p ${outDir}

(
job=1
while [ "$job" -le "$Njob" ]; do
   ((iii=iii%Ncpu)); ((iii++==0)) && wait
   CurrChr=${arrCHR[$job-1]}
   echo "Job ${job} / ${Njob} RUNNING : "
   echo "     - running Domainogram for chr=${CurrChr}"
   #R --vanilla --slave -f runDomainogram.R --args ${TEMP_RAW_ALL} ${outNamePrefix}_${CurrChr} ${outNamePrefix}_${CurrChr} ${CurrChr} 500 50 1 &> ${outNamePrefix}_${CurrChr}_Process.log &
   callDomainogram.sh ${TEMP_RAW_ALL} ${outNamePrefix} ${CurrChr} ${outDir} &> ${outDir}/${outNamePrefix}_${CurrChr}_Process.log &
   #./TestSleep.sh &>> ${outDir}/${outNamePrefix}_${CurrChr}_Process.log &
   sleep 1
   job=$(($job+1))
done
wait
)

rm ${TEMP_RAW_ALL}

echo "*****************************"
echo "done `date +%Y-%m-%d:%H:%M:%S`"

