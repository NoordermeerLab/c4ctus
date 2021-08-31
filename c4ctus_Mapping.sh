#!/bin/sh


# To map the demultiplexed FastQ file using Bowtie2
# and process the SAM file in a similar manner than done by HtsStation
#
# CAREFUL: UNTIL KNOW, ONLY WORKS FOR DM6 GENOME (because only dm6 index path has been setup !)
#
#   author=Benoit Moindrot
#   credit goes to bbcf/bbcflib 
#
#
# List of changes
#     12/02/2019 : Add bam Index for the filtered bam file
#     17/02/2019 : Add Bowtie2 indexes for mm10 assembly
#     26/02/2019 : Change file naming (to trim fastq.gz from the file name)
#     04/03/2019 : Change usage and argument parsing
#     04/03/2019 : maxhits variable can now be set. Default value is 5




usage ()
{
	echo 'c4ctus_Mapping.sh  -f <FastQReadFile> -g <Genome Assembly>'
	echo '                   [--maxhits <VALUE>]'
	echo ''
	echo '   -f/--fastq          De-multiplexexed fastq file'
	echo '                       Note: can be compressed'
	echo '   -g/--genome         Genome Assembly'
	echo '                       ex: dm6, mm10, hg19, hg38...'
	echo '                       Note: Only works for dm6 and mm10 until now'
	echo ''
	echo '  [-m/--maxhits        Optional: Maximum amounts of hits tolerated for'
	echo '                       the reads => criteria to filter the bam file'
	echo '                       default = 5]'
	exit
}



#################################
#--    Parse Arguments        --#
#################################
#FastQReads=$1
#GenomeAssembly=$2

while test -n "$1"; do
    case "$1" in
      -f|--fastq)
          FastQReads=$2
          shift 2
          ;;  
      -g|--genome)
          GenomeAssembly=$2
          [[ "$GenomeAssembly" = "dm6" || "$GenomeAssembly" = "mm10" ]] || usage
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


# for mandatory arguments: check that have been defined
if [ -z "$FastQReads" ]; then echo "FastQReads is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$GenomeAssembly" ]; then echo "GenomeAssembly is not set => abort"; echo "####"; usage; exit 1;fi

# for optional arguments: use default value
if [ -z "$maxhits" ]; then echo "'  maxhits number' not set => Using default value (5)"; maxhits=5;fi

# Check for file existence. Die if not found
if [ ! -e "$FastQReads" ]; then echo "FastQ Reads: ${FastQReads} => FILE NOT FOUND"; exit 1;fi



#Trim the Extension of input fastq file (<SampleName>.fastq or <SampleName>.fastq.gz)
dbname=$(basename "$FastQReads")
case "$FastQReads" in
    *.fastq.gz)
         SampleName="${dbname%.*}"
         SampleName="${SampleName%.*}"
         ;;
    *.fastq)
         SampleName="${dbname%.*}"
         ;;
    *)
         echo "Extension not recognized... Files may not be named correctly"
         echo "this will affect the downstream file recongnition for the multi4C script"
         ;;
esac




####################################
#-- Bowtie2 indexes paths        --#
#   UPDATE ACORDING TO YOUR NEEDS  #
####################################

if [ $GenomeAssembly = 'dm6' ]
then
    BowtieIndexPath="/store/EQUIPES/CHRODY/COMMON/_GENOME_INDEX/Drosophila_melanogaster/UCSC/dm6/Bowtie2Index/genome"
elif [ $GenomeAssembly = 'mm10' ]
then
    BowtieIndexPath="/store/EQUIPES/CHRODY/COMMON/_GENOME_INDEX/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
else
    echo "!!! Genome not known or not available yet !!!"
    echo " "
    usage
    exit 
fi



####################################################################
#-- STEP 1: 
#   Map Using bowtie2
####################################################################
#     k = 20 => HtsStation paramterer = str(max(20,maxhits))
#     --end-to-end --sensitive => #default parameters in htsstation, default paramters in bowtie2 as well
echo "############ Processing $FastQReads ###############"
echo " "
echo "--------- Start Mapping reads : $(date +'%T') ---------"
echo "Bowtie2 Index: ${BowtieIndexPath}"

bowtie2 -k 20 \
        --end-to-end --sensitive \
        --un "${SampleName}_Unmapped.txt" \
        -x $BowtieIndexPath \
        -U $FastQReads \
        -S "${SampleName}_${GenomeAssembly}.sam"
echo "    done"

#Get Number of Mapped reads
nReads="$( samtools view -F 0x04 ${SampleName}_${GenomeAssembly}.sam  | cut -f1 | sort | uniq | wc -l )"
echo "    Number of mapped reads = " $nReads




####################################################################
#-- STEP 2: 
#   Python - To add NH_tag to bowtie output
#   And restrict bam file to reads having NH tag < maxhits
####################################################################
echo " "
echo "--------- Start Filtering reads : $(date +'%T') ---------"

#it needs a file sorted by READ name, which seems to be the case with bowtie2 output with the -k option
## !!!! BEFORE, CHECK THAT THE FILE IS SORTED BY READ NAME !!!! ##
AddNHtag.py ${SampleName}_${GenomeAssembly}.sam

#Sort BAM file by chromosome and coordinates
samtools sort -o ${SampleName}_${GenomeAssembly}_NHflaged_sorted.bam ${SampleName}_${GenomeAssembly}_NHflaged.bam

#Restrict bam file to reads having NH tag < maxhits
FilterBamFileBasedOnNH.py ${SampleName}_${GenomeAssembly}_NHflaged_sorted.bam $maxhits ${SampleName}_${GenomeAssembly}_filtered.bam
samtools index ${SampleName}_${GenomeAssembly}_filtered.bam

echo "    done"





####################################################################
#-- STEP 3: 
#   CREATE DENSITY FILES from filtered bam files
####################################################################
echo " "
echo "--------- Start Creating Density files : $(date +'%T') ---------"

#### Convert bam to bed, with score = 1/NH, in a strand specific manner
bedtools bamtobed -i ${SampleName}_${GenomeAssembly}_filtered.bam -tag NH > ${SampleName}_${GenomeAssembly}_filtered_toBed.bed
ParseAndScoreBed.pl ${SampleName}_${GenomeAssembly}_filtered_toBed.bed ${SampleName}_${GenomeAssembly}_fwd_temp.bed ${SampleName}_${GenomeAssembly}_rev_temp.bed

#### Convert bam to bedgraph using bedtools in a strand specific manner
bedtools genomecov -ibam ${SampleName}_${GenomeAssembly}_filtered.bam -bg -strand + > ${SampleName}_${GenomeAssembly}_fwd_temp.bedgraph
bedtools genomecov -ibam ${SampleName}_${GenomeAssembly}_filtered.bam -bg -strand - > ${SampleName}_${GenomeAssembly}_rev_temp.bedgraph

#### Bedtools Map
bedtools map -a ${SampleName}_${GenomeAssembly}_fwd_temp.bedgraph \
             -b ${SampleName}_${GenomeAssembly}_fwd_temp.bed \
             -c 4 -o sum | \
             awk -v nReads="$nReads" 'BEGIN {OFS="\t"} {print $1,$2,$3,$5/nReads*1e7}' \
             > ${SampleName}_${GenomeAssembly}_Density_fwd.txt

bedtools map -a ${SampleName}_${GenomeAssembly}_rev_temp.bedgraph \
             -b ${SampleName}_${GenomeAssembly}_rev_temp.bed \
             -c 4 -o sum | \
             awk -v nReads="$nReads" 'BEGIN {OFS="\t"} {print $1,$2,$3,$5/nReads*1e7}' \
             > ${SampleName}_${GenomeAssembly}_Density_rev.txt

echo "    done"



####################################################################
#-- STEP 3: 
#   Clean-up
####################################################################
rm ${SampleName}_${GenomeAssembly}_filtered_toBed.bed
rm ${SampleName}_${GenomeAssembly}_fwd_temp.bed
rm ${SampleName}_${GenomeAssembly}_rev_temp.bed
rm ${SampleName}_${GenomeAssembly}_fwd_temp.bedgraph
rm ${SampleName}_${GenomeAssembly}_rev_temp.bedgraph

mkdir -p Mapping/debug
mv ${SampleName}_Unmapped.txt Mapping/debug/${SampleName}_Unmapped.txt
mv ${SampleName}_${GenomeAssembly}.sam Mapping/debug/${SampleName}_${GenomeAssembly}.sam
mv ${SampleName}_${GenomeAssembly}_NHflaged.bam Mapping/debug/${SampleName}_${GenomeAssembly}_NHflaged.bam
mv ${SampleName}_${GenomeAssembly}_NHflaged_sorted.bam Mapping/debug/${SampleName}_${GenomeAssembly}_NHflaged_sorted.bam

mv ${SampleName}_${GenomeAssembly}_Density_fwd.txt Mapping/${SampleName}_${GenomeAssembly}_Density_fwd.txt
mv ${SampleName}_${GenomeAssembly}_Density_rev.txt Mapping/${SampleName}_${GenomeAssembly}_Density_rev.txt
mv ${SampleName}_${GenomeAssembly}_filtered.bam Mapping/${SampleName}_${GenomeAssembly}_filtered.bam
mv ${SampleName}_${GenomeAssembly}_filtered.bam.bai Mapping/${SampleName}_${GenomeAssembly}_filtered.bam.bai

echo " "
echo "############ ALL DONE : $(date +'%T') #############"

