#!/bin/sh


# Script to create the 4C libraries
# Require: Bedtools version > 2.24
# Two python scripts: ParseFragFile.py   ParseRmsk.py
# One perl script:    getRestEnzymeOccAndSeq.pl
#
#   author=Benoit Moindrot
#   credit goes to bbcf/bbcflib
#
#
# !!! disclaimer: only works for TypeI Libraries !!!
#
# V2 - 27/02/2019
#    * change of the argument parser and usage
#    * Now use BioPython to get the restriction from the restriction enzyme name




usage ()
{
	echo 'c4ctus_4CLib.sh --fasta <GenomeFasta>'
	echo '                --genome <GenomeNickName>'
	echo '                --PrimaryEnzyme <PrimaryRestrEnzyme>'
	echo '                --SecondaryEnzyme <SecondaryRestrEnzyme>'
	echo '                --length <SegmentLength>'
	echo '               [--output <OutputFolder]'
	echo '               [--rmsk GenomeRmskTab]'
	echo ''
	echo '!!! disclaimer: only works for TypeI Libraries !!!'
	echo ''
	echo '    -f/--fasta              Genome Sequence in fasta format'
	echo ''
	echo '    -g/--genome             Genome Assembly (dm6, mm9, mm10, hg19, hg38)'
	echo ''
	echo '    -p/--PrimaryEnzyme      Name of the 1st cutter enzyme. Ex DpnII, NlaIII...'
	echo ''
	echo '    -s/--SecondaryEnzyme    Name of the 2nd cutter enzyme. Ex DpnII, NlaIII...'
	echo ''
	echo '    -l/--length             Fragment length for generating the library'
	echo '                            30 is the classical value'
	echo ''
	echo '    [-o/--output]           Output folder where final file will be moved'
	echo '                            default= current folder'
	echo ''
	echo '    [-r/--rmsk]             Optional. Rmsk table (UCSC) to repeat mask the 4C lib'
	echo "                            Leave empty if you don't want to Repeat Mask it"
	echo '                            [chr = 6th column, START = 7th column, END = 8th column]'
	exit
}


# -------------------------------
#   Parse Arguments
# -------------------------------
while test -n "$1"; do
    case "$1" in
      -f|--fasta)
          GenomeFastaFile=$2
          shift 2
          ;;  
      -g|--genome)
          GenomeNickName=$2
          shift 2
          ;;  
      -p|--PrimaryEnzyme)
          PrimaryRestrEnz=$2
          shift 2
          ;;
      -s|--SecondaryEnzyme)
          SecondaryRestrEnz=$2
          shift 2
          ;;
      -l|--length)
          SegmentLength=$2
          shift 2
          ;;
      -r|--rmsk)
          GenomeRmskTab=$2
          shift 2
          ;;
      -o|--output)
          OutputDestination=$2
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
if [ -z "$GenomeFastaFile" ]; then echo "GenomeFastaFile is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$GenomeNickName" ]; then echo "GenomeNickName is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$PrimaryRestrEnz" ]; then echo "PrimaryRestrEnz is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$SecondaryRestrEnz" ]; then echo "SecondaryRestrEnz is not set => abort"; echo "####"; usage; exit 1;fi
if [ -z "$SegmentLength" ]; then echo "SegmentLength is not set => abort"; echo "####"; usage; exit 1;fi


# for optional aguments
if [ -z "$GenomeRmskTab" ]; then GenomeRmskTab="0";fi
if [ -z "$OutputDestination" ]; then OutputDestination=".";fi


# get Restriction site from Restriction enzyme Name
PrimaryRestrSite="$( Restr.py ${PrimaryRestrEnz})"
SecondaryRestrSite="$( Restr.py ${SecondaryRestrEnz})"



# -------------------------------
#   Test For file existence
# -------------------------------
if [ ! -e "$GenomeFastaFile" ]; then echo "genome file not found"; exit; fi

if [ ! "$GenomeRmskTab" = "0" ]
then
	if [ ! -e "$GenomeRmskTab" ]; then echo "RMSK file not found"; exit; fi
fi



# -------------------------------
#   Sort summary
# -------------------------------

echo "-----------"
echo "Genome Used:            $GenomeFastaFile"			#dm6.fa
echo "Genome Nickname:        $GenomeNickName"			#dm6
echo "Primary Restr Enzyme:   $PrimaryRestrEnz"		#DpnII
echo "Primary Restr Site:     $PrimaryRestrSite"			#GATC
echo "Secondary Restr Enzyme: $SecondaryRestrEnz"		#NlaIII
echo "Secondary Restr Site:   $SecondaryRestrSite"		#CATG
echo "Segment Length:         $SegmentLength"				#30
echo "Rmsk file used:         $GenomeRmskTab"				#dm6_rmsk.txt
echo "-----------"









# -------------------------------
#   Step 1:
#   Convert Genome to Uppercase
#   The following steps require an uppercase genome
# -------------------------------

echo "############ Converting Genome to Uppercase : $(date +'%T') #############"

#Variables Definition
GenomeFastaFile2=$(basename "$GenomeFastaFile")
GenomeFastaName="${GenomeFastaFile2%.*}"

#Actual Work
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' $GenomeFastaFile > ${GenomeFastaName}_uppercase.fa
echo "   > done"; echo ""; echo ""



# -------------------------------------
#   Step 2:
#   Get List of Restriction Fragments
# -------------------------------------

echo "############ Getting list of Restriction Fragments : $(date +'%T') #############"

getRestEnzymeOccAndSeq.pl \
         -i ${GenomeFastaName}_uppercase.fa \
         -m $PrimaryRestrSite \
         -s $SecondaryRestrSite \
         -l $SegmentLength \
         -o segFile.txt \
         -f FragFile.txt \
         -x logFile.txt \




# -------------------------------------
#   Step 3:
#   Parse the restriction File and keep only ValidFragments
#      (+ some cosmetic changes)
# -------------------------------------

var=$( ParseFragFile.py FragFile.txt)
echo "    > Generating two files: ${var}"
SeqInfoBedFile=$(echo $var | perl -lne 'print $& if /seqInfo\S*/')

echo "    > Renaming SegInfoBedFile"
C4LibUnmaskedFile="Library_${GenomeNickName}_${PrimaryRestrEnz}_${SecondaryRestrEnz}_${SegmentLength}bps_segmentsInfos.bed"
mv $SeqInfoBedFile $C4LibUnmaskedFile
echo "    > 4Cseq Library file (ValidFragments, but not repeat masked) = ${C4LibUnmaskedFile}"; echo ""; echo ""






if [ ! "$GenomeRmskTab" = "0" ]
then
	# -------------------------------------
	#   Step 4:
	#   Masking of Valid restriction fragments if covered by repeats
	# -------------------------------------
	
	echo "############ Computing overlap with rmsk : $(date +'%T') #############"
	echo "   > Repeat Masked table from UCSC = " $(basename "$GenomeRmskTab")
	
	# Step 4a: Convert the RMSK to bed file [chr = 6th column, START = 7th column, END = "8th column"]
	echo "   > Convert the UCSC table to bed file"
	awk 'BEGIN {OFS="\t"} {print $6,$7,$8}' $GenomeRmskTab > $(basename "$GenomeRmskTab")_Bed.bed
	
	# Step 4b: Compute overlap (using Bedtools)"]
	echo "   > Compute overlap between 4C-Lib and Rsmk (using CoverageBed from bedtools)"
	coverageBed -a $C4LibUnmaskedFile -b $(basename "$GenomeRmskTab")_Bed.bed > OverlapRMSK.txt
	
	# Step 4c: Cosmetic changes. Generating finale file
	echo "   > Generating Final file : 4C-seq library with repeats information"
	C4LibRmskFile="Library_${GenomeNickName}_${PrimaryRestrEnz}_${SecondaryRestrEnz}_${SegmentLength}bps_Rmsk_segmentsInfos.bed"
	echo "     finalFile = " $C4LibRmskFile
	ParseRmsk.py OverlapRMSK.txt $C4LibRmskFile
	echo ""; echo ""
	
fi


# -------------------------------------
#   Step 5:
#   clean-up
# -------------------------------------

mkdir -p ${OutputDestination}/4CLib/debug

mv $C4LibUnmaskedFile ${OutputDestination}/4CLib/$C4LibUnmaskedFile

mv ${GenomeFastaName}_uppercase.fa ${OutputDestination}/4CLib/debug/${GenomeFastaName}_uppercase.fa
mv segFile.txt ${OutputDestination}/4CLib/debug/segFile.txt
mv FragFile.txt ${OutputDestination}/4CLib/debug/FragFile.txt
mv logFile.txt ${OutputDestination}/4CLib/debug/logFile.txt
mv $(echo $var | perl -lne 'print $& if /fragmentBed\S*/') ${OutputDestination}/4CLib/debug/$(echo $var | perl -lne 'print $& if /fragmentBed\S*/')

if [ ! "$GenomeRmskTab" = "0" ]
then
	mv $C4LibRmskFile ${OutputDestination}/4CLib/$C4LibRmskFile
	mv $(basename "$GenomeRmskTab")_Bed.bed ${OutputDestination}/4CLib/debug/$(basename "$GenomeRmskTab")_Bed.bed
	mv OverlapRMSK.txt ${OutputDestination}/4CLib/debug/OverlapRMSK.txt
fi



if [ ! "$GenomeRmskTab" = "0" ]
then
	echo "############ ALL DONE : $(date +'%T') #############"
	echo "### the following files are the ${OutputDestination}/4CLib folder (TO KEEP)"
	echo "       - $C4LibUnmaskedFile   (Valid Fragments, no repeat information)"
	echo "       - $C4LibRmskFile   (Valid Fragments, with repeats information)"
	echo ""
	echo "### the following files are the in ${OutputDestination}/4CLib/debug folder [can be deleted (if not debugging)]"
	echo "       - ${GenomeFastaName}_uppercase.fa    (Uppercase Genome)"
	echo "       - segFile.txt / FragFile.txt / logFile.txt   (all restriction fragments, including not valid ones)"
	echo "       -" $(echo $var | perl -lne 'print $& if /fragmentBed\S*/') "      (List of Valid restriction Fragment in a Bed format)"
	echo "       - $(basename "$GenomeRmskTab")_Bed.bed    (bed file of RMSK table from UCSC)"
	echo "       - OverlapRMSK.txt    (Overlap between valid Restr fragments and RMSK bed file)"
else
	echo "############ ALL DONE : $(date +'%T') #############"
	echo "### the following files are the ${OutputDestination}/4CLib folder (TO KEEP)"
	echo "       - $C4LibUnmaskedFile   (Valid Fragments, no repeat information)"
	echo ""
	echo "### the following files are the in ${OutputDestination}/4CLib/debug folder [can be deleted (if not debugging)]"
	echo "       - ${GenomeFastaName}_uppercase.fa    (Uppercase Genome)"
	echo "       - segFile.txt / FragFile.txt / logFile.txt   (all restriction fragments, including not valid ones)"
	echo "       -" $(echo $var | perl -lne 'print $& if /fragmentBed\S*/') "      (List of Valid restriction Fragment in a Bed format)"
fi




