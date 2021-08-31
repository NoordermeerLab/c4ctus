#!/bin/sh


##################################
#-- Parse Arguments
##################################
FileToProcess=$1
outNamePrefix=$2
CurrChr=$3
outDir=$4

mkdir -p ${outDir}
cd ${outDir}

##################################
#-- Run Instruction
##################################
echo "Domainogram: Processing ${FileToProcess} for ${CurrChr}"
echo "Output files will have the following prefix: ${outNamePrefix}"
echo "      -- Start: $(date +'%T') --"
echo "***************";echo""


runDomainogram.R ../${FileToProcess} \
                 ${outNamePrefix}_${CurrChr} \
                 ${outNamePrefix}_${CurrChr} \
                 ${CurrChr} \
                 500 \
                 50 \
                 0 \
                 /store/EQUIPES/CHRODY/Benoit/Lab/c4ctus/c4ctus_bin



##################################
#-- Clean-up
##################################

#Bring Log file and RData in the Domainogram folder
mv ../${outNamePrefix}_${CurrChr}_domainograms.RData ${outNamePrefix}_${CurrChr}_domainograms.RData

echo ""
echo "############ ALL DONE : $(date +'%T') #############"