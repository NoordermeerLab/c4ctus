#! /usr/bin/env Rscript
Args <- commandArgs(TRUE)
#print("PlotDemultStat.R --args")
#print(Args)
infile <- Args[1]
NbrTotal <- as.numeric(Args[2])
Group <- Args[3]





library(ggplot2)
library(grid)




#NbrTotal = 5e7

#WordCountFile <- read.table(file="/Volumes/EQUIPES/CHRODY/Benoit/ExchangeWithCluster/MyAttemptHTSstation/Demultiplexing/ReadCount.mDxmyU06")
WordCountFile <- read.table(file=infile)
colnames(WordCountFile) <- c("file", "Viewpoint","ReadCount")

#Give NickName based on file name
WordCountFile$Name <- sapply(WordCountFile$file, function(x) {sub("Myunaligned_", "LowConfidence", x)})
WordCountFile$Viewpoint <- sapply(WordCountFile$Viewpoint, function(x) {sub("unaligned", "LowConfidence", x)})


WordCountFile$flag=1	#To keep or not for pie chart
WordCountFile$flag[WordCountFile$Viewpoint %in% c("discarded","LowConfidence","ambiguous_fastq")]=0
Unclassified = NbrTotal - sum(WordCountFile$ReadCount) + sum(WordCountFile$ReadCount[WordCountFile$flag == 0])
Unaligned    = NbrTotal - sum(WordCountFile$ReadCount)

WordCountFile <- rbind(WordCountFile, data.frame(file = ".", "Viewpoint" = "UNCLASSIFIED", "ReadCount" = Unclassified, "Name" = "UNCLASSIFIED", "flag" = 1))
WordCountFile <- rbind(WordCountFile, data.frame(file = ".", "Viewpoint" = "Unaligned", "ReadCount" = Unaligned, "Name" = paste0(Group,"_Unaligned"), "flag" = 0))




pdf(paste0(Group, "_DemutliplexingReport.pdf"), width = 8.26, height = 11.7, paper = "a4")

######################
# Generate Pie Graph
######################
PieData = subset(WordCountFile, flag == 1)
Lab <- paste(PieData$Name, paste(round(PieData$ReadCount/sum(PieData$ReadCount)*100,2), '%', sep=''), sep='\n')

layout(matrix(c(1,2), nrow = 2))
par(mar=c(0.5, 0.5, 1, 0.5))
pie(PieData$ReadCount, label = Lab, main = paste(Group, "   [Total = ", sum(PieData$ReadCount), "reads]"), cex=0.7)




######################
# Generate Bar Graph
######################
# Exclude the unclassified category
WordCountFile$flag[WordCountFile$flag == 0] = -1
WordCountFile$flag[WordCountFile$Name == "UNCLASSIFIED"] = 0

# sort
WordCountFile$Name <- factor(WordCountFile$Name,
       levels = WordCountFile$Name[order(WordCountFile$flag,WordCountFile$ReadCount)])

# create an apporpriate viewport.  Modify the dimensions and coordinates as needed
vp.Bottom <- viewport(height=unit(.5, "npc"), width=unit(1, "npc"), 
                           just=c("left","top"), 
                           y=0.5, x=0)

#ggplot
PP <- ggplot(subset(WordCountFile, flag != 0), aes(x=Name,y=ReadCount, fill = as.factor(flag))) + 
			geom_bar(stat="identity") + 
			coord_flip() +
			geom_text(aes(y = 0, label=paste("  ",format(ReadCount, big.mark=","),"reads")), vjust=0.5, hjust=0, color="black", size=3) + 
			scale_y_continuous(expand=c(0,0)) +
			theme_bw() +
			theme(axis.text.y = element_text(face = "bold", size = 8),
			      axis.title.y=element_blank(),
			      legend.position = "none") +
			ylab("Read count")
print(PP, vp=vp.Bottom)


garbage <- dev.off()
