#! /usr/bin/env Rscript
Args <- commandArgs(TRUE)

infile       <- Args[1]
Exclude      <- Args[2]
RollWindow   <- as.numeric(Args[3])
NormLeft     <- as.numeric(Args[4])
NormRight    <- as.numeric(Args[5])
OutFilePrefix <- Args[6]

#
# USAGE: Norm4C.r segToFrag_MyS2_50M_Abd-A_dm6_rep57795_all.sql chr3R:16827223-16833379 11 16.7e6 17e6 MyS2_50M_Abd-A_dm6_rep57795
#    Output file will be <OutFilePrefix>_<CHR>_smoothed_11FragPerWin_Norm.bedGraph



#####  CHANGE IF NECESSARY  #####
NormTo = 1e6
#################################


library(stringr)
library(RSQLite)
suppressMessages(library(zoo))




###########################
#  FUNCTION: Load Raw 4C
###########################
LoadRaw4C <- function(CurrFile, chr) {
	print(paste("Loading File 1 : ", CurrFile))
	
	#Connect to database
	data.Raw4C = data.frame(chr=character(),start=integer(),end=integer(),score=double())
	con = dbConnect(SQLite(), dbname=CurrFile)
	
	#Retrieve data
	myQuery <- dbSendQuery(con, paste("SELECT * FROM ",chr," ORDER BY start"))
	data.Raw4C = dbFetch(myQuery, n = -1)
	
	#disconnect
	bool <- dbClearResult(myQuery)
	dbDisconnect(con)
	
	#polish
	data.Raw4C$chr = chr
	data.Raw4C = data.Raw4C[order(data.Raw4C$start),]
	col_order <- c("chr", "start", "end", "score")
	data.Raw4C <- data.Raw4C[, col_order]
}


###########################
# FUNCTION: Normalize and 
#      smooth Raw 4C data
###########################
NormSmoothRaw4C <- function(df, ExcludeLeft, ExcludeRight, NormLeft, NormRight, RollMeanSize) {
	
	# Compute Rolling mean
	df$scoreRollMean = round(rollmean(df$score, RollMeanSize, fill = 0), digits = 2)
	
	# Exclude region close to the VP
	Idx = df$end > ExcludeLeft & df$start < ExcludeRight
	df$score[Idx] = 0
	df$scoreRollMean[Idx] = 0
	
	# Normalize FILE 1 data
	R = NormTo / sum(subset(df, end > NormLeft & start < NormRight)$score)
	df$scoreNormSmooth = round(df$scoreRollMean * R, digits = 2)
	
	DD = subset(df, end > NormLeft & start < NormRight)
	return(DD)
}




#Parse Exclude (CHR):(ExcludeLeft)-(ExcludeRight)
res <- str_match(Exclude, "^(chr.+):([:digit:]+)-([:digit:]+)$")
CHR          <- res[,2]
ExcludeLeft  <- as.numeric(res[,3])
ExcludeRight <- as.numeric(res[,4])

#Load Raw 4C of CHR, Normalize and Smooth it
Data <- LoadRaw4C(infile, CHR)
Data <- NormSmoothRaw4C(Data, ExcludeLeft, ExcludeRight, NormLeft, NormRight, RollWindow)

#write to file
write.table(Data[,c("chr","start","end","scoreNormSmooth")],
            file = paste0(OutFilePrefix,"_",CHR,"_smoothed_",RollWindow,"FragPerWin_Norm.bedGraph"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE, col.names = FALSE)