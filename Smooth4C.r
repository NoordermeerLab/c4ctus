#! /usr/bin/env Rscript
Args <- commandArgs(TRUE)
infile      <- Args[1]
Exclude     <- Args[2]
RollWindow  <- as.numeric(Args[3])
outPrefix   <- Args[4]



library(stringr)
library(RSQLite)
suppressMessages(library(zoo))




###########################
#  FUNCTION: Load Raw 4C
###########################
LoadRaw4C <- function(CurrFile, chr) {
	#print(paste("Loading File 1 : ", CurrFile))
	
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


###############################
# FUNCTION: smooth Raw 4C data
###############################

SmoothRaw4C <- function(df, RollMeanSize, ExcludeLeft, ExcludeRight) {
	
	# Compute RollMean
	df$scoreRollMean = round(rollmean(df$score, RollMeanSize, fill = 0), digits = 3)
	
	# Exclude region close to the VP
	Idx = df$end > ExcludeLeft & df$start < ExcludeRight
	df$score[Idx] = 0
	df$scoreRollMean[Idx] = 0
	
	return(df)
}



#infile = "4Cseq/segToFrag_MyS2_50M_Ubx_dm6_rep57751_all.sql"
#Exclude <- "chr3R:16730708-16737272"
#RollWindow = 11

#Parse Exclude (CHR):(ExcludeLeft)-(ExcludeRight)
res <- str_match(Exclude, "^(chr.+):([:digit:]+)-([:digit:]+)$")
CHR          <- res[,2]
ExcludeLeft  <- as.numeric(res[,3])
ExcludeRight <- as.numeric(res[,4])

#Load Raw 4C of CHR and smooth it
Data <- LoadRaw4C(infile, CHR)
Data <- SmoothRaw4C(Data, RollWindow, ExcludeLeft, ExcludeRight)

#write to file
write.table(Data[,c("chr","start","end","scoreRollMean")],
            file = paste0(outPrefix,"_",CHR,"_smoothed_",RollWindow,"FragPerWin.bedGraph"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE, col.names = FALSE)