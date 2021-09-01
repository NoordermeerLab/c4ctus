#! /usr/bin/env Rscript
Args <- commandArgs(TRUE)

if (length(Args) < 2) {
	cat(paste("RePlotDomainogram.R <domainograms.RData file> <OutputFile>",
	            "The output file type depends in the extension of OutputFile",
	            "it can be jpg / png / pdf",
	            "\n",
	            sep="\n"))
	quit()
}

inRDataFile <- Args[1]
OutFile     <- Args[2]

#inRDataFile = "MyS2_50M_Abd-B-P1a_dm6_rep57798_chr3R_domainograms.RData"
#OutFile = "ttututu.jpg"

library(plotrix)
library(reshape2)
library(ggplot2)
suppressMessages(library(scales))
suppressMessages(library(fields))



#####   README   ######
# HTSstation used to plot the full domainogram in an extremely large PDF
# This script aim to replot the domainogram using ggplot and ggsave to export
#
# This is very experimental for two raisons
#  - 1: Plotting the full resDomainogram image using ggplot works on my MAC
#       but it crashes on the server, for a totally unclear raison
#  - 2: For this reason, I chose to resize the resDomainogram image (bilinear interpolation)
#       I haven't checked if everything is correct
#  - 3: The x-axis is tricky. My FEELING is that HTSstation didn't bother and used regularly 
#       spaced coordinates despite the fact that each bin correspond to a valid restriction
#       fragment, and thus bins are not linear with genomic distances.
#       I therefore tried to also do a linear interpolation of the x-axis.
#       Keep in mind that I could be wrong there.
#       Benoit Moindrot (08/03/2019)
#
# USAGE: RePlotDomainogram.R <domainograms.RData file> <OutputFile>
#        The output file type depends in the extension of OutputFile
#        it can be jpg / png / pdf




#############################################
###  MATRIX RESIZING FUNCTION
#############################################
#source: https://stackoverflow.com/questions/11123152/function-for-resizing-matrices-in-r


rescale <- function(x, newrange=range(x)){
  xrange <- range(x)
  mfac <- (newrange[2]-newrange[1])/(xrange[2]-xrange[1])
  newrange[1]+(x-xrange[1])*mfac
}

ResizeMat <- function(mat, ndim=dim(mat)){

  # input object
  odim <- dim(mat)
  obj <- list(x= 1:odim[1], y=1:odim[2], z= mat)

  # output object
  ans <- matrix(NA, nrow=ndim[1], ncol=ndim[2])
  ndim <- dim(ans)

  # rescaling
  ncord <- as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
  loc <- ncord
  loc[,1] = rescale(ncord[,1], c(1,odim[1]))
  loc[,2] = rescale(ncord[,2], c(1,odim[2]))

  # interpolation
  ans[ncord] <- interp.surface(obj, loc)

  ans
}



#############################################
###  ACUTAL WORK
#############################################

print(paste("Loading data from :", inRDataFile))
load(inRDataFile)

#clipping data
resDomainogram[resDomainogram<1e-6]=1e-6



#Define Colormap and transformation function
nGrad = 10
myCols <- smoothColors("black",(nGrad-2)/2,"purple",(nGrad-2)/2,"red")
MinusLog10 <- trans_new("-log10", function(x) -log10(x), function(x) 10^(-x))


print("Converting data")


###  -- If you want to replot the full domainogram, uncomment the lines below --- ####

#tdm <- melt(resDomainogram)
##tdm <- tdm[c(1:1e6),]		#To troubleshoot ggplot

##HTSstation stategy for the x-axis (no interpolation)
#lensc = length(dataToTreat$V1)
#I=seq(1,lensc,by=2000)
#ticks=as.integer(0.5+(as.numeric(dataToTreat[I,2])+as.numeric(dataToTreat[I,3]))*5e-7)

#ll <- ggplot(tdm, aes(x = Var2, y = Var1, fill = value)) + geom_raster() +
#		scale_fill_gradientn(colours=myCols,
#		                     na.value="white",
#		                     trans=MinusLog10,
#		                     breaks=c(1 %o% 10^seq(-8,0,by=2))) +
#		scale_x_continuous(expand=c(0,0),
#		                   limits=c(0,lensc),
#		                   breaks=I,
#		                   label =ticks) +
#		scale_y_continuous(expand=c(0,0)) +
#		xlab("Position in Mb") + ylab("size") +
#		theme_light() +
#		theme(legend.position="top",legend.key.width = unit(0.1, units="npc"))



#############################################
# bilinear interpolation of resDomainogram
#############################################
PixelWidth = 5000
ss <- ResizeMat(resDomainogram, c(500,PixelWidth))
tdm2 <- melt(ss)

CoordMax <- round(max(dataToTreat$V3)/1e6)
CoordMin <- round(min(dataToTreat$V2)/1e6)
I = seq(CoordMin*1e6,CoordMax*1e6,by=2e6)

dataToTreat$mid = (dataToTreat$V3 + dataToTreat$V3)/2
#Linear Approx
RealBreaks = approx(dataToTreat$mid, 1:nrow(dataToTreat), I)$y
#Rescale RealBreaks between 0 and PixelWidth (instead of between 0 and length(dataToTreat$V3))
RescaledBreaks = RealBreaks * PixelWidth/length(dataToTreat$V3)


ll <- ggplot(tdm2, aes(x = Var2, y = Var1, fill = value)) + geom_raster() +
		scale_fill_gradientn(colours=myCols,
		                     na.value="white",
		                     trans=MinusLog10,
		                     breaks=c(1 %o% 10^seq(-8,0,by=2))) +
		scale_x_continuous(expand=c(0,0),
		                   limits=c(0,5000),
		                   breaks=RescaledBreaks,
		                   label =I / 1e6) +
		scale_y_continuous(expand=c(0,0)) +
		xlab("Position in Mb") + ylab("size") +
		theme_light() +
		theme(legend.position="top",legend.key.width = unit(0.1, units="npc"),
		      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


print(paste("Saving data plot in :",OutFile))
print("   This may take some times...")
ggsave(OutFile, plot = ll, width = 29.7, height = 21, units ="cm", dpi = 300)
print("All done")
