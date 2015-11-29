#This script has been used to generate the tSNE-plots for figure 8g in the paper "The
#heterogeneity of human CD127+ innate lymphoid cells revealed by single-cell RNA
#sequencing" by Björklund, Forkel et al. 

# All data needed for the analysis is included in FACS_raw_protein_data_ILC.zip

# Copyright (c) 2015, Jakob Theorell <jakob.theorell@ki.se> for this script. For citation
#purposes, the user is referred to The Core Team of R, as well as the authors for every
#included package. 
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

#Packages needed for preprocessing
library(plyr)

#Packages needed for Barnes-Hut SNE
library(Rtsne)

#Packages needed for graphic display/density estimates
library(gplots) 
library(ggplot2)
library(RColorBrewer)
library(grid)
library(rgl)
library(MASS)

#Create a new folder if it is not already present and change the working directory there (this will be iterated multiple times below)
setwd("~/Desktop/")
subDir <- "SNE_mapper"
if (! file.exists(subDir)){
  dir.create(subDir)
}
setwd("~/Desktop/SNE_mapper/")

subDir <- "Data"
if (! file.exists(subDir)){
  dir.create(subDir)
}
setwd("~/Desktop/SNE_mapper/Data/")

#Move raw files into this directory
file.rename("~/Desktop/Raw_data_ILC", "~/Desktop/SNE_mapper/Data/Raw_data_ILC")

##########################

#Import the raw data. For three of the files, a warning is generated, but it does not have any impact on the data. 
setwd("~/Desktop/SNE_mapper/Data/Raw_data_ILC/")
filenames <- list.files(path = ".", pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
 read_csv_filename <- function(filename){
    ret <- read.csv(filename)
    id <- gsub(pattern="ILC_..._......._|\\.csv", "\\1", filename)
    ret$id <- id
    label <- gsub(pattern="ILC_..._|_.......\\.csv", "\\1", filename)
    ret$label <- label
        
    ret
 }  
raw_data <- ldply(filenames, read_csv_filename)

setwd("~/Desktop/SNE_mapper/Data/") 
write.csv(raw_data, "all_data_before_raw.csv", row.names=FALSE)

colnames(raw_data)

colnames(raw_data) <- c("CD45RA", "HLA.DR", "CD62L", "NKp44", "ids", "label")

#############################################################


#Creates a data frame with the full-range indata-parameters.
indata_full <- subset(raw_data, select=c("CD45RA", "HLA.DR", "CD62L", "NKp44"))

#Display the column names
colnames(indata_full)


#Turn each column of the data frame into percentages of the range between the highest and the lowest values
indata_percent <- 100*((indata_full-apply(indata_full, 2, min, na.rm=TRUE)[col(indata_full)])/(apply(indata_full, 2, max, na.rm=TRUE)[col(indata_full)]-apply(indata_full, 2, min, na.rm=TRUE)[col(indata_full)]))

setwd("~/Desktop/SNE_mapper/Data")
write.csv(indata_percent, "All_indata.csv", row.names=FALSE)

#############################################################
#Generate the tSNE parameters

twodimSNE <- Rtsne(as.matrix(indata_percent), perplexity=30, theta=0.5)

write.csv(twodimSNE, "Result.csv", row.names=FALSE)

twodimSNEY <- twodimSNE$Y

colnames(twodimSNEY) <- c("V1", "V2")

#Turn each column of the data frame into percentages of the range between the highest and the lowest values
percent_result <- 100*((twodimSNEY-apply(twodimSNEY, 2, min, na.rm=TRUE)[col(twodimSNEY)])/(apply(twodimSNEY, 2, max, na.rm=TRUE)[col(twodimSNEY)]-apply(twodimSNEY, 2, min, na.rm=TRUE)[col(twodimSNEY)]))

#############################################################

#Sort the variables into categorical and quantitative variables, i e
quantitative <- raw_data[ , -which(colnames(raw_data) %in% c("ids", "label"))]
categorical <- raw_data[ , which(colnames(raw_data) %in% c("ids", "label"))]


#4 Normalize the most extreme low and high per mille of the events to their lower end, so they will not influence the graphical display. This does not have any effect of the SNE parameter generation. 
for (i in 1:length(quantitative)){    
    cat(i,"\n")    

    x <- quantile(quantitative[[i]], c(.999))
    y <- quantile(quantitative[[i]], c(.001))    

    print(x)
    print(y)    

    quantitative[[i]][quantitative[[i]] > x] <- x
    quantitative[[i]][quantitative[[i]] < y] <- y    
}



#Turn each column of the data frame into percentages of the range between the highest and the lowest values
percent_quantitative <- 100*((quantitative-apply(quantitative, 2, min, na.rm=TRUE)[col(quantitative)])/(apply(quantitative, 2, max, na.rm=TRUE)[col(quantitative)]-apply(quantitative, 2, min, na.rm=TRUE)[col(quantitative)]))

Percent_all_data <- cbind(percent_quantitative, percent_result, categorical)

#Export this
setwd("~/Desktop/SNE_mapper/Data")
write.csv(Percent_all_data, "Percent_all_data.csv", row.names=FALSE)

#Stupid to solve problem with id not being values, but a string of some other kind. 
Percent_all_data <- read.csv("Percent_all_data.csv")
#############################################################

#Plotting of all the flow cytometry variables

setwd("~/Desktop/SNE_mapper/")
subDir <- "Graphics"
if (! file.exists(subDir)){
  dir.create(subDir)
}
setwd("~/Desktop/SNE_mapper/Graphics/")


#Save the palette. This is necessary for a correct color representation below. 
pdf("palette.pdf")
palette(rev(rich.colors(100, plot=TRUE)))
dev.off()

#The label parameter has been given values that corresponds to nice shapes in ggplot2. Therefore: 
cols <- Percent_all_data$ids

#Make contour plot of the subpopulations
fname3 = paste("Subpop_combined",'.pdf')
p <- ggplot(Percent_all_data,aes(x=V1,y=V2)) +
    	geom_point(col=cols, size=2.2, alpha=0.7) 	+  xlim(-10, 110) + ylim(-10, 110) +
    	theme (line = element_blank(),
        		text = element_blank(),
        		line = element_blank(),
        		title = element_blank(),
				panel.background = element_rect(fill = "white"))
gt <- ggplot_gtable(ggplot_build(p))
ge <- subset(gt$layout, name == "panel")
grid.draw(gt[ge$t:ge$b, ge$l:ge$r])
ggsave(filename = fname3, dpi=300)


##################################################################################

#Figures for all the fluorescent parameters are sorted into a specific folder
setwd("~/Desktop/SNE_mapper/Graphics")
subDir <- "Indata_params"
if (! file.exists(subDir)){
  dir.create(subDir)
}
setwd("~/Desktop/SNE_mapper/Graphics/Indata_params")

Indata_plotting <- Percent_all_data[ , which(colnames(Percent_all_data) %in% c(colnames(indata_percent)))]

#Make contour plot of all the populations with color representing each csv column. 
for (i in 1:length(Indata_plotting)){   
    cat(i,"\n")
    graphics.off()
    
    #graphs
    fname3 = paste(colnames(Indata_plotting[i]),'.pdf')
    
    p <- ggplot(Percent_all_data,aes(x=V1,y=V2)) +
    	geom_point(col=(1 + 0.98*(102-Indata_plotting[[i]])), size=2.5, alpha=0.7) + xlim(-10, 110) + ylim(-10, 110) +
    	theme (line = element_blank(),
        		text = element_blank(),
        		line = element_blank(),
        		title = element_blank(),
				panel.background = element_rect(fill = "white"))
	gt <- ggplot_gtable(ggplot_build(p))
	ge <- subset(gt$layout, name == "panel")

	grid.draw(gt[ge$t:ge$b, ge$l:ge$r])

    #Save file
    ggsave(filename = fname3, dpi=300)
}

###########################################

#Figures for the label plot is sorted separately
setwd("~/Desktop/SNE_mapper/Graphics")
subDir <- "Donor_distribution"
if (! file.exists(subDir)){
  dir.create(subDir)
}
setwd("~/Desktop/SNE_mapper/Graphics/Donor_distribution")

cols <- Percent_all_data$label

#Make contour plot of the label parameter
fname3 = paste("Donors_adults",'.pdf')
p <- ggplot(Percent_all_data,aes(x=V1,y=V2)) +
    	geom_point(col=cols, size=2.5, alpha=0.5) 	+  xlim(-10, 110) + ylim(-10, 110) +
    	theme (line = element_blank(),
        		text = element_blank(),
        		line = element_blank(),
        		title = element_blank(),
				panel.background = element_rect(fill = "white"))
gt <- ggplot_gtable(ggplot_build(p))
ge <- subset(gt$layout, name == "panel")
grid.draw(gt[ge$t:ge$b, ge$l:ge$r])
ggsave(filename = fname3, dpi=300)
