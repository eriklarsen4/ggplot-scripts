

## Script for obtaining Enrichr-based bioinformatics data
## Developed by Erik Larsen


##### Environment Prep #####
## Load previously installed packages from the library (if not previously installed, install!)

library(ggplot2)
library(tidyverse)
library(stringr)
library(plyr)
library(ggrepel)
#library(calibrate)
library(reshape2)
#library(gridExtra)
library(stats)
library(dplyr)
library(readr)
library(readxl)
library(tibble)

## Install packages to use enrichr, panther packages
## At least for the first installation, use the following code; THIS REQUIRES R version 4.0.0+ !!!

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")

#install.packages("BiocGenerics")
library(AnnotationDbi)
library(stats4)
library(BiocGenerics)
library(parallel)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("PANTHER.db")

## Load the Panther 2016 Pathway Database package and use "mouse" as reference genome/specified species for downstream pathway analysis
library(PANTHER.db)
pthOrganisms(PANTHER.db) = "MOUSE"

## Load the enrichr package
library(enrichR)
## Connect live to the enrichr server/master database (website) and store the "space" as a variable. This contains a vast array of bioinformatic databases/websites
## Store it as a variable for better visualization and eventual subsetting
DBs = listEnrichrDbs()

## View the variable in the editor to find the relevant databases and their indeces for subsetting (click on the dataframe in the global environment and manually peruse)

  ## In this case, "GO_Biological_Process_2018" (index # 130), and "Panther_2016" (index # 102)
    ## Concatenate them in a list for subsetting, or slice directly; we'll slice directly in the next line
DBs = DBs$libraryName[c(130,102)]

##### Data prep #####
## Import the dataset you want to analyze
DESeq2_Expression_Results = read.csv("M:/Erik/Data/Omics/TimeCourse/Processed Galaxy Output/Test Results to Upload/DESeq2_result_file_on_GT_e13_PARA.csv")
## Filter (subset) genes that went undetected or were outliers in terms of counts; new dataframe should not contain any NAs in p-value columns
DESeq2_Expression_Results3 = subset(DESeq2_Expression_Results, (!is.na(DESeq2_Expression_Results[,"AdjP"])))


## Perform relevant analysis through enrichr; in this case, FDR <= 0.01
  ## first: GO Bio Processes; this analyzes genes from your dataset in enrichr, which draws terms from geneontology.org
GOs = as.data.frame(
  enrichr(
    c(DESeq2_Expression_Results3$GeneID[which(DESeq2_Expression_Results3$AdjP <= 0.01)]), DBs[1])
  )

  ## second: Panther Pathways; this analyzes genes from your dataset in Enrichr, which draws terms (pathways) from pantherdb.org
PATHWAYS = as.data.frame(
  enrichr(
    c(DESeq2_Expression_Results3$GeneID[which(DESeq2_Expression_Results3$AdjP <= 0.01)]), DBs[2])
  )

## Clean up both dataframes and export both dataframes.
  ## The pathway analysis dataframe needs additional cleaning in Python;
    ## Following that manipulation, data will be re-imported in other scripts for visualization using ggplot

## Remove unnecessary columns
GOs = GOs[,-c(5,6,7)]
PATHWAYS = PATHWAYS[,-c(5,6,7)]

## Export data as .CSVs to the server
#write.table(GOs, "M:\\PAPER ASSEMBLY\\Itch paper\\Timecourse RNAseq\\e13 Mut GO Bio Processes.csv", sep = ",", col.names = T)
#write.table(PATHWAYS, "M:\\PAPER ASSEMBLY\\Itch paper\\Timecourse RNAseq\\e13 Mut Panther Pathways.csv", sep = ",", col.names = T)

##### Volcano Plot of e13 GT DRG #####

## Create a column in the DESeq2 dataframe that scales the Adjusted P-value by log-base 10
DESeq2_Expression_Results3$log10Pval = -log10(DESeq2_Expression_Results3$AdjP)


## Create the data points of interest for functional display.. ggplot struggles to plot subsets of dataframes, so creating a variable that directly labels all the points accordingly helps.
  ## Create a new column (points of interest); Using "GeneID" is irrelevant; the values will be replaced in the next step
DESeq2_Expression_Results3$p.o.i. = DESeq2_Expression_Results3$GeneID

## Fill the points with appropriately indexed data;
  ## "Differentially Expressed Genes" defined as being "points of interest" with Adjusted P-Values of < 0.05
DESeq2_Expression_Results3$p.o.i.[which(DESeq2_Expression_Results3$AdjP < 0.01)]= "Differentially Expressed Genes"

  ## "Non-differentially Expressed Genes" defined as being "points of interest" with Adjusted P-values > 0.05
DESeq2_Expression_Results3$p.o.i.[which(DESeq2_Expression_Results3$AdjP > 0.01)]= "Non-differentially Expressed Genes"


## Create labels on points we're interested in labeling in volcano plots
## Create the column, "labs"; hide all of the text labels with: ""
DESeq2_Expression_Results3$labs = ""


## Label only these selected items; adjust for each experiment; for the current example and a more global volcano plot, highlight the most extreme data points on the plot
ix_label1 = c(DESeq2_Expression_Results3$GeneID[DESeq2_Expression_Results3$log10Pval > 28])

  ## Fill the appropriate rows of the "labs" column with the corresponding Gene names for labeling those points
DESeq2_Expression_Results3$labs[DESeq2_Expression_Results3$log10Pval > 28] = DESeq2_Expression_Results3$GeneID[DESeq2_Expression_Results3$log10Pval > 28]

    ## Notes on ggplot commands:
    ## geom-point = alpha is a ggplot add-on function that controls point color transparency
    ## coord_cartesian is a ggplot add-on function that determines axes sizes
    ## labs is a ggplot add-on function that determines the plot labels; x-axis, y-axis, title, legend colors
    ## scale_color_manual is a ggplot add-on function that sets the colors of the plot, corresponding to the groups determined in the the base function ("col" in aes)
    ## Theme is used to manipulate the physical location of the plot title

    ## Sometimes, a column from a dataframe will contain apostrophes around the column name. Pay attention to this, you may need to revise depending on what R version you are running..
vol1 = ggplot(data = DESeq2_Expression_Results3, aes(x = DESeq2_Expression_Results3$`log2(FC)`, y = DESeq2_Expression_Results3$log10Pval, col = DESeq2_Expression_Results3$p.o.i.)) + geom_point(alpha = 0.2) + coord_cartesian(xlim = c(-7.9,7.6), ylim = c(0, 123)) + labs(x = "log2(FC) in mRNA Exp. (Relative to WT)", y = "-log10(Adj. P-value)", title = "e13 Tmem184b-GT DRG Gene Exp. Changes", col = "Gene Data Type") + scale_color_manual(values = c("blue", "gray")) + theme(plot.title = element_text(hjust = 0.5)) + theme_light()


    ## geom_text_repel labels the selected points in the color commanded
## Re-save the variable with the added functionality; for unlabeled plot
vol1 = vol1  + geom_text_repel(color = "black", aes(label=DESeq2_Expression_Results3$labs)) 

## Plot it
vol1



##### Bar Plots of GO Terms Using ggplot2 #####
## Filter some of the Processes out, including the top 19 Processes
GO_Bar_ADJP = GOs[1:19,]

## Plot the bar plot
## Store as a variable
## Order the Processes by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
## scale_fill_continuous() creates the color gradient and legend
## coord_flip() switches axes for better visualization of the long GO terms
BAR_ADJP = ggplot(GO_Bar_ADJP, aes(x = reorder(GO_Bar_ADJP$GO_Biological_Process_2018.Term, -GO_Bar_ADJP$GO_Biological_Process_2018.Adjusted.P.value), y = GO_Bar_ADJP$GO_Biological_Process_2018.Adjusted.P.value)) + geom_col(stat = "identity" , aes(fill = GO_Bar_ADJP$GO_Biological_Process_2018.Adjusted.P.value)) + scale_fill_continuous(name = "Enrichr Adj.P-val") + labs(x = "GO Biological Process", y = "Enrichr (BHM) Adjusted P-value", title = "e13 Tmem184b-GT GO Analysis") + theme_light() + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5))  + coord_flip()


## Plot it.
BAR_ADJP




##### Bar Plots of Panther Pathways Using ggplot2 #####

## Import data you want to visualize
Pathway_Bar = read.csv("M:/PAPER ASSEMBLY/Itch paper/Timecourse RNAseq/e13 Mut Panther Pathways.csv", as.is = TRUE)

## Sort and Filter some of the Pathways out by AdjP val
Pathway_Bar = arrange(df = Pathway_Bar, Pathway_Bar$Panther_2016.Adjusted.P.value)
## Remove all but the top 19 pathways
Pathway_Bar = Pathway_Bar[1:19,]
## Remove the " Homo sapiens.." pathway name strings in the Panther_2016.Term column
Pathway_Bar$Panther_2016.Term = Pathway_Bar$Panther_2016.Term %>% gsub(x = Pathway_Bar$Panther_2016.Term, pattern = " Homo sapiens .+.?$", replacement = "")

## Plot the bar plot
## Store as a variable
## Order the pathways by smallest Adjusted P-value to largest, simultaneously "mapping" the values of each term to a color (re-order the x variable, y variable within the ggplot generic function).
## Adding the geom_col() enables the plot to look like a bar plot, and filling the aesthetics of this graphing function.. completes the color mapping
## scale_fill_continuous() creates the color gradient and legend
## coord_flip() switches axes for better visualization of the long pathway names
PATHWAY_BAR_ADJP = ggplot(Pathway_Bar, aes(x = reorder(Panther_2016.Term, -Panther_2016.Adjusted.P.value), y = Panther_2016.Adjusted.P.value)) + geom_col(stat = "identity" , aes(fill = Panther_2016.Adjusted.P.value)) + scale_fill_continuous(name = "Enrichr Adj.P-val") + labs(x = "Panther 2016 Pathways", y = "Enrichr (BHM) Adjusted P-value", title = "e13 Tmem184b-GT Pathway Analysis") + theme_light() + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), plot.title = element_text(hjust = 0.5))  + coord_flip()

## Plot it.
PATHWAY_BAR_ADJP



##### e13 Panther Pathways Heatmap with labels using ggplot #####

## Import the relevant data
e13_DEGs_and_Pathways =  read_csv("M:/Omics/Custom Python Enrichr Clustergram e13 GT.csv")

## Export the gene list for photoshopping
#write.csv(e13_DEGs_and_Pathways$DEGs, "M:\\PAPER ASSEMBLY\\Itch paper\\Figure Drafts\\e13 FDR 01 Heatmap DEGs.csv", sep = ",")


## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
e13_Heatmap = melt(e13_DEGs_and_Pathways, id = "DEGs")
colnames(e13_Heatmap) = c("DEGs", "Pathways", "Presence In Pathway")

## Remove excessive pathway name strings
Pathway_Names = e13_Heatmap$Pathways %>% gsub(x = e13_Heatmap$Pathways, pattern = " Homo sapiens .+.?", replacement = "")

## Put the new strings back into the dataframe
e13_Heatmap[,2] = Pathway_Names

## Make the digital values factors so it's compatible to graph
e13_Heatmap$`Presence In Pathway` = as.factor(e13_Heatmap$`Presence In Pathway`)



## Create and store the heatmap core as a variable
## Use geom_tile to map the binary values
## Customize by removing the filler surrounding the graph and the tickmarks
## Center the Title

HEATlabs = ggplot(e13_Heatmap, aes(e13_Heatmap$Pathways, e13_Heatmap$DEGs)) + geom_tile(aes(fill = e13_Heatmap$`Presence In Pathway`), color = "black") + scale_fill_manual(values = c("dark blue", "gold")) + labs(x = "Panther 2016 Pathways", y = "Differentially Expressed Genes (FDR < 0.01)", title = "Heatmap of Tmem184b-GT-affected Pathways") + guides(fill = guide_legend(title = "Presence in Pathway")) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0))


## Plot the heatmap
HEATlabs

##### e13 Panther Pathways Heatmap without labels using ggplot #####

## Without labels
HEATnolabs = ggplot(e13_Heatmap, aes(e13_Heatmap$Pathways, e13_Heatmap$DEGs)) + geom_tile(aes(fill = e13_Heatmap$`Presence In Pathway`), color = "black") + scale_fill_manual(values = c("dark blue", "gold")) + labs(x = "Panther 2016 Pathways", y = "Differentially Expressed Genes (FDR < 0.01)", title = "Heatmap of Tmem184b-GT-affected Pathways") + guides(fill = guide_legend(title = "Presence in Pathway")) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), panel.grid.major = element_line(colour = "black", size = 0.1), panel.grid.minor = element_line(colour = "black", size = 0.1), plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0), axis.text= element_blank())

## Plot the heatmap
HEATnolabs