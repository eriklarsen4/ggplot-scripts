

## Developed by Erik Larsen #

## Upload the packages that include ggplot2, which is a set of functions and commands that produce images, "Grammar of Graphics", by mapping data to aesthetic properties;

## plyr and dplyr help with some functions (copied and pasted from online; can't remember which);
## ggrepel is necessary for neatly labeling data points;
## reshape2 is necessary for re-formatting data ("melt" function that stacks groups of variables)
## stats enables statistical functionality
## gridExtra doesn't do shit yet
## readr enables importing excel and text files into your Environment to then manipulate; also enables exporting

## Install the tidyverse and reshape2 packages onto your drive if not already done
#install.packages("tidyverse")
#install.packages("reshape2")

## Upload ("call") the above-mentioned functions into your R session. Once installed, they are stored in your R "library" of packages
library(ggplot2)
library(plyr)
library(ggrepel)
library(calibrate)
library(reshape2)
library(gridExtra)
library(stats)
library(dplyr)
library(readr)


## Import the data if not already done ("Import Dataset" tab in the Environment Window); it is directly called via the script below. These are the DESeq2 results from Galaxy, computing which genes are differentially expressed

## If working with your own data files from Galaxy, save the DESeq2 Output (downloaded from Galaxy as a .txt file) as a CSV (with headings!) to import into the environment (top right window panel) to then manipulate and eventually graph with the "ggplot" function

## Also, in Excel, add a column that converts log2 FC to % WT Expression if not already done (copy and paste the formula "POWER(2,"_value_from_log2FC_cell")" for the entire column. Save and import.


#DESeq2_Expression_Results = read_csv("M:/Omics/Hippocampus/Age Comparisons/Processed Galaxy Output/Test Results to Upload/DESeq2 Expression Results for Age Mutant Comparison without omu1849 and ymu4141.csv")

DESeq2_Expression_Results = read_csv("M:/Omics/TimeCourse/Processed Galaxy Output/Test Results to Upload/DESeq2_result_file_on_GT_e13_PARA.csv")


## Process the relevant data, removing the transcripts that lack an Adjusted P-value

## Use this new data frame for manipulation throughout the script for graphics and statistics
DESeq2_Expression_Results3 = subset(DESeq2_Expression_Results, (!is.na(DESeq2_Expression_Results[,"AdjP"])))


#### How many genes are significantly different from WT expression? ####
padj.10 = which(DESeq2_Expression_Results3$AdjP <= 0.10)
padj.05 = which(DESeq2_Expression_Results3$AdjP <= 0.05)
padj.01 = which(DESeq2_Expression_Results3$AdjP <= 0.01)


##### Take a subset the dataframe to select differentially expressed genes with Adjusted P-values < 0.05 for use in pathway analysis to find potential mechanistic pathways or other genes/interactions of interest #####
HITS.01 = subset(DESeq2_Expression_Results3, DESeq2_Expression_Results3$AdjP <= 0.01)
HITS.05 = subset(DESeq2_Expression_Results3, DESeq2_Expression_Results3$AdjP <= 0.05)
HITS.10 = subset(DESeq2_Expression_Results3, DESeq2_Expression_Results3$AdjP <= 0.10)

## You can save new dataframes/variables/vectors as CSVs and export them to the server (edit path if the server is not accessible!).
## This is useful for analysis in other bioinformatics programs online that will require you to attach a list or copy and paste an entire dataset.

## To do this, you'll need to know where you're telling R to send the exported file.

## As an example, I've included a variable I created in another script. The variable is called, "Candidate Screen" and it contained just a list of     differentially expressed genes.
##      Use these following lines as a template, editing the destination (i.e. change the entire destination; the following destination is to my folder in the server.  Server -> "Erik" folder -> "Data" folder within "Erik" folder, etc. Then include the title of the file with the file extension (".csv")

## Create a vector of the DEGs (FDR < 0.01) for input into PantherDatabase.org for Pathway Analysis
#e13_Candidate_Screen = HITS.01$GeneID
#write.table(e13_Candidate_Screen, "M:\\Omics\\TimeCourse\\Processed Results\\e13 Candidate Screen.csv", sep=",",row.names=F, col.names = F)


##### Volcano Plot Pre-processing #####


## Create a column in the dataframe you want to plot as a volcano. Scale the Adjusted P-value by log-base 10.
## Open the variable to view the dataframe and confirm
DESeq2_Expression_Results3$log10Pval = -log10(DESeq2_Expression_Results3$AdjP)


## Create the data points of interest for functional display.. ggplot struggles to color specific, different groups of data in a dataframe.
##  My way around this is to create a new column, entirely, and assign certain names ("strings") to the genes I want color differently.


## Create a new column (points of interest, "p.o.i"); assigning Genes with the values from their "GeneID" column is irrelevant-- we just need to put something in the cells; the values will be replaced in the next step
DESeq2_Expression_Results3$p.o.i. = DESeq2_Expression_Results3$GeneID


## Fill the points with appropriately indexed data; this will vary depending on what you want to display. For this example, we'll stick to two colors:
##  genes differentially expressed with a significance below 0.05, and all others

  ## "Differentially Expressed Genes" defined as being "points of interest" with Adjusted P-Values of < 0.05
DESeq2_Expression_Results3$p.o.i.[which(DESeq2_Expression_Results3$AdjP < 0.05)]= "Differentially Expressed Genes"


  ## "Non-differentially Expressed Genes" defined as being "points of interest" with Adjusted P-values > 0.05
DESeq2_Expression_Results3$p.o.i.[which(DESeq2_Expression_Results3$AdjP > 0.05)]= "Non-differentially Expressed Genes"




##### Create the first volcano plot and save the base of the plot settings as a variable ####
  

  ## Use the DESeq2....3 data frame as the data to plot
  ## Set the x-axis variable as log-base 2 fold change
  ## Set the y-axis as log-base 10 Adjusted P-value
  ## Set the points to be colored by the p.o.i. column previously created!

  ## Create labels on points we're interested in labeling in volcano plots
    ## Create the column, "labs"; hide all of the text labels with: ""
DESeq2_Expression_Results3$labs = ""


  ## Label only these selected items; adjust for each experiment; for the current example and a more global volcano plot, highlight the most extreme data points on the plot
ix_label1 = c(DESeq2_Expression_Results3$GeneID[DESeq2_Expression_Results3$log10Pval > 11])

  ## Fill the appropriate rows of the "labs" column with the corresponding Gene names for labeling those points
DESeq2_Expression_Results3$labs[DESeq2_Expression_Results3$log10Pval > 11] = DESeq2_Expression_Results3$GeneID[DESeq2_Expression_Results3$log10Pval > 11]

  ## Store the first volcano plot as a variable with the necessary components (ggplot function)

  ## Notes on ggplot commands:

  ## geom-point = alpha is a ggplot add-on function that controls point color transparency
  ## coord_cartesian is a ggplot add-on function that determines axes sizes
  ## labs is a ggplot add-on function that determines the plot labels; x-axis, y-axis, title, legend colors
  ## scale_color_manual is a ggplot add-on function that sets the colors of the plot, corresponding to the groups determined in the the base function ("col" in aes)
  ## Theme is used to manipulate the physical location of the plot title

  ## Sometimes, a column from a dataframe will contain apostrophes around the column name. Pay attention to this, you may need to revise depending on how R version you are running..
vol1 = ggplot(data = DESeq2_Expression_Results3, aes(x = DESeq2_Expression_Results3$`Log2(FC)`, y = DESeq2_Expression_Results3$log10Pval, col = DESeq2_Expression_Results3$p.o.i.)) + geom_point(alpha = 0.2) + coord_cartesian(xlim = c(-1.22,1.0), ylim = c(0, 13)) + labs(x = "log2(FC) in mRNA Exp. (Relative to WT)", y = "-log10(Adj. P-value)", title = "Tmem184b-GT Gene Exp. Changes (P 180)", col = "Gene Data Type") + scale_color_manual(values = c("blue", "gray")) + theme(plot.title = element_text(hjust = 0.5)) + theme_light()

  ## Plot it
vol1

  ## The parameters above will need to be adjusted by dataset. You have to pay attention to the dataframe's extremes (what's the highest/lowest log2(FC), the largest -log10Pval number?)

  ## Save it


  ## Add in more specificity
  ##   geom_text_repel is a ggplot add-on function that labels the selected points in the colors commanded

  ## Re-save the variable with the added functionality; for unlabeled plot
vol1 = vol1  + geom_text_repel(color = "black", aes(label=DESeq2_Expression_Results3$labs)) 

  ## Plot it with labels for the top two genes
vol1

## Save it



##### Create a second volcano plot AFTER PATHWAY ANALYSIS to create a more specific graph; save the base of the plot settings as a variable ####

## Re-write the p.o.i. column (or create a new column)
DESeq2_Expression_Results3$p.o.i. = ""

## Fill the points with appropriately indexed data;

## "Notable Growth Genes" defined as being "points of interest" with Adjusted P-Values of < 0.05
DESeq2_Expression_Results3$p.o.i.[c(10,17,35,77,218,288)]= "Transport Channels"

## "Differentially Expressed Genes" defined as being "points of interest"
DESeq2_Expression_Results3$p.o.i.[c(1:297)][-c(10,17,35,77,218,288)] = "D.E.G.s"


## "Non-differentially Expressed Genes" defined as being "points of interest" with Adjusted P-values > 0.05
DESeq2_Expression_Results3$p.o.i.[which(DESeq2_Expression_Results3$AdjP > 0.05)]= "Non-D.E.G.s"


## Create labels on points we're interested in labeling on the volcano plot

  ## Create the column, "labs2"; hide all of the text labels.
DESeq2_Expression_Results3$labs2 = ""
  ## Fill the appropriate rows of the "labs2" column with the corresponding Gene names used to label those points
DESeq2_Expression_Results3$labs2[c(c(10,17,35,77,218,288))] = DESeq2_Expression_Results3$GeneID[c(c(10,17,35,77,218,288))]



vol2 = ggplot(data = DESeq2_Expression_Results3, aes(x = DESeq2_Expression_Results3$`Log2(FC)`, y = DESeq2_Expression_Results3$log10Pval, col = DESeq2_Expression_Results3$p.o.i.)) + geom_point(alpha = 0.6) + coord_cartesian(xlim = c(-1.22,1.0), ylim = c(0, 13)) + labs(col = "Gene Data Type", title = "Tmem184b GT/GT Mouse HC Trx Profile", x = "log2(FC) mRNA Rel. to WT", y = "-log10 (Adj.P-val)") + scale_color_manual(values = c("sky blue", "gray", "navy")) + theme(plot.title = element_text(hjust = 0.5)) + theme_bw() + geom_hline(linetype = "dashed", yintercept = 1.30, color = "red")

  ## Plot it
vol2

  ## Save it


  ## Add in specificity, labeling all the most interesting differentially expressed genes, overwriting the old stored      plot variable
vol2 = vol2 + geom_text_repel(color = "black", aes(label=DESeq2_Expression_Results3$labs2)) 

  ## Plot it
vol2

  ## Save it



