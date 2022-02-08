Binomial Heatmap Script
================
Erik Larsen
7/17/2020

The following code was developed to display which differentially
expressed genes are involved in a given number of `PANTHER` pathways
using `ggplot2`.  
The data was processed using `Salmon`, `DESeq2`, and run through another
script in `Python` to acquire the proper data structure for `ggplot2` to
handle below in a `geom_tile` call.  
Other scripts use the `pheatmap` package for more traditional
transcriptional profiling.

## Environment Prep

Note that code can be sectioned and condensed with the ‘Alt + O’
command.

List of packages for this script:
[tidyverse](https://cran.r-project.org/package=tidyverse),
[stringr](https://cran.r-project.org/package=stringr),
[reshape2](https://cran.r-project.org/package=reshape2),[readr](https://cran.r-project.org/package=readr),
[ggplot2](https://cran.r-project.org/package=ggplot2)

``` r
library(readr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(stringr)
```

``` r
  ## Import the relevant data
e13_DEGs_and_Pathways =  read_csv("C:/Users/Erik/Desktop/BoxCopy/Lab/Omics/RNAseq/Embryonic DRG/Analysis/Heatmaps/Custom Python Enrichr Pathway Clustergram e13 GT.csv")

  ## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
e13_Heatmap = melt(e13_DEGs_and_Pathways, id = "DEGs")
colnames(e13_Heatmap) = c("DEGs", "Pathways", "Presence In Pathway")

  ## Remove excessive pathway name strings
Pathway_Names = e13_Heatmap$Pathways %>%
  gsub(x = e13_Heatmap$Pathways, pattern = " Homo sapiens .+.?", replacement = "")

  ## Re-shape the pathways for aesthetic viewing
Term_list = c("Axon guidance mediated by Slit/Robo", "De novo pyrimidine\ndeoxribonucleotide biosynthesis",
              "Dopamine receptor mediated\nsignaling pathway", "FAS signaling pathway", "Heme biosynthesis",
              "Integrin signaling pathway", "Muscarinic acetylcholine\nreceptor 2 and 4 signaling pathway",
              "Oxidative stress response", "Synaptic vesicle trafficking", "Wnt signaling pathway")

  ## Create a for loop to replace the strings
for (i in 1:length(Term_list)){
  Pathway_Names[ which(Pathway_Names == unique(Pathway_Names)[i]) ] = Term_list[i]
}
  ## Put the new strings back into the dataframe
e13_Heatmap[,2] = Pathway_Names

  ## Make the digital values factors so it's compatible to graph
e13_Heatmap$`Presence In Pathway` = as.factor(e13_Heatmap$`Presence In Pathway`)
```

## Create and store the heatmap

(There are 84 genes involved, so for aesthetic readability, they have
been omitted from the plot)

    ## Warning: Use of `e13_Heatmap$`Presence In Pathway`` is discouraged. Use
    ## `Presence In Pathway` instead.

    ## Warning: Use of `e13_Heatmap$Pathways` is discouraged. Use `Pathways` instead.

    ## Warning: Use of `e13_Heatmap$DEGs` is discouraged. Use `DEGs` instead.

![](Binomial-Heatmap-Script_files/figure-gfm/Heatmap%20of%20DEGs%20Identified%20Across%20PANTHER%20Pathways%20wLabs-1.png)<!-- -->
