Bioinformatics
================
Erik Larsen
7/17/2021

I developed this Markdown as a means to document some of bioinformatics
analyses of data published by M. Bhattacharya Lab in the publication,
[Transmembrane protein TMEM184B is necessary for interleukin-31–induced
itch](https://journals.lww.com/pain/Abstract/9000/Transmembrane_protein_TMEM184B_is_necessary_for.97918.aspx).

These analyses include MA & volcano plots, heatmaps, and bar graphs on a
number of different RNA-sequencing datasets. The RNAseq data was
pseudoaligned with the `Salmon` algorithm, then run through
`Usegalaxy.org`’s `DESeq2` wrapper algorithm to obtain DGEA files.
`Galaxy` also returned TPM data for use in transcriptional profiling
(heatmaps). To obtain Pathway and GO info, the
[Enrichr](https://maayanlab.cloud/Enrichr/) database was queried from
the console using the
[enrichR](https://github.com/guokai8/Enrichr#:~:text=Description%20EnrichR%20is%20a%20package%20can%20be%20used,species%20pubished%20by%20ENSEMBL%20and%20included%20with%20Bioconductor.)
package.

The source script of this file is the [Bioinformatics R
Script](C:/Users/Erik/Desktop/BoxCopy/Programming%20Scripts%20and%20Data/Bio/Scripts/R/Broad/Bioinformatics%20Script.R).

Standalone analyses or tutorials, including volcano plotting, and gene
ontology and pathway analysis can be found in other R Markdown files and
scripts: [Gene Ontology and Pathway Analysis R
Markdown](C:/Users/Erik/Desktop/BoxCopy/Lab/Omics/Gene-Ontology-and-Pathway-Analysis-for-Lab.html),
companion [GO Analysis R
Script](C:/Users/Erik/Desktop/BoxCopy/Programming%20Scripts%20and%20Data/Bio/Scripts/R/Specific/Searches/GO%20Analysis.R),
[Volcano plot tutorial R
Markdown](C:/Users/Erik/Desktop/BoxCopy/Lab/Omics/Volcano-plot-tutorial.html),
and companion [Volcano R
Script](C:/Users/Erik/Desktop/BoxCopy/Programming%20Scripts%20and%20Data/Bio/Scripts/R/Specific/Graphs/Volcano%20Script.R).

## Environment Prep

Note that code can be sectioned and condensed with the `Alt + O`
command.

List of packages for this script:
[tidyverse](https://cran.r-project.org/package=tidyverse),
[stringr](https://cran.r-project.org/package=stringr),
[plyr](https://cran.r-project.org/package=plyr),
[reshape2](https://cran.r-project.org/package=reshape2),
[dplyr](https://cran.r-project.org/package=dplyr),
[readr](https://cran.r-project.org/package=readr),
[biomaRt](https://bioconductor.org/packages/biomaRt/),
[GO.db](http://bioconductor.riken.jp/packages/3.0/data/annotation/html/GO.db.html),
[PANTHER.db](https://bioconductor.org/packages/PANTHER.db/),
[BiocGenerics](https://bioconductor.org/packages/BiocGenerics),
[pheatmap](https://cran.r-project.org/package=pheatmap)

Install biology-based packages with `BiocManager`. Load the
[BiocGenerics](https://bioconductor.org/packages/BiocGenerics) package
for Bioconductor-relevant functionality (installing packages from
[Bioconductor](https://www.bioconductor.org/))

``` r
library(BiocGenerics) ## Needed to install and/or load Bioconductor packages
```

Load the [PANTHER Pathway html
metadatabase](https://bioconductor.org/packages/PANTHER.db/) package and
specify the species of interest. Ignore warnings

-   Don’t update the PANTHER database if the option pops up. Let the
    package load before continuing.

``` r
library(PANTHER.db) ## Needed in the Enrichr package when querying pathways
pthOrganisms(PANTHER.db) = "MOUSE"
```

Load the remaining packages

``` r
library(enrichR) ## Taps the Enrichr database; much better utility than the website

  ## For data wrangling (slicing/adding/removing/melting/rearranging dataframes and their columns and rows):
library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)

library(readr) ## For importing data and files

library(stringr) ## Awesome for manipulating strings

  ## Databases for downstream analyses of lists via gene ontology:
library(biomaRt)
#library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(GO.db)

library(pheatmap) ## For creating awesome heatmaps

library(httpuv) ## For including plots in Markdown file output
```

Obtain databases from enrichr by creating a variable and viewing it.
Subset from the variable a few databases of interest.

``` r
DBs = listEnrichrDbs() ## Store it as a variable for better visualization and eventual subsetting
diff_DBs = DBs$libraryName[c(50,127,91,92,105,162)]
DBs = DBs$libraryName[c(130,131,132,102,14,110)]

  ## Peek at the selections
diff_DBs
```

    ## [1] "ENCODE_TF_ChIP-seq_2015"                 
    ## [2] "Enrichr_Submissions_TF-Gene_Coocurrence" 
    ## [3] "L1000_Kinase_and_GPCR_Perturbations_down"
    ## [4] "L1000_Kinase_and_GPCR_Perturbations_up"  
    ## [5] "huMAP"                                   
    ## [6] "ProteomicsDB_2020"

``` r
DBs
```

    ## [1] "GO_Biological_Process_2018" "GO_Cellular_Component_2018"
    ## [3] "GO_Molecular_Function_2018" "Panther_2016"              
    ## [5] "PPI_Hub_Proteins"           "Jensen_DISEASES"

Create dataframes to house bioinformatics analyses returned from
enrichr.

``` r
  ## Create dataframes to store returned analyses as variables
GO_Processes = data.frame()
GO_Cell_Comps = data.frame()
GO_Mol_Funcs = data.frame()
PATHWAYS = data.frame()
PPIs = data.frame()
DISEASES = data.frame()
ENCODE_TFs = data.frame()
Enrichr_TFs = data.frame()
GPCR_Down_Perturbs = data.frame()
GPCR_Up_Perturbs = data.frame()
GO_PROCESS_GENES = data.frame()
GO_MOL_FUNC_GENES = data.frame()
GO_CELL_COMP_GENES = data.frame()
PATHWAY_GENES = data.frame()
PPI_GENES = data.frame()
DISEASE_GENES = data.frame()
ENCODE_TF_GENES = data.frame()
ENRICHR_TF_GENES = data.frame()
GPCR_DOWN_PERTURB_GENES = data.frame()
GPCR_UP_PERTURB_GENES = data.frame()
```

Create a data frame to collect data returned by enrichr; then create the
function that will tap enrichr’s database and populate all the
dataframes.

These dataframes are created verbosely, so only one (**GO\_Processes**)
is included as an example for brevity.

``` r
  ## Create a function that will return Enrichr-based functional analysis on a given gene set
    ## First, create a data frame that will house combined enrichr output
      ## (this is for exporting collapsed, collective data into a CSV)
Enrichr_Functional_Output = data.frame()

  ## Create the analysis function
Enrichr_Analysis = function(Genes){
  ## first: Gene ontology biological processes; this identifies cell/biological processes with which the genes in the dataset are associated
  GO_Processes = as.data.frame(
    enrichr(
      c(Genes), DBs[1])
  )
  
  
  ## Remove unnecessary columns
  GO_Processes = GO_Processes[,-c(5,6,7)]
   
  
  ## Clean up the dataframes
  GO_Processes = GO_Processes %>%
    separate(GO_Biological_Process_2018.Overlap, c("Overlap", "Process Size"), "/")
  ## Convert the values into integers for computation
  GO_Processes$Overlap = as.numeric(GO_Processes$Overlap)
  GO_Processes$`Process Size` = as.numeric(GO_Processes$`Process Size`)
  
  
  ## Add Enrichr's Z-score. This determines a term/pathway's rank change across
  ## all terms/pathways when compared to the term/pathway's "pre-rank" based on
  ## multiple simulations of the same number of enriched genes.
  GO_Processes = add_column(GO_Processes, EnrichrZscore = GO_Processes$GO_Biological_Process_2018.Combined.Score/log(GO_Processes$GO_Biological_Process_2018.P.value), .before = 6)
  ## Add the "weighted overlap ratio" of a term/pathway's number of enriched genes. 
    ## The fraction of differentially expressed genes belonging to a term/pathway
    ## ("overlap"; numerator) out of the total genes in that term/pathway/hub (denominator)
    ## is scaled by the overlap. This weights differentially expressed genes by pathway size.
  GO_Processes = add_column(GO_Processes, Weighted_Overlap_Ratio = GO_Processes$Overlap*(GO_Processes$Overlap/GO_Processes$`Process Size`), .before = 4)
  ## Add the modifcation of Enrichr's "Combined Score". This scales the term/pathway's Adjusted P-value by the Enrichr Z-score. 
  GO_Processes = add_column(GO_Processes, Modified_Combined_Score = (abs(GO_Processes$EnrichrZscore))*(-log(GO_Processes$GO_Biological_Process_2018.Adjusted.P.value)), .before = 6)
  ## Remove GO IDs
  GO_Processes$GO_Biological_Process_2018.Term = GO_Processes$GO_Biological_Process_2018.Term %>% gsub(x = GO_Processes$GO_Biological_Process_2018.Term, pattern = " \\(.*\\)$", replacement = "")
  
 
  ## Extract a concatenated string containing genes involved in each GO Biological Process
  GO_PROCESS_GENES = str_split(GO_Processes$GO_Biological_Process_2018.Genes, pattern = ";", simplify = TRUE)
  for (i in 1:nrow(GO_PROCESS_GENES)){
    GO_PROCESS_GENES[i,] = paste(substr(GO_PROCESS_GENES[i,], 1, 1),
                                 tolower(substr(GO_PROCESS_GENES[i,], 2, 7)), sep = "") 
  }
  
  
  ## Make sure the Enrichr reference results are correctly ordered
  GO_Processes = GO_Processes[order(GO_Processes$Modified_Combined_Score, decreasing = TRUE), ]
  ## Re-set the index to make sure indeces are correct
  row.names(GO_Processes) = NULL
  
  
  ## Export these data frames to the global environment
  GO_Processes <<- GO_Processes
  
  return(c("GO BIOLOGICAL PROCESSES", GO_Processes[1:10,1]))
}
```

## Data import

Import the adult DRG DESeq2 file

``` r
aDRG = read.csv("C:/Users/Erik/Desktop/BoxCopy/Lab/Omics/RNAseq/Adult DRG/Processed Galaxy Output/Test Results to Upload/DESeq2 Expression Results.csv")
  
  ## Filter (subset) genes that went undetected or were outliers in terms of counts;
  ## new dataframe should not contain any NAs in p-value columns
aDRG3 = subset(aDRG, (!is.na(aDRG[,"AdjP"])))
 ## Filter the DEGs by removing rRNAs and mitochondrial tRNAs, along with pseudogenes, etc.
aDRG9 = aDRG3 %>% filter(!grepl(aDRG3$GeneID, pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))
#mt.+.?$|  <-- string identifier for mitochondrial tRNAs

  ## Create a column in the DESeq2 dataframe that scales the Adjusted P-value by log-base 10
aDRG9$log10ADJP = -log10(aDRG9$AdjP)

  ## Add a column to the DEG dataset that contains a string, describing whether the gene is differentially expressed
    ## First create the column and use Gene IDs as place-holders
aDRG9$g.o.i. = aDRG9$GeneID
  ## Replace DEGs with the string, "DEGs"
aDRG9$g.o.i.[which(aDRG9$AdjP <= 0.05)]= "DEGs"
  ## Replace the remaining genes with "Non-DEGs"
aDRG9$g.o.i.[which(aDRG9$AdjP >= 0.05)]= "Non-DEGs"

aDRG_DEG_list = c(aDRG9$GeneID[which(aDRG9$g.o.i. == "DEGs")])
```

Import the e13 DRG DESeq2 file as above (not shown).

Import the P0 DRG DESeq2 file and proceed as before (not shown).

Import the P10 DRG DESeq2 file and proceed as before (not shown).

Import the in vitro DRG data.

``` r
DIV14 = read.csv("C:/Users/Erik/Desktop/BoxCopy/Lab/Omics/RNAseq/Embryonic DRG/Processed Galaxy Output/Test Results to Upload/Cultured Embryos for Genotype at DIV14.csv")

colnames(DIV14) = colnames(aDRG)[c(1:7)]

DIV143 = subset(DIV14, (!is.na(DIV14[,"AdjP"])))
DIV149 = DIV143 %>% filter(!grepl(DIV143$GeneID, pattern = "Rps.+.?$|RP.+.?$|Rpl.+.?$|MRPL.+.?$|Mrpl.+.?$|MRPS.+.?$|Mrps.+.?$|.*Rik.+$|.*Rik$|Gm.+.?$|^[A-Z]+[A-Z].+.?$|^[0-9]+.+.?$"))

DIV149 = add_column(DIV149, X..WT = 2^(DIV149$log2.FC.[]))
colnames(DIV149)[8] = colnames(aDRG)[8]

DIV149$log10ADJP = -log10(DIV149$AdjP)

DIV149$g.o.i. = DIV149$GeneID
DIV149$g.o.i.[which(c(DIV149$AdjP < 0.01 & DIV149$log2.FC. > 0))]= "Up DEGs"
DIV149$g.o.i.[which(c(DIV149$AdjP < 0.01 & DIV149$log2.FC. < 0))]= "Down DEGs"
```

## Query the enrichr database

Search the `Enrichr` database using the function created above with any
of the above gene lists.

``` r
Enrichr_Analysis(Genes = aDRG_DEG_list)
```

    ## Uploading data to Enrichr... Done.
    ##   Querying GO_Biological_Process_2018... Done.
    ## Parsing results... Done.

    ##  [1] "GO BIOLOGICAL PROCESSES"                                                  
    ##  [2] "negative regulation of calcium ion-dependent exocytosis"                  
    ##  [3] "positive regulation of CD4-positive, alpha-beta T cell activation"        
    ##  [4] "positive regulation of interferon-gamma secretion"                        
    ##  [5] "vocal learning"                                                           
    ##  [6] "peripheral nervous system neuron development"                             
    ##  [7] "positive regulation of fibroblast migration"                              
    ##  [8] "regulation of interleukin-12 secretion"                                   
    ##  [9] "positive regulation of microtubule motor activity"                        
    ## [10] "heparan sulfate proteoglycan biosynthetic process, enzymatic modification"
    ## [11] "imitative learning"

## Prep data for heatmaps

Import the adult transcriptional profile (normalized TPM).

``` r
  ## Import the adult normalized counts file.
    ## These are expression estimates for each gene, for each sample/replicate,
    ## where each gene's value is normalized to its sample's effect size
aTPM = read_csv("C:/Users/Erik/Desktop/BoxCopy/Lab/Omics/RNAseq/Adult DRG/Processed Galaxy Output/Counts files to Upload/RNASeqRepResults.csv", col_names = TRUE)
  ## Rename columns
colnames(aTPM) = c("GeneID", "WT1", "WT2", "WT3", "WT4", "Mut1", "Mut2", "Mut3", "Mut4")

  ## Subset by only genes in the filtered aDRG DESeq2 file
ExpProfile = aTPM[aTPM$GeneID %in% c(aDRG9$GeneID),]
```

Create a function that makes a dataframe housing the Z-scored
transcriptional profile to then graph using `pheatmap`.

Run the Z-score function.

``` r
Find_Row_Z(Expression_Profile = ExpProfile)
```

Arrange the data for generating a full transcriptional profile.

``` r
  ## Don't remove non-DEGs
Z = Z %>%
  filter(GeneID %in% aDRG9$GeneID[])

  ## Create a list of gene names ordered by euclidean distance by which to re-order the dataframe
Euclid_dist_order = hclust(dist(Z[,c(2:9)], method = "euclidean"))$order
  ## The names (not the numbers)
Euclid_dist_ord_Genes = c(Z$GeneID[Euclid_dist_order])

  ## Transform again to order by clusters (Euclidean distances)
Za = Z %>%
  mutate(GeneID =  factor(GeneID, levels = Euclid_dist_ord_Genes)) %>%
  arrange(GeneID)

  ## Find where the "Itch-related DEGs" are in the clustered matrix subsetted to DEGs for the heatmap
itch_DEGs = c("Il31ra", "Cysltr2", "Npy2r", "Sst", "Htr1a", "P2rx3", "Lpar3", "Lpar5", "Scn11a", "Scn10a", "Mrgprd", "Trpc6", "Trpc3", "F2rl2", "Htr1f", "Osmr", "Fxyd2", "Htr4", "Mrgprx1", "Ptgdr", "Trpa1", "Trpm6", "Hrh1", "Mrgpra3", "Tmem184b")

temp = vector()
for( i in 1:length(itch_DEGs)){
  temp[i] = which(Euclid_dist_ord_Genes == itch_DEGs[i])
}
itch_index = temp
#itch_index

  ## Confirm by indexing; copy and paste to customize the clustered heatmap
#Euclid_dist_ord_Genes[c(itch_index)]
  ## Confirm the location/order of the transformed matrix/dataframe
  ## Should be "Il31ra"
Za[itch_index[1],]
```

    ## # A tibble: 1 x 9
    ##   GeneID   WT1   WT2   WT3   WT4   Mut1   Mut2   Mut3   Mut4
    ##   <fct>  <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
    ## 1 Il31ra  1.19  1.18 0.441 0.830 -0.943 -0.859 -0.927 -0.909

## Plot the heatmaps

View the full transcriptional profile of **Adult Tmem184b-mutant DRG
neurons**.

``` r
 ## Visualize the Z-scored, Euclidean-clustered and ordered gene TPMs across replicates with "pheatmap"
pheatmap(mat = Za[,2:9],
         color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)),
         clustsering_distance_rows = "euclidean",
         angle_col = 0,
         cutree_rows = 2,
         treeheight_row = 35,
         treeheight_col = 7,
         show_rownames = F)
```

![](Bioinformatics_files/figure-gfm/Full%20Z-scored,%20Euclidean-clustered%20and%20ordered%20pheatmap-1.png)<!-- -->

**View the cluster near Tmem184b**.

``` r
  ## Visualize around the Tmem184b cluster
pheatmap(mat = Za[1597:1637,2:9],
         color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)),
         clustsering_distance_rows = "euclidean",
         angle_col = 0,
         treeheight_row = 35,
         treeheight_col = 9,
         labels_row = Euclid_dist_ord_Genes[1597:1637])
```

![](Bioinformatics_files/figure-gfm/Tmem184b%20cluster%20Zoom%20on%20Z-scored%20Euclidean-clustered%20and%20ordered%20pheatmap-1.png)<!-- -->

**View the cluster including select itch transcripts**.

``` r
  ## Visualize the "Itch-related DEGs"
pheatmap(mat = Za[c(itch_index),2:9],
         color = colorRampPalette(c("navy", "white", "darkgoldenrod4"))((50)),
         clustsering_distance_rows = "euclidean",
         angle_col = 0,
         treeheight_row = 35,
         treeheight_col = 9,
         cutree_rows = 4,
         labels_row = Euclid_dist_ord_Genes[c(itch_index)])
```

![](Bioinformatics_files/figure-gfm/Itch%20trx%20cluster%20Zoom%20on%20Z-scored%20Euclidean-clustered%20and%20ordered%20pheatmap-1.png)<!-- -->

Arrange the data to include only DEGs.

``` r
  ## Transform the profile to include only DEGs
Zb = Z %>%
  filter(GeneID %in% aDRG9$GeneID[c(1:376)])
```

Prep it for a heatmap as before (not shown).

**Plot the results of the profile of only DEGs in Tmem184b-mutant DRG
neurons.**

![](Bioinformatics_files/figure-gfm/Z-scored%20Euclidean-clustered%20and%20ordered%20pheatmap%20of%20DEGS-1.png)<!-- -->

**Zoom to Tmem184b**.

![](Bioinformatics_files/figure-gfm/Tmem184b%20Cluster%20Zoom%20Z-scored%20Euclidean-clustered%20ordered%20pheatmap%20of%20DEGs-1.png)<!-- -->

**Zoom to the area comprising many itch genes**.

![](Bioinformatics_files/figure-gfm/Itch%20DEG%20Cluster%20Zoom%20Z-scored%20Euclidean-clustered%20ordered%20pheatmap%20of%20DEGs-1.png)<!-- -->

**View select itch genes**.  
The difference between the first three heatmaps and the last four
heatmaps is hierarchically clustering only on DEGs in the last four.
Clustering was performed on all genes in the first three.

![](Bioinformatics_files/figure-gfm/Itch%20DEG%20only%20of%20Z-scored%20Euclidean-clustered%20ordered%20pheatmap%20of%20DEGs-1.png)<!-- -->
