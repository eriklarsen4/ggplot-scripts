

##### Heatmaps ####


library(readr)
library(tidyverse)
library(ggplot)
library(reshape2)
library(dplyr)
library(plyr)

##### e13 Panther GO Bio Processes Heatmap #####

## Import the relevant data
e13_DEGs_and_Pathways =  read_csv("M:/Omics/Custom Python Enrichr Clustergram e13 GT.csv")

## Export the gene list for photoshopping
#write.csv(e13_DEGs_and_Pathways$DEGs, "M:\\PAPER ASSEMBLY\\Itch paper\\Figure Drafts\\e13 FDR 01 Heatmap DEGs.csv", sep = ",")


## Re-arrange the data so that columns and rows can be run appropriately in a heatmap
e13_Heatmap = melt(e13_DEGs_and_Pathways, id = "DEGs")
colnames(e13_Heatmap) = c("DEGs", "Pathways", "Presence In Pathway")

## Remove excessive pathway name strings
Pathway_Names = e13_Heatmap$Pathways %>%
  gsub(x = e13_Heatmap$Pathways, pattern = " Homo sapiens .+.?", replacement = "")

## Put the new strings back into the dataframe
e13_Heatmap[,2] = Pathway_Names

## Make the digital values factors so it's compatible to graph
e13_Heatmap$`Presence In Pathway` = as.factor(e13_Heatmap$`Presence In Pathway`)

## Create and store the heatmap core as a variable
## Use geom_tile to map the binary values

## Customize by removing the filler surrounding the graph and the tickmarks
## Center the Title

HEATER = ggplot(e13_Heatmap, aes(e13_Heatmap$Pathways, e13_Heatmap$DEGs)) +
  geom_tile(aes(fill = e13_Heatmap$`Presence In Pathway`), color = "black") +
  scale_fill_manual(values = c("dark blue", "gold")) +
  labs(x = "Panther GO Biological Processes",
       y = "Differentially Expressed Genes (FDR < 0.01)",
       title = "Heatmap of Tmem184b-GT-affected Processes and Pathways") +
  guides(fill = guide_legend(title = "Presence in Pathway")) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(colour = "black", size = 0.1),
        panel.grid.minor = element_line(colour = "black", size = 0.1),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0))

# axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), , axis.text= element_blank()

## Plot the heatmap
HEATER


## Without labels
HEATER = ggplot(e13_Heatmap, aes(e13_Heatmap$Pathways, e13_Heatmap$DEGs)) +
  geom_tile(aes(fill = e13_Heatmap$`Presence In Pathway`), color = "black") +
  scale_fill_manual(values = c("dark blue", "gold")) +
  labs(x = "Panther GO Biological Processes",
       y = "Differentially Expressed Genes (FDR < 0.01)",
       title = "Heatmap of Tmem184b-GT-affected Processes and Pathways") +
  guides(fill = guide_legend(title = "Presence in Pathway")) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 20, l = 0)),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major = element_line(colour = "black", size = 0.1),
        panel.grid.minor = element_line(colour = "black", size = 0.1),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0),
        axis.text= element_blank())

## Plot the heatmap
HEATER
