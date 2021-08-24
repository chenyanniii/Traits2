library(ape)
library(nlme)
library(geiger)
library(corpcor)
library(nloptr)
library(RColorBrewer)
library(OUwie)
library(readxl)
library(plyr)
library(tidyverse)
library(phytools)
library(MASS)
library(ggplot2)
library(ggrepel)
library(AICcmodavg)
library(MASS)
library(broom)
library(patchwork)
library(ggtree)


##########################################################################################
################################## 0. Import Data ########################################
##########################################################################################

## Trait Value Data
InputMatrix <- read.csv("../Input/RawData.csv") # Import Raw Measurement Data

Germination = read_xlsx("../Input/Germination(sd_se).xlsx") %>%
  mutate(mean_Germination = mean) %>%
  mutate(se_Germination = se) %>%
  dplyr::select(Species, mean_Germination, se_Germination)  # Import Germination Data

Rep = 10 # Number of replicates in measurements

dat = InputMatrix %>% group_by(Family, Species, Classification) %>% 
  dplyr::summarise (mean_Area = mean(Area), se_Area = sd(Area)/sqrt(Rep), 
                    mean_Height = mean(Height), se_Height = sd(Height)/sqrt(Rep),
                    mean_Mass = mean(Mass), 
                    se_Mass = sd(Mass/sqrt(Rep))) %>% 
  as.data.frame() 

## Input Phylogenetic Tree
## Some prep work
treeVascularPlants <- read.tree("../Input/Vascular_Plants_rooted.tre")
tips <-treeVascularPlants$tip.label
Genera<-unique(sapply(strsplit(tips,"_"),function(x) x[1]))

## PruneTree Function
## Keeptip + add tips under the same genus  
PruneTree <- function(x){
  ## x, a list of desired species in "xx_xx" format
  DesiredSpecies <- unique(x)
  ## Desired Genera
  DesiredGenera <- sapply(strsplit(DesiredSpecies, "_"),function(x) x[1])
  
  ## TreeSpecies Data Frame
  SpeciesGeneraSpecies <- data_frame(TreeSpecies = treeVascularPlants$tip.label,
                                     TreeGenera = sapply(strsplit(tips,"_"),function(x) x[1]),
                                     DesiringGenera = TreeGenera %in% intersect(Genera,DesiredGenera) ,
                                     DesiringSpecies = TreeSpecies %in% DesiredSpecies) 
  
  SpeciesListSpecies <- filter(SpeciesGeneraSpecies, DesiringSpecies == "TRUE")
  
  GeneraSpecies <- filter(SpeciesGeneraSpecies, SpeciesGeneraSpecies$TreeGenera %in% 
                            setdiff(DesiredGenera, SpeciesListSpecies$TreeGenera))
  SpeciesListGenera <- group_by(GeneraSpecies, TreeGenera) %>% group_modify(~ head(.x, 1L))
  
  LISTALLSPECIES <- rbind.data.frame(SpeciesListSpecies, SpeciesListGenera)
  
  treeTestedSpecies <- keep.tip(treeVascularPlants, LISTALLSPECIES$TreeSpecies)
  
  ## Replacing tip label
  
  aaa <- data_frame(DesiredSpecies = DesiredSpecies, 
                    TreeGenera = sapply(strsplit(DesiredSpecies,"_"),function(x) x[1]))
  bbb <- merge(aaa, SpeciesListGenera)
  
  treeTestedSpecies$tip.label <- mapvalues(treeTestedSpecies$tip.label, c(bbb$TreeSpecies), c(bbb$DesiredSpecies))
  tree_x <- treeTestedSpecies
  plotTree(tree_x, ftype="i")
  return(tree_x)
}
