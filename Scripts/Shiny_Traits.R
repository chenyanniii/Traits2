library(shiny)
library(tidyverse)

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

#setwd("~/Desktop/Traits_Validated") 

##########################################################################################
################################## 0. Import Data ########################################
##########################################################################################

## Trait Value Data
InputMatrix <- read.csv("Input/RawData.csv")  # Import Raw Measurement Data

Rep = 10 # Number of replicates in measurements

dat = InputMatrix %>% group_by(Family, Species, Classification) %>% 
  dplyr::summarise (mean_Area = mean(Area), se_Area = sd(Area)/sqrt(Rep), 
             mean_Height = mean(Height), se_Height = sd(Height)/sqrt(Rep),
             mean_Mass = mean(Mass), 
             se_Mass = sd(Mass/sqrt(Rep))) %>% 
  as.data.frame() 

## Input Phylogenetic Tree
## Some prep work
treeVascularPlants <- read.tree("Input/Vascular_Plants_rooted.tre")
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


##########################################################################################
############################### 1. Define UI for application #############################
##########################################################################################

ui <- fluidPage(
  #### sidebarLayout(sidebarPanel(),mainPanel())
  pageWithSidebar(
    headerPanel("Phylogenetic Comparative Methods Interactive Learning"),
    sidebarPanel(
      checkboxGroupInput(
        inputId = "checked_species",
        label = "Which species do you want to using for building phylogeny?",
        choices = c("Eryngium_leavenworthii","Polytaenia_nuttalli","Asclepias_asperula",
                    "Centaurea_americana","Coreopsis_lanceolata","Coreopsis_tinctoria",
                    "Echinacea_angustifolia","Gutierrezia_sarothrae","Helianthus_annuus",
                    "Liatris_mucronata","Ratibida_columnifera","Tradescantia_occidentalis",
                    "Astragalus_crassicarpus","Desmanthus_illinoensis","Corydalis_curvisiliqua",
                    "Phacelia_congesta","Herbertia_lahue","Monarda_citriodora","Salvia_azurea",
                    "Salvia_coccinea","Salvia_farinacea","Salvia_lyrata","Linum_rigidum",
                    "Callirhoe_involucrata","Pavonia_lasiopetala","Oenothera_rhombipetala",
                    "Argemone_albiflora","Phytolacca_americana","Rivina_humilis","Andropogon_gerardii",
                    "Aristida_purpurea","Bouteloua_curtipendula","Bouteloua_gracilis", 
                    "Chasmanthium_atifolium","Chloris_cucullata",
                    "Digitaria_ciliaris","Eragrostis_trichodes","Schizachyrium_scoparium",
                    "Sorghastrum_nutans","Sporobolus_airoides","Sporobolus_cryptandrus","Ipomopsis_rubra"),
        selected = c("Eryngium_leavenworthii","Polytaenia_nuttalli","Asclepias_asperula",
                     "Centaurea_americana","Coreopsis_lanceolata","Coreopsis_tinctoria",
                     "Echinacea_angustifolia","Gutierrezia_sarothrae","Helianthus_annuus",
                     "Liatris_mucronata","Ratibida_columnifera","Tradescantia_occidentalis",
                     "Astragalus_crassicarpus","Desmanthus_illinoensis","Corydalis_curvisiliqua",
                     "Phacelia_congesta","Herbertia_lahue","Monarda_citriodora","Salvia_azurea",
                     "Salvia_coccinea","Salvia_farinacea","Salvia_lyrata","Linum_rigidum",
                     "Callirhoe_involucrata","Pavonia_lasiopetala","Oenothera_rhombipetala",
                     "Argemone_albiflora","Phytolacca_americana","Rivina_humilis","Andropogon_gerardii",
                     "Aristida_purpurea","Bouteloua_curtipendula","Bouteloua_gracilis", 
                     "Chasmanthium_atifolium","Chloris_cucullata",
                     "Digitaria_ciliaris","Eragrostis_trichodes","Schizachyrium_scoparium",
                     "Sorghastrum_nutans","Sporobolus_airoides","Sporobolus_cryptandrus","Ipomopsis_rubra")
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Phylogenetic_Tree", plotOutput("Phylogenetic_Tree")),
        tabPanel("Phylogenetic_Position", plotOutput("Phylogenetic_Position"))
  )
)))

##########################################################################################
################################## 2. Define Server Logic ################################
##########################################################################################

server <- function(input, output, session) {

  ## Selected Species As Input  
  plot_dat = reactive({filter(dat, Species %in% input$checked_species)})
  tree_species = PruneTree(plot_dat$Species)
  
  ## Display Phylogenetic Tree
  output$Phylogenetic_Tree <- renderPlot({
    
    plotTree(tree_species, ftype="i")
    
  })
  
  ## Display Phylogenetic Position
  output$Phylogenetic_Position <- renderPlot({
    
    ## Phyloposit of input species
    dist = reactive({cophenetic.phylo(tree_species)})
    phyloposi = reactive({isoMDS(dist) %>% as.data.frame()})
    
    phyloposi_species = reactive({phyloposi %>% 
        mutate(Species = row.names(phyloposi)) %>% 
        mutate(phy1 = round(scale(points.1),digits = 1)) %>% 
        mutate(phy2 = round(scale(points.2), digits = 1))})
    
    plot_phyloposit = reactive({merge(plot_dat, phyloposi_species)})
    plot_phyloposit_family = reactive({plot_phyloposit %>% group_by(Family) %>% 
        summarise(Family, Classification, phy1, phy2) %>%
        distinct()})
    ## ggplot: phyloposit_family
    phyloposit_family = ggplot(plot_phyloposit_family, aes(x=phy1, y=phy2, color=Classification)) +
      geom_point() + labs(x="phy1", y="phy2") + 
      geom_text_repel(aes(label = Family), size =3.5) + 
      theme_bw()
    
     })
}

##########################################################################################
################################## 3. Run the Application ################################
##########################################################################################

shinyApp(ui = ui, server = server)
