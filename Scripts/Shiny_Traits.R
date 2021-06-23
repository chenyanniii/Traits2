library(shiny)


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

setwd("~/Desktop/Traits_Validated") 

##########################################################################################
################################## 0. Import Data ########################################
##########################################################################################

## Trait Value Data
InputMatrix <- read.csv("Input/RawData.csv") # Import Raw Measurement Data

Germination = read_xlsx("Input/Germination(sd_se).xlsx") %>%
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
        tabPanel("Phylogenetic_Position", plotOutput("Phylogenetic_Position")),
        tabPanel("Phylogenetic_Signal", plotOutput("Phylogenetic_Signal")),
        tabPanel("Model_Selection"), plotOutput("Model_Selection")
  )
)))

##########################################################################################
################################## 2. Define Server Logic ################################
##########################################################################################

server <- function(input, output, session) {
  
  #################### Selected Species As Input #################### 
  plot_react = reactive({filter(dat, Species %in% input$checked_species)})
  
  
  ## Display Phylogenetic Tree
  output$Phylogenetic_Tree <- renderPlot({
    plot_dat = plot_react()
    tree_species = PruneTree(plot_dat$Species)
    plotTree(tree_species, ftype="i")
    
  })
  
  #################### Display Phylogenetic Position ####################
  output$Phylogenetic_Position <- renderPlot({
  
  plot_phyloposit_fam = reactive({
    
    plot_dat = plot_react()
    
    tree_species = PruneTree(plot_dat$Species)
    dist = cophenetic.phylo(tree_species)
    phyloposi = isoMDS(dist) %>% as.data.frame()
    
    phyloposi_species = phyloposi %>% mutate(Species = row.names(phyloposi)) %>% 
      mutate(phy1 = round(scale(points.1),digits = 1)) %>% 
      mutate(phy2 = round(scale(points.2), digits = 1))
    
    merge(plot_dat, phyloposi_species) %>% 
      group_by(Family) %>% 
      summarise(Family, Classification, phy1, phy2) %>%
      distinct()
    
  })
  
  plot_phyloposit_family = plot_phyloposit_fam()
  
  phyloposit_family = ggplot(plot_phyloposit_family, aes(x=phy1, y=phy2, color=Classification)) +
    geom_point() + labs(x="phy1", y="phy2") + 
    geom_text_repel(aes(label = Family), size =3.5) + 
    theme_bw()
  
  phyloposit_family
  
  })
  
  #################### Display Phylogenetic Signal ####################
  ## need the variable from previous section: phyloposi_species
  ## have all the variable store globally, then display in reactive function
  
  output$Phylogenetic_Signal <- renderText({
    
    PhyloSig = reactive({
      
      plot_dat = plot_react()
      
      tree_species = PruneTree(plot_dat$Species)
      dist = cophenetic.phylo(tree_species)
      phyloposi = isoMDS(dist) %>% as.data.frame()
      
      phyloposi_species = phyloposi %>% mutate(Species = row.names(phyloposi)) %>% 
        mutate(phy1 = round(scale(points.1),digits = 1)) %>% 
        mutate(phy2 = round(scale(points.2), digits = 1))
      
      AllMatrix = merge(plot_dat, phyloposi_species)
      
      AllMatrix_G = filter(merge(AllMatrix, Germination), Species %in% input$checked_species)
      
      ### Format data before test phylogenetic signal
      Mass.All <- AllMatrix_G$mean_Mass
      names(Mass.All) <- AllMatrix$Species
      
      Height.All <- AllMatrix_G$mean_Height
      names(Height.All) <- AllMatrix$Species
      
      Area.All <- AllMatrix_G$mean_Area
      names(Area.All) <- AllMatrix$Species
      
      Germination.All <- AllMatrix_G$mean_Germination
      names(Germination.All) <- AllMatrix_G$Species

      PhyloSigMass.k <- phylosig(tree_species, Mass.All, method="K", test=TRUE, nsim=999)
      
      PhyloSigHeight.k <- phylosig(tree_species, Height.All, method="K", test=TRUE, nsim=999)
      
      phylosigArea.k <- phylosig(tree_species, Area.All, method="K", test=TRUE, nsim=999)
      
      phylosigGermination.k <- phylosig(tree_species, Germination.All, method="K", test=TRUE, nsim=999)
    
      paste0("The phylogenetic signal for seed mass is", PhyloSigMass.k, ".")
      
      paste0("The phylogenetic signal for seed height is", PhyloSigHeight.k, ".")
      
      paste0("The phylogenetic signal for seed area is", phylosigArea.k, ".")
      
      paste0("The phylogenetic signal for seed area is", phylosigArea.k, ".")
      
      paste0("The phylogenetic signal for seed germiantion rate is", phylosigGermination.k, ".")
      
    })
    
    PhyloSig()
    
  })
  
  #################### Display Model Selection ####################
  output$Model_Selection <- renderPlot({
    
    
    
    
  })
}
  

##########################################################################################
################################## 3. Run the Application ################################
##########################################################################################

shinyApp(ui = ui, server = server)
  
  
  
  
