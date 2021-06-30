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