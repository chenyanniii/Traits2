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
  
  output$Phylogenetic_Signal <- renderPlot({
    
    ## prepdata
    plot_dat = plot_react()
    # prep species list
    tree_species = PruneTree(plot_dat$Species)
    
    dist = cophenetic.phylo(tree_species)
    phyloposi = isoMDS(dist) %>% as.data.frame()
    
    phyloposi_species = phyloposi %>% mutate(Species = row.names(phyloposi)) %>% 
      mutate(phy1 = round(scale(points.1),digits = 1)) %>% 
      mutate(phy2 = round(scale(points.2), digits = 1))
    
    AllMatrix = merge(plot_dat, phyloposi_species)
    
    AllMatrix_G = filter(merge(AllMatrix, Germination), Species %in% input$checked_species)
    
    Mass.All <- AllMatrix_G$mean_Mass
    names(Mass.All) <- AllMatrix_G$Species
    
    cal_phyloSigMass.k = reactive({phylosig(tree_species, Mass.All, method="K", test=TRUE, nsim=999)})
    
    cal_phyloSigHeight.k = reactive({phylosig(tree_species, Height.All, method="K", test=TRUE, nsim=999)})
    
    cal_phylosigArea.k = reactive({phylosig(tree_species, Area.All, method="K", test=TRUE, nsim=999)})
    
    cal_phylosigGermination.k = reactive({phylosig(tree_species, Germination.All, method="K", test=TRUE, nsim=999)})
    
    phyloSigMass.k = cal_phyloSigMass.k()
    
    phyloSigHeight.k = cal_phyloSigHeight.k()
    
    phylosigArea.k = cal_phylosigArea.k()
    
    #phylosigGermination.k = cal_phylosigGermination.k()
    
    
    pMass = plot(phyloSigMass.k)
    
    pHeight = plot(phyloSigHeight.k)
    
    pArea = plot(phylosigArea.k)
    
    pGermination = plot(phylosigGermination.k)
    
    p1 = qplot(1:10,rnorm(10))
    p2 = qplot(1:10,rnorm(10))
    
    #ptlist = list(cal_phyloSigMass.k(),cal_phyloSigHeight.k())
    
    #pMass|pHeight
    #pMass|pArea
    pMass|pGermination
    
    p1+p2+plot_layout(widths = c(2,1))
    
  }) 
  #################### Display Model Selection ####################
  output$Model_Selection <- renderPlot({
    
    
    
    
  })
}


