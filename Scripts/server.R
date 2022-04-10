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
  

  output$Phylogenetic_Signal <- renderTable({

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
    
    Height.All <- AllMatrix_G$mean_Height
    names(Height.All) <- AllMatrix$Species
    
    Area.All <- AllMatrix_G$mean_Area
    names(Area.All) <- AllMatrix$Species
    
    Germination.All <- AllMatrix_G$mean_Germination
    names(Germination.All) <- AllMatrix_G$Species
    
    # calculate the phylogenetic signals
    
    cal_phyloSigMass.k = reactive({phylosig(tree_species, Mass.All, method="K", test=TRUE, nsim=999)})
    
    cal_phyloSigHeight.k = reactive({phylosig(tree_species, Height.All, method="K", test=TRUE, nsim=999)})
    
    cal_phyloSigArea.k = reactive({phylosig(tree_species, Area.All, method="K", test=TRUE, nsim=999)})
    
    cal_phyloSigGermination.k = reactive({phylosig(tree_species, Germination.All, method="K", test=TRUE, nsim=999)})
    
    phyloSigMass.k = cal_phyloSigMass.k()
    
    phyloSigHeight.k = cal_phyloSigHeight.k()
    
    phyloSigArea.k = cal_phyloSigArea.k()
    
    phyloSigGermination.k = cal_phyloSigGermination.k()
    

    # text display results
    
    #paste0("The phylogenetic signal of seed area is ", round(phyloSigArea.k$K, digits = 2), " (p = ", round(phyloSigArea.k$P, digits = 2) , ").",
    #       " The phylogenetic signal of seed germination is ", round(phyloSigGermination.k$K, digits = 2), " (p = ", round(phyloSigGermination.k$P, digits = 2) , ").",
    #       " The phylogenetic signal of seed mass is ", round(phyloSigMass.k$K, digits = 2), " (p = ", round(phyloSigMass.k$P, digits = 2) , ").",
    #       " The phylogenetic signal of seed height is ", round(phyloSigHeight.k$K, digits = 2), " (p = ", round(phyloSigHeight.k$P, digits = 2), ").")
 
  
    # table display results
    tab_phyloSig = matrix(c(phyloSigMass.k$K, phyloSigHeight.k$K, phyloSigArea.k$K, 
                            phyloSigGermination.k$K, phyloSigMass.k$P, phyloSigHeight.k$P, phyloSigArea.k$P, 
                            phyloSigGermination.k$P), ncol = 2)
    
    colnames(tab_phyloSig) = c("Blomberg's K", "P-value")
    
    Variables = c("Mass", "Height", "Area", "Germination")
    tab_phyloSig1 = cbind(Variables,tab_phyloSig)
    tab_phyloSig1
  
  })
  
  #################### Display Model Selection ####################

  output$Model_Selection <- renderTable({
    
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
    
    plot_react_scaled = reactive({filter(AllMatrix_G, Species %in% input$checked_species) %>%
        mutate(scaled_Area = scale(mean_Area)) %>%
        mutate(scaled_Height = scale(mean_Height)) %>%
        mutate(scaled_Mass = scale(mean_Mass)) %>%
        mutate(scaled_Germination = scale(mean_Germination)) %>%
        dplyr::select(scaled_Area, scaled_Height, scaled_Mass, scaled_Germination, phy1, phy2)
    })
    
    AllMatrix_G_scaled = plot_react_scaled ()
   
    vars = names(AllMatrix_G_scaled[-4]) 
    models = list()
    for (i in seq_along(vars)){
      vc = combn(vars, i)
      for (j in 1:ncol(vc)){
        model = as.formula(paste0('scaled_Germination ~', paste0(vc[, j], collapse = '+')))
        models = c(models, model)
      }
    }
    
    glmmodels = lapply(models, function(x) glm(x, data = AllMatrix_G_scaled))
    
    AICc = sapply(models, function(x) AICc(glm(x, data = AllMatrix_G_scaled), 
                                           return.K = FALSE))
    
    names(AICc) = sapply(glmmodels, function(x) x$formula)
    
    ## textPlot
    #paste0("The minimun AICc vallue is equal to ", round(min(AICc), digits = 2), "; 
             #The maximum AICc is equal to ", round(max(AICc), digits = 2), ".")

    ## tablePlot
    tab_models = data_frame(names(AICc),AICc)
    colnames(tab_models) = c("Model", "AICc")
    
    tab_models
    ## tablePlot1
    #tab_model = matrix(sapply(glmmodels, function(x) x$formula))
    #tab_model1 = cbind(tab_model, AICc)
    
    #rownames(tab_model1) = NULL
    #colnames(tab_model1) = c("Model", "AICc")

    
  }) 
  

 
}
      



