##########################################################################################
############################### 1. Define UI for application #############################
##########################################################################################
library(shiny)

source("./helper.R")

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
        tabPanel("Phylogenetic_Signal", tableOutput("Phylogenetic_Signal")),
        tabPanel("Model_Selection", tableOutput("Model_Selection"))
      )
    )#mainPanel
    )#pageWithsidebar
  )#fluidPage
