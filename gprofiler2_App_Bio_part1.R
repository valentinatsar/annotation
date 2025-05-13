#search the available ensembl datasets-genome assemblies per species and select which one to implement for annotation
library("BiocManager")
BiocManager::install("biomaRt")

# Collection of Ensembl databases to extract and analyse large amounts of data 
library(biomaRt)

mart = useEnsembl("genes")
listDatasets(mart)

########################################################################################################
#### gprofiler2 package - Functional Enrichment Analysis ####

install.packages(gprofiler2)
library(gprofiler2)

# annotate one gene list - gost() -> functional enrichment analysis -> chr3, locations 145489075 + 131527909
gostres = gost(query = c("ENSOARG00020000431","ENSOARG00020000548","ENSOARG00020006159","ENSOARG00020006505",
                         "ENSOARG00020026711","ENSOARG00020006808","ENSOARG00020036268","ENSOARG00020036542",
                         "ENSOARG00020014046","ENSOARG00020012556","ENSOARG00020013661","ENSOARG00020013827",
                         "ENSOARG00020013894","ENSOARG00020014249","ENSOARG00020014442","ENSOARG00020014571",
                         "ENSOARG00020014635","ENSOARG00020015356","ENSOARG00020015408","ENSOARG00020015432",
                         "ENSOARG00020015494","ENSOARG00020015560","ENSOARG00020016080","ENSOARG00020016299",
                         "ENSOARG00020016394","ENSOARG00020016620","ENSOARG00020016959","ENSOARG00020017012",
                         "ENSOARG00020017475","ENSOARG00020017499","ENSOARG00020017550","ENSOARG00020017769",
                         "ENSOARG00020038789","ENSOARG00020017880","ENSOARG00020017906","ENSOARG00020017951",
                         "ENSOARG00020017976","ENSOARG00020018004","ENSOARG00020018036","ENSOARG00020018067",
                         "ENSOARG00020018091","ENSOARG00020036004","ENSOARG00020029218","ENSOARG00020037444",
                         "ENSOARG00020038302","ENSOARG00020036697","ENSOARG00020031790","ENSOARG00020032407",
                         "ENSOARG00020038393","ENSOARG00020029262","ENSOARG00020035683","ENSOARG00020033603","ENSOARG00020027934"),
               organism = "oarambouillet", # for ARS-UI_Ramb_v2.0 sheep genome assembly, "chircus" for ARS1 - goats  
               significant = TRUE, user_threshold = 0.05, correction_method = "fdr", # filter for significance and adjust for multiple testing
               evcodes = TRUE, highlight = TRUE) 

#Interactive Manhattan plot - uses the output object from gost() as an input 
gostplot(gostres, capped = TRUE, interactive = TRUE)

#write specific terms to highlight on the manhattan plot
p1 = gostplot(gostres, capped = TRUE, interactive = FALSE)
publish_gostplot(p1, highlight_terms = c("GO:0001501", "GO:0009952")) 

#significant_terms contains only the enriched terms across all your regions or gene lists.
gostres_filtered <- gostres$result[gostres$result$highlight == TRUE, ]




##########################################################################################################################
##Annotate multiple genomic regions together

#Create a list with two genomic regions
gene_lists <- list(
  "Region1" = c("ENSOARG00020000431","ENSOARG00020000548","ENSOARG00020006159","ENSOARG00020006505",
                "ENSOARG00020026711","ENSOARG00020006808","ENSOARG00020036268","ENSOARG00020036542",
                "ENSOARG00020014046","ENSOARG00020012556","ENSOARG00020013661","ENSOARG00020013827",
                "ENSOARG00020013894","ENSOARG00020014249","ENSOARG00020014442","ENSOARG00020014571",
                "ENSOARG00020014635","ENSOARG00020015356","ENSOARG00020015408","ENSOARG00020015432",
                "ENSOARG00020015494","ENSOARG00020015560","ENSOARG00020016080","ENSOARG00020016299",
                "ENSOARG00020016394","ENSOARG00020016620","ENSOARG00020016959","ENSOARG00020017012",
                "ENSOARG00020017475","ENSOARG00020017499","ENSOARG00020017550","ENSOARG00020017769",
                "ENSOARG00020038789","ENSOARG00020017880","ENSOARG00020017906","ENSOARG00020017951",
                "ENSOARG00020017976","ENSOARG00020018004","ENSOARG00020018036","ENSOARG00020018067",
                "ENSOARG00020018091","ENSOARG00020036004","ENSOARG00020029218","ENSOARG00020037444",
                "ENSOARG00020038302","ENSOARG00020036697","ENSOARG00020031790","ENSOARG00020032407",
                "ENSOARG00020038393","ENSOARG00020029262","ENSOARG00020035683","ENSOARG00020033603","ENSOARG00020027934"),
  "Region2" = c("ENSOARG00020003528","ENSOARG00020039732","ENSOARG00020003767","ENSOARG00020004466",
                "ENSOARG00020004542","ENSOARG00020004965","ENSOARG00020034002","ENSOARG00020040152",
                "ENSOARG00020030797","ENSOARG00020034889","ENSOARG00020040754","ENSOARG00020036306",
                "ENSOARG00020036459","ENSOARG00020000900","ENSOARG00020000924","ENSOARG00020001084",
                "ENSOARG00020001354","ENSOARG00020001535","ENSOARG00020001619","ENSOARG00020001651",
                "ENSOARG00020001686","ENSOARG00020028642","ENSOARG00020040416","ENSOARG00020038003",
                "ENSOARG00020033036","ENSOARG00020037313","ENSOARG00020029157")
)



