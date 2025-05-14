#search the available ensembl datasets-genome assemblies per species and select which one to implement for annotation
library("BiocManager")
BiocManager::install("biomaRt")

# Collection of Ensembl databases to extract and process large amounts of data 
library(biomaRt)

mart = useEnsembl("genes")
listDatasets(mart)

########################################################################################################
#### gprofiler2 package - Functional Enrichment Analysis ####

install.packages(gprofiler2)
library(gprofiler2)

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


# Annotate regions together using `multi_query = TRUE`
gostres <- gost(query = gene_lists,
                organism = "oarambouillet",  # for ARS-UI_Ramb_v2.0 sheep genome assembly
                multi_query = TRUE,   # this allows parallel querying of all regions
                significant = TRUE, user_threshold = 0.05, correction_method = "fdr",
                highlight = TRUE)
                

# View the results
head(gostres$result)

#Interactive Manhattan plot - uses the output object from gost() as an input 
gostplot(gostres, capped = TRUE, interactive = TRUE)

# plot the results
p2 = gostplot(gostres, capped = TRUE, interactive = FALSE)
publish_gostplot(p2, highlight_terms = c("GO:0001501", "GO:0009952")) #or gostres$result[c(1, 5, 10),] #NULL


##########################################################################################################################                                        
##visualization of gene lists using the g:Profiler - gost() output 
# install.packages("BiocManager")
library("BiocManager")
BiocManager::install(c("clusterProfiler", "enrichplot", "DOSE"), force = T)
library(clusterProfiler)
library(enrichplot)
library(DOSE) # needed to convert to enrichResult object
library(ggplot2)

#run again gost() function with multi_query=FALSE to visualize with packages requiring single-enrichment results
gostres2 <- gost(query = gene_lists,
                organism = "oarambouillet",  # for ARS-UI_Ramb_v2.0 sheep genome assembly
                multi_query = FALSE,   # this allows parallel querying of all regions
                significant = TRUE, user_threshold = 0.05, correction_method = "fdr", # adjust for multiple testing
                evcodes = TRUE, highlight = TRUE) # filter for significant results  

head(gostres2$result)


#significant_terms contains only the enriched terms across all your regions or gene lists.
gostres2_filtered <- gostres2$result[gostres2$result$highlight == TRUE, ]

# modify the g:Profiler data frame
#Plotting functions in clusterProfiler/enrichplot) expect certain columns like: GeneRatio (how many genes map to a term), 
#BgRatio (how common the term is in the annotation background generally), p.adjust, geneID
gostres2_mod <- gostres2_filtered[,c("query", "source", "term_id",
                            "term_name", "p_value", "query_size",
                            "intersection_size", "term_size",
                            "effective_domain_size", "intersection")]

gostres2_mod$GeneRatio = paste0(gostres2_mod$intersection_size, "/", gostres2_mod$query_size)

#term_size=total number of genes annotated to a given term, effective_domain_size=total number of genes annotated of the studied species
gostres2_mod$BgRatio = paste0(gostres2_mod$term_size, "/", gostres2_mod$effective_domain_size)

names(gostres2_mod) = c("Cluster", "Category", "ID", "Description", "p.adjust",
                  "query_size", "Count", "term_size", "effective_domain_size",
                  "geneID", "GeneRatio", "BgRatio")

gostres2_mod$geneID = gsub(",", "/", gostres2_mod$geneID)
row.names(gostres2_mod) = gostres2_mod$ID

# define as compareClusterResult object
gostres2_mod_cluster = new("compareClusterResult", compareClusterResult = gostres2_mod)


# define as enrichResult object
gostres2_mod_enrich = new("enrichResult", result = gostres2_mod)

#dotplot - gene ratio
enrichplot::dotplot(gostres2_mod_cluster)

#barplot - gene count and abundance (Highlighting terms with many genes)
barplot(gostres2_mod_enrich, showCategory = 80, font.size = 8) +
  ggplot2::facet_grid(~Cluster) +
  ggplot2::xlab("Intersection size")

