################################################################################
### Network analysis with NetCoMi
# with R version 4.3.0 (2023-04-21) -- "Already Tomorrow"
## network on control subset using spearman and all subsets on genus level

## bacterial analysis of Inoculation concentration experiment
## Sarah Benning

################################################################################
### 1.- Working Directory
setwd("/home/snst1/ZMHZ/pespind1/2_Sarah") # in Linux
getwd()   #  returns the current working direct

###################################################
### 2.- required libraries. 
#      .libPaths()  ##    Check where are your packages were/will be installed

library(NetCoMi)


###################################################
#### 3.-Load/Read data
load("./inputs/data_networksSarah_IV.RData")  


###################################################
#### 4.-subset phyloseq object for all treatments, 
#### 	plus new phyloseq objets agglomerated on genus level 

library(phyloseq)   
library(microbiome)

control_subset <- subset_samples(uBiome, Treatment == "control")
summarize_phyloseq(control_subset)
control_subset <- prune_taxa(taxa_sums(control_subset) > 0, control_subset)

control_genus <- tax_glom(control_subset, taxrank = "Genus")



###################################################
# III.- network on ASV level for control samples
#   i) netConstruct
netConstructed_control <- netConstruct(control_subset,
                                       filtTax = "numbSamp",
                                       filtTaxPar = list(numbSamp = 3),
                                       # jointPrepro = FALSE,
                                       measure = "spearman",
                                       # measurePar = list(pulsar.params = list(rep.num = 10), nlambda = 10),
                                       dissFunc = "unsigned",
                                       normMethod = "clr",
                                       zeroMethod = "pseudo",
                                       sparsMethod = "threshold",
                                       thresh = 0.3, 
                                       verbose = 3)

saveRDS(netConstructed_control, file= "./outputs/netConstructed_control_spearman.rds" )   #Pame# This is not really necesary


#    ii) netAnalyze
netAnalyzed_control <- netAnalyze(netConstructed_control,
                                  clustMethod = "cluster_fast_greedy",
                                  hubPar = "degree", 
                                  centrLCC = FALSE, # if TRUE, then compute centralities only for the largest connected component
                                  hubQuant = 0.95,  # a node is identified as hub if for degree", the node’s centrality value is above the 90% quantile of the fitted log-normal distribution.
                                  lnormFit = TRUE,  # a log-normal distribution is fitted to the centrality values to identify nodes with “highest” centrality values.
                                  weightDeg = FALSE, normDeg = FALSE)     
                        # weightDeg and normDeg are set to FALSE so that the degree of
                        # a node is simply defined as number of nodes that are adjacent to the node.

summary(netAnalyzed_control, numbNodes = 22L)

############ Save networkfile
save.image("./outputs/output_networksSarah_CONTROL.RData")   #Pame# intermediate point, save until here


#############################################
#############################################
## genus level

###################################################
# IV.- network on genus level for control samples
#   i) netConstruct
netConstructed_control_gen <- netConstruct(control_genus,
                                           filtTax = "numbSamp",
                                           filtTaxPar = list(numbSamp = 3),
                                           # jointPrepro = FALSE,
                                           measure = "spearman",
                                           # measurePar = list(pulsar.params = list(rep.num = 10), nlambda = 10),
                                           dissFunc = "unsigned",
                                           normMethod = "clr",
                                           zeroMethod = "pseudo",
                                           sparsMethod = "threshold",
                                           thresh = 0.3, 
                                           verbose = 3)

saveRDS(netConstructed_control_gen, file= "./outputs/netConstructed_control_genus_spearman.rds" )   #Pame# This is not really necesary


#    ii) netAnalyze
netAnalyzed_control_gen <- netAnalyze(netConstructed_control_gen,
                                      clustMethod = "cluster_fast_greedy",
                                      hubPar = "degree", 
                                      centrLCC = FALSE, # if TRUE, then compute centralities only for the largest connected component
                                      hubQuant = 0.95,  # a node is identified as hub if for degree", the node’s centrality value is above the 90% quantile of the fitted log-normal distribution.
                                      lnormFit = TRUE,  # a log-normal distribution is fitted to the centrality values to identify nodes with “highest” centrality values.
                                      weightDeg = FALSE, normDeg = FALSE)     
                                    # weightDeg and normDeg are set to FALSE so that the degree of
                                    # a node is simply defined as number of nodes that are adjacent to the node.

summary(netAnalyzed_control_gen, numbNodes = 22L)

############ Save networkfile
save.image("./outputs/output_networksSarah_CONTROL.RData")   #Pame# intermediate point, save until here


