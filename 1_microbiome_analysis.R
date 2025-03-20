## Analysis of Inoculation concentration experiment
## Using SRS data

## 05.02.2024

## Sarah Benning
## R Version 4.3.1

###############

## 16S amplicon sequence analysis of Sarah's samples
## Inoculation concentration experiment (2021)

###############


# pacman to install and load libraries

if (!require("pacman")) install.packages(
  "pacman",
  verbose = F)

# BiocManager for Bioconductor libraries

if (!require("BiocManager")) install.packages(
  "BiocManager",
  verbose = F)

# devtools for GitHub repositories 

if (!require("devtools")) install.packages(
  "devtools",
  verbose = F)

# GitHub libraries

pacman::p_load_gh(
  "jbisanz/qiime2R",
  "benjjneb/decontam",
  "mikemc/speedyseq",
  "adw96/breakaway",
  "adw96/DivNet")

# install/load the remainder of libraries

pacman::p_load(
  tidyverse,
  vegan,
  ape,
  microbiome,
  DESeq2,
  ggpubr,
  patchwork,
  hrbrthemes,
  SRS,
  phangorn,
  DECIPHER,
  MicrobiotaProcess,
  compositions,
  hues,
  zCompositions)


library("dplyr")


select = dplyr::select
transform = microbiome::transform


# select path to files (easier if the script is always in the same folder as your files)

setwd("D:/Users/sarah.benning/OneDrive - Helmholtz Zentrum M?nchen/R79_inoculation/R_data_analysis/amplicon_new/")

metadata_table = "mapping_subset.csv"

taxonomic_assignment = "taxa_table_ASV.txt"

counts_table = "count_table_ASV_subset.txt"

tree = "Galaxy4-[FASTTREE_on_data_3_tree.nhx].nhx"



control_samples = c("sarahextneg1", "sarahextneg2", "sarahpcrneg1", "sarahextneg3", "sarahpcrneg3")


# import dada2 results, clean them and integrate in a phyloseq object

# import metadata (if not working check separators of your csv)

metaTab = 
  read_csv2(metadata_table) %>% 
  mutate_if(is_character, as.factor)

print(metaTab)

# Import the phylogenetic tree

phytree = treeio::read.nhx(tree)

# import and clean taxonomy (check in taxa:table if capital letters, add "ASV" column name if neccessary)

pre_taxTab =
  read_tsv(taxonomic_assignment) %>% 
  map_df(~gsub(
    pattern = "metagenome|uncultured|unidentified|Unknown",
    replacement = NA,
    .x)) %>%
  mutate_if(is_character, str_trim) %>% 
  mutate(Domain = ifelse(is.na(Domain),
                         "U. Domain",
                         Domain),
         Phylum = coalesce(Phylum,
                           ifelse(grepl("^U.", Domain),
                                  Domain,
                                  paste("U.", Domain))),
         Class = coalesce(Class,
                          ifelse(grepl("^U.", Phylum),
                                 Phylum,
                                 paste("U.", Phylum))),
         Order = coalesce(Order,
                          ifelse(grepl("^U.", Class),
                                 Class,
                                 paste("U.", Class))),
         Family = coalesce(Family,
                           ifelse(grepl("^U.", Order),
                                  Order,
                                  paste("U.", Order))),
         Genus = coalesce(Genus,
                          ifelse(grepl("^U.", Family),
                                 Family,
                                 paste("U.", Family))),
         Species = coalesce(Species,
                            ifelse(grepl("^U.", Genus),
                                   Genus,
                                   paste("U.", Genus)))) %>% 
  column_to_rownames("ASV") %>% 
  filter(Domain %in% "Bacteria" &
           !Order %in% "Chloroplast" &
           !Family %in% "Mitochondria")

print(head(pre_taxTab %>% 
             as_tibble()))




# if you have batches with different negative controls

batch_nr  = c(rep("A", 21), rep("B", 7), rep("A", 3)) 


# import ASVs counts

raw_abuTab = 
  read_tsv(counts_table) %>% 
  column_to_rownames("...1") %>% 
  `colnames<-`(gsub(pattern = "_.*",replacement = "", x = colnames(.)))


pre_abuTab = 
  read_tsv(counts_table) %>% 
  column_to_rownames("...1") %>% 
  `colnames<-`(gsub(pattern = "_.*",replacement = "",x = colnames(.))) %>% 
  filter(!rownames(.) %in%
           isContaminant(
             t(.),
             neg = c(colnames(.) %in%
                       control_samples),
             batch = batch_nr,
             threshold = .1,
             normalize = T,
             detailed = F)) %>% 
  select_if(!names(.) %in% control_samples) %>%
  filter(rowSums(.) > 0 &
           rownames(.) %in% rownames(pre_taxTab))

colSums(pre_abuTab)


# final filtering of taxa table

taxTab = 
  pre_taxTab %>%
  filter(rownames(.) %in%
           rownames(pre_abuTab)) %>% 
  as.matrix()



# import into phyloseq-object

## if you got an error Error in (function (classes, fdef, mtable): 
## unable to find an inherited method for function 'tax_table' for signature '"matrix"', 
## change  tax_table(taxTab) to phyloseq::tax_table(taxTab)


uBiome =  
  phyloseq(
    otu_table(pre_abuTab,
              taxa_are_rows = T),
    phyloseq::tax_table(taxTab),
    phy_tree(phytree@phylo),
    sample_data(metaTab %>% 
                  column_to_rownames("id")))

uBiome


# save phyloseq-class object

save(uBiome,
     file = paste("uBiome_large_metadata_subset",
                  format(Sys.time(), "_%d%m%y_%H%M"),
                  ".RData",
                  sep = ""))



# Congratulations, now you have your phyloseq object prepared
# you can finally start to explore the microbial wonderland in your samples


######

## normalization with SRS

#####

srs_abuTab_16S = pre_abuTab %>% 
  SRS(Cmin = min(colSums(.))) %>% 
  `rownames<-`(rownames(pre_abuTab))

colSums(srs_abuTab_16S)



###import into phyloseq-object (srs_abuTab)

uBiome_16S_srs =  phyloseq(
  otu_table(srs_abuTab_16S,
            taxa_are_rows = T),
  phyloseq::tax_table(taxTab),
  phy_tree(phytree@phylo), 
  sample_data(metaTab %>% 
                column_to_rownames("id")))

###save phyloseq-class object

save(uBiome_16S_srs,
     file = paste("uBiome_large_metadata_subset_SRS",
                  format(Sys.time(), "_%d%m%y_%H%M"),
                  ".RData",
                  sep = ""))



# without normalitzation
load("uBiome_large_metadata_subset_050224_1330.RData")

phylo_subset <- uBiome

summarize_phyloseq(phylo_subset)

# load uBiome scaled with SRS
load("uBiome_large_metadata_subset_SRS_050224_1331.RData")

phylo_SRS <- uBiome_16S_srs

summarize_phyloseq(phylo_SRS)


##filtering singletons
#Filter  singletons (OTU with only one sequence)from the data sets, meaning all ASV that are present in only one sample (in case you have and want to filter them)(Benoit's script)

phylo_SRS_filtered <- prune_taxa(taxa_sums(phylo_SRS) > 1, 
                                 phylo_SRS)

summarize_phyloseq(phylo_SRS_filtered)


save(phylo_SRS_filtered,
     file = paste("uBiome_large_metadata_subset_SRS_filtered_singletons",
                  format(Sys.time(), "_%d%m%y_%H%M"),
                  ".RData",
                  sep = ""))


#variables_of_interest

variables_of_interest = c("Treatment", "Phytoalexins", "Shoot_growth", "Root_dry_weight")


#extracting the asv table

OTU1 = as(otu_table(phylo_SRS_filtered), "matrix")


#if(taxa_are_rows(uBiome_16S_pre_abutab)){OTU1 <- t(OTU1)}
# Coerce to data.frame

OTUdf = as.data.frame(OTU1)

tax1 = as(phyloseq::tax_table(phylo_SRS_filtered), "matrix")
taxdf = as.data.frame(tax1)


meta1 = as(sample_data(phylo_SRS_filtered), "matrix")
meta1df = as.data.frame(meta1)


# setwd("D:/Users/sarah.benning/OneDrive - Helmholtz Zentrum M?nchen/R79_inoculation/R_data_analysis/amplicon_new/")

write.csv2(OTUdf, "asvtable_removed_contaminants_SRS_singletons_removed.csv", row.names = TRUE)

write.csv2(taxdf, "taxtable_removed_contaminants_SRS_singletons_removed.csv", row.names = TRUE)

write.csv2(meta1df, "metadata_table_SRS_singletons_removed.csv", row.names = TRUE)



#####################################################################

## Microbiome analysis

#####################################################################



library(knitr)
library(tidyr)
library(rstatix)
library(svglite)
library(phyloseq)
library(ggplot2)    
library(readxl)
library(reshape2)
library(dplyr)        
library(microbiome)
library(DESeq2)
library(writexl)
library(microeco)
library(vegan)
library(microbiomeutilities)
library(lme4)
library(lmerTest)
library(nlme)
library(metagMisc)
library(gridExtra)
library(file2meco)
library(ggalluvial)
library(ggcorrplot)
library(ggpubr)
library(car)
library(DT)
library(tidyverse)
library(microViz)


setwd("C:/Users/SarahBenning/OneDrive - Helmholtz Zentrum München/R79_inoculation/R_data_analysis/amplicon_new/")





### I will use the SRS normalization for diversity metrices, 
### for other things you can also use the original one and to other transformations on it

#########
# Rarefaction curve and normalization


mat <- t(otu_table(phylo_subset))
class(mat) <- "matrix"
rarecurve(mat, step=50, cex=0.5, label = TRUE, 
          xlab = "Reads", ylab = "ASVs")


paste("max read count: ",max(sample_sums(phylo_subset)))

paste("min read count: ",min(sample_sums(phylo_subset)))

sample_sums(phylo_subset)


# after SRS scaling, obviously rarefaction curves look nice

mat_SRS <- t(otu_table(phylo_SRS_filtered))
class(mat_SRS) <- "matrix"
rarecurve(mat_SRS, step=50, cex=0.5, label = TRUE, 
          xlab = "Reads", ylab = "ASVs")

paste("max read count: ",max(sample_sums(phylo_SRS_filtered)))

paste("min read count: ",min(sample_sums(phylo_SRS_filtered)))

sample_sums(phylo_SRS_filtered)




## normalize to relative read counts

comp_ps <- microbiome::transform(phylo_subset, "compositional")

paste("max read count: ",max(sample_sums(comp_ps)))

paste("min read count: ",min(sample_sums(comp_ps)))



#####################################################################################

## Alpha Diversity

#####################################################################################


###### Estimating alpha diversity using normalized data

alpha_diversity_frame_SRS = 
  phylo_SRS_filtered %>%
  microbiome::alpha(index = c("Shannon", "inverse_simpson", "observed", "pielou")) %>%
  merge(sample_data(phylo_SRS_filtered),
        by = 0) 

write.csv2(alpha_diversity_frame_SRS, "alpha_diversity_indices_SRS.csv", row.names = TRUE)


#############

## Statistics on alpha diversity

############

data_stat <- alpha_diversity_frame_SRS


# Density plot

ggdensity(data_stat$observed, fill = "lightgray")+
  stat_overlay_normal_density(color = "red", linetype = "dashed")  #normal

ggdensity(data_stat$diversity_shannon, fill = "lightgray")+
  stat_overlay_normal_density(color = "red", linetype = "dashed")     #normal

ggdensity(data_stat$diversity_inverse_simpson, fill = "lightgray")+
  stat_overlay_normal_density(color = "red", linetype = "dashed")   #normal

ggdensity(data_stat$evenness_pielou, fill = "lightgray")+
  stat_overlay_normal_density(color = "red", linetype = "dashed")         #non-normal ?


# QQ plot

ggqqplot(data_stat$observed)     #normal

ggqqplot(data_stat$diversity_shannon)        #normal

ggqqplot(data_stat$diversity_inverse_simpson)      #normal

ggqqplot(data_stat$evenness_pielou)            #non-normal



# Shapiro test (if significant, data is non-normally distributed)

data_stat %>% shapiro_test(
  observed, diversity_shannon, diversity_inverse_simpson, evenness_pielou)


#only significant for evenness => non-normal, the rest is normally distributed

######################################################

## Check for Homogeneity of Variance
#####################################################

library(tidyverse)
library(car)


# Levene's test with one independent variable

# not significant, groups are not significantly different, have equal variances

leveneTest(observed ~ Treatment, data = data_stat)

leveneTest(diversity_shannon ~ Treatment, data = data_stat)

leveneTest(diversity_inverse_simpson ~ Treatment, data = data_stat)

leveneTest(evenness_pielou ~ Treatment, data = data_stat)

# all have equal variances


############################################################################

# observed (normal, equal_variances) --> t.test or anova

# diversity_shannon (normal, equal_variances) --> t.test or anova

# evenness_pielou (non-normal, equal_variances) --> wilcoxon

# diversity_inverse_simpson (normal, equal_variances) --> t.test or anova

############################################################################



#plotting shannon index

alpha_shannon_plot = ggplot(alpha_diversity_frame_SRS,
                            aes(x = Treatment,
                                y = diversity_shannon,
                                fill = Treatment)) + 
  stat_boxplot(geom ='errorbar', width=0.2, lwd=0.5, position = position_dodge(width = 0.8)) +
  geom_boxplot(
    linewidth = 1/4,
    alpha = 3/4,
    outlier.shape = NA) +
  geom_point(color = "black", size = 2, alpha = 0.9) +
  scale_fill_manual(values = c("#90EE90", "#228822", "#5F7FC7", "orange","#DA5724")) +
  #theme(axis.title.x = element_text(color="black", size=14, face="bold"),
  #      axis.title.y = element_text(color="black", size=14, face="bold")) +
  labs(x = "Treatment", y = "Shannon Index") +
  theme_classic(base_size = 20)

alpha_shannon_plot

ggsave(paste0("alpha_shannon_plot_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = alpha_shannon_plot)

ggsave(paste0("alpha_shannon_plot_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = alpha_shannon_plot)


#plotting observed diversity

alpha_observed_plot = ggplot(alpha_diversity_frame_SRS,
                             aes(x = Treatment,
                                 y = observed,
                                 fill = Treatment)) + 
  stat_boxplot(geom ='errorbar', width=0.2, lwd=0.5, position = position_dodge(width = 0.8)) +
  geom_boxplot(
    linewidth = 1/4,
    alpha = 3/4,
    outlier.shape = NA) +
  geom_point(color = "black", size = 2, alpha = 0.9) +
  scale_fill_manual(values = c(control="cyan4", R6="deeppink4", R7="deeppink4", R8="deeppink4", R9="deeppink4")) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  labs(x = "Treatment", y = "Observed Diversity") +
  theme_classic(base_size = 20)

alpha_observed_plot

ggsave(paste0("alpha_observed_plot_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = alpha_observed_plot)

ggsave(paste0("alpha_observed_plot_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = alpha_observed_plot)


#plotting inverse simpson diversity

alpha_inv_simpson_plot = ggplot(alpha_diversity_frame_SRS,
                                aes(x = Treatment,
                                    y = diversity_inverse_simpson,
                                    fill = Treatment)) + 
  stat_boxplot(geom ='errorbar', width=0.2, lwd=0.5, position = position_dodge(width = 0.8)) +
  geom_boxplot(
    linewidth = 1/4,
    alpha = 3/4,
    outlier.shape = NA) +
  geom_point(color = "black", size = 2, alpha = 0.9) +
  scale_fill_manual(values=c(control="cyan4", R6="deeppink4", R7="deeppink4", R8="deeppink4", R9="deeppink4")) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  labs(x = "Treatment", y = "Inverse Simpson Index") +
  theme_classic(base_size = 20)

alpha_inv_simpson_plot


ggsave(paste0("alpha_inv_simpson_plot_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = alpha_inv_simpson_plot)

ggsave(paste0("alpha_inv_simpson_plot_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = alpha_inv_simpson_plot)


# plotting evenness pielou index

alpha_evenness_pielou_plot = ggplot(alpha_diversity_frame_SRS,
                                    aes(x = Treatment,
                                        y = evenness_pielou,
                                        fill = Treatment)) + 
  stat_boxplot(geom ='errorbar', width=0.2, lwd=0.5, position = position_dodge(width = 0.8)) +
  geom_boxplot(
    linewidth = 1/4,
    alpha = 3/4,
    outlier.shape = NA) +
  geom_point(color = "black", size = 2, alpha = 0.9) +
  scale_fill_manual(values = c("#90EE90", "#228822", "#5F7FC7", "orange","#DA5724")) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  labs(x = "Treatment", y = "Evenness Pielou") +
  theme_classic(base_size = 20)

alpha_evenness_pielou_plot

ggsave(paste0("alpha_evenness_pielou_plot_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = alpha_evenness_pielou_plot)

ggsave(paste0("alpha_evenness_pielou_plot_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = alpha_evenness_pielou_plot)



library(gridExtra)

# arrange alpha diversity plots

grid.arrange(alpha_shannon_plot, 
             alpha_observed_plot, 
             alpha_inv_simpson_plot, 
             alpha_evenness_pielou_plot, 
             ncol = 2)



## Calculation diversity metrics 

diversity = estimate_richness(pseq.subset.rarefied)
shannon = diversity$Shannon
observed = diversity$Observed
observed.log <- log(observed)
evenness=shannon/observed.log
diversity$evennes=evenness
richness=as.data.frame(diversity)


DT::datatable(richness, class = "cell-border stripe", rownames = T,
              filter = "top", editable = TRUE, extensions = 'Buttons', options = list(
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                scrollY = '300px',
                scrollCollapse = TRUE
              ))


######################################################

## Beta Diversity

#####################################################

theme_set(theme_classic())


## PCoA

########################

PCoA_bray_ord <- ordinate(physeq = phylo_SRS_filtered,  
                          method = "PCoA", 
                          distance = "bray")

plot_scree(PCoA_bray_ord, "Scree Plot")


plot_PCoA_bray_ord <- plot_ordination(physeq = phylo_SRS_filtered, 
                                      ordination = PCoA_bray_ord, 
                                      color = "Treatment",
                                      axes = c(1,2)) + 
  geom_point(aes(colour = Treatment), size = 3) + 
  geom_point(size = 3) +
  scale_color_manual(
    values = c( "#90EE90", "#228822", "#5F7FC7", "orange","#DA5724")) +
  stat_ellipse(type = "t") + 
  ggtitle("PCoA - Bray Curtis Distance") 

plot_PCoA_bray_ord

ggsave(paste0("PCoA_Bray_ord_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = plot_PCoA_bray_ord)

ggsave(paste0("PCoA_Bray_ord_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = plot_PCoA_bray_ord)

###

set.seed(123123)

PCoA_wunifrac_ord <- ordinate(physeq = phylo_SRS_filtered,  
                              method = "PCoA", 
                              distance = "wunifrac")

plot_scree(PCoA_wunifrac_ord, "Scree Plot")


plot_PCoA_wunifrac_ord <- plot_ordination(physeq = phylo_SRS_filtered, 
                                          ordination = PCoA_wunifrac_ord, 
                                          color = "Treatment",
                                          axes = c(1,2)) + 
  geom_point(aes(colour = Treatment), size = 3) + 
  geom_point(size = 3) +
  scale_color_manual(
    values = c( "#90EE90", "#228822", "#5F7FC7", "orange","#DA5724")) +
  stat_ellipse(type = "t") + 
  ggtitle("PCoA - weighted Unifrac") 

## does not look very nice

###

set.seed(123123)

PCoA_jaccard_ord <- ordinate(physeq = phylo_SRS_filtered,  
                             method = "PCoA", 
                             distance = "jaccard")

plot_scree(PCoA_jaccard_ord, "Scree Plot")


plot_PCoA_jaccardc_ord <- plot_ordination(physeq = phylo_SRS_filtered, 
                                          ordination = PCoA_jaccard_ord, 
                                          color = "Treatment",
                                          axes = c(1,2)) + 
  geom_point(aes(colour = Treatment), size = 3) + 
  geom_point(size = 3) +
  scale_color_manual(
    values = c( "#90EE90", "#228822", "#5F7FC7", "orange","#DA5724")) +
  stat_ellipse(type = "t") + 
  ggtitle("PCoA - Jaccard") 



grid.arrange(
             plot_PCoA_bray_ord, 
             plot_PCoA_wunifrac_ord, 
             plot_PCoA_jaccardc_ord, 
             ncol = 2)


########################

## NMDS

########################

NMDS_bray_ord <- ordinate(physeq = phylo_SRS_filtered,  
                          method = "NMDS", 
                          distance = "bray")

plot_scree(PCoA_bray_ord, "Scree Plot")


plot_NMDS_bray_ord <- plot_ordination(physeq = phylo_SRS_filtered, 
                                      ordination = NMDS_bray_ord, 
                                      color = "Treatment",
                                      axes = c(1,2)) + 
  geom_point(aes(colour = Treatment), size = 3) + 
  geom_point(size = 3) +
  scale_color_manual(
    values = c( "#90EE90", "#228822", "#5F7FC7", "orange","#DA5724")) +
  stat_ellipse(type = "t") + 
  ggtitle("NMDS - Bray Curtis Distance") 

plot_NMDS_bray_ord

ggsave(paste0("NMDS_Bray_ord_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = plot_NMDS_bray_ord)

ggsave(paste0("NMDS_Bray_ord_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = plot_NMDS_bray_ord)

###

set.seed(123123)

NMDS_wunifrac_ord <- ordinate(physeq = phylo_SRS_filtered,  
                              method = "NMDS", 
                              distance = "wunifrac")

plot_scree(PCoA_wunifrac_ord, "Scree Plot")


plot_NMDS_wunifrac_ord <- plot_ordination(physeq = phylo_SRS_filtered, 
                                          ordination = NMDS_wunifrac_ord, 
                                          color = "Treatment",
                                          axes = c(1,2)) + 
  geom_point(aes(colour = Treatment), size = 3) + 
  geom_point(size = 3) +
  scale_color_manual(
    values = c( "#90EE90", "#228822", "#5F7FC7", "orange","#DA5724")) +
  stat_ellipse(type = "t") + 
  ggtitle("NMDS - weighted Unifrac") 

plot_NMDS_wunifrac_ord

###

set.seed(123123)

NMDS_jaccard_ord <- ordinate(physeq = phylo_SRS_filtered,  
                             method = "NMDS", 
                             distance = "jaccard")

plot_scree(PCoA_jaccard_ord, "Scree Plot")


plot_NMDS_jaccardc_ord <- plot_ordination(physeq = phylo_SRS_filtered, 
                                          ordination = NMDS_jaccard_ord, 
                                          color = "Treatment",
                                          axes = c(1,2)) + 
  geom_point(aes(colour = Treatment), size = 3) + 
  geom_point(size = 3) +
  scale_color_manual(
    values = c( "#90EE90", "#228822", "#5F7FC7", "orange","#DA5724")) +
  stat_ellipse(type = "t") + 
  ggtitle("NMDS - Jaccard") 

plot_NMDS_jaccardc_ord



grid.arrange(
             plot_NMDS_bray_ord, 
             plot_NMDS_wunifrac_ord, 
             plot_NMDS_jaccardc_ord, 
             ncol = 2)


### PCoA and Redundancy Analysis (RDA)

sarah1 <- tax_fix(phylo_SRS_filtered) # try tax_fix_interactive if you have problems with your own data
ibd <- phyloseq_validate(sarah1, remove_undetected = TRUE)


plot = ibd %>% 
  tax_transform("clr", rank = "unique") %>%
  ord_calc(constraints = c("Phytoalexins", "Shoot_growth", "Root_dry_weight")) %>%
  ord_plot( colour = "Treatment", size = 2, alpha = 0.9) 

### Code - https://david-barnett.github.io/microViz/reference/ord_plot.html
# https://david-barnett.github.io/microViz/reference/ord_plot.html
#https://david-barnett.github.io/microViz/articles/web-only/ordination.html?q=rda#rda

plot_ellipses <- plot + 
#  lims(x = c(-10, 10), y = c(-9, 9)) +
  stat_ellipse(aes(colour = Treatment)) +
  scale_colour_brewer(palette = "Paired")



ggsave(paste0("RDS_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = plot_ellipses)

ggsave(paste0("RDS_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = plot_ellipses)




### PCoA and Constrained Analysis of Principal Coordinates (CAP)

cap_ord <- ordinate(physeq = phylo_SRS_filtered,  
                    method = "CAP", 
                    distance = "bray", 
                    formula = ~  Shoot_growth + Root_dry_weight)

scree.cap <- plot_scree(cap_ord, "Scree Plot for MCs in Constrained Analysis of Principal Coordinates - CAPSCALE")
print(scree.cap)


cap_plot_bray <- plot_ordination(physeq = phylo_SRS_filtered, 
                                 ordination = cap_ord, 
                                 color = "Treatment", 
                                 # shape = "Fungicide_Treatment",
                                 axes = c(1,2)) + 
  geom_point(aes(colour = Treatment), size = 3) + 
  geom_point(size = 3) +
  scale_color_manual(
    values = c( "#90EE90", "#228822", "#5F7FC7", "orange","#DA5724"))+ 
  stat_ellipse(type = "t") + 
  ggtitle("CAP_Plot - Bray Curtis")


cap_plot_bray

ggsave(paste0("cap_plot_bray_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = cap_plot_bray)

ggsave(paste0("cap_plot_bray_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = cap_plot_bray)


arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP1, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP1, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))


cap_plot_bray_arrows <-  cap_plot_bray +
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "black", 
    arrow = arrowhead) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE)


ggsave(paste0("cap_plot_bray_SRS_arrows", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = cap_plot_bray_arrows)

ggsave(paste0("cap_plot_bray_SRS_arrows", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = cap_plot_bray_arrows)



#####################################################

# Statistical analyzes

#####################################################


# Permanova

bray <- phyloseq::distance(phylo_SRS_filtered, method = "bray")
braydf <- data.frame(sample_data(phylo_SRS_filtered))

adonis2(formula = bray ~ Treatment, data = braydf)


#           Df   SumOfSqs      R2       F       Pr(>F)    
# Treatment  4  0.65428       0.28731   2.1164  0.001 ***
# Residual  21  1.62300       0.71269                  
# Total     25  2.27728       1.00000 



# also shoot growth makes significant clusters

adonis2(formula = bray ~ Shoot_growth, data = braydf)


#               Df SumOfSqs   R2        F       Pr(>F)  
# Shoot_growth  1  0.18864    0.08284   2.1676  0.011 *
# Residual     24  2.08864    0.91716                
# Total        25  2.27728    1.00000 




## permanova of only inoculated samples

## all inoculated samples (n =21)
 inoc_subset <- subset_samples(phylo_SRS_filtered, Treatment != "control")
 
 inoc_subset <- prune_taxa(taxa_sums(inoc_subset) > 0, inoc_subset) # remove singletons
 summarize_phyloseq(inoc_subset)

bray <- phyloseq::distance(inoc_subset, method = "bray")
braydf <- data.frame(sample_data(inoc_subset))

adonis2(formula = bray ~ Treatment, data = braydf)


#                 Df     SumOfSqs      R2    F      Pr(>F)   
# Treatment        3    0.41957   0.26107 2.0021  0.003 **
# Residual         17   1.18753   0.73893                 
# Total            20   1.60711   1.00000  




## PCoA

########################

PCoA_bray_ord_inoc <- ordinate(physeq = inoc_subset,  
                          method = "PCoA", 
                          distance = "bray")

plot_scree(PCoA_bray_ord_inoc, "Scree Plot")


plot_PCoA_bray_ord_inoc <- plot_ordination(physeq = inoc_subset, 
                                      ordination = PCoA_bray_ord_inoc, 
                                      color = "Treatment",
                                      axes = c(1,2)) + 
  geom_point(aes(colour = Treatment), size = 3) + 
  geom_point(size = 3) +
  scale_color_manual(
    values = c( "#90EE90", "#228822", "#5F7FC7", "orange","#DA5724")) +
  stat_ellipse(type = "t") + 
  ggtitle("PCoA - Bray Curtis Distance") 

plot_PCoA_bray_ord_inoc

ggsave(paste0("PCoA_Bray_ord_SRS_inoc", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = plot_PCoA_bray_ord_inoc)

ggsave(paste0("PCoA_Bray_ord_SRS_inoc", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = plot_PCoA_bray_ord_inoc)



############

# comparison of means - alpha diversity
# everything is normally distributed, except evenness


data_stat <- alpha_diversity_frame_SRS


# shannon

stat_shannon_wilcox <- data.frame(
  compare_means
  (formula = diversity_shannon ~ Treatment, data = data_stat,
    method = "wilcox.test", paired = FALSE, group.by = NULL, p.adjust.method = "bonferroni"))

View(stat_shannon_wilcox)

stat_shannon_anova <- data.frame(
  compare_means
  (formula = diversity_shannon ~ Treatment, data = data_stat,
    method = "anova", paired = FALSE, group.by = NULL, p.adjust.method = "bonferroni"))

View(stat_shannon_anova)

write.csv2(stat_shannon_wilcox, "stat_shannon_wilcox_SRS.csv", row.names = FALSE)
write.csv2(stat_shannon_anova, "stat_shannon_anova_SRS.csv", row.names = FALSE)



# observed

stat_observed_wilcox <- data.frame(
  compare_means
  (formula = observed ~ Treatment, data = data_stat,
    method = "wilcox.test", paired = FALSE, group.by = NULL, p.adjust.method = "bonferroni"))

View(stat_observed_wilcox)

stat_observed_anova <- data.frame(
  compare_means
  (formula = observed ~ Treatment, data = data_stat,
    method = "anova", paired = FALSE, group.by = NULL, p.adjust.method = "bonferroni"))

View(stat_shannon_anova)

write.csv2(stat_observed_wilcox, "stat_observed_wilcox_SRS.csv", row.names = FALSE)
write.csv2(stat_observed_anova, "stat_observed_anova_SRS.csv", row.names = FALSE)



# inverse_simpson

stat_inverse_simpson_wilcox <- data.frame(
  compare_means
  (formula = diversity_inverse_simpson ~ Treatment, data = data_stat,
    method = "wilcox.test", paired = FALSE, group.by = NULL, p.adjust.method = "bonferroni"))

View(stat_inverse_simpson_wilcox)

stat_inverse_simpson_anova <- data.frame(
  compare_means
  (formula = diversity_inverse_simpson ~ Treatment, data = data_stat,
    method = "anova", paired = FALSE, group.by = NULL, p.adjust.method = "bonferroni"))

View(stat_inverse_simpson_anova)

write.csv2(stat_inverse_simpson_wilcox, "stat_inverse_simpson_wilcox_SRS.csv", row.names = FALSE)
write.csv2(stat_inverse_simpson_anova, "stat_inverse_simpson_anova_SRS.csv", row.names = FALSE)



# evenness

stat_evenness_pielou_wilcox <- data.frame(
  compare_means
  (formula = evenness_pielou ~ Treatment, data = data_stat,
    method = "wilcox.test", paired = FALSE, group.by = NULL, p.adjust.method = "bonferroni"))

View(stat_evenness_pielou_wilcox)

stat_evenness_pielou_kruskal <- data.frame(
  compare_means
  (formula = evenness_pielou ~ Treatment, data = data_stat,
    method = "kruskal.test", paired = FALSE, group.by = NULL, p.adjust.method = "bonferroni"))

View(stat_evenness_pielou_kruskal)

write.csv2(stat_evenness_pielou_wilcox, "stat_evenness_pielou_wilcox_SRS.csv", row.names = FALSE)
write.csv2(stat_evenness_pielou_kruskal, "stat_evenness_pielou_anova_SRS.csv", row.names = FALSE)



##################################################

# Taxonomic composition

# Box plot


mycols <- c( "#90EE90", "#228822", "#5F7FC7", "orange","#DA5724")


pn <- plot_taxa_boxplot(phylo_SRS_filtered,
                        taxonomic.level = "Family",
                        top.otu = 9,
                        group = "Treatment",
                        title = "Top 9 taxa at family level - Bacteria",
                        keep.other = FALSE,
                        add.violin= FALSE,
                        #  group.order = c("ACCO Drechs_NO", "ACCO Drechs_YES",
                        #                 "KLAR Drechs_NO", "KLAR Drechs_YES" ,
                        #                  "PLAN Drechs_YES", "PLAN Drechs_NO"),
                        group.colors = mycols)


box_plot_bacteria = pn + theme_biome_utils() +  theme(axis.title.x=element_blank(),
                                                      axis.text.x=element_blank(),
                                                      axis.ticks.x=element_blank()) +
  theme(legend.position = "bottom")

box_plot_bacteria


ggsave(paste0("Top 9 taxa at family level SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = box_plot_bacteria)

ggsave(paste0("Top 9 taxa at family level SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = box_plot_bacteria)


### Exporting data most abundant taxa abundances for linear models ###

top_nine_bacteria <- subset_taxa(pseq.subset.rarefied, Family %in% c("Bacillaceae",
                                                                     "Gemmatimonadaceae",
                                                                     "Itrasporangiaceae",
                                                                     "Rhodanobacteraceae",
                                                                     "Sphingomonadaceae",
                                                                     "U. Gaiellales",
                                                                     "U. KD4-96",
                                                                     "U. Vicinamibacteriales",
                                                                     "Xanthobacteraceae"))

top_nine_bacteria


top_nine_melted=psmelt(top_nine_bacteria)


write_xlsx(top_nine_melted,"top_nine_bacteria.xlsx")



## top taxa genus level

pnf <- plot_taxa_boxplot(phylo_SRS_filtered,
                        taxonomic.level = "Genus",
                        top.otu = 9,
                        group = "Treatment",
                        title = "Top 9 taxa at genus level - Bacteria",
                        keep.other = FALSE,
                        add.violin= FALSE,
                        #  group.order = c("ACCO Drechs_NO", "ACCO Drechs_YES",
                        #                 "KLAR Drechs_NO", "KLAR Drechs_YES" ,
                        #                  "PLAN Drechs_YES", "PLAN Drechs_NO"),
                        group.colors = mycols)


box_plot_bacteria_genus = pnf + theme_biome_utils() +  theme(axis.title.x=element_blank(),
                                                      axis.text.x=element_blank(),
                                                      axis.ticks.x=element_blank()) +
  theme(legend.position = "bottom")

box_plot_bacteria_genus


ggsave(paste0("Top 9 taxa at genus level SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = box_plot_bacteria_genus)

ggsave(paste0("Top 9 taxa at genus level SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = box_plot_bacteria_genus)



######################

# Heatmap

# this looks very strange. why is abundance of top families not much more than 1 % ??

meco_dataset <- phyloseq2meco(phylo_SRS_filtered)

t1 <- trans_abund$new(dataset = meco_dataset, taxrank = "Family", ntaxa = 20)
t1$plot_heatmap(facet = "Treatment", xtext_keep = FALSE, withmargin = FALSE)


###############################################

# Relative abundance using non normalised phyloseq object 
# and use function compositional for normalization


#Family_level

pseq.subset_family <- phylo_subset %>% aggregate_taxa(level = "Family") %>% 
  microbiome::transform(transform = "compositional") 

pseq.subset_family

# extracting the OTU table at the family level
OTU4 = as(otu_table(pseq.subset_family), "matrix")  

#if(taxa_are_rows(uBiome)){OTU1 <- t(OTU1)}

# Coerce to data.frame
OTUdf4 = as.data.frame(OTU4)

#exporting to csv
write.csv2(OTUdf4, "relative_abundance_16S_family.csv", row.names = TRUE) 


#for aggregate top taxa function I had to download and load microbiomeutilities package
pseq.subset_family_top20 <- phylo_subset %>% 
  microbiomeutilities::aggregate_top_taxa2(
    level = "Family", top = 20) %>%
  microbiome::transform(transform = "compositional")     

#for having a look on each individual sample in each soil
pseq.subset_family_top20 %>%  plot_composition(group_by = "Treatment")


#colour palette
mypalette = rev(c('blue', 'orange', 'green', 'darkgreen', 'navyblue', 'cyan', 
                  'hotpink','purple', 'burlywood1', 'skyblue', 'gray', 'red', 
                  'yellow', 'darkolivegreen3', 'firebrick4',"#5B2C6F","#000000",
                  "#808000","#D4AC0D","#28B463","#A2D9CE","#FADBD8","#E8DAEF",
                  "#D4E6F1","#A3E4D7","#F9E79F","#F5B041","#FAE5D3","#AEB6BF",
                  "#FEF5E7","#D35400","#797D7F"))

color_21 <- c("#FEF5E7", "#6b8e23", "#191970", "#008080", "#b03060", "#ff0000", "#ff8c00", "#ffff00", 
              "#0000cd", "#7cfc00", "#00fa9a", "#00ffff", "#d8bfd8", "#ff00ff", "#1e90ff", "#eee8aa", "#ee82ee",
              'darkolivegreen3', 'firebrick4',"#5B2C6F","#28B463")

#ploting 10 higly abundand families by soil
relative_abund_family_16S_top20 = pseq.subset_family_top20 %>%  
  plot_composition(average_by = "Treatment",
                   otu.sort = "abundance")+ 
  scale_color_manual(values = color_21, aesthetics = "fill")+
  labs(x = "Treatment", y = "Relative Abundance %") +
  theme_bw(base_size = 16) +
  theme(axis.title = element_text(colour = "black"), 
        axis.text = element_text(colour = "black"), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank")) +
  theme(axis.title.x = element_text(size = 16,
                                    color = "black", 
                                    face = "bold", 
                                    angle = 0)) + 
  theme(axis.text.x = element_text(size = 14, 
                                   color = "black", 
                                   #face = "bold",
                                   angle = 0)) + 
  theme(axis.title.y = element_text(size = 16, 
                                    face = "bold", 
                                    color = "black", 
                                    angle = 90)) + 
  theme(axis.text.y = element_text(size = 14, 
                                   color = "black",
                                   #face = "bold",
                                   angle = 0)) +
  theme(
    legend.title = element_text(size = 10,
                                #face = "bold"
    ), 
    legend.text = element_text(size = 10,
                               #face = "bold"
    ))

relative_abund_family_16S_top20 


ggsave(paste0("relative_abund_family_16S_top20_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = relative_abund_family_16S_top20)

ggsave(paste0("relative_abund_family_16S_top20_SRS", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = relative_abund_family_16S_top20)




#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

## Plotting top 30 abundant taxa with ampvis


#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------


## heatmap using ampvis2 package

library(dplyr)
library(ampvis2)

# Upload OTU table, Abundance table & meta data files

setwd("C:/Users/SarahBenning/OneDrive - Helmholtz Zentrum M?nchen/R79_inoculation/R_data_analysis/amplicon_new/SRS_data/")

data1 = read.csv2("asvtable_removed_contaminants_SRS_singletons_removed.csv")
data2 = read.csv2("ampvis_metadata_table_SRS_singletons_removed.csv")
data3 = read.csv2("taxtable_removed_contaminants_SRS_singletons_removed.csv")

data2$Treatment<-factor(data2$Treatment,levels = c("control", "R6", "R7", "R8", "R9"))


# Create the ampvis2 file object

data_heatmap = amp_load(
  otutable = data1,
  metadata = data2,
  taxonomy = data3,
  pruneSingletons = TRUE,
  removeAbsentOTUs = TRUE)


facet_label_theme <- theme(
  strip.text = element_text(size = 16, 
                            color = "black",
                            face = "bold"),  # Adjust size, color, and style as needed
  strip.background = element_rect(color = "black",
                                  fill = "white", 
                                  linewidth = 1),  # Adjust border color, fill, and linewidth
  strip.placement = "outside"  # Place the facet labels outside the plotting area
)





####subsetting ampvis object if necessary like Fatma did it

#amp_16S_ARD_6d <- amp_subset_samples(
#  amp_16S,
#  Soil == "ARD" & Time_Point == "6",
#  removeAbsentOTUs = FALSE ##when it is enabled it removes OTUs that are absent
#)


#########

## genus level

#########

heatmap_relative_abundance_top_30_genera <-amp_heatmap(
  data_heatmap,
  group_by = "Treatment",
#  facet_by = "Time_Point",
  normalise = TRUE,
  tax_aggregate = "Genus",
  tax_add = "Class",
  tax_show = 30,
  showRemainingTaxa = TRUE,
  tax_class = TRUE,
  tax_empty = "best",
  order_x_by = NULL,
  order_y_by = NULL,
  plot_values = TRUE,
  plot_values_size = 5,
  #plot_legendbreaks = c(1, 5, 15),
  plot_colorscale = "log10",
  plot_na = TRUE,
  measure = "mean",
  min_abundance = 0.1,
  max_abundance = NULL,
  sort_by = NULL,
  normalise_by = NULL,
  scale_by = NULL,
  color_vector = c("royalblue4",
                   "white",
                   "darkred"),
  round = 1,
  textmap = FALSE,
  plot_functions = FALSE,
  function_data = NULL,
  #functions = c("MiDAS", "Filamentous", "AOB", "NOB", "PAO", "GAO"),
  rel_widths = c(0.75, 0.25)
)

heatmap_relative_abundance_top_30_genera




heatmap_relative_abundance_top_30_genera_pretty <- heatmap_relative_abundance_top_30_genera +
  theme_bw(base_size = 13,
  ) +
#  facet_wrap(~ Time_Point, 
 #            labeller = labeller(Time_Point = c("6" = "6 days post inoculation"))) +
  theme(axis.title = element_text(colour = "black"), 
        axis.text = element_text(colour = "black"), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank")) +
  theme(axis.title.x = element_text(size = 16,
                                    color = "black", 
                                    face = "bold", 
                                    family = "Arial",
                                    angle = 0)) + 
  theme(axis.text.x = element_text(size = 14, 
                                   color = "black", 
                                   #face = "bold",
                                   family = "Arial",
                                   angle = 0)) + 
  theme(axis.title.y = element_text(size = 16, 
                                    face = "bold", 
                                    color = "black",
                                    family = "Arial",
                                    angle = 90)) + 
  theme(axis.text.y = element_text(size = 14, 
                                   color = "black",
                                   #face = "bold",
                                   family = "Arial",
                                   angle = 0)) +
  theme(
    legend.title = element_text(size = 16,
                                face = "bold",
                                family = "Arial"
    ), 
    legend.text = element_text(size = 14,
                               # family = "Arial",
                               #face = "bold"
    )
  )  +
  facet_label_theme 

heatmap_relative_abundance_top_30_genera_pretty


ggsave(paste0("heatmap_relative_abundance_top_30_genera_ampvis", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = heatmap_relative_abundance_top_30_genera_pretty)
ggsave(paste0("heatmap_relative_abundance_top_30_genera_ampvis", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = heatmap_relative_abundance_top_30_genera_pretty)



### textmap dataframe for tidyheatmap


textmap_relative_abundance_top_30_genera <-amp_heatmap(
  data_heatmap,
  group_by = "Treatment",
  #  facet_by = "Time_Point",
  normalise = TRUE,
  tax_aggregate = "Genus",
  tax_add = "Class",
  tax_show = 30,
  showRemainingTaxa = TRUE,
  tax_class = TRUE,
  tax_empty = "best",
  order_x_by = NULL,
  order_y_by = NULL,
  plot_values = TRUE,
  plot_values_size = 5,
  #plot_legendbreaks = c(1, 5, 15),
  plot_colorscale = "log10",
  plot_na = TRUE,
  measure = "mean",
  min_abundance = 0.1,
  max_abundance = NULL,
  sort_by = NULL,
  normalise_by = NULL,
  scale_by = NULL,
  color_vector = c("royalblue4",
                   "white",
                   "darkred"),
  round = 1,
  textmap = TRUE,
  plot_functions = FALSE,
  function_data = NULL,
  #functions = c("MiDAS", "Filamentous", "AOB", "NOB", "PAO", "GAO"),
  rel_widths = c(0.75, 0.25)
)

textmap_relative_abundance_top_30_genera


write.csv2(textmap_relative_abundance_top_30_genera, file = "rel_abundance_top_30_genera_ampvis.csv")



########

## family level

########



heatmap_relative_abundance_top_30_families <-amp_heatmap(
  data_heatmap,
  group_by = "Treatment",
  #  facet_by = "Time_Point",
  normalise = TRUE,
  tax_aggregate = "Family",
  tax_add = "Class",
  tax_show = 30,
  showRemainingTaxa = TRUE,
  tax_class = TRUE,
  tax_empty = "best",
  order_x_by = NULL,
  order_y_by = NULL,
  plot_values = TRUE,
  plot_values_size = 5,
  #plot_legendbreaks = c(1, 5, 15),
  plot_colorscale = "log10",
  plot_na = TRUE,
  measure = "mean",
  min_abundance = 0.1,
  max_abundance = NULL,
  sort_by = NULL,
  normalise_by = NULL,
  scale_by = NULL,
  color_vector = c("royalblue4",
                   "white",
                   "darkred"),
  round = 1,
  textmap = FALSE,
  plot_functions = FALSE,
  function_data = NULL,
  #functions = c("MiDAS", "Filamentous", "AOB", "NOB", "PAO", "GAO"),
  rel_widths = c(0.75, 0.25)
)

heatmap_relative_abundance_top_30_families




heatmap_relative_abundance_top_30_families_pretty <- heatmap_relative_abundance_top_30_families +
  theme_bw(base_size = 13,
  ) +
  #  facet_wrap(~ Time_Point, 
  #            labeller = labeller(Time_Point = c("6" = "6 days post inoculation"))) +
  theme(axis.title = element_text(colour = "black"), 
        axis.text = element_text(colour = "black"), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank")) +
  theme(axis.title.x = element_text(size = 16,
                                    color = "black", 
                                    face = "bold", 
                                    family = "Arial",
                                    angle = 0)) + 
  theme(axis.text.x = element_text(size = 14, 
                                   color = "black", 
                                   #face = "bold",
                                   family = "Arial",
                                   angle = 0)) + 
  theme(axis.title.y = element_text(size = 16, 
                                    face = "bold", 
                                    color = "black",
                                    family = "Arial",
                                    angle = 90)) + 
  theme(axis.text.y = element_text(size = 14, 
                                   color = "black",
                                   #face = "bold",
                                   family = "Arial",
                                   angle = 0)) +
  theme(
    legend.title = element_text(size = 16,
                                face = "bold",
                                family = "Arial"
    ), 
    legend.text = element_text(size = 14,
                               # family = "Arial",
                               #face = "bold"
    )
  )  +
  facet_label_theme 

heatmap_relative_abundance_top_30_families_pretty


ggsave(paste0("heatmap_relative_abundance_top_30_families_ampvis", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = heatmap_relative_abundance_top_30_families_pretty)
ggsave(paste0("heatmap_relative_abundance_top_30_families_ampvis", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = heatmap_relative_abundance_top_30_families_pretty)



### textmap dataframe for tidyheatmap


textmap_relative_abundance_top_30_families <-amp_heatmap(
  data_heatmap,
  group_by = "Treatment",
  #  facet_by = "Time_Point",
  normalise = TRUE,
  tax_aggregate = "Family",
  tax_add = "Class",
  tax_show = 30,
  showRemainingTaxa = TRUE,
  tax_class = TRUE,
  tax_empty = "best",
  order_x_by = NULL,
  order_y_by = NULL,
  plot_values = TRUE,
  plot_values_size = 5,
  #plot_legendbreaks = c(1, 5, 15),
  plot_colorscale = "log10",
  plot_na = TRUE,
  measure = "mean",
  min_abundance = 0.1,
  max_abundance = NULL,
  sort_by = NULL,
  normalise_by = NULL,
  scale_by = NULL,
  color_vector = c("royalblue4",
                   "white",
                   "darkred"),
  round = 1,
  textmap = TRUE,
  plot_functions = FALSE,
  function_data = NULL,
  #functions = c("MiDAS", "Filamentous", "AOB", "NOB", "PAO", "GAO"),
  rel_widths = c(0.75, 0.25)
)

textmap_relative_abundance_top_30_families


write.csv2(textmap_relative_abundance_top_30_families, file = "rel_abundance_top_30_families_ampvis.csv")


#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------


#  trying tidy heatmap on taxonomy data, to get clustering of treatments


library(dbplyr)
library(tidyr)
library(tidyverse)
library(reshape2)
library(tidyHeatmap)


# use exported csv from textmap step, rename 1. column as "Taxon", capitalize control to Control


data <- read.csv2("rel_abundance_top_30_families_ampvis.csv", check.names = FALSE)
#data <- read.csv2("rel_abundance_top_30_genera_ampvis.csv", check.names = FALSE)

data$Taxon <- as.factor(data$Taxon)

#data_1 <- data[data$Taxon != "Remaining taxa (683)", ]
data_1 <- data[data$Taxon != "Remaining taxa (355)", ]

data_2 <- melt(data_1)

data_2 <- as_tibble(data_2)




my_palette =  grDevices::colorRampPalette(c("lightgray", "#255668", "#0e7279", "#0e8f81", "#3aab7e", "#6ec574", "#a9dc67", "#edef5d"))

data_2 |> 
  heatmap(Taxon, variable, value,  scale = "none" ,
          palette_value = circlize::colorRamp2(c(0, seq(1, 7, length.out = 6)), my_palette(7)),
          rect_gp = grid::gpar(col = "#161616", lwd = 0.5))|>
  as_ComplexHeatmap() |>
  ComplexHeatmap::draw(
    column_title = "Relative abundance top 30 Families"#, 
    #  column_title_gp = gpar(fontsize = 16)
  )




#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

## CORE MICROBIOME ANALYSIS


#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------

library(dplyr)
library(ampvis2)

# 16S Bacterial core microbiome

# Upload OTU table, Abundance table & meta data files

setwd("C:/Users/SarahBenning/OneDrive - Helmholtz Zentrum M?nchen/R79_inoculation/R_data_analysis/amplicon_new/SRS_data/")

data1 = read.csv2("asvtable_removed_contaminants_SRS_singletons_removed.csv")
data2 = read.csv2("ampvis_metadata_table_SRS_singletons_removed.csv")
data3 = read.csv2("taxtable_removed_contaminants_SRS_singletons_removed.csv")

#data2$Soil<-factor(data2$Soil,levels = c("M-U", "C-U", "Ha-R-U", "L-U", "ST-U"))


# Create the ampvis2 file object

data = amp_load(
  otutable = data1,
  metadata = data2,
  taxonomy = data3,
  pruneSingletons = TRUE,
  removeAbsentOTUs = TRUE)


# Compute the venn diagram using reproducibly ASV in 5 out 6 replicates (in 80% of samples = cut_f)

venn_all_80_0.1 <- amp_venn(
  data = data,
  # group_by = "genotypes", # you can group with this venn diagram in max 3 levels of your variable, so in this case we don't group in the diagram, but still get the core microbiome
  cut_a = 0.1, # be aware of abundance cutoff 0.1 or 0.01, depends on what you want to look at
  cut_f = 80,
  text_size = 5,
  normalise = TRUE,
  detailed_output = TRUE
)
venn_all_80_0.1

# Export results as data frame table, this step is not necessary unless you wish it

write.csv2(data.frame(venn_all_80_0.1$Otutable), "core_16S_all_80_0.1.csv")

# write.csv2(data.frame(venn_all_80_0.1$Otutable), "core_ITS_all_80_0.1.csv") # only do in ITS step


#Decipher the core microbiome composition
#Create a data frame of the OTUtable

df_all_80_0.1 <- venn_all_80_0.1$Otutable

#Keep rows that are only the core (filter only core microbiome),so it will discard the unique ASVs. Unique ASVs are the one influenced by your treatments
df_core_all_80_0.1 <- df_all_80_0.1[df_all_80_0.1$Shared == "Core", ]

#Remove the last column == Shared
df_core_all_80_0.1 <- select(df_core_all_80_0.1, -Shared)


#rename row headers
df_core_all_80_0.1 = df_core_all_80_0.1 %>%
  rename(OTU2 = OTU
  )

# To merge the core microbiome table with the original abundance table, 
# the header ASV should be replaced by OTU. If, in the abundance table used to create the ampvis2 file object, you already named the first row "OTU", then no need to run this line code.
data1 = data1 %>%
  rename(
    OTU = ASV,
  )

# This code allow you to filter the samples used for the core microbiome with their original values 
# of the abundance as in the core venn diagram, the abundance values were normalized 
# (transformed in percentage)
# grep is used to select the samples ID, I named my samples like this:
# soil_B1, soil_B2,etc, so all samples have in common soil_.


heat_all_80_0.1 = 
  data1 %>%
  filter(OTU %in% df_core_all_80_0.1$OTU2) %>%
  select(c(OTU, names(df_core_all_80_0.1) %>% grep("sarah", ., value = T))) %>%
  right_join(df_core_all_80_0.1[, 1:8], by = c("OTU" = "OTU2"))

write.csv2(heat_all_80_0.1, "core_80_0.1_abundance.csv")

# Once you have the core microbiome samples with the original abundance values, you can create a second ampvis2 file object.
# heat is the OTU abundance table of samples that are only found in the core microbiome
# metadata does not change (so you have the same samples ID, treatments etc)

data_heat_all_80_0.1 <- amp_load(
  heat_all_80_0.1,
  metadata = data2,
  taxonomy = NULL,
  fasta = NULL,
  tree = NULL,
  pruneSingletons = TRUE)

# Plot the core microbiome composition and their relative abundance, you can choose which level (Phylum, Family, Genus etc) you want. 
# Here I plotted at the genus level, and I added which class those genera belonged to, using the argmuent tax_add
# You can read more about this amp_heatmap in the help desk of R studio. 


######## heatmap all 80%, cut_a = 0.1

## plot the heatmap

plot_all_80_0.1 <-amp_heatmap(
  data_heat_all_80_0.1,
  group_by = "Treatment",
  normalise = TRUE,
  tax_aggregate = "Family",
  tax_add = "Class",
  tax_show = 28,
  showRemainingTaxa = TRUE,
  tax_class = TRUE,
  tax_empty = "best",
  order_x_by = NULL,
  order_y_by = NULL,
  plot_values = TRUE,
  plot_values_size = 5,
  # plot_legendbreaks = c(1, 7.5, 15),
  plot_colorscale = "log10",
  plot_na = TRUE,
  measure = "mean",
  min_abundance = 0.1,
  max_abundance = NULL,
  sort_by = NULL,
  normalise_by = NULL,
  scale_by = NULL,
  color_vector = c("royalblue4",
                   "whitesmoke",
                   "darkred"),
  round = 1,
  textmap = FALSE,
  plot_functions = FALSE,
  function_data = NULL,
  functions = c("MiDAS", "Filamentous", "AOB", "NOB", "PAO", "GAO"),
  rel_widths = c(0.75, 0.25)
)

plot_all_80_0.1

ggsave(paste0("heat_core_genus_all_80_0.1", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = plot_all_80_0.1)
ggsave(paste0("heat_core_genus_all_80_0.1", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = plot_all_80_0.1)

# ggsave(paste0("heat_core_family_all_80_0.1", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = plot_all_80_0.1)
# ggsave(paste0("heat_core_family_all_80_0.1", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = plot_all_80_0.1)

## give heatmap values as dataframe (eg. for making other heatmaps with tidyheatmap)

textmap_all_80_0.1 <-amp_heatmap(
  data_heat_all_80_0.1,
  group_by = "Treatment",
  normalise = TRUE,
  tax_aggregate = "Genus",
  tax_add = "Class",
  tax_show = 70,
  showRemainingTaxa = TRUE,
  tax_class = TRUE,
  tax_empty = "best",
  order_x_by = NULL,
  order_y_by = NULL,
  plot_values = TRUE,
  plot_values_size = 5,
  # plot_legendbreaks = c(1, 7.5, 15),
  plot_colorscale = "log10",
  plot_na = TRUE,
  measure = "mean",
  min_abundance = 0.1,
  max_abundance = NULL,
  sort_by = NULL,
  normalise_by = NULL,
  scale_by = NULL,
  color_vector = c("royalblue4",
                   "whitesmoke",
                   "darkred"),
  round = 1,
  textmap = TRUE,
  plot_functions = FALSE,
  function_data = NULL,
  functions = c("MiDAS", "Filamentous", "AOB", "NOB", "PAO", "GAO"),
  rel_widths = c(0.75, 0.25)
)
textmap_all_80_0.1

write.csv2(textmap_all_80_0.1, file = "genus_heat_core_all_80_0.1.csv")
# write.csv2(textmap_all_80_0.1, file = "ITS_heat_core_all_80_0.1.csv")


## run everything again with ITS as data (make sure to also export the correct filenames)


##-------------------------------------------------------------------------------------------------------------------------------

##-------------------------------------------------------------------------------------------------------------------------------



#  trying tidy heatmap on  core microbiome data, to get clustering of treatments

#  17.01.2023, Sarah Benning

library(dbplyr)
library(tidyr)
library(tidyverse)
library(reshape2)
library(tidyHeatmap)


#16S

# use exported csv from last core microbiome step

# beware abundance cotoff of 0.01 or 0.1

#data <- read.csv2("16S_heat_core_all_80_0.1_SRS.csv", check.names = FALSE)
data <- read.csv2("genus_heat_core_all_80_0.1.csv", check.names = FALSE)

data$Genus <- as.factor(data$Genus)

data_1 <- melt(data)

data_1 <- as_tibble(data_1)




my_palette =  grDevices::colorRampPalette(c("lightgray", "#255668", "#0e7279", "#0e8f81", "#3aab7e", "#6ec574", "#a9dc67", "#edef5d"))

data_1 |> 
  heatmap(Genus, variable, value,  scale = "none" ,
          palette_value = circlize::colorRamp2(c(0, seq(1, 15, length.out = 14)), my_palette(15)),
          rect_gp = grid::gpar(col = "#161616", lwd = 0.5))|>
  as_ComplexHeatmap() |>
  ComplexHeatmap::draw(
    column_title = "Bacterial Core Microbiome 16S"#, 
    #  column_title_gp = gpar(fontsize = 16)
  )

