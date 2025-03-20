## Differential abundance Analysis using DeSeq2
## Inoculation concentration Experiment
## Sarah Benning
## R Version 4.3.1


##Load libraries needed
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(microbiome)
packageVersion("DESeq2") #1.42.0

#set seed for reproducibility
set.seed(713)


load("C:/Users/SarahBenning/OneDrive - Helmholtz Zentrum München/R79_inoculation/R_data_analysis/amplicon_new/Diff_abundance/uBiome_for_DA_200324_1335.RData")

meta <- meta(uBiome_DA)

taxonomy <- data.frame(tax_table(uBiome_DA))

sample_sums(uBiome_DA)

#transform to relative abundance for the boxplots

ps_rel <- transform_sample_counts(uBiome_DA, function (x) 100*x/sum(x))


# Differential Abundance - genus level
##convert the phyloseq object into a DESeq2 object
##we want to include the major factor for the dataset: Treatment 

##note: variable of interest should be at the last position in the design formula if you have more than one

#do on ASV level

dds_ASV <- phyloseq_to_deseq2(uBiome_DA, ~ Treatment)
dds_ASV$Treatment <- relevel(dds_ASV$Treatment, ref="control")


#run DESeq
#i followed the recommendations of the DESeq2 vignette -> LRT is recommended for testing multiple terms with 3 or more levels of a factor at once
# wald for a single parameter.
#all other settings were kept as default expect for the fitType, since the curve doesn't fit the observed dispersion. 

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds_ASV), 1, gm_mean)


dds_ASV <- estimateSizeFactors(dds_ASV, type="poscounts")
dds_ASV <- DESeq(dds_ASV, test = "Wald", full = ~ Treatment,  fitType = "local")

plotDispEsts(dds_ASV)

resultsNames(dds_ASV)
alpha = 0.05

#extract results from full model (seems to be the more reliable option)
#start with R6

res_R6_ASV <- results(dds_gen, contrast = list(c("Treatment_R6_vs_control")), alpha = alpha)
summary(res_R6_ASV)
res_R6_ASV <- subset(res_R6_ASV, padj < 0.05)
res_R6_ASV <- data.frame(res_R6_ASV)
res_R6_ASV <- merge(res_R6_ASV, taxonomy, by=0) #add taxa names
R6_ASV <- res_R6_ASV[(res_R6_ASV$log2FoldChange <= -1) | (res_R6_ASV$log2FoldChange >=1),]

#to control if there are some false positives
R6_ASV_da_tax <- R6_ASV$Row.names

ps_R6_ASV <- prune_taxa(taxa_names(ps_rel) %in% R6_ASV_da_tax, ps_rel)

ps_R6_ASV <- psmelt(ps_R6_ASV)

DA_R6_ASV <- ggplot(ps_R6_ASV, aes(x = ps_R6_ASV$Treatment, y = ps_R6_ASV$Abundance, col = ps_R6_ASV$Treatment)) + 
  geom_boxplot() + 
  facet_wrap(~ OTU, scales="free") 


#extract results from full model (seems to be the more reliable option)
#R7

res_R7_ASV <- results(dds_gen, contrast = list(c("Treatment_R7_vs_control")), alpha = alpha)
summary(res_R7_ASV)
res_R7_ASV <- subset(res_R7_ASV, padj < 0.05)
res_R7_ASV <- data.frame(res_R7_ASV)
res_R7_ASV <- merge(res_R7_ASV, taxonomy, by=0) #add taxa names
R7_ASV <- res_R7_ASV[(res_R7_ASV$log2FoldChange <= -1) | (res_R7_ASV$log2FoldChange >=1),]

#to control if there are some false positives
R7_ASV_da_tax <- R7_ASV$Row.names

ps_R7_ASV <- prune_taxa(taxa_names(ps_rel) %in% R7_ASV_da_tax, ps_rel)

ps_R7_ASV <- psmelt(ps_R7_ASV)

DA_R7_ASV <- ggplot(ps_R7_ASV, aes(x = ps_R7_ASV$Treatment, y = ps_R7_ASV$Abundance, col = ps_R7_ASV$Treatment)) + 
  geom_boxplot() + 
  facet_wrap(~ OTU, scales="free") 


#extract results from full model (seems to be the more reliable option)
#R8

res_R8_ASV <- results(dds_gen, contrast = list(c("Treatment_R8_vs_control")), alpha = alpha)
summary(res_R8_ASV)
res_R8_ASV <- subset(res_R8_ASV, padj < 0.05)
res_R8_ASV <- data.frame(res_R8_ASV)
res_R8_ASV <- merge(res_R8_ASV, taxonomy, by=0) #add taxa names
R8_ASV <- res_R8_ASV[(res_R8_ASV$log2FoldChange <= -1) | (res_R8_ASV$log2FoldChange >=1),]

#to control if there are some false positives
R8_ASV_da_tax <- R8_ASV$Row.names

ps_R8_ASV <- prune_taxa(taxa_names(ps_rel) %in% R8_ASV_da_tax, ps_rel)

ps_R8_ASV <- psmelt(ps_R8_ASV)

DA_R8_ASV <- ggplot(ps_R8_ASV, aes(x = ps_R8_ASV$Treatment, y = ps_R8_ASV$Abundance, col = ps_R8_ASV$Treatment)) + 
  geom_boxplot() + 
  facet_wrap(~ OTU, scales="free") 



#extract results from full model (seems to be the more reliable option)
#R9

res_R9_ASV <- results(dds_gen, contrast = list(c("Treatment_R9_vs_control")), alpha = alpha)
summary(res_R9_ASV)
res_R9_ASV <- subset(res_R9_ASV, padj < 0.05)
res_R9_ASV <- data.frame(res_R9_ASV)
res_R9_ASV <- merge(res_R9_ASV, taxonomy, by=0) #add taxa names
R9_ASV <- res_R9_ASV[(res_R9_ASV$log2FoldChange <= -1) | (res_R9_ASV$log2FoldChange >=1),]

#to control if there are some false positives
R9_ASV_da_tax <- R9_ASV$Row.names

ps_R9_ASV <- prune_taxa(taxa_names(ps_rel) %in% R9_ASV_da_tax, ps_rel)

ps_R9_ASV <- psmelt(ps_R9_ASV)

DA_R9_ASV <- ggplot(ps_R9_ASV, aes(x = ps_R9_ASV$Treatment, y = ps_R9_ASV$Abundance, col = ps_R9_ASV$Treatment)) + 
  geom_boxplot() + 
  facet_wrap(~ OTU, scales="free") 




#for plotting, all diff. abundant asv of all three treatments must be combined into one data.frame with an additional column called treatment.
#add treatment column

Treatment <- rep("R6", 19)
R6_ASV$Treatment <- Treatment

Treatment <- rep("R7", 7)
R7_ASV$Treatment <- Treatment

Treatment <- rep("R8", 13)
R8_ASV$Treatment <- Treatment

Treatment <- rep("R9", 25)
R9_ASV$Treatment <- Treatment

#modify label column -> combine family & genus
cols <- c("Phylum", "Family", "Genus", "Row.names")

R6_ASV$Label <- apply(R6_ASV[, cols], 1, paste, collapse=";")
R7_ASV$Label <- apply(R7_ASV[, cols], 1, paste, collapse=";")
R8_ASV$Label <- apply(R8_ASV[, cols], 1, paste, collapse=";")
R9_ASV$Label <- apply(R9_ASV[, cols], 1, paste, collapse=";")

#save full dataframe as excel

diff_abu_gen_deseq_ASV <- rbind(R6_ASV, R7_ASV, R8_ASV, R9_ASV)

write.csv2(diff_abu_gen_deseq_ASV, file ="DiffAbundance_ASV.csv", row.names = F)


############################################################################

#bar chart
#remove unclassified on order level to reduce dimensions

difftax_ASV <- read.csv2("DiffAbundance_ASV.csv")

#exclude FP Taxa

difftax_ASV <- subset(difftax_ASV, Row.names != "ASV353")


#difftax_ad_filt <- subset(difftax_ad, !Order == "unclassified")
#difftax_juv_filt <- subset(difftax_juv, !Order == "unclassified")

#adjust Timepoint labeling

#difftax_juv_filt$Timepoint <- gsub("t10", "Treatment", difftax_juv_filt$Timepoint)
#difftax_juv_filt$Timepoint <- gsub("t14", "4d post-treatment", difftax_juv_filt$Timepoint)
#difftax_juv_filt$Timepoint <- gsub("t21", "11d post-treatment", difftax_juv_filt$Timepoint)
#difftax_juv_filt$Timepoint <- gsub("t28", "18d post-treatment", difftax_juv_filt$Timepoint)

#set order for labeling
#difftax_ad_filt$Timepoint <- factor(difftax_ad_filt$Timepoint, levels=c("t0", "Treatment", "4d post-treatment", "11d post-treatment", "18d post-treatment"))

#adjust labeling, make shorter
label_short <- c("Genus", "Row.names")
difftax_ASV$Label_short <- apply(difftax_ASV[, label_short], 1, paste, collapse="; ")

# Ändere die Reihenfolge der Faktorstufen von A nach Z
difftax_ASV$Label_short <- factor(difftax_ASV$Label_short, levels = unique(difftax_ASV$Label_short))

difftax_ASV_ord <- difftax_ASV %>%
  mutate(Label_short = factor(Label_short, levels = Label_short[order(log2FoldChange)]))  # Order Label by logchange

difftax_ASV_bar <- ggplot(data = difftax_ASV, 
                          aes(x = factor(Label_short, levels = rev(unique(Label_short))),
                              y = log2FoldChange, fill = Treatment)) 




difftax_ASV_bar <- ggplot(data = difftax_ASV, 
                          aes(x = reorder(Label_short, -log2FoldChange),  # Reorder based on log2FoldChange
                              y = log2FoldChange, fill = Treatment))+
  geom_bar(position = "dodge", stat = "identity", width = 0.7, color = "black", alpha = 0.8) +
  scale_fill_manual(values = c(#"#90EE90", 
                               "#228822", "#5F7FC7", "orange","#DA5724"))+
  theme_classic() +
  labs(x = "Taxa", y = "log2FoldChange") +
  theme(
    axis.text.y = element_text(face = "italic", family = "Arial", color = "black", size=12),
    axis.text.x = element_text(family = "Arial", color = "black"), 
    axis.title = element_text(family = "Arial", color = "black"),    
    legend.text = element_text(family = "Arial", color = "black"),
    legend.title = element_text(family = "Arial", color = "black", face="bold"),
    strip.text = element_text(family = "Arial", color = "black", face="bold")     
  ) +
  guides(fill = guide_legend(title = "Treatment")) +
#  scale_y_continuous(limits = c(-30, 30), breaks = seq(-30, 30, by = 10)) +
  coord_flip() +
#  facet_grid(~ Treatment) +
  geom_vline(xintercept = seq(0.5, length(unique(difftax_ASV$Label)) - 0.5, by = 1), color = "gray", linetype = "dotted")

difftax_ASV_bar


ggsave(paste0("Diff_abundance_barplot_ASV_sorted", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = difftax_ASV_bar)

ggsave(paste0("Diff_abundance_barplot_ASV_sorted", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = difftax_ASV_bar)




difftax_ASV_bar_purple <- ggplot(data = difftax_ASV, 
                          aes(x = reorder(Label_short, -log2FoldChange),  # Reorder based on log2FoldChange
                              y = log2FoldChange, fill = Treatment))+
  geom_bar(position = "dodge", stat = "identity", width = 0.7, color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("#b3cde3", "#8c96c6", "#8856a7","#810f7c"))+
  theme_classic() +
  labs(x = "Taxa", y = "log2FoldChange") +
  theme(
    axis.text.y = element_text(face = "italic", family = "Arial", color = "black", size=12),
    axis.text.x = element_text(family = "Arial", color = "black"), 
    axis.title = element_text(family = "Arial", color = "black"),    
    legend.text = element_text(family = "Arial", color = "black"),
    legend.title = element_text(family = "Arial", color = "black", face="bold"),
    strip.text = element_text(family = "Arial", color = "black", face="bold")     
  ) +
  guides(fill = guide_legend(title = "Treatment")) +
  #  scale_y_continuous(limits = c(-30, 30), breaks = seq(-30, 30, by = 10)) +
  coord_flip() +
  #  facet_grid(~ Treatment) +
  geom_vline(xintercept = seq(0.5, length(unique(difftax_ASV$Label)) - 0.5, by = 1), color = "gray", linetype = "dotted")

difftax_ASV_bar_purple

difftax_ASV_bar_rgsafe <- ggplot(data = difftax_ASV, 
                                 aes(x = reorder(Label_short, -log2FoldChange),  # Reorder based on log2FoldChange
                                     y = log2FoldChange, fill = Treatment))+
  geom_bar(position = "dodge", stat = "identity", width = 0.7, color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("#fdae61", "#ffffbf", "#abd9e9","#2c7bb6"))+
  theme_classic() +
  labs(x = "Taxa", y = "log2FoldChange") +
  theme(
    axis.text.y = element_text(face = "italic", family = "Arial", color = "black", size=12),
    axis.text.x = element_text(family = "Arial", color = "black"), 
    axis.title = element_text(family = "Arial", color = "black"),    
    legend.text = element_text(family = "Arial", color = "black"),
    legend.title = element_text(family = "Arial", color = "black", face="bold"),
    strip.text = element_text(family = "Arial", color = "black", face="bold")     
  ) +
  guides(fill = guide_legend(title = "Treatment")) +
  #  scale_y_continuous(limits = c(-30, 30), breaks = seq(-30, 30, by = 10)) +
  coord_flip() +
  #  facet_grid(~ Treatment) +
  geom_vline(xintercept = seq(0.5, length(unique(difftax_ASV$Label)) - 0.5, by = 1), color = "gray", linetype = "dotted")

difftax_ASV_bar_rgsafe



ggsave(paste0("Diff_abundance_barplot_ASV_sorted_blind", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = difftax_ASV_bar_rgsafe)

ggsave(paste0("Diff_abundance_barplot_ASV_sorted_blind", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = difftax_ASV_bar_rgsafe)


