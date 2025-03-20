## correlation of inoculation concentration level with ASVs

## Inoculation concentration Experiment
## Sarah Benning
## R Version 4.3.1

## first subset SRS data according to treatment

## filter for ASVs present in ALL samples

## export on long format with psmelt

## calculate an average abundance per replicate

## calculate correlation with cor and p values with cor.mtest


library(phyloseq)   
library(microbiome)
library(metagMisc)
library(dplyr)



all_long_filtered_100 <- all_long %>%
  group_by(OTU) %>%  # OTU
  filter(all(Abundance != 0)) %>%  # filter out everything that is 0
  ungroup()


all_long_filtered_100 <- all_long_filtered_100 %>%
  select(-Phytoalexins, -Shoot_growth, -Root_dry_weight)



all_long_filtered_100$Load <- with(all_long_filtered_100, ifelse(Treatment == "control", 0,
                                                                 ifelse(Treatment == "R6", 1000000,
                                                                        ifelse(Treatment == "R7", 10000000,
                                                                               ifelse(Treatment == "R8", 100000000, 
                                                                                      ifelse(Treatment == "R9", 1000000000, NA))))))




all_long_minimal_100 <- all_long_filtered_100 %>%
  select(-sample_Sample, -Domain, -Phylum, -Class, -Order, -Family, -Genus, -Species)


all_flipped_100 <- all_long_minimal_100 %>%
  pivot_wider(names_from = OTU, values_from = Abundance)


all_for_cor_100 <- all_flipped_100 %>%
  select(-Sample, -Treatment)


cor_matrix_replicates_100 <- cor(all_for_cor_100, method = "spearman")


cor_p_matrix_replicates_100 <- cor.mtest(all_for_cor_100, conf.level = 0.95)


setwd("C:/Users/SarahBenning/OneDrive - Helmholtz Zentrum München/R79_inoculation/R_data_analysis/amplicon_new/SRS_data/correlation/")

write.csv2(cor_matrix_replicates_100, "correlation_matrix_ASV_in_all_samples.csv", row.names = TRUE)

write.csv2(cor_p_matrix_replicates_100, "correlation_pvalue_matrix_ASV_in_all_samples.csv", row.names = TRUE)


################################

## visualisation


# Step 1: Extract the taxonomy table
taxonomy_table <- as.data.frame(tax_table(phylo_SRS_filtered))

# Step 2: Convert row names to a column for merging
taxonomy_table$ASV <- rownames(taxonomy_table)
genus_table <- taxonomy_table %>%
  select(ASV, Genus)

genus_table <- genus_table %>%
  mutate(name = paste(Genus, ASV, sep = "_"))

cor_mat_all <- read.csv2(
  "C:/Users/SarahBenning/OneDrive - Helmholtz Zentrum München/R79_inoculation/R_data_analysis/amplicon_new/SRS_data/correlation/correlation_matrix_ASV_in_all_samples.csv", 
  row.names = NULL)

colnames(cor_mat_all)[colnames(cor_mat_all) == "X"] <- "ASV"  # Rename "X" to "ASV"


# Step 3: Create a data frame from cor_mat_all
#cor_mat_all_df <- as.data.frame(cor_mat_all)
#cor_mat_all_df$ASV <- rownames(cor_mat_all_df)

# Step 4: Merge the data frames
merged_data <- merge(cor_mat_all, genus_table[, c("ASV", "name")], by = "ASV", all.x = TRUE)


write.csv2(merged_data, "correlation_matrix_genus_name_all.csv", row.names = TRUE)


# Step 5: Set row names back if needed
#rownames(merged_data) <- merged_data$ASV
#merged_data$ASV <- NULL  # Remove ASV column if you want to keep only the data

# View the merged data
#print(merged_data)

## load dataframe of ASVs were p < 0.05 of cor_mat_all_samples


cor_sig_all <- read.csv2(
  "C:/Users/SarahBenning/OneDrive - Helmholtz Zentrum München/R79_inoculation/R_data_analysis/amplicon_new/SRS_data/correlation/cor_mat_all_genus_signif.csv")


library(ggplot2)

cor_sig_all <- cor_sig_all %>%
  mutate(Name = factor(Name, levels = Name[order(Load)]))  # Order Name by Load

plot <- ggplot(cor_sig_all, aes(x = Name, y = Load)) +
        geom_bar(stat = "identity", width = 0.3, color = "black", fill = "gray90") +
        coord_flip() +
        theme_minimal(base_size = 14) +
        labs(x = "Name", y = "Correlation") #+
       # scale_y_continuous(breaks = seq(-0.6, 0.6, by = 0.2))#+
       # theme(axis.text.x = element_text(angle = 90, hjust = 1)
              


ggsave(paste0("cor_sig_all", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = plot)

ggsave(paste0("cor_sig_all", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = plot)

###############################

#partial correlation


library(corpcor)


 


all_for_cor_trans =
  all_for_cor_100[,2:147] %>% 
  log() %>% 
  sweep(1, apply(., 1, mean)) %>% 
  cbind(all_for_cor_100[,1])


all_for_cor_100 %>% 
  cor() %>% 
  corpcor::cor2pcor() %>% 
  as.data.frame() %>%
  `colnames<-`(names(all_for_cor_100)) %>% 
  `rownames<-`(names(all_for_cor_100)) %>% 
  rownames_to_column("from") %>%
  pivot_longer(cols = -from, names_to = "to",
               values_to = "pcor") %>%
  filter(str_detect(from, "Load"),
         abs(pcor) > 0.3) %>%
  left_join(tax_table(phylo_SRS_filtered) %>% 
              as.data.frame() %>% 
              rownames_to_column("to")) %>% 
  ggplot() +
  geom_col(aes(to,  pcor, fill = Phylum))




partial_spearman <- pcor(all_for_cor_trans)


partial_cor_100 <- all_for_cor_100 %>% 
  cor() %>% 
  corpcor::cor2pcor() %>% 
  as.data.frame() %>%
  `colnames<-`(names(all_for_cor_100)) %>% 
  `rownames<-`(names(all_for_cor_100)) %>% 
  rownames_to_column("from")


setwd("C:/Users/SarahBenning/OneDrive - Helmholtz Zentrum München/R79_inoculation/R_data_analysis/amplicon_new/SRS_data/correlation/")

write.csv2(partial_cor_100, "partial_correlation_matrix_ASV_in_all_samples.csv", row.names = TRUE)


#######

## visualisation
library(dplyr)

# Step 1: Extract the taxonomy table
taxonomy_table <- as.data.frame(tax_table(phylo_SRS_filtered))

# Step 2: Convert row names to a column for merging
taxonomy_table$ASV <- rownames(taxonomy_table)

# from here its not working anymore ...
genus_table <- taxonomy_table %>%
  select(ASV, Genus, Phylum)

genus_table <- genus_table %>%
  mutate(Name = paste(Genus, ASV, sep = "_"))

par_cor_mat_all <- read.csv2(
  "C:/Users/SarahBenning/OneDrive - Helmholtz Zentrum München/R79_inoculation/R_data_analysis/amplicon_new/SRS_data/correlation/partial_cor_long_greater0.3.csv", 
  row.names = NULL)

#colnames(cor_mat_all)[colnames(cor_mat_all) == "X"] <- "ASV"  # Rename "X" to "ASV"



# Step 4: Merge the data frames
merged_data <- merge(par_cor_mat_all, genus_table[, c("ASV", "Name", "Phylum")], by = "ASV", all.x = TRUE)


write.csv2(merged_data, "partial_correlation_matrix_genus_name_all_0.3.csv", row.names = TRUE)


library(ggplot2)

merged_data <- merged_data %>%
  mutate(Name = factor(Name, levels = Name[order(Load)]))  # Order Name by Load

plot <- ggplot(merged_data, aes(x = Name, y = Load, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.3, color = "black") +
  scale_fill_brewer(palette = "Paired") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(x = "Name", y = "Correlation") #+
# scale_y_continuous(breaks = seq(-0.6, 0.6, by = 0.2))#+
# theme(axis.text.x = element_text(angle = 90, hjust = 1)



ggsave(paste0("partial_cor_0.3", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = plot)

ggsave(paste0("partial_cor_0.3", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = plot)

###############################