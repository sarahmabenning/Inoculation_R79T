##################################################################

## 16S qPCR on rhizosphere microbiome of Inoculation concentration Experiment 2021

##################################################################


library("ggplot2")
library("reshape2")
library("ggpubr")
library("gridExtra")
library("grid")
library("tidyverse")
library("RColorBrewer")
library("hues")
library("ggsci")
library("svglite")


setwd("C:/Users/SarahBenning/OneDrive - Helmholtz Zentrum München/R79_inoculation/cDNA_qPCR/checked/")

#data_all <- read.csv2("all_values_qPCR.csv")
data_means <- read.csv2("mean_qPCR.csv")


#str(data_all)

#data_all$Sample <- as.factor(data_all$Sample)
#data_all$DNA <- as.factor(data_all$DNA)
#data_all$Treatment <- as.factor(data_all$Treatment)
#data_all$Gene <- as.factor(data_all$Gene)

str(data_means)

data_means$DNA <- as.factor(data_means$DNA)
data_means$Treatment <- as.factor(data_means$Treatment)
data_means$Gene <- as.factor(data_means$Gene)
data_means$qPCR <- as.factor(data_means$qPCR)


subset_16S <- subset(data_means, Gene == "16S")
subset_bphd <- subset(data_means, Gene == "bphd")
subset_29_11 <- subset(data_means, qPCR == "29_11_23")
subset_28_11 <- subset(data_means, qPCR == "28_11_23")

subset_29_11_DNA <- subset(subset_29_11, DNA == "DNA")
subset_28_11_DNA <- subset(subset_28_11, DNA == "DNA")

## check the different bphd runs

data_plot <- subset_bphd

qPCR_bphd_ng_DNA <- ggplot(data_plot, 
                      aes( x = Treatment, 
                           y = mean_copy_ng_DNA,
                           fill = qPCR)) +
  geom_bar(position = position_dodge(preserve = "single"), 
           stat = "identity",
           color = "black")+
  geom_errorbar(aes(ymin = mean_copy_ng_DNA, 
                    ymax = mean_copy_ng_DNA + stdev_copy_ng_DNA), #keep only upper error bars
                width = 0.2,
                position = position_dodge(0.9, preserve = "single")) +
  scale_fill_manual(values = c("darkcyan", "rosybrown3", "darkblue"),
                    breaks = c("28_11_23", "29_11_23", "05_03_24")) +
  #  coord_cartesian( ylim = c(0, 50000000) )+
  labs(x = "Strain", y = "copies/ ng DNA") +
  theme_classic(base_size = 16) +
  theme(axis.title = element_text(colour = "black"), 
        axis.text.x =  element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black")) +
  scale_y_continuous(labels = scales::scientific) + #scientific notation
  facet_wrap(~ Gene, scales = "free")


qPCR_bphd_ng_DNA


setwd("C:/Users/SarahBenning/OneDrive - Helmholtz Zentrum München/R79_inoculation/cDNA_qPCR")

ggsave(paste0("qPCR_bphd_ng_DNA", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = qPCR_bphd_ng_DNA)

ggsave(paste0("qPCR_bphd_ng_DNA", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = qPCR_bphd_ng_DNA)


qPCR_bphd_g_soil <- ggplot(data_plot, 
                      aes( x = Treatment, 
                           y = mean_copy_g_soil,
                           fill = qPCR)) +
  geom_bar(position = position_dodge(preserve = "single"), 
           stat = "identity",
           color = "black")+
  geom_errorbar(aes(ymin = mean_copy_g_soil, 
                    ymax = mean_copy_g_soil + stdev_copy_g_soil), #keep only upper error bars
                width = 0.2,
                position = position_dodge(0.9, preserve = "single")) +
  scale_fill_manual(values = c("darkcyan", "rosybrown3", "darkblue"),
                    breaks = c("28_11_23", "29_11_23", "05_03_24")) +
  #  coord_cartesian( ylim = c(0, 50000000) )+
  labs(x = "Strain", y = "copies/ g soil") +
  theme_classic(base_size = 16) +
  theme(axis.title = element_text(colour = "black"), 
        axis.text.x =  element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black")) +
  scale_y_continuous(labels = scales::scientific) + #scientific notation
  facet_wrap(~ Gene, scales = "free")


qPCR_bphd_g_soil

ggsave(paste0("qPCR_bphd_g_soil", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = qPCR_bphd_g_soil)

ggsave(paste0("qPCR_bphd_g_soil", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = qPCR_bphd_g_soil)




## check the different 16S runs

data_plot <- subset_16S

qPCR_16S_ng_DNA <- ggplot(data_plot, 
                           aes( x = Treatment, 
                                y = mean_copy_ng_DNA,
                                fill = DNA)) +
  geom_bar(position = position_dodge(preserve = "single"), 
           stat = "identity",
           color = "black")+
  geom_errorbar(aes(ymin = mean_copy_ng_DNA, 
                    ymax = mean_copy_ng_DNA + stdev_copy_ng_DNA), #keep only upper error bars
                width = 0.2,
                position = position_dodge(0.9, preserve = "single")) +
  scale_fill_manual(values = c("darkcyan", "rosybrown3"),
                    breaks = c("cDNA", "DNA")) +
  #  coord_cartesian( ylim = c(0, 50000000) )+
  labs(x = "Strain", y = "copies/ ng DNA") +
  theme_classic(base_size = 16) +
  theme(axis.title = element_text(colour = "black"), 
        axis.text.x =  element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black")) +
  scale_y_continuous(labels = scales::scientific) + #scientific notation
  facet_wrap(~ qPCR, scales = "free")


qPCR_16S_ng_DNA


setwd("C:/Users/SarahBenning/OneDrive - Helmholtz Zentrum München/R79_inoculation/cDNA_qPCR")

ggsave(paste0("qPCR_16S_ng_DNA", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = qPCR_16S_ng_DNA)

ggsave(paste0("qPCR_16S_ng_DNA", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = qPCR_16S_ng_DNA)


qPCR_16S_g_soil <- ggplot(data_plot, 
                           aes( x = Treatment, 
                                y = mean_copy_g_soil,
                                fill = DNA)) +
  geom_bar(position = position_dodge(preserve = "single"), 
           stat = "identity",
           color = "black")+
  geom_errorbar(aes(ymin = mean_copy_g_soil, 
                    ymax = mean_copy_g_soil + stdev_copy_g_soil), #keep only upper error bars
                width = 0.2,
                position = position_dodge(0.9, preserve = "single")) +
  scale_fill_manual(values = c("darkcyan", "rosybrown3"),
                    breaks = c("cDNA", "DNA")) +
  #  coord_cartesian( ylim = c(0, 50000000) )+
  labs(x = "Strain", y = "copies/ g soil") +
  theme_classic(base_size = 16) +
  theme(axis.title = element_text(colour = "black"), 
        axis.text.x =  element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black")) +
  scale_y_continuous(labels = scales::scientific) + #scientific notation
  facet_wrap(~ qPCR, scales = "free")


qPCR_16S_g_soil

ggsave(paste0("qPCR_16S_g_soil", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = qPCR_16S_g_soil)

ggsave(paste0("qPCR_16S_g_soil", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = qPCR_16S_g_soil)



########################


## compare qPCR of 29.11.23

#data_plot <- subset_29_11
data_plot <- subset_28_11

qPCR_ng_DNA <- ggplot(data_plot, 
                   aes( x = Treatment, 
                        y = mean_copy_ng_DNA,
                        fill = DNA)) +
  geom_bar(position = position_dodge(preserve = "single"), 
           stat = "identity",
           color = "black")+
  geom_errorbar(aes(ymin = mean_copy_ng_DNA, 
                    ymax = mean_copy_ng_DNA + stdev_copy_ng_DNA), #keep only upper error bars
                width = 0.2,
                position = position_dodge(0.9, preserve = "single")) +
  scale_fill_manual(values = c("darkcyan", "rosybrown3"),
                    breaks = c("cDNA", "DNA")) +
#  coord_cartesian( ylim = c(0, 100000000) )+
  labs(x = "Strain", y = "copies/ ng DNA") +
  theme_classic(base_size = 16) +
  theme(axis.title = element_text(colour = "black"), 
        axis.text.x =  element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black")) +
  scale_y_continuous(labels = scales::scientific) + #scientific notation
  facet_wrap(~ Gene, scales = "free")


qPCR_ng_DNA


setwd("C:/Users/SarahBenning/OneDrive - Helmholtz Zentrum München/R79_inoculation/cDNA_qPCR/checked/")

ggsave(paste0("qPCR_28_11_23_copies_per_ng_DNA", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = qPCR_ng_DNA)

ggsave(paste0("qPCR_28_11_23_copies_per_ng_DNA", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = qPCR_ng_DNA)



qPCR_g_soil <- ggplot(data_plot, 
                   aes( x = Treatment, 
                        y = mean_copy_g_soil,
                        fill = DNA)) +
  geom_bar(position = position_dodge(preserve = "single"), 
           stat = "identity",
           color = "black")+
  geom_errorbar(aes(ymin = mean_copy_g_soil, 
                    ymax = mean_copy_g_soil + stdev_copy_g_soil), #keep only upper error bars
                width = 0.2,
                position = position_dodge(0.9, preserve = "single")) +
  scale_fill_manual(values = c("darkcyan", "rosybrown3"),
                    breaks = c("cDNA", "DNA")) +
  #  coord_cartesian( ylim = c(0, 50000000) )+
  labs(x = "Strain", y = "copies/ g soil") +
  theme_classic(base_size = 16) +
  theme(axis.title = element_text(colour = "black"), 
        axis.text.x =  element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black")) +
  scale_y_continuous(labels = scales::scientific) + #scientific notation
  facet_wrap(~ Gene, scales = "free")


qPCR_g_soil


setwd("C:/Users/SarahBenning/OneDrive - Helmholtz Zentrum München/R79_inoculation/cDNA_qPCR")

ggsave(paste0("qPCR_28_11_23_copies_per_g_soil", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = qPCR_g_soil)

ggsave(paste0("qPCR_28_11_23_copies_per_g_soil", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = qPCR_g_soil)


#############

data_plot <- subset_29_11_DNA
data_plot <- subset_28_11_DNA

qPCR_g_soil_DNA <- ggplot(data_plot, 
                      aes( x = Treatment, 
                           y = mean_copy_g_soil,
                           fill = DNA)) +
  geom_bar(position = position_dodge(preserve = "single"), 
           stat = "identity",
           color = "black")+
  geom_errorbar(aes(ymin = mean_copy_g_soil, 
                    ymax = mean_copy_g_soil + stdev_copy_g_soil), #keep only upper error bars
                width = 0.2,
                position = position_dodge(0.9, preserve = "single")) +
  scale_fill_manual(values = c("rosybrown3")) +
  #  coord_cartesian( ylim = c(0, 50000000) )+
  labs(x = "Strain", y = "copies/ g soil") +
  theme_classic(base_size = 16) +
  theme(axis.title = element_text(colour = "black"), 
        axis.text.x =  element_text(colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(colour = "black")) +
  scale_y_continuous(labels = scales::scientific) + #scientific notation
  facet_wrap(~ Gene, scales = "free")


qPCR_g_soil_DNA


setwd("C:/Users/SarahBenning/OneDrive - Helmholtz Zentrum München/R79_inoculation/cDNA_qPCR/checked/")

ggsave(paste0("qPCR_only_DNA_28_11_23_copies_per_g_soil", format(Sys.time(), "_%d%m%y_%H%M"), ".svg"), plot = qPCR_g_soil_DNA)

ggsave(paste0("qPCR_only_DNA_28_11_23_copies_per_g_soil", format(Sys.time(), "_%d%m%y_%H%M"), ".png"), plot = qPCR_g_soil_DNA)



###########################################################

## statistics

##########################################################


setwd("C:/Users/SarahBenning/OneDrive - Helmholtz Zentrum München/R79_inoculation/cDNA_qPCR/checked/")

#data_all <- read.csv2("all_values_qPCR_28.11.23_checked.csv")

data_all <- read.csv2("all_values_DNA_qPCR_29.11.23_checked.csv")

str(data_all)

data_all$Sample <- as.factor(data_all$Sample)
data_all$DNA <- as.factor(data_all$DNA)
data_all$Treatment <- as.factor(data_all$Treatment)
data_all$Detector <- as.factor(data_all$Detector)
data_all$qPCR <- as.factor(data_all$qPCR)


subset_16S_all <- subset(data_all, Detector == "16S")
subset_bphd_all <- subset(data_all, Detector == "bphd")

#subset_16S_DNA <- subset(subset_16S, DNA == "DNA")
#subset_16S_cDNA <- subset(subset_16S, DNA == "cDNA")


data_stat <- subset_16S_DNA

stat_DNA_16S_all <- data.frame(compare_means(formula = copies.g.soil ~ Treatment, data = data_stat, 
                                                method = "wilcox.test", paired = FALSE, group.by = NULL, ref.group = NULL))

View(stat_DNA_16S)


write.csv2(stat_DNA_16S, file = "stats_16S_DNA_28.11.checked.csv")


data_stat <- subset_16S_cDNA

stat_cDNA_16S <- data.frame(compare_means(formula = copies.g.soil ~ Treatment, data = data_stat, 
                                         method = "wilcox.test", paired = FALSE, group.by = NULL, ref.group = NULL))

View(stat_cDNA_16S)


write.csv2(stat_cDNA_16S, file = "stats_16S_cDNA_28.11.checked.csv")



data_stat <- subset_bphd_all

stat_bphd <- data.frame(compare_means(formula = copies.g.soil ~ Treatment, data = data_stat, 
                                          method = "wilcox.test", paired = FALSE, group.by = NULL, ref.group = NULL))

View(stat_bphd)


write.csv2(stat_bphd, file = "stats_bphd_28.11.checked.csv")

######


data_stat <- subset_16S_all

stat_DNA_16S_all <- data.frame(compare_means(formula = copies.g.soil ~ Treatment, data = data_stat, 
                                             method = "wilcox.test", paired = FALSE, group.by = NULL, ref.group = NULL))

View(stat_DNA_16S_all)


write.csv2(stat_DNA_16S_all, file = "stats_16S_DNA_29.11.checked.csv")




data_stat <- subset_bphd_all

stat_bphd <- data.frame(compare_means(formula = copies.g.soil ~ Treatment, data = data_stat, 
                                      method = "wilcox.test", paired = FALSE, group.by = NULL, ref.group = NULL))

View(stat_bphd)


write.csv2(stat_bphd, file = "stats_bphd_29.11.checked.csv")
