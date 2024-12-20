
#Flume-mesocosms-biofilms

###Analysis amplicon sequences and biolog statistics

##Load libraries
library(reshape2)
library(ggpubr)
library(gridExtra)
library(grid)
library(patchwork)
library(vegan)
library(cowplot)
library(compositions)
library(data.table)
library(pairwiseAdonis)
library(phyloseq)
library(tidyverse)
library(viridis)
library(microViz)
library(scales)
library(ggtext) 
library(ggrepel)
library(metR)

options(scipen=999)

### Preparation -------------------------------------------------------------

##Average BA and ChlA
sample_data(Sentinel16S_2022_bacteria_filt) <- data.frame(sample_data(Sentinel16S_2022_bacteria_filt)) %>% 
  rownames_to_column() %>% 
  group_by(Sample) %>% 
  mutate(BA = mean(BA), ChlA = mean(ChlA)) %>% 
  column_to_rownames() %>% 
  sample_data()
sample_data(Sentinel18S_2022_eukaryote_filt) <- data.frame(sample_data(Sentinel18S_2022_eukaryote_filt)) %>% 
  rownames_to_column() %>% 
  group_by(Sample) %>% 
  mutate(BA = mean(BA), ChlA = mean(ChlA)) %>% 
  column_to_rownames() %>% 
  sample_data()
sample_data(Sentinel18S_2022_photo_filt) <- data.frame(sample_data(Sentinel18S_2022_photo_filt)) %>% 
  rownames_to_column() %>% 
  group_by(Sample) %>% 
  mutate(BA = mean(BA), ChlA = mean(ChlA)) %>% 
  column_to_rownames() %>% 
  sample_data()
sample_data(Sentinel18S_2022_others_filt) <- data.frame(sample_data(Sentinel18S_2022_others_filt)) %>% 
  rownames_to_column() %>% 
  group_by(Sample) %>% 
  mutate(BA = mean(BA), ChlA = mean(ChlA)) %>% 
  column_to_rownames() %>% 
  sample_data()

#Taxonomy files
taxonomy_16S <- data.frame(tax_table(Sentinel16S_2022_bacteria_filt)) %>% 
  select(-Species) %>% 
  distinct(Genus, .keep_all = TRUE)

taxonomy_18S <- data.frame(tax_table(Sentinel18S_2022_eukaryote_filt)) %>% 
  select(-Species) %>% 
  distinct(Genus, .keep_all = TRUE)


### Taxonomy ----------------------------------------------------------------

## 16S Figure 2A + Figure S3

#Transform in relative abundance
Sentinel16S_2022_bacteria_filt_merged_relative <- transform_sample_counts(Sentinel16S_2022_bacteria_filt_merged, function(x) x / sum(x))
Sentinel16S_2022_bacteria_filt_merged_relative

#Export ASV_TAX table
ASV16S <- data.frame(otu_table(Sentinel16S_2022_bacteria_filt_merged_relative)) %>% rownames_to_column()
TAX16S <- data.frame(tax_table(Sentinel16S_2022_bacteria_filt_merged_relative)) %>% rownames_to_column()
ASVTAX16S <- left_join(ASV16S, TAX16S)
write.csv(ASVTAX16S, "/Users/touchett/Documents/Éducation/EPFL/Sentinel/Chapter1/ASVTAX16S.csv", row.names=FALSE)

#Table ASV percent per samples 
write.table(Sentinel16S_2022_bacteria_filt_merged_relative %>% tax_glom(taxrank = "Phylum") %>% 
              transform_sample_counts(function(x) {x*100}) %>% psmelt() %>%
              dplyr::select(Phylum, Sample, Abundance) %>% spread(Sample, Abundance),
            file = "16S.relative_abundance.phylum.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(Sentinel16S_2022_bacteria_filt_merged_relative %>% tax_glom(taxrank = "Family") %>% 
              transform_sample_counts(function(x) {x*100}) %>% psmelt() %>%
              dplyr::select(Family, Sample, Abundance) %>% spread(Sample, Abundance),
            file = "16S.relative_abundance.family.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(Sentinel16S_2022_bacteria_filt_merged_relative %>% tax_glom(taxrank = "Genus") %>% 
              transform_sample_counts(function(x) {x*100}) %>% psmelt() %>%
              dplyr::select(Genus, Sample, Abundance) %>% spread(Sample, Abundance),
            file = "16S.relative_abundance.genus.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

#Get taxa summary
Phylum_16S <- Sentinel16S_2022_bacteria_filt_merged_relative %>% 
  tax_glom(taxrank = "Phylum") %>% 
  transform_sample_counts(function(x) {x*100}) %>% 
  psmelt() %>% 
  dplyr::select(Phylum, Abundance) %>% 
  group_by(Phylum) %>% 
  mutate(Average = mean(Abundance), SD = sd(Abundance)) %>% 
  select(-Abundance) %>% 
  distinct(Phylum, .keep_all = TRUE)

Family_16S <- Sentinel16S_2022_bacteria_filt_merged_relative %>%
  tax_glom(taxrank = "Family") %>% 
  transform_sample_counts(function(x) {x*100}) %>% 
  psmelt() %>% 
  dplyr::select(Family, Abundance) %>% 
  group_by(Family) %>% 
  mutate(Average = mean(Abundance), SD = sd(Abundance)) %>% 
  select(-Abundance) %>% 
  distinct(Family, .keep_all = TRUE)

Genus_16S <- Sentinel16S_2022_bacteria_filt_merged_relative %>% 
  tax_glom(taxrank = "Genus") %>% 
  transform_sample_counts(function(x) {x*100}) %>% 
  psmelt() %>% 
  dplyr::select(Genus, Abundance) %>% 
  group_by(Genus) %>% 
  mutate(Average = mean(Abundance), SD = sd(Abundance)) %>% 
  select(-Abundance) %>% 
  distinct(Genus, .keep_all = TRUE)

#Plot Graph-Bar at Genus level (Figure 2A)
Graphbar_genus_16S <- Sentinel16S_2022_bacteria_filt_merged %>%
  merge_samples(group = "Succession") %>% 
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 14,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = TRUE,
    sample_order = "default",
    bar_width = 0.9) + 
  labs(x = "Days of growth") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(size=16, face="bold"),
        axis.title.x = element_text(size=16, face="bold"),
        strip.text = element_text(size = 16), 
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  guides(fill = guide_legend(ncol = 3))
Graphbar_genus_16S
ggsave("F2A.png", Graphbar_genus_16S, width = 30,
       height = 20,
       units = "cm")

#Figure S5A
Graphbar_phylum_16S_all <- Sentinel16S_2022_bacteria_filt_merged %>%
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 19,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    bar_width = 0.9) + 
  labs(x = "Days of growth") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  scale_y_continuous(expand = c(0, 0.01)) +
  facet_grid(cols = vars(factor(Flow, levels = c("Natural", "Intermittent", "Stochastic", "Constant"))) ,rows = vars(Temperature), scales = "free_y", space = "free_y") +
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 14, face = "bold"), 
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        strip.text = element_text(size = 16))
Graphbar_phylum_16S_all
ggsave("FS5A.png", Graphbar_phylum_16S_all, width = 50,
       height = 20,
       units = "cm")

#Figure S5B
Graphbar_family_16S_all <- Sentinel16S_2022_bacteria_filt_merged %>%
  tax_fix(unknowns = c("uncultured", "Unknown_Family")) %>% 
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Family", n_taxa = 19,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    bar_width = 0.9) + 
  labs(x = "Days of growth") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  scale_y_continuous(expand = c(0, 0.01)) +
  facet_grid(cols = vars(factor(Flow, levels = c("Natural", "Intermittent", "Stochastic", "Constant"))) ,rows = vars(Temperature), scales = "free_y", space = "free_y") +
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 14, face = "bold"), 
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        strip.text = element_text(size = 16))
Graphbar_family_16S_all
ggsave("FS5B.png", Graphbar_family_16S_all, width = 50,
       height = 20,
       units = "cm")

#Figure S5C
Graphbar_genus_16S_all <- Sentinel16S_2022_bacteria_filt_merged %>%
  tax_fix(unknowns = c("uncultured", "Unknown_Family")) %>% 
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 19,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    bar_width = 0.9) + 
  labs(x = "Days of growth") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  scale_y_continuous(expand = c(0, 0.01)) +
  facet_grid(cols = vars(factor(Flow, levels = c("Natural", "Intermittent", "Stochastic", "Constant"))) ,rows = vars(Temperature), scales = "free_y", space = "free_y") +
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 14, face = "bold"), 
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        strip.text = element_text(size = 16))
Graphbar_genus_16S_all
ggsave("FS5C.png", Graphbar_genus_16S_all, width = 50,
       height = 20,
       units = "cm")

#Merging Figure S5
FigureS5 <- ggarrange(Graphbar_phylum_16S_all, Graphbar_family_16S_all, Graphbar_genus_16S_all, ncol = 1, nrow = 3, labels = c("A", "B", "C"), font.label = list(size = 18))
FigureS5 <- annotate_figure(FigureS5, right = text_grob("Temperature regimes", size = 16, face = "bold", rot = 270), top = text_grob("Flow regimes", size = 16, face = "bold"))
FigureS5
ggsave("FigureS5.pdf", FigureS5, width = 50,
       height = 70,
       units = "cm")

##18S Figure 2B + Figure S6

#Transform in relative abundance
Sentinel18S_2022_eukaryote_filt_merged_relative <- transform_sample_counts(Sentinel18S_2022_eukaryote_filt_merged, function(x) x / sum(x))
Sentinel18S_2022_eukaryote_filt_merged_relative

Sentinel18S_2022_photo_filt_merged_relative <- transform_sample_counts(Sentinel18S_2022_photo_filt_merged, function(x) x / sum(x))
Sentinel18S_2022_photo_filt_merged_relative

Sentinel18S_2022_others_filt_merged_relative <- transform_sample_counts(Sentinel18S_2022_others_filt_merged, function(x) x / sum(x))
Sentinel18S_2022_others_filt_merged_relative

#Export ASV_TAX table
ASV18S <- data.frame(otu_table(Sentinel18S_2022_eukaryote_filt_merged_relative)) %>% rownames_to_column()
TAX18S <- data.frame(tax_table(Sentinel18S_2022_eukaryote_filt_merged_relative)) %>% rownames_to_column()
ASVTAX18S <- left_join(ASV18S, TAX18S)
write.csv(ASVTAX18S, "/Users/touchett/Documents/Éducation/EPFL/Sentinel/Chapter1/ASVTAX18S.csv", row.names=FALSE)

#Table ASV percent per samples 
write.table(Sentinel18S_2022_eukaryote_filt_merged_relative %>% tax_glom(taxrank = "Phylum") %>% 
              transform_sample_counts(function(x) {x*100}) %>% psmelt() %>%
              dplyr::select(Phylum, Sample, Abundance) %>% spread(Sample, Abundance),
            file = "18S.relative_abundance.phylum.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(Sentinel18S_2022_eukaryote_filt_merged_relative %>% tax_glom(taxrank = "Order") %>% 
              transform_sample_counts(function(x) {x*100}) %>% psmelt() %>%
              dplyr::select(Order, Sample, Abundance) %>% spread(Sample, Abundance),
            file = "18S.relative_abundance.order.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(Sentinel18S_2022_eukaryote_filt_merged_relative %>% tax_glom(taxrank = "Genus") %>% 
              transform_sample_counts(function(x) {x*100}) %>% psmelt() %>%
              dplyr::select(Genus, Sample, Abundance) %>% spread(Sample, Abundance),
            file = "18S.relative_abundance.genus.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

#Get taxa summary 
Phylum_18S_photo <- Sentinel18S_2022_photo_filt_merged_relative %>% 
  tax_glom(taxrank = "Phylum") %>% 
  transform_sample_counts(function(x) {x*100}) %>% 
  psmelt() %>% 
  dplyr::select(Phylum, Abundance) %>% 
  dplyr::group_by(Phylum) %>% 
  mutate(Average = mean(Abundance), SD = sd(Abundance)) %>% 
  select(-Abundance) %>% 
  distinct(Phylum, .keep_all = TRUE)

Class_18S_photo <- Sentinel18S_2022_photo_filt_merged_relative %>% 
  tax_glom(taxrank = "Class") %>% 
  transform_sample_counts(function(x) {x*100}) %>% 
  psmelt() %>% 
  dplyr::select(Class, Abundance) %>% 
  group_by(Class) %>% 
  mutate(Average = mean(Abundance), SD = sd(Abundance)) %>% 
  select(-Abundance) %>% 
  distinct(Class, .keep_all = TRUE)

Genus_18S_photo <- Sentinel18S_2022_photo_filt_merged_relative %>%
  tax_glom(taxrank = "Genus") %>% 
  transform_sample_counts(function(x) {x*100}) %>% 
  psmelt() %>% 
  dplyr::select(Genus, Abundance) %>% 
  group_by(Genus) %>% 
  mutate(Average = mean(Abundance), SD = sd(Abundance)) %>% 
  select(-Abundance) %>% 
  distinct(Genus, .keep_all = TRUE)

Phylum_18S_others <- Sentinel18S_2022_others_filt_merged_relative %>% 
  tax_glom(taxrank = "Phylum") %>% 
  transform_sample_counts(function(x) {x*100}) %>% 
  psmelt() %>% 
  dplyr::select(Phylum, Abundance) %>% 
  dplyr::group_by(Phylum) %>% 
  mutate(Average = mean(Abundance), SD = sd(Abundance)) %>% 
  select(-Abundance) %>% 
  distinct(Phylum, .keep_all = TRUE)

Class_18S_others <- Sentinel18S_2022_others_filt_merged_relative %>% 
  tax_glom(taxrank = "Class") %>% 
  transform_sample_counts(function(x) {x*100}) %>% 
  psmelt() %>% 
  dplyr::select(Class, Abundance) %>% 
  group_by(Class) %>% 
  mutate(Average = mean(Abundance), SD = sd(Abundance)) %>% 
  select(-Abundance) %>% 
  distinct(Class, .keep_all = TRUE)

Genus_18S_others <- Sentinel18S_2022_others_filt_merged_relative %>%
  tax_glom(taxrank = "Genus") %>% 
  transform_sample_counts(function(x) {x*100}) %>% 
  psmelt() %>% 
  dplyr::select(Genus, Abundance) %>% 
  group_by(Genus) %>% 
  mutate(Average = mean(Abundance), SD = sd(Abundance)) %>% 
  select(-Abundance) %>% 
  distinct(Genus, .keep_all = TRUE)

#Plot Graph-Bar at Class level photo only (Figure 2B)
Graphbar_class_18S_photo <- Sentinel18S_2022_photo_filt_merged %>%
  merge_samples(group = "Succession") %>%
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Class", n_taxa = 14,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = TRUE,
    sample_order = "default",
    palette = distinct_palette(14, pal = "kelly"),
    bar_width = 0.9) + 
  labs(x = "Days of growth") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14), 
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(size=16, face="bold"),
        axis.title.x = element_text(size=16, face="bold"),
        strip.text = element_text(size = 16), 
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  guides(fill = guide_legend(ncol = 3))
Graphbar_class_18S_photo
ggsave("F2B.png", Graphbar_class_18S_photo, width = 30,
       height = 20,
       units = "cm")

#Figure S7A
Graphbar_phylum_18S_photo_all <- Sentinel18S_2022_photo_filt_merged %>%
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 14,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    palette = distinct_palette(14, pal = "kelly"),
    bar_width = 0.9) + 
  labs(x = "Days of growth") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  scale_y_continuous(expand = c(0, 0.01)) +
  facet_grid(cols = vars(factor(Flow, levels = c("Natural", "Intermittent", "Stochastic", "Constant"))) ,rows = vars(Temperature), scales = "free_y", space = "free_y") +
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 14, face = "bold"), 
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        strip.text = element_text(size = 16))
Graphbar_phylum_18S_photo_all
ggsave("FS6A.png", Graphbar_phylum_18S_photo_all, width = 50,
       height = 20,
       units = "cm")

#Figure S6B
Graphbar_class_18S_photo_all <- Sentinel18S_2022_photo_filt_merged %>%
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Class", n_taxa = 14,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    palette = distinct_palette(14, pal = "kelly"),
    bar_width = 0.9) + 
  labs(x = "Days of growth") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  scale_y_continuous(expand = c(0, 0.01)) +
  facet_grid(cols = vars(factor(Flow, levels = c("Natural", "Intermittent", "Stochastic", "Constant"))) ,rows = vars(Temperature), scales = "free_y", space = "free_y") +
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 14, face = "bold"), 
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        strip.text = element_text(size = 16))
Graphbar_class_18S_photo_all
ggsave("FS6B.png", Graphbar_class_18S_photo_all, width = 50,
       height = 20,
       units = "cm")

#Figure S6C
Graphbar_genus_18S_photo_all <- Sentinel18S_2022_photo_filt_merged %>%
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 14,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    palette = distinct_palette(14, pal = "kelly"),
    bar_width = 0.9) + 
  labs(x = "Days of growth") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  scale_y_continuous(expand = c(0, 0.01)) +
  facet_grid(cols = vars(factor(Flow, levels = c("Natural", "Intermittent", "Stochastic", "Constant"))) ,rows = vars(Temperature), scales = "free_y", space = "free_y") +
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 14, face = "bold"), 
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        strip.text = element_text(size = 16))
Graphbar_genus_18S_photo_all
ggsave("FS6C.png", Graphbar_genus_18S_photo_all, width = 50,
       height = 20,
       units = "cm")

#Merging Figure S6
FigureS8 <- ggarrange(Graphbar_phylum_18S_photo_all, Graphbar_class_18S_photo_all, Graphbar_genus_18S_photo_all, ncol = 1, nrow = 3, labels = c("A", "B", "C"), font.label = list(size = 18))
FigureS8 <- annotate_figure(FigureS8, right = text_grob("Temperature regimes", size = 16, face = "bold", rot = 270), top = text_grob("Flow regimes", size = 16, face = "bold"))
FigureS8
ggsave("FigureS8.pdf", FigureS8, width = 50,
       height = 70,
       units = "cm")

#New palettes
myPal2_phylum <- tax_palette(data = Sentinel18S_2022_others_filt_merged, rank = "Phylum", n = 20, pal = "greenArmytage", add = c(Other = "lightgray"))
myPal2_phylum[1] <- "darkorange"
myPal2_class <- tax_palette(data = Sentinel18S_2022_others_filt_merged, rank = "Class", n = 20, pal = "greenArmytage", add = c(Other = "lightgray"))
myPal2_class[1] <- "darkorange"
myPal2_genus <- tax_palette(data = Sentinel18S_2022_others_filt_merged, rank = "Genus", n = 20, pal = "greenArmytage", add = c(Other = "lightgray"))
myPal2_genus[1] <- "darkorange"


#Figure S8A
Graphbar_phylum_18S_others_all <- Sentinel18S_2022_others_filt_merged %>%
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 14,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    palette = myPal2_phylum,
    bar_width = 0.9) + 
  labs(x = "Days of growth") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  scale_y_continuous(expand = c(0, 0.01)) +
  facet_grid(cols = vars(factor(Flow, levels = c("Natural", "Intermittent", "Stochastic", "Constant"))) ,rows = vars(Temperature), scales = "free_y", space = "free_y") +
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 14, face = "bold"), 
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        strip.text = element_text(size = 16))
Graphbar_phylum_18S_others_all
ggsave("FS8A.png", Graphbar_phylum_18S_others, width = 50,
       height = 20,
       units = "cm")

#Figure S8B
Graphbar_class_18S_others_all <- Sentinel18S_2022_others_filt_merged %>%
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Class", n_taxa = 14,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    palette = myPal2_class,
    bar_width = 0.9) + 
  labs(x = "Days of growth") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  scale_y_continuous(expand = c(0, 0.01)) +
  facet_grid(cols = vars(factor(Flow, levels = c("Natural", "Intermittent", "Stochastic", "Constant"))) ,rows = vars(Temperature), scales = "free_y", space = "free_y") +
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 14, face = "bold"), 
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        strip.text = element_text(size = 16))
Graphbar_class_18S_others_all
ggsave("FS8B.png", Graphbar_class_18S_others, width = 50,
       height = 20,
       units = "cm")

#Figure S8C
Graphbar_genus_18S_others_all <- Sentinel18S_2022_others_filt_merged %>%
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 14,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    palette = myPal2_genus,
    bar_width = 0.9) + 
  labs(x = "Days of growth") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  scale_y_continuous(expand = c(0, 0.01)) +
  facet_grid(cols = vars(factor(Flow, levels = c("Natural", "Intermittent", "Stochastic", "Constant"))) ,rows = vars(Temperature), scales = "free_y", space = "free_y") +
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 14, face = "bold"), 
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"),
        strip.text = element_text(size = 16))
Graphbar_genus_18S_others_all
ggsave("FS8C.png", Graphbar_genus_18S_others, width = 50,
       height = 20,
       units = "cm")

#Merging Figure S7
FigureS9 <- ggarrange(Graphbar_phylum_18S_others_all, Graphbar_class_18S_others_all, Graphbar_genus_18S_others_all, ncol = 1, nrow = 3, labels = c("A", "B", "C"), font.label = list(size = 18))
FigureS9 <- annotate_figure(FigureS9, right = text_grob("Temperature regimes", size = 16, face = "bold", rot = 270), top = text_grob("Flow regimes", size = 16, face = "bold"))
FigureS9
ggsave("FigureS9.pdf", FigureS8, width = 50,
       height = 70,
       units = "cm")


### Pattern of succession ---------------------------------------------------

##16S - Figure 2 

# Calculate alpha diversity metrics
alphadiv_16S <- estimate_richness(Sentinel16S_2022_bacteria_filt, measures = c("Observed", "Shannon"))
alphadiv_16S_df <- cbind(sample_data(Sentinel16S_2022_bacteria_filt), alphadiv_16S) %>% mutate(Evenness = Shannon/log(Observed))
sample_data(Sentinel16S_2022_bacteria_filt) <- alphadiv_16S_df

#Extract metadata
metadata_16S <- data.frame(sample_data(Sentinel16S_2022_bacteria_filt))
metadata_16S$Treatment <- factor(metadata_16S$Treatment, levels = c("Nc", "Nw", "Ic", "Iw", "Sc", "Sw", "Cc", "Cw"))

#Plot Bacterial abundance over time (Figure 2C)
BA_plot <- ggplot(metadata_16S, aes(x = Succession, group=Treatment, y = BA)) +
  geom_point(aes(shape = Temperature), color= "gray", size = 3, alpha = 0.8) +
  geom_smooth(aes(color = Treatment, group = Treatment),linewidth = 1.25, se = FALSE, span = 0.5, lineend = "round") +
  scale_color_manual(values = c("#283747","#85929e","#b03a2e","#ec7063","#1e8449","#52be80","#2874a6","#5dade2"),
                     labels= c(bquote("N"[flow]*"C"[temp]), bquote("N"[flow]*"W"[temp]), bquote("I"[flow]*"C"[temp]), bquote("I"[flow]*"W"[temp]), bquote("S"[flow]*"C"[temp]), bquote("S"[flow]*"W"[temp]), bquote("C"[flow]*"C"[temp]), bquote("C"[flow]*"W"[temp]))) +
  scale_y_continuous(n.breaks=5, trans = "log10", labels = scientific) + 
  theme_classic() +
  labs(x = "Days of growth", y = bquote(bold('Bacterial abundance'~(cells/cm^2)))) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14), 
        axis.text.y =  element_text(size = 14, angle=90, hjust = 0.5),
        legend.position = "none") +
  guides(shape = guide_legend(override.aes = list(size = 3.5)),
         colour = guide_legend(override.aes = list(linewidth= 3, linetype=2)))
BA_plot


#Plot Shannon diversity index over time (Figure 2E)
Shannon_16S_plot <- ggplot(metadata_16S, aes(x = Succession, group=Treatment, y = Shannon)) +
  geom_point(aes(shape = Temperature), color= "gray", size = 3, alpha = 0.8) +
  geom_smooth(aes(color = Treatment, group = Treatment),linewidth = 1.25, se = FALSE, span = 0.5, lineend = "round") +
  scale_color_manual(values = c("#283747","#85929e","#b03a2e","#ec7063","#1e8449","#52be80","#2874a6","#5dade2"),
                     labels= c(bquote("N"[flow]*"C"[temp]), bquote("N"[flow]*"W"[temp]), bquote("I"[flow]*"C"[temp]), bquote("I"[flow]*"W"[temp]), bquote("S"[flow]*"C"[temp]), bquote("S"[flow]*"W"[temp]), bquote("C"[flow]*"C"[temp]), bquote("C"[flow]*"W"[temp]))) +
  scale_y_continuous(n.breaks=6) + 
  theme_classic() +
  ylab("Shannon index") +
  xlab("Days of growth") +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14), 
        axis.text.y =  element_text(size = 14, angle=90, hjust = 0.5),
        legend.position = "none",
        plot.margin = margin(0, 0, 1, 0.45, "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 3.5)),
         colour = guide_legend(override.aes = list(linewidth= 3, linetype=2)))
Shannon_16S_plot

##18S - Figure 2 

# Calculate alpha diversity metrics
alphadiv_18S_photo <- estimate_richness(Sentinel18S_2022_photo_filt, measures = c("Observed", "Shannon"))
alphadiv_18S_df_photo <- cbind(sample_data(Sentinel18S_2022_photo_filt), alphadiv_18S_photo) %>% mutate(Evenness = Shannon/log(Observed))
sample_data(Sentinel18S_2022_photo_filt) <- alphadiv_18S_df_photo

alphadiv_18S_others <- estimate_richness(Sentinel18S_2022_others_filt, measures = c("Observed", "Shannon"))
alphadiv_18S_df_others <- cbind(sample_data(Sentinel18S_2022_others_filt), alphadiv_18S_others) %>% mutate(Evenness = Shannon/log(Observed))
sample_data(Sentinel18S_2022_others_filt) <- alphadiv_18S_df_others

#Extract metadata
metadata_18S_photo <- data.frame(sample_data(Sentinel18S_2022_photo_filt))
metadata_18S_photo$Treatment <- factor(metadata_18S_photo$Treatment, levels = c("Nc", "Nw", "Ic", "Iw", "Sc", "Sw", "Cc", "Cw"))

metadata_18S_others <- data.frame(sample_data(Sentinel18S_2022_others_filt))
metadata_18S_others$Treatment <- factor(metadata_18S_others$Treatment, levels = c("Nc", "Nw", "Ic", "Iw", "Sc", "Sw", "Cc", "Cw"))

#Plot Bacterial abundance over time (Figure 2D)
ChlA_plot_photo <- ggplot(metadata_18S_others, aes(x = Succession, group=Treatment, y = ChlA)) +
  geom_point(aes(shape = Temperature), color= "gray", size = 3, alpha = 0.8) +
  geom_smooth(aes(color = Treatment, group = Treatment),linewidth = 1.25, se = FALSE, span = 0.5, lineend = "round") +
  scale_color_manual(values = c("#283747","#85929e","#b03a2e","#ec7063","#1e8449","#52be80", "#2874a6","#5dade2"),
                     labels= c(bquote("N"[flow]*"C"[temp]), bquote("N"[flow]*"W"[temp]), bquote("I"[flow]*"C"[temp]), bquote("I"[flow]*"W"[temp]), bquote("S"[flow]*"C"[temp]), bquote("S"[flow]*"W"[temp]), bquote("C"[flow]*"C"[temp]), bquote("C"[flow]*"W"[temp]))) +
  scale_y_continuous(n.breaks=5) + 
  theme_classic() +
  labs(x = "Days of growth", y = bquote(bold('Chlorophyll-a concentration'~(µg/cm^2)))) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14), 
        axis.text.y =  element_text(size = 14, angle=90, hjust = 0.5),
        legend.position = "none") +
  guides(shape = guide_legend(override.aes = list(size = 3.5)),
         colour = guide_legend(override.aes = list(linewidth= 3, linetype=2)))
ChlA_plot_photo

#Plot Shannon diversity index phototrophs over time (Figure 2F)
Shannon_18S_plot_photo <- ggplot(metadata_18S_photo, aes(x = Succession, group=Treatment, y = Shannon)) +
  geom_point(aes(shape = Temperature), color= "gray", size = 3, alpha = 0.8) +
  geom_smooth(aes(color = Treatment, group = Treatment),linewidth = 1.25, se = FALSE, span = 0.5, lineend = "round") +
  scale_color_manual(values = c("#283747","#85929e","#b03a2e","#ec7063","#1e8449","#52be80", "#2874a6","#5dade2"),
                     labels= c(bquote("N"[flow]*"C"[temp]), bquote("N"[flow]*"W"[temp]), bquote("I"[flow]*"C"[temp]), bquote("I"[flow]*"W"[temp]), bquote("S"[flow]*"C"[temp]), bquote("S"[flow]*"W"[temp]), bquote("C"[flow]*"C"[temp]), bquote("C"[flow]*"W"[temp]))) +
  scale_y_continuous(n.breaks=6) + 
  theme_classic() +
  ylab("Shannon index") +
  xlab("Days of growth") +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14), 
        axis.text.y =  element_text(size = 14, angle=90, hjust = 0.5),
        legend.position = "none",
        plot.margin = margin(0, 0, 1, 0.4, "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 3.5)),
         colour = guide_legend(override.aes = list(linewidth= 3, linetype=2)))
Shannon_18S_plot_photo

#Plot Shannon diversity index phototrophs over time (Figure S10A)
Shannon_18S_plot_others <- ggplot(metadata_18S_others, aes(x = Succession, group=Treatment, y = Shannon)) +
  geom_point(aes(shape = Temperature), color= "gray", size = 3, alpha = 0.8) +
  geom_smooth(aes(color = Treatment, group = Treatment),linewidth = 1.25, se = FALSE, span = 0.5, lineend = "round") +
  scale_color_manual(values = c("#283747","#85929e","#b03a2e","#ec7063","#1e8449","#52be80", "#2874a6","#5dade2"),
                     labels= c(bquote("N"[flow]*"C"[temp]), bquote("N"[flow]*"W"[temp]), bquote("I"[flow]*"C"[temp]), bquote("I"[flow]*"W"[temp]), bquote("S"[flow]*"C"[temp]), bquote("S"[flow]*"W"[temp]), bquote("C"[flow]*"C"[temp]), bquote("C"[flow]*"W"[temp]))) +
  scale_y_continuous(n.breaks=6) + 
  theme_classic() +
  ylab("Shannon index") +
  xlab("Days of growth") +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14), 
        axis.text.y =  element_text(size = 14, angle=90, hjust = 0.5),
        legend.position = "none",
        plot.margin = margin(0, 0, 1, 0.4, "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 3.5)),
         colour = guide_legend(override.aes = list(linewidth= 3, linetype=2)))
Shannon_18S_plot_others

#Merging Figure 2
top_title1 <- text_grob("", size = 18, hjust = 0.5, vjust = 0.5, face ="bold")
top_title2 <- text_grob("Bacteria", size = 18, hjust = 0.5, vjust = 0.5, face ="bold")
top_title3 <- text_grob("Phototrophic eukaryotes", size = 18, hjust = 0.5, vjust = 0.5, face ="bold")
left_title1 <- text_grob("Taxonomy", size = 18, rot = 90, hjust = 0.5, vjust = 0.5, face ="bold")
left_title2 <- text_grob("Biomass", size = 18, rot = 90, hjust = 0.5, vjust = 0.5, face ="bold")
left_title3 <- text_grob("Diversity", size = 18, rot = 90, hjust = 0.5, vjust = 0.5, face ="bold")

Figure2.1 <- ggarrange(top_title1, top_title2, top_title3, left_title1, Graphbar_genus_16S, Graphbar_class_18S_photo, ncol = 3, nrow = 2, labels = c("", "", "", "", "A", "B"), font.label = list(size = 18), common.legend = F, align = "h", heights = c(0.05, 1), widths = c(0.05, 1, 1))
Figure2.1
Figure2.2 <- ggarrange(left_title2, BA_plot, ChlA_plot_photo, ncol = 3, nrow = 1, labels = c("", "C", "D"), font.label = list(size = 18), align = "h", widths = c(0.05, 1, 1))
Figure2.2
Figure2.3a <- ggarrange(Shannon_16S_plot, Shannon_18S_plot_photo, ncol = 2, nrow = 1, labels = c("E", "F"), font.label = list(size = 18), align = "h", common.legend = T, legend ="bottom")
Figure2.3a
Figure2.3 <- ggarrange(left_title3, Figure2.3a, ncol =, nrow = 1, widths = c(0.05, 2))  
Figure2.3
Figure2 <- ggarrange(Figure2.1, Figure2.2, Figure2.3, ncol = 1, nrow = 3, align = "v", heights = c(1, 0.70, 0.8))
Figure2
ggsave("Figure2.pdf", Figure2, width = 40,
       height = 50,
       units = "cm")


### Environmental colinearity -----------------------------------------------

library(car)
df_colin <- metadata_16S[9:47]
df_colin <- df_colin %>% select(-Primers, -Type, -TAD, -Drought, -PO4, -NO3, -NH4, -NO2, -AWCD, -E, -R, -H, -Shannon, -Evenness, -Observed, -ChlA, -BA)  
model <- lm(Succession ~ ., data=df_colin)
vif(model)
#There is a severe correlations if VIF > 5. It is the case of C25, O2, Temp, Turbidity, Fluoride, Chloride, Nitrite, Bromide, Sulfate, Sodium, Magnesium, Potassium, Calcium, Strontium.


### Compositional patterns --------------------------------------------------

##16S - Figure 3 
#Plot NMDS of Bray-Curtis dissimilarity, on Hellinger-transformed merged data (Figure 3A)
factor(sample_data(Sentinel16S_2022_bacteria_filt_merged_hell)$Treatment, levels = c("Nc", "Nw", "Ic", "Iw", "Sc", "Sw", "Cc", "Cw"))

matrix_dist_16S_merged_2022 <- as.matrix(t(data.frame(otu_table(Sentinel16S_2022_bacteria_filt_merged_hell))))
metadata_dist_16S_merged_2022 <- data.frame(sample_data(Sentinel16S_2022_bacteria_filt_merged_hell)) %>% rownames_to_column("Samples")
set.seed(19950930)
dist_16S_merged_2022_bray <- vegdist(matrix_dist_16S_merged_2022, method ="bray") %>% 
  as.matrix() %>% 
  data.frame() %>% 
  rownames_to_column("Samples")
dist_16S_merged_2022 <- dist_16S_merged_2022_bray %>% 
  dplyr::select(all_of(.[["Samples"]])) %>% 
  as.dist()
nmds_16S_merged_2022 <- metaMDS(dist_16S_merged_2022, trymax=100)
stress <- nmds_16S_merged_2022$stress
stress

scores_nmds_16S_merged_2022 <- scores(nmds_16S_merged_2022) %>% 
  as_tibble(rownames = "Samples") %>% 
  inner_join(., metadata_dist_16S_merged_2022, by="Samples")

scores_nmds_16S_merged_2022  <- scores_nmds_16S_merged_2022[order(scores_nmds_16S_merged_2022$Treatment, scores_nmds_16S_merged_2022$Time),]

ordi_16S <- ordisurf(nmds_16S_merged_2022, metadata_dist_16S_merged_2022$Succession, plot = FALSE, bs="ds")
ordi_16S.grid <- ordi_16S$grid #extracts the ordisurf object
ordi_16S.mite <- expand.grid(x = ordi_16S.grid$x, y = ordi_16S.grid$y) #get x and ys
ordi_16S.mite$z <- as.vector(ordi_16S.grid$z) #unravel the matrix for the z scores
ordi_16S.mite.na <- data.frame(na.omit(ordi_16S.mite)) #gets rid of the nas

NMDS_16S_merged_trajectories_2022_NI <- scores_nmds_16S_merged_2022 %>%
  filter(Flow %in% c("Natural", "Intermittent")) %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape = Temperature), color = "darkgray", size = 3) +
  theme_classic() +
  geom_contour(data = ordi_16S.mite.na, aes(x = x, y = y, z = z, colour=after_stat(level)), linewidth = 0.75, linetype = "dashed", color = "darkgray") +
  geom_label_contour(data = ordi_16S.mite.na, aes(x = x, y = y, z = z), skip = 0, size = 4) +
  geom_segment(data=scores_nmds_16S_merged_2022[scores_nmds_16S_merged_2022$Treatment == "Nc",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#283747") +
  geom_segment(data=scores_nmds_16S_merged_2022[scores_nmds_16S_merged_2022$Treatment == "Nw",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#85929e") +
  geom_segment(data=scores_nmds_16S_merged_2022[scores_nmds_16S_merged_2022$Treatment == "Ic",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#b03a2e") +
  geom_segment(data=scores_nmds_16S_merged_2022[scores_nmds_16S_merged_2022$Treatment == "Iw",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#ec7063") +
  labs(color = "Days of growth") +
  ylim(-0.2,0.5) +
  xlim(-0.5,0.4) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text =  element_text(size = 14),
        strip.text = element_text(size=16, face="bold"),
        strip.background = element_blank(),
        plot.caption = element_text(size=14, hjust=0.95, vjust = 25),
        axis.line=element_line(),
        legend.direction= "horizontal",
        legend.position = "none",
        legend.background = element_rect(fill = "transparent"))
NMDS_16S_merged_trajectories_2022_NI

NMDS_16S_merged_trajectories_2022_NS <- scores_nmds_16S_merged_2022 %>%
  filter(Flow %in% c("Natural", "Stochastic")) %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape = Temperature), color = "darkgray", size = 3) +
  theme_classic() +
  geom_contour(data = ordi_16S.mite.na, aes(x = x, y = y, z = z, colour=after_stat(level)), linewidth = 0.75, linetype = "dashed", color = "darkgray") +
  geom_label_contour(data = ordi_16S.mite.na, aes(x = x, y = y, z = z), skip = 0, size = 4) +
  geom_segment(data=scores_nmds_16S_merged_2022[scores_nmds_16S_merged_2022$Treatment == "Nc",], aes(xend = lead(NMDS1), yend = lead(NMDS2)),
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#283747") +
  geom_segment(data=scores_nmds_16S_merged_2022[scores_nmds_16S_merged_2022$Treatment == "Nw",], aes(xend = lead(NMDS1), yend = lead(NMDS2)),
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#85929e") +
  geom_segment(data=scores_nmds_16S_merged_2022[scores_nmds_16S_merged_2022$Treatment == "Sc",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#1e8449") +
  geom_segment(data=scores_nmds_16S_merged_2022[scores_nmds_16S_merged_2022$Treatment == "Sw",], aes(xend = lead(NMDS1), yend = lead(NMDS2)),
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#52be80") +
  labs(color = "Days of growth") +
  ylim(-0.2,0.5) +
  xlim(-0.5,0.4) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text =  element_text(size = 14),
        strip.text = element_text(size=16, face="bold"),
        strip.background = element_blank(),
        plot.caption = element_text(size=14, hjust=0.95, vjust = 25),
        axis.line=element_line(),
        legend.direction= "horizontal",
        legend.position = "none",
        legend.background = element_rect(fill = "transparent"))
NMDS_16S_merged_trajectories_2022_NS

NMDS_16S_merged_trajectories_2022_NC <- scores_nmds_16S_merged_2022 %>%
  filter(Flow %in% c("Natural", "Constant")) %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape = Temperature), color = "darkgray", size = 3) +
  theme_classic() +
  geom_contour(data = ordi_16S.mite.na, aes(x = x, y = y, z = z, colour=after_stat(level)), linewidth = 0.75, linetype = "dashed", color = "darkgray") +
  geom_label_contour(data = ordi_16S.mite.na, aes(x = x, y = y, z = z), skip = 0, size = 4) +
  geom_segment(data=scores_nmds_16S_merged_2022[scores_nmds_16S_merged_2022$Treatment == "Nc",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#283747") +
  geom_segment(data=scores_nmds_16S_merged_2022[scores_nmds_16S_merged_2022$Treatment == "Nw",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#85929e") +
  geom_segment(data=scores_nmds_16S_merged_2022[scores_nmds_16S_merged_2022$Treatment == "Cc",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#2874a6") +
  geom_segment(data=scores_nmds_16S_merged_2022[scores_nmds_16S_merged_2022$Treatment == "Cw",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#5dade2") +
  labs(caption="Stress 0.075") +
  labs(color = "Days of growth") +
  ylim(-0.2,0.5) +
  xlim(-0.5,0.4) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text =  element_text(size = 14),
        strip.text = element_text(size=16, face="bold"),
        strip.background = element_blank(),
        plot.caption = element_text(size=14, hjust=1, vjust = 15),
        axis.line=element_line(),
        legend.direction= "horizontal",
        legend.position = "none",
        legend.background = element_rect(fill = "transparent"))
NMDS_16S_merged_trajectories_2022_NC

Figure3A <- ggarrange(top_title2, NMDS_16S_merged_trajectories_2022_NI + rremove("ylab") + rremove("xlab"), NMDS_16S_merged_trajectories_2022_NS + rremove("ylab") + rremove("xlab"), NMDS_16S_merged_trajectories_2022_NC + rremove("ylab") + rremove("xlab"), ncol = 1, nrow = 4, heights = c(0.1, 1, 1, 1))
Figure3A  <- annotate_figure(Figure3A, left = text_grob("NMDS2", rot = 90, vjust = 1, size = 16, face = "bold"),
                             bottom = text_grob("NMDS1", size = 16, face = "bold"))
Figure3A 


# Figure S12
matrix_dist_16S_2022 <- as.matrix(t(data.frame(otu_table(Sentinel16S_2022_bacteria_filt_hell))))
metadata_dist_16S_2022 <- data.frame(sample_data(Sentinel16S_2022_bacteria_filt_hell)) %>% rownames_to_column("Samples")
set.seed(19950930)
dist_16S_2022_bray <- vegdist(matrix_dist_16S_2022, method ="bray") %>% 
  as.matrix() %>% 
  data.frame() %>% 
  rownames_to_column("Samples")
dist_16S_2022 <- dist_16S_2022_bray %>% 
  dplyr::select(all_of(.[["Samples"]])) %>% 
  as.dist()
nmds_16S_2022 <- metaMDS(dist_16S_2022, trymax=100)
stress_all <- nmds_16S_2022$stress
stress_all

scores_nmds_16S_2022 <- scores(nmds_16S_2022) %>% 
  as_tibble(rownames = "Samples") %>% 
  inner_join(., metadata_dist_16S_2022, by="Samples")

scores_nmds_16S_2022  <- scores_nmds_16S_2022[order(scores_nmds_16S_2022$Treatment, scores_nmds_16S_2022$Time),]

env_16S <- data.frame(sample_data(Sentinel16S_2022_bacteria_filt_hell)) %>% 
  rownames_to_column("Samples") %>% 
  select(-c("Samples", "Sample", "Run", "Flow", "FlumeID", "Type", "Primers", "Temperature", "Time","Year","Treatment","NO3", "NO2", "NH4", "PO4", "R", "E", "AWCD", "H")) %>% 
  mutate(TAD = coalesce(TAD, 0))
set.seed(19950930)
fit_16S <- envfit(nmds_16S_2022, env_16S, permutations = 999, na.rm = TRUE)

fit_16S_vec <- data.frame(fit_16S$vectors$arrows * sqrt(fit_16S$vectors$r), P = fit_16S$vectors$pvals)
fit_16S_vec$Environment <- rownames(fit_16S_vec)
fit_16S_fac <- data.frame(fit_16S$factors$centroids * sqrt(fit_16S$factors$r), P = fit_16S$factors$pvals)
fit_16S_fac$Environment <- rownames(fit_16S_fac) 
fit_16S_fac <- fit_16S_fac %>% slice(2L) %>% mutate(Environment = "Drought")
fit_16S_env <- rbind(fit_16S_vec,fit_16S_fac)
fit_16S_env_pv <- subset(fit_16S_env, P<=0.05)

fit_16S_env_pv_star <- fit_16S_env_pv %>%
  group_by(Environment) %>%
  mutate(EnvironmentSign = case_when(
    Environment == "O2" ~ "O[2]*\"*\"",
    Environment %in% c("C25", "Turbidity", "Fluoride", "Chloride", "Nitrite", "Bromide",
                       "Sulfate", "Sodium", "Magnesium", "Potassium", "Calcium",
                       "Strontium") ~ paste0("\"", Environment, "*\""),
    TRUE ~ paste0("\"", Environment, "\"")
  )) %>%
  mutate(EnvironmentSign = as.character(EnvironmentSign))

NMDS_16S_2022 <- ggplot(scores_nmds_16S_2022, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = Succession, shape = Flow), size = 2, alpha = 0.8) +
  scale_color_viridis(discrete = FALSE, option = "viridis", breaks = c(0, 20, 40, 60, 80, 100)) +
  theme_classic() +
  guides(colour = guide_colourbar(barheight = 10)) +
  labs(color = "Days of growth", caption = "Stress 0.101") +
  scale_x_continuous(limits = c(-0.6, 0.5)) +
  theme(legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12),
        strip.background = element_blank(),
        plot.caption = element_text(size = 16, hjust = 1, vjust = 20),
        strip.text = element_text(size = 16, face = "bold"),
        axis.line = element_line()) +
  geom_segment(data = fit_16S_env_pv_star, aes(x = 0, xend = NMDS1 * 0.5, y = 0, yend = NMDS2 * 0.5),
               arrow = arrow(length = unit(0.2, "cm")), colour = "darkgrey", lwd = 0.5) +
  ggrepel::geom_text_repel(data = fit_16S_env_pv_star,
                           aes(x = NMDS1 * 0.5, y = NMDS2 * 0.5,
                               label = EnvironmentSign),
                           parse = TRUE,
                           cex = 5, direction = "x", colour = "black",
                           segment.size = 0.25,
                           nudge_x = c(-0.105, 0, 0, -0.02, -0.01, -0.08, 0.02, -0.020, -0.02, -0.12, -0.105, -0.01, -0.04, 0.10, -0.09, -0.01, -0.10, -0.05, -0.115, -0.01, -0.11, -0.125, -0.005, -0.005, 0.07),
                           nudge_y = c(-0.010, -0.01, 0.01, 0.01, -0.03, 0.07, -0.07, -0.012, 0, -0.04, -0.030, -0.05, -0.01, 0.08, 0.05, 0.05, 0.02, 0.11, 0.030, -0.02, 0.10, 0.010, -0.010, 0.010, -0.05))
NMDS_16S_2022

ggsave("FigureS12.pdf", NMDS_16S_2022, width = 40,
       height = 35,
       units = "cm")

##18S - Figure 3
#Plot NMDS of Bray-Curtis dissimilarity, on Hellinger-transformed merged data (Figure 3B)
factor(sample_data(Sentinel18S_2022_photo_filt_merged_hell)$Treatment, levels = c("Nc", "Nw", "Ic", "Iw", "Sc", "Sw", "Cc", "Cw"))

matrix_dist_18S_photo_merged_2022 <- as.matrix(t(data.frame(otu_table(Sentinel18S_2022_photo_filt_merged_hell))))
metadata_dist_18S_photo_merged_2022 <- data.frame(sample_data(Sentinel18S_2022_photo_filt_merged_hell)) %>% rownames_to_column("Samples")
set.seed(19950930)
dist_18S_photo_merged_2022_bray <- vegdist(matrix_dist_18S_photo_merged_2022, method ="bray") %>% 
  as.matrix() %>% 
  data.frame() %>% 
  rownames_to_column("Samples")
dist_18S_photo_merged_2022 <- dist_18S_photo_merged_2022_bray %>% 
  dplyr::select(all_of(.[["Samples"]])) %>% 
  as.dist()
nmds_18S_photo_merged_2022 <- metaMDS(dist_18S_photo_merged_2022, trymax=100)
stress <- nmds_18S_photo_merged_2022$stress
stress

scores_nmds_18S_photo_merged_2022 <- scores(nmds_18S_photo_merged_2022) %>% 
  as_tibble(rownames = "Samples") %>% 
  inner_join(., metadata_dist_18S_photo_merged_2022, by="Samples")

scores_nmds_18S_photo_merged_2022  <- scores_nmds_18S_photo_merged_2022[order(scores_nmds_18S_photo_merged_2022$Treatment, scores_nmds_18S_photo_merged_2022$Time),]

ordi_18S_photo <- ordisurf(nmds_18S_photo_merged_2022, metadata_dist_18S_photo_merged_2022$Succession, plot = FALSE, bs="ds")
ordi_18S_photo.grid <- ordi_18S_photo$grid #extracts the ordisurf object
ordi_18S_photo.mite <- expand.grid(x = ordi_18S_photo.grid$x, y = ordi_18S_photo.grid$y) #get x and ys
ordi_18S_photo.mite$z <- as.vector(ordi_18S_photo.grid$z) #unravel the matrix for the z scores
ordi_18S_photo.mite.na <- data.frame(na.omit(ordi_18S_photo.mite)) #gets rid of the nas

NMDS_18S_photo_merged_trajectories_2022_NI <- scores_nmds_18S_photo_merged_2022 %>%
  filter(Flow %in% c("Natural", "Intermittent")) %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape = Temperature), color = "darkgray", size = 3) +
  theme_classic() +
  geom_contour(data = ordi_18S_photo.mite.na, aes(x = x, y = y, z = z, colour=after_stat(level)), linewidth = 0.75, linetype = "dashed", color = "darkgray") +
  geom_label_contour(data = ordi_18S_photo.mite.na, aes(x = x, y = y, z = z), skip = 0, size = 4) +
  geom_segment(data=scores_nmds_18S_photo_merged_2022[scores_nmds_18S_photo_merged_2022$Treatment == "Nc",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#283747") +
  geom_segment(data=scores_nmds_18S_photo_merged_2022[scores_nmds_18S_photo_merged_2022$Treatment == "Nw",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#85929e") +
  geom_segment(data=scores_nmds_18S_photo_merged_2022[scores_nmds_18S_photo_merged_2022$Treatment == "Ic",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#b03a2e") +
  geom_segment(data=scores_nmds_18S_photo_merged_2022[scores_nmds_18S_photo_merged_2022$Treatment == "Iw",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#ec7063") +
  labs(color = "Days of growth") +
  scale_y_continuous(breaks=c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3), limit=c(-0.33, 0.30)) +
  scale_x_continuous(breaks=c(-0.2, 0.0, 0.2, 0.4), limit=c(-0.25, 0.46)) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text =  element_text(size = 14),
        strip.text = element_text(size=16, face="bold"),
        strip.background = element_blank(),
        plot.caption = element_text(size=14, hjust=0.95, vjust = 25),
        axis.line=element_line(),
        legend.direction= "horizontal",
        legend.position = "none",
        legend.background = element_rect(fill = "transparent"))
NMDS_18S_photo_merged_trajectories_2022_NI

NMDS_18S_photo_merged_trajectories_2022_NS <- scores_nmds_18S_photo_merged_2022 %>%
  filter(Flow %in% c("Natural", "Stochastic")) %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape = Temperature), color = "darkgray", size = 3) +
  theme_classic() +
  geom_contour(data = ordi_18S_photo.mite.na, aes(x = x, y = y, z = z, colour=after_stat(level)), linewidth = 0.75, linetype = "dashed", color = "darkgray") +
  geom_label_contour(data = ordi_18S_photo.mite.na, aes(x = x, y = y, z = z), skip = 0, size = 4) +
  geom_segment(data=scores_nmds_18S_photo_merged_2022[scores_nmds_18S_photo_merged_2022$Treatment == "Nc",], aes(xend = lead(NMDS1), yend = lead(NMDS2)),
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#283747") +
  geom_segment(data=scores_nmds_18S_photo_merged_2022[scores_nmds_18S_photo_merged_2022$Treatment == "Nw",], aes(xend = lead(NMDS1), yend = lead(NMDS2)),
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#85929e") +
  geom_segment(data=scores_nmds_18S_photo_merged_2022[scores_nmds_18S_photo_merged_2022$Treatment == "Sc",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#1e8449") +
  geom_segment(data=scores_nmds_18S_photo_merged_2022[scores_nmds_18S_photo_merged_2022$Treatment == "Sw",], aes(xend = lead(NMDS1), yend = lead(NMDS2)),
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#52be80") +
  scale_y_continuous(breaks=c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3), limit=c(-0.33, 0.30)) +
  scale_x_continuous(breaks=c(-0.2, 0.0, 0.2, 0.4), limit=c(-0.25, 0.46)) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text =  element_text(size = 14),
        strip.text = element_text(size=16, face="bold"),
        strip.background = element_blank(),
        plot.caption = element_text(size=14, hjust=0.95, vjust = 25),
        axis.line=element_line(),
        legend.direction= "horizontal",
        legend.position = "none",
        legend.background = element_rect(fill = "transparent"))
NMDS_18S_photo_merged_trajectories_2022_NS

NMDS_18S_photo_merged_trajectories_2022_NC <- scores_nmds_18S_photo_merged_2022 %>%
  filter(Flow %in% c("Natural", "Constant")) %>% 
  ggplot(aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape = Temperature), color = "darkgray", size = 3) +
  theme_classic() +
  geom_contour(data = ordi_18S_photo.mite.na, aes(x = x, y = y, z = z, colour=after_stat(level)), linewidth = 0.75, linetype = "dashed", color = "darkgray") +
  geom_label_contour(data = ordi_18S_photo.mite.na, aes(x = x, y = y, z = z), skip = 0, size = 4) +
  geom_segment(data=scores_nmds_18S_photo_merged_2022[scores_nmds_18S_photo_merged_2022$Treatment == "Nc",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#283747") +
  geom_segment(data=scores_nmds_18S_photo_merged_2022[scores_nmds_18S_photo_merged_2022$Treatment == "Nw",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#85929e") +
  geom_segment(data=scores_nmds_18S_photo_merged_2022[scores_nmds_18S_photo_merged_2022$Treatment == "Cc",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#2874a6") +
  geom_segment(data=scores_nmds_18S_photo_merged_2022[scores_nmds_18S_photo_merged_2022$Treatment == "Cw",], aes(xend = lead(NMDS1), yend = lead(NMDS2)), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#5dade2") +
  labs(caption="Stress 0.121") +
  labs(color = "Days of growth") +
  scale_y_continuous(breaks=c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3), limit=c(-0.33, 0.30)) +
  scale_x_continuous(breaks=c(-0.2, 0.0, 0.2, 0.4), limit=c(-0.25, 0.46)) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text =  element_text(size = 14),
        strip.text = element_text(size=16, face="bold"),
        strip.background = element_blank(),
        plot.caption = element_text(size=14, hjust=1, vjust = 15),
        axis.line=element_line(),
        legend.direction= "horizontal",
        legend.position = "none",
        legend.background = element_rect(fill = "transparent"))
NMDS_18S_photo_merged_trajectories_2022_NC


Figure3B <- ggarrange(top_title3, NMDS_18S_photo_merged_trajectories_2022_NI + rremove("ylab") + rremove("xlab"), NMDS_18S_photo_merged_trajectories_2022_NS + rremove("ylab") + rremove("xlab"), NMDS_18S_photo_merged_trajectories_2022_NC + rremove("ylab") + rremove("xlab"), ncol = 1, nrow = 4, heights = c(0.1, 1, 1, 1))
Figure3B  <- annotate_figure(Figure3B, left = text_grob("NMDS2", rot = 90, vjust = 1, size = 16, face = "bold"),
                             bottom = text_grob("NMDS1", size = 16, face = "bold"))
Figure3B

#Merging Figure 3
Figure3 <- ggarrange(Figure3A, Figure3B, ncol = 2, nrow = 1, labels = c("A", "B"), font.label = list(size = 18), common.legend = F)
Figure3

ggsave("Figure3.pdf", Figure3, width = 35,
       height = 45,
       units = "cm")


#Plot NMDS of Bray-Curtis dissimilarity of non-phototrophic eukaryotes, on Hellinger-transformed merged data (Figure S11B)
factor(sample_data(Sentinel18S_2022_others_filt_merged_hell)$Treatment, levels = c("Nc", "Nw", "Ic", "Iw", "Sc", "Sw", "Cc", "Cw"))

matrix_dist_18S_others_merged_2022 <- as.matrix(t(data.frame(otu_table(Sentinel18S_2022_others_filt_merged_hell))))
metadata_dist_18S_others_merged_2022 <- data.frame(sample_data(Sentinel18S_2022_others_filt_merged_hell)) %>% rownames_to_column("Samples")
set.seed(19950930)
dist_18S_others_merged_2022_bray <- vegdist(matrix_dist_18S_others_merged_2022, method ="bray") %>% 
  as.matrix() %>% 
  data.frame() %>% 
  rownames_to_column("Samples")
dist_18S_others_merged_2022 <- dist_18S_others_merged_2022_bray %>% 
  dplyr::select(all_of(.[["Samples"]])) %>% 
  as.dist()
nmds_18S_others_merged_2022 <- metaMDS(dist_18S_others_merged_2022, trymax=100)
stress <- nmds_18S_others_merged_2022$stress
stress

scores_nmds_18S_others_merged_2022 <- scores(nmds_18S_others_merged_2022) %>% 
  as_tibble(rownames = "Samples") %>% 
  inner_join(., metadata_dist_18S_others_merged_2022, by="Samples")

scores_nmds_18S_others_merged_2022  <- scores_nmds_18S_others_merged_2022[order(scores_nmds_18S_others_merged_2022$Treatment, scores_nmds_18S_others_merged_2022$Time),]

ordi_18S_others <- ordisurf(nmds_18S_others_merged_2022, metadata_dist_18S_others_merged_2022$Succession, plot = FALSE, bs="ds")
ordi_18S_others.grid <- ordi_18S_others$grid #extracts the ordisurf object
ordi_18S_others.mite <- expand.grid(x = ordi_18S_others.grid$x, y = ordi_18S_others.grid$y) #get x and ys
ordi_18S_others.mite$z <- as.vector(ordi_18S_others.grid$z) #unravel the matrix for the z scores
ordi_18S_others.mite.na <- data.frame(na.omit(ordi_18S_others.mite)) #gets rid of the nas

NMDS_18S_others_merged_trajectories_2022 <- scores_nmds_18S_others_merged_2022 %>%
  ggplot(aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(shape = Temperature), color = "darkgray", size = 3) +
  theme_classic() +
  geom_contour(data = ordi_18S_others.mite.na, aes(x = x, y = y, z = z, colour=after_stat(level)), linewidth = 0.75, linetype = "dashed", color = "darkgray") +
  geom_label_contour(data = ordi_18S_others.mite.na, aes(x = x, y = y, z = z), skip = 0, size = 4) +
  geom_segment(data=scores_nmds_18S_others_merged_2022[scores_nmds_18S_others_merged_2022$Treatment == "Nc",], aes(xend = after_stat(lead(x)), yend = after_stat(lead(y))), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#283747") +
  geom_segment(data=scores_nmds_18S_others_merged_2022[scores_nmds_18S_others_merged_2022$Treatment == "Nw",], aes(xend = after_stat(lead(x)), yend = after_stat(lead(y))), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#85929e") +
  geom_segment(data=scores_nmds_18S_others_merged_2022[scores_nmds_18S_others_merged_2022$Treatment == "Ic",], aes(xend = after_stat(lead(x)), yend = after_stat(lead(y))), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#b03a2e") +
  geom_segment(data=scores_nmds_18S_others_merged_2022[scores_nmds_18S_others_merged_2022$Treatment == "Iw",], aes(xend = after_stat(lead(x)), yend = after_stat(lead(y))), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#ec7063") +
  geom_segment(data=scores_nmds_18S_others_merged_2022[scores_nmds_18S_others_merged_2022$Treatment == "Sc",], aes(xend = after_stat(lead(x)), yend = after_stat(lead(y))), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#1e8449") +
  geom_segment(data=scores_nmds_18S_others_merged_2022[scores_nmds_18S_others_merged_2022$Treatment == "Sw",], aes(xend = after_stat(lead(x)), yend = after_stat(lead(y))), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#52be80") +
  geom_segment(data=scores_nmds_18S_others_merged_2022[scores_nmds_18S_others_merged_2022$Treatment == "Cc",], aes(xend = after_stat(lead(x)), yend = after_stat(lead(y))), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#2874a6") +
  geom_segment(data=scores_nmds_18S_others_merged_2022[scores_nmds_18S_others_merged_2022$Treatment == "Cw",], aes(xend = after_stat(lead(x)), yend = after_stat(lead(y))), 
               arrow = arrow(length = unit(3, "mm")), linewidth = 0.75, color = "#5dade2") +
  labs(caption="Stress 0.089") +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text =  element_text(size = 14),
        strip.text = element_text(size=16, face="bold"),
        strip.background = element_blank(),
        plot.caption = element_text(size=14, hjust=0.95, vjust = 25),
        axis.line=element_line(),
        legend.direction= "horizontal",
        legend.position = "none",
        legend.background = element_rect(fill = "transparent"))
NMDS_18S_others_merged_trajectories_2022


FigureS11 <- ggarrange(Shannon_18S_plot_others, NMDS_18S_others_merged_trajectories_2022, ncol = 1, nrow = 2, labels = c("A", "B"), font.label = list(size = 18), common.legend = F)
FigureS11

ggsave("FigureS11.pdf", FigureS11, width = 20,
       height = 25,
       units = "cm")

#Figure S13
matrix_dist_18S_photo_2022 <- as.matrix(t(data.frame(otu_table(Sentinel18S_2022_photo_filt_hell))))
metadata_dist_18S_photo_2022 <- data.frame(sample_data(Sentinel18S_2022_photo_filt_hell)) %>% rownames_to_column("Samples")
set.seed(19950930)
dist_18S_photo_2022_bray <- vegdist(matrix_dist_18S_photo_2022, method ="bray") %>% 
  as.matrix() %>% 
  data.frame() %>% 
  rownames_to_column("Samples")
dist_18S_photo_2022 <- dist_18S_photo_2022_bray %>% 
  dplyr::select(all_of(.[["Samples"]])) %>% 
  as.dist()
nmds_18S_photo_2022 <- metaMDS(dist_18S_photo_2022, trymax=100)
stress_18S_photo_all <- nmds_18S_photo_2022$stress
stress_18S_photo_all

scores_nmds_18S_photo_2022 <- scores(nmds_18S_photo_2022) %>% 
  as_tibble(rownames = "Samples") %>% 
  inner_join(., metadata_dist_18S_photo_2022, by="Samples")

scores_nmds_18S_photo_2022  <- scores_nmds_18S_photo_2022[order(scores_nmds_18S_photo_2022$Treatment, scores_nmds_18S_photo_2022$Time),]

env_18S_photo <- data.frame(sample_data(Sentinel18S_2022_photo_filt_hell)) %>% 
  rownames_to_column("Samples") %>% 
  select(-c("Samples", "Sample", "Run", "Flow", "FlumeID", "Type", "Primers", "Temperature", "Time","Year","Treatment","NO3", "NO2", "NH4", "PO4")) %>% 
  mutate(TAD = coalesce(TAD, 0))
set.seed(19950930)
fit_18S_photo <- envfit(nmds_18S_photo_2022, env_18S_photo, permutations = 999, na.rm = TRUE)

fit_18S_photo_vec <- data.frame(fit_18S_photo$vectors$arrows * sqrt(fit_18S_photo$vectors$r), P = fit_18S_photo$vectors$pvals)
fit_18S_photo_vec$Environment <- rownames(fit_18S_photo_vec)
fit_18S_photo_fac <- data.frame(fit_18S_photo$factors$centroids * sqrt(fit_18S_photo$factors$r), P = fit_18S_photo$factors$pvals)
fit_18S_photo_fac$Environment <- rownames(fit_18S_photo_fac) 
fit_18S_photo_fac <- fit_18S_photo_fac %>% slice(2L) %>% mutate(Environment = "Drought")
fit_18S_photo_env <- rbind(fit_18S_photo_vec,fit_18S_photo_fac)
fit_18S_photo_env_pv <- subset(fit_18S_photo_env, P<=0.05)

fit_18S_photo_env_pv_star <- fit_18S_photo_env_pv %>%
  group_by(Environment) %>%
  mutate(EnvironmentSign = case_when(
    Environment == "O2" ~ "O[2]*\"*\"",
    Environment %in% c("C25", "Temp", "Turbidity", "Fluoride", "Chloride", "Nitrite", "Bromide",
                       "Sulfate", "Sodium", "Magnesium", "Potassium", "Calcium",
                       "Strontium") ~ paste0("\"", Environment, "*\""),
    TRUE ~ paste0("\"", Environment, "\"")
  )) %>%
  mutate(EnvironmentSign = as.character(EnvironmentSign))

NMDS_18S_photo_2022 <- ggplot(scores_nmds_18S_photo_2022 , aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = Succession, shape = Flow), size = 2, alpha = 0.8) +
  scale_color_viridis(discrete=FALSE, option="magma", breaks= c(0, 20, 40, 60, 80, 100)) +
  theme_classic() +
  guides(colour = guide_colourbar(barheight = 10)) +
  labs(color="Days of growth", caption="Stress 0.144") +
  scale_x_continuous(limits=c(-0.2,0.3)) +
  scale_y_continuous(limits=c(-0.25,0.2)) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 16), 
        axis.title.y = element_text(size=16, face="bold"), 
        axis.title.x = element_text(size=16, face="bold"),
        axis.text =  element_text(size = 12),
        strip.background = element_blank(),
        plot.caption = element_text(size=16, hjust=1, vjust = 20),
        strip.text = element_text(size=16, face="bold"),
        axis.line=element_line()) +
  geom_segment(data = fit_18S_photo_env_pv_star, aes(x = 0, xend=NMDS1*0.3, y=0, yend=NMDS2*0.3), arrow=arrow(length=unit(0.2, "cm")), colour = "darkgrey", lwd=0.5) +
  ggrepel::geom_text_repel(data = fit_18S_photo_env_pv_star, aes(x=NMDS1*0.3, y=NMDS2*0.3, label = EnvironmentSign), parse = TRUE, cex = 5, direction = "x", colour = "black", segment.size = 0.25)
NMDS_18S_photo_2022

ggsave("FigureS13.pdf", NMDS_18S_photo_2022, width = 40,
       height = 35,
       units = "cm")

#Figure S14
matrix_dist_18S_others_2022 <- as.matrix(t(data.frame(otu_table(Sentinel18S_2022_others_filt_hell))))
metadata_dist_18S_others_2022 <- data.frame(sample_data(Sentinel18S_2022_others_filt_hell)) %>% rownames_to_column("Samples")
set.seed(19950930)
dist_18S_others_2022_bray <- vegdist(matrix_dist_18S_others_2022, method ="bray") %>% 
  as.matrix() %>% 
  data.frame() %>% 
  rownames_to_column("Samples")
dist_18S_others_2022 <- dist_18S_others_2022_bray %>% 
  dplyr::select(all_of(.[["Samples"]])) %>% 
  as.dist()
nmds_18S_others_2022 <- metaMDS(dist_18S_others_2022, trymax=100)
stress_18S_others_all <- nmds_18S_others_2022$stress
stress_18S_others_all

scores_nmds_18S_others_2022 <- scores(nmds_18S_others_2022) %>% 
  as_tibble(rownames = "Samples") %>% 
  inner_join(., metadata_dist_18S_others_2022, by="Samples")

scores_nmds_18S_others_2022  <- scores_nmds_18S_others_2022[order(scores_nmds_18S_others_2022$Treatment, scores_nmds_18S_others_2022$Time),]

env_18S_others <- data.frame(sample_data(Sentinel18S_2022_others_filt_hell)) %>% 
  rownames_to_column("Samples") %>% 
  select(-c("Samples", "Sample", "Run", "Flow", "FlumeID", "Type", "Primers", "Temperature", "Time","Year","Treatment","NO3", "NO2", "NH4", "PO4")) %>% 
  mutate(TAD = coalesce(TAD, 0))
set.seed(19950930)
fit_18S_others <- envfit(nmds_18S_others_2022, env_18S_others, permutations = 999, na.rm = TRUE)

fit_18S_others_vec <- data.frame(fit_18S_others$vectors$arrows * sqrt(fit_18S_others$vectors$r), P = fit_18S_others$vectors$pvals)
fit_18S_others_vec$Environment <- rownames(fit_18S_others_vec)
fit_18S_others_fac <- data.frame(fit_18S_others$factors$centroids * sqrt(fit_18S_others$factors$r), P = fit_18S_others$factors$pvals)
fit_18S_others_fac$Environment <- rownames(fit_18S_others_fac) 
fit_18S_others_fac <- fit_18S_others_fac %>% slice(2L) %>% mutate(Environment = "Drought")
fit_18S_others_env <- rbind(fit_18S_others_vec,fit_18S_photo_fac)
fit_18S_others_env_pv <- subset(fit_18S_others_env, P<=0.05)

fit_18S_others_env_pv_star <- fit_18S_others_env_pv %>%
  group_by(Environment) %>%
  mutate(EnvironmentSign = case_when(
    Environment == "O2" ~ "O[2]*\"*\"",
    Environment %in% c("C25", "Temp", "Turbidity", "Fluoride", "Chloride", "Nitrite", "Bromide",
                       "Sulfate", "Sodium", "Magnesium", "Potassium", "Calcium",
                       "Strontium") ~ paste0("\"", Environment, "*\""),
    TRUE ~ paste0("\"", Environment, "\"")
  )) %>%
  mutate(EnvironmentSign = as.character(EnvironmentSign))

NMDS_18S_others_2022 <- ggplot(scores_nmds_18S_others_2022 , aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = Succession, shape = Flow), size = 2, alpha = 0.8) +
  scale_color_viridis(discrete=FALSE, option="mako", breaks= c(0, 20, 40, 60, 80, 100)) +
  theme_classic() +
  guides(colour = guide_colourbar(barheight = 10)) +
  labs(color="Days of growth", caption="Stress 0.115") +
  scale_x_continuous(limits=c(-0.4,0.5)) +
  scale_y_continuous(limits=c(-0.2,0.3)) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 16), 
        axis.title.y = element_text(size=16, face="bold"), 
        axis.title.x = element_text(size=16, face="bold"),
        axis.text =  element_text(size = 12),
        strip.background = element_blank(),
        plot.caption = element_text(size=16, hjust=1, vjust = 20),
        strip.text = element_text(size=16, face="bold"),
        axis.line=element_line()) +
  geom_segment(data = fit_18S_others_env_pv_star, aes(x = 0, xend=NMDS1*0.3, y=0, yend=NMDS2*0.3), arrow=arrow(length=unit(0.2, "cm")), colour = "darkgrey", lwd=0.5) +
  ggrepel::geom_text_repel(data = fit_18S_others_env_pv_star, aes(x=NMDS1*0.3, y=NMDS2*0.3, label = EnvironmentSign), parse = TRUE, cex = 5, direction = "x", colour = "black", segment.size = 0.25)
NMDS_18S_others_2022

ggsave("FigureS14.pdf", NMDS_18S_others_2022, width = 40,
       height = 35,
       units = "cm")



### Statistics on patterns of succession ------------------------------------

##BA
#Is the data normally distributed? (yes: p > 0.05) If yes: ANOVA. If no: Kruskal-Wallis
shapiro.test(metadata_16S$BA) #Not normal
#Is their a significant group? (yes: p < 0.05)
#anova <- aov(Shannon ~ Temperature*Flow*Succession, data = metadata_16S)
#summary(anova)
#Is their a significant group? (yes: p < 0.05)
kruskal.test(BA ~ Succession, data = metadata_16S) #Yes 
kruskal.test(BA ~ Flow, data = metadata_16S) #No
kruskal.test(BA ~ Temperature, data = metadata_16S) #No
kruskal.test(BA ~ Treatment, data = metadata_16S) #No interaction

##16S alpha diversity
shapiro.test(alphadiv_16S_df$Shannon) #Not normal
kruskal.test(Shannon ~ Succession, data = alphadiv_16S_df) #Yes
kruskal.test(Shannon ~ Flow, data = alphadiv_16S_df) #No
kruskal.test(Shannon ~ Temperature, data = alphadiv_16S_df) #No
kruskal.test(Shannon ~ Treatment, data = alphadiv_16S_df) #No interaction

##16S beta diversity
adonis_bc_16S <- adonis2(matrix_dist_16S_2022 ~ Succession+Temperature+Flow+Treatment, data=metadata_dist_16S_2022, permutations=999, strata=metadata_dist_16S_2022$FlumeID, method="bray")
adonis_bc_16S #Succession, Flow, Temperature and Treatment are significant 
set.seed(19950930)
pairwise.adonis_16S <- pairwise.adonis(matrix_dist_16S_2022, metadata_dist_16S_2022$Treatment, sim.function='vegdist', sim.method='bray',p.adjust.m='BH', perm = 999)
pairwise.adonis_16S

##Betadisper
dist_16S <- vegdist(matrix_dist_16S_2022, method ="bray")
#Succession
disp_16S_Succession <- betadisper(dist_16S, metadata_dist_16S_2022$Succession, type=c("centroid"))
anov_disp_16S_Succession <- anova(disp_16S_Succession)
anov_disp_16S_Succession #Significant, so homogeneity of dispersion is not assumed. Should I trust adonis?
tukey_bdisp_mr<- TukeyHSD(disp_16S_Succession)
#Flow
disp_16S_Flow <- betadisper(dist_16S, metadata_dist_16S_2022$Flow, type=c("centroid"))
anov_disp_16S_Flow <- anova(disp_16S_Flow)
anov_disp_16S_Flow #Not significant, so homogeneity of dispersion is assumed. I can trust adonis.
#Temperature
disp_16S_Temperature <- betadisper(dist_16S, metadata_dist_16S_2022$Temperature, type=c("centroid"))
anov_disp_16S_Temperature <- anova(disp_16S_Temperature)
anov_disp_16S_Temperature #Not significant, so homogeneity of dispersion is assumed. I can trust adonis.
#Treatment
disp_16S_Treatment <- betadisper(dist_16S, metadata_dist_16S_2022$Treatment, type=c("centroid"))
anov_disp_16S_Treatment <- anova(disp_16S_Treatment)
anov_disp_16S_Treatment #Not significant, so homogeneity of dispersion is assumed. I can trust adonis.

###Phototrophs
##ChlA
shapiro.test(metadata_18S_photo$ChlA) #Not normal
kruskal.test(ChlA ~ Succession, data = metadata_18S_photo) #Yes 
kruskal.test(ChlA ~ Flow, data = metadata_18S_photo) #Yes
kruskal.test(ChlA ~ Temperature, data = metadata_18S_photo) #Yes
kruskal.test(ChlA ~ Treatment, data = metadata_18S_photo)

#Which group are different (p < 0.05)?
wilcox_observed_ChlA_Flow <- pairwise.wilcox.test(metadata_18S_photo$ChlA, metadata_18S_photo$Treatment, p.adjust.method = "BH")
wilcox_observed_ChlA_Flow #Natural vs Intermittent

##18S alpha diversity
shapiro.test(alphadiv_18S_df_photo$Shannon) #Not normal
kruskal.test(Shannon ~ Succession, data = alphadiv_18S_df_photo) #Yes
kruskal.test(Shannon ~ Flow, data = alphadiv_18S_df_photo) #No
kruskal.test(Shannon ~ Temperature, data = alphadiv_18S_df_photo) #No
kruskal.test(Shannon ~ Treatment, data = alphadiv_18S_df_photo) #No interaction

##18S beta diversity
adonis_bc_18S_photo <- adonis2(matrix_dist_18S_photo_2022 ~ Succession+Temperature+Flow+Treatment, data=metadata_dist_18S_photo_2022, permutations=999, strata=metadata_dist_18S_photo_2022$FlumeID, method="bray")
adonis_bc_18S_photo #Succession, Flow, Temperature and Treatment are significant 
set.seed(19950930)
pairwise.adonis_18S_photo <- pairwise.adonis(matrix_dist_18S_photo_2022, metadata_dist_18S_photo_2022$Treatment, sim.function='vegdist', sim.method='bray',p.adjust.m='BH', perm = 999)
pairwise.adonis_18S_photo

##Betadisper
dist_18S_photo <- vegdist(matrix_dist_18S_photo_2022, method ="bray")
#Succession
disp_18S_photo_Succession <- betadisper(dist_18S_photo, metadata_dist_18S_photo_2022$Succession, type=c("centroid"))
anov_disp_18S_photo_Succession <- anova(disp_18S_photo_Succession)
anov_disp_18S_photo_Succession #Significant, so homogeneity of dispersion is not assumed. Should I trust adonis?
#Flow
disp_18S_photo_Flow <- betadisper(dist_18S_photo, metadata_dist_18S_photo_2022$Flow, type=c("centroid"))
anov_disp_18S_photo_Flow <- anova(disp_18S_photo_Flow)
anov_disp_18S_photo_Flow #Not significant, so homogeneity of dispersion is assumed. I can trust adonis.
#Temperature
disp_18S_photo_Temperature <- betadisper(dist_18S_photo, metadata_dist_18S_photo_2022$Temperature, type=c("centroid"))
anov_disp_18S_photo_Temperature <- anova(disp_18S_photo_Temperature)
anov_disp_18S_photo_Temperature #Significant, so homogeneity of dispersion is not assumed. Should I trust adonis?


###Non-phototrophs

##18S alpha diversity
shapiro.test(alphadiv_18S_df_others$Shannon) #Not normal
kruskal.test(Shannon ~ Succession, data = alphadiv_18S_df_others) #Yes
kruskal.test(Shannon ~ Flow, data = alphadiv_18S_df_others) #No
kruskal.test(Shannon ~ Temperature, data = alphadiv_18S_df_others) #No
kruskal.test(Shannon ~ Treatment, data = alphadiv_18S_df_others) #No interactions

##18S beta diversity
adonis_bc_18S_others <- adonis2(matrix_dist_18S_others_2022 ~ Succession+Temperature+Flow+Treatment, data=metadata_dist_18S_others_2022, permutations=999, strata=metadata_dist_18S_others_2022$FlumeID, method="bray")
adonis_bc_18S_others #Succession, Flow, Temperature and Treatment are significant 
set.seed(19950930)
pairwise.adonis_18S_others <- pairwise.adonis(matrix_dist_18S_others_2022, metadata_dist_18S_others_2022$Treatment, sim.function='vegdist', sim.method='bray',p.adjust.m='BH', perm = 999)
pairwise.adonis_18S_others

##Betadisper
dist_18S_others <- vegdist(matrix_dist_18S_others_2022, method ="bray")
#Succession
disp_18S_others_Succession <- betadisper(dist_18S_others, metadata_dist_18S_others_2022$Succession, type=c("centroid"))
anov_disp_18S_others_Succession <- anova(disp_18S_others_Succession)
anov_disp_18S_others_Succession #Significant, so homogeneity of dispersion is not assumed. Should I trust adonis?
#Flow
disp_18S_others_Flow <- betadisper(dist_18S_others, metadata_dist_18S_others_2022$Flow, type=c("centroid"))
anov_disp_18S_others_Flow <- anova(disp_18S_others_Flow)
anov_disp_18S_others_Flow #Not significant, so homogeneity of dispersion is assumed. I can trust adonis.
#Temperature
disp_18S_others_Temperature <- betadisper(dist_18S_others, metadata_dist_18S_others_2022$Temperature, type=c("centroid"))
anov_disp_18S_others_Temperature <- anova(disp_18S_others_Temperature)
anov_disp_18S_others_Temperature #Not significant, so homogeneity of dispersion is assumed. I can trust adonis.

  
### Differential abundance --------------------------------------------------
library(ANCOMBC)
library(DT)

##16S

#Glom at the genus level
Sentinel16S_2022_bacteria_filt_genus <- tax_glom(Sentinel16S_2022_bacteria_filt, "Genus")
Sentinel16S_2022_bacteria_filt_genus <- Sentinel16S_2022_bacteria_filt_genus %>% 
  ps_mutate(Incubation = as.numeric(Succession), FlumeID = as.numeric(FlumeID))
Sentinel16S_2022_bacteria_filt_genus

#Flow DA
#Filter per temperature
Sentinel16S_2022_bacteria_filt_genus_flowincontrol <- ps_filter(Sentinel16S_2022_bacteria_filt_genus, Temperature == "control")
Sentinel16S_2022_bacteria_filt_genus_flowinwarm <- ps_filter(Sentinel16S_2022_bacteria_filt_genus, Temperature == "warm")

set.seed(19950930)
TSE_bacteria_2022_genus_flowincontrol <- mia::makeTreeSummarizedExperimentFromPhyloseq(Sentinel16S_2022_bacteria_filt_genus_flowincontrol)
TSE_bacteria_2022_genus_flowincontrol$Flow <- factor(TSE_bacteria_2022_genus_flowincontrol$Flow, levels = c("Natural", "Intermittent", "Stochastic", "Constant"))
TSE_bacteria_2022_genus_flowincontrol$FlumeID <- factor(TSE_bacteria_2022_genus_flowincontrol$FlumeID, levels = c("1", "2", "3", "4", "5", "6","7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))

TSE_bacteria_2022_genus_flowinwarm <- mia::makeTreeSummarizedExperimentFromPhyloseq(Sentinel16S_2022_bacteria_filt_genus_flowinwarm)
TSE_bacteria_2022_genus_flowinwarm$Flow <- factor(TSE_bacteria_2022_genus_flowinwarm$Flow, levels = c("Natural", "Intermittent", "Stochastic", "Constant"))
TSE_bacteria_2022_genus_flowinwarmFlumeID <- factor(TSE_bacteria_2022_genus_flowinwarm$FlumeID, levels = c("1", "2", "3", "4", "5", "6","7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))

ancombc_16S_flowincontrol_genus <- ancombc2(data = TSE_bacteria_2022_genus_flowincontrol, assay_name = "counts", tax_level = "Genus",
                                  fix_formula = "Flow", rand_formula = "((1|Succession) + (1|FlumeID))",
                                  p_adj_method = "BH",
                                  prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                  group = "Flow", struc_zero = TRUE, neg_lb = TRUE,
                                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                                  iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                                  em_control = list(tol = 1e-5, max_iter = 100),
                                  lme_control = lme4::lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                  mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
                                  trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100))
res_dunn_16S_flowincontrol_genus <- ancombc_16S_flowincontrol_genus$res_dunn

ancombc_16S_flowinwarm_genus <- ancombc2(data = TSE_bacteria_2022_genus_flowinwarm, assay_name = "counts", tax_level = "Genus",
                                  fix_formula = "Flow", rand_formula = "((1|Succession) + (1|FlumeID))",
                                  p_adj_method = "BH",
                                  prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                  group = "Flow", struc_zero = TRUE, neg_lb = TRUE,
                                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                                  iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                                  em_control = list(tol = 1e-5, max_iter = 100),
                                  lme_control = lme4::lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                  mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
                                  trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100))
res_dunn_16S_flowinwarm_genus <- ancombc_16S_flowinwarm_genus$res_dunn

#Temperature DA

#Filter per Flow
Sentinel16S_2022_bacteria_filt_genus_tempinnatural <- ps_filter(Sentinel16S_2022_bacteria_filt_genus, Flow == "Natural")
Sentinel16S_2022_bacteria_filt_genus_tempinintermittent <- ps_filter(Sentinel16S_2022_bacteria_filt_genus, Flow == "Intermittent")
Sentinel16S_2022_bacteria_filt_genus_tempinstochastic <- ps_filter(Sentinel16S_2022_bacteria_filt_genus, Flow == "Stochastic")
Sentinel16S_2022_bacteria_filt_genus_tempinconstant <- ps_filter(Sentinel16S_2022_bacteria_filt_genus, Flow == "Constant")

set.seed(19950930)
TSE_bacteria_2022_genus_tempinnatural <- mia::makeTreeSummarizedExperimentFromPhyloseq(Sentinel16S_2022_bacteria_filt_genus_tempinnatural)
TSE_bacteria_2022_genus_tempinnatural$Temperature <- factor(TSE_bacteria_2022_genus_tempinnatural$Temperature, levels = c("control", "warm"))
TSE_bacteria_2022_genus_tempinnatural$FlumeID <- factor(TSE_bacteria_2022_genus_tempinnatural$FlumeID, levels = c("1", "2", "3", "4", "5", "6","7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))

TSE_bacteria_2022_genus_tempinintermittent <- mia::makeTreeSummarizedExperimentFromPhyloseq(Sentinel16S_2022_bacteria_filt_genus_tempinintermittent)
TSE_bacteria_2022_genus_tempinintermittent$Temperature <- factor(TSE_bacteria_2022_genus_tempinintermittent$Temperature, levels = c("control", "warm"))
TSE_bacteria_2022_genus_tempinintermittent$FlumeID <- factor(TSE_bacteria_2022_genus_tempinintermittent$FlumeID, levels = c("1", "2", "3", "4", "5", "6","7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))

TSE_bacteria_2022_genus_tempinstochastic <- mia::makeTreeSummarizedExperimentFromPhyloseq(Sentinel16S_2022_bacteria_filt_genus_tempinstochastic)
TSE_bacteria_2022_genus_tempinstochastic$Temperature <- factor(TSE_bacteria_2022_genus_tempinstochastic$Temperature, levels = c("control", "warm"))
TSE_bacteria_2022_genus_tempinstochastic$FlumeID <- factor(TSE_bacteria_2022_genus_tempinstochastic$FlumeID, levels = c("1", "2", "3", "4", "5", "6","7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))

TSE_bacteria_2022_genus_tempinconstant <- mia::makeTreeSummarizedExperimentFromPhyloseq(Sentinel16S_2022_bacteria_filt_genus_tempinconstant)
TSE_bacteria_2022_genus_tempinconstant$Temperature <- factor(TSE_bacteria_2022_genus_tempinconstant$Temperature, levels = c("control", "warm"))
TSE_bacteria_2022_genus_tempinconstant$FlumeID <- factor(TSE_bacteria_2022_genus_tempinconstant$FlumeID, levels = c("1", "2", "3", "4", "5", "6","7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))

ancombc_16S_tempinnatural_genus <- ancombc2(data = TSE_bacteria_2022_genus_tempinnatural, assay_name = "counts", tax_level = "Genus",
                                   fix_formula = "Temperature", rand_formula = "((1|Succession) + (1|FlumeID))",
                                   p_adj_method = "BH",
                                   prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                   group = "Temperature", struc_zero = TRUE, neg_lb = TRUE,
                                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                                   global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                                   iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                                   em_control = list(tol = 1e-5, max_iter = 100),
                                   lme_control = lme4::lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                   mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
                                   trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100))
res_dunn_16S_tempinnatural_genus <- ancombc_16S_tempinnatural_genus$res

ancombc_16S_tempinintermittent_genus <- ancombc2(data = TSE_bacteria_2022_genus_tempinintermittent, assay_name = "counts", tax_level = "Genus",
                                            fix_formula = "Temperature", rand_formula = "((1|Succession) + (1|FlumeID))",
                                            p_adj_method = "BH",
                                            prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                            group = "Temperature", struc_zero = TRUE, neg_lb = TRUE,
                                            alpha = 0.05, n_cl = 2, verbose = TRUE,
                                            global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                                            iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                                            em_control = list(tol = 1e-5, max_iter = 100),
                                            lme_control = lme4::lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                            mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
                                            trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100))
res_dunn_16S_tempinintermittent_genus <- ancombc_16S_tempinintermittent_genus$res

ancombc_16S_tempinstochastic_genus <- ancombc2(data = TSE_bacteria_2022_genus_tempinstochastic, assay_name = "counts", tax_level = "Genus",
                                                 fix_formula = "Temperature", rand_formula = "((1|Succession) + (1|FlumeID))",
                                                 p_adj_method = "BH",
                                                 prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                                 group = "Temperature", struc_zero = TRUE, neg_lb = TRUE,
                                                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                                                 global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                                                 iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                                                 em_control = list(tol = 1e-5, max_iter = 100),
                                                 lme_control = lme4::lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                                 mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
                                                 trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100))
res_dunn_16S_tempinstochastic_genus <- ancombc_16S_tempinstochastic_genus$res

ancombc_16S_tempinconstant_genus <- ancombc2(data = TSE_bacteria_2022_genus_tempinconstant, assay_name = "counts", tax_level = "Genus",
                                               fix_formula = "Temperature", rand_formula = "((1|Succession) + (1|FlumeID))",
                                               p_adj_method = "BH",
                                               prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                               group = "Temperature", struc_zero = TRUE, neg_lb = TRUE,
                                               alpha = 0.05, n_cl = 2, verbose = TRUE,
                                               global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                                               iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                                               em_control = list(tol = 1e-5, max_iter = 100),
                                               lme_control = lme4::lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                               mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
                                               trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100))
res_dunn_16S_tempinconstant_genus <- ancombc_16S_tempinconstant_genus$res

df_volcano_16S_flow <- left_join(res_dunn_16S_flowincontrol_genus, res_dunn_16S_flowinwarm_genus, by = "taxon") %>% 
  mutate(Intermittent.c = round(lfc_FlowIntermittent.x, 2),
         Stochastic.c = round(lfc_FlowStochastic.x, 2),
         Constant.c = round(lfc_FlowConstant.x, 2), 
         Intermittent.w = round(lfc_FlowIntermittent.y, 2),
         Stochastic.w = round(lfc_FlowStochastic.y, 2),
         Constant.w = round(lfc_FlowConstant.y, 2)) %>% 
  pivot_longer(cols = "Intermittent.c":"Constant.w", names_to = "Effect", values_to = "lfc") %>% 
  mutate(p_value = case_when(Effect == "Intermittent.c" ~ q_FlowIntermittent.x,
                             Effect == "Intermittent.w" ~ q_FlowIntermittent.y,
                             Effect == "Stochastic.c" ~ q_FlowStochastic.x,
                             Effect == "Stochastic.w" ~ q_FlowStochastic.y,
                             Effect == "Constant.c" ~ q_FlowConstant.x,
                             Effect == "Constant.w" ~ q_FlowConstant.y,
                             TRUE ~ 0),
         DA = case_when((Effect == "Intermittent.c" & diff_FlowIntermittent.x == 1 & passed_ss_FlowIntermittent.x == 1 & lfc >= 0) ~ "More",
                        (Effect == "Intermittent.w" & diff_FlowIntermittent.y == 1 & passed_ss_FlowIntermittent.y == 1 & lfc >= 0) ~ "More",
                        (Effect == "Stochastic.c" & diff_FlowStochastic.x == 1 & passed_ss_FlowStochastic.x == 1 & lfc >= 0) ~ "More",
                        (Effect == "Stochastic.w" & diff_FlowStochastic.y == 1 & passed_ss_FlowStochastic.y == 1 & lfc >= 0) ~ "More",
                        (Effect == "Constant.c" & diff_FlowConstant.x == 1 & passed_ss_FlowConstant.x == 1 & lfc >= 0) ~ "More",
                        (Effect == "Constant.w" & diff_FlowConstant.y == 1 & passed_ss_FlowConstant.y == 1 & lfc >= 0) ~ "More",
                        (Effect == "Intermittent.c" & diff_FlowIntermittent.x == 1 & passed_ss_FlowIntermittent.x == 1 & lfc <= 0) ~ "Less",
                        (Effect == "Intermittent.w" & diff_FlowIntermittent.y == 1 & passed_ss_FlowIntermittent.y == 1 & lfc <= 0) ~ "Less",
                        (Effect == "Stochastic.c" & diff_FlowStochastic.x == 1 & passed_ss_FlowStochastic.x == 1 & lfc <= 0) ~ "Less",
                        (Effect == "Stochastic.w" & diff_FlowStochastic.y == 1 & passed_ss_FlowStochastic.y == 1 & lfc <= 0) ~ "Less",
                        (Effect == "Constant.c" & diff_FlowConstant.x == 1 & passed_ss_FlowConstant.x == 1 & lfc <= 0) ~ "Less",
                        (Effect == "Constant.w" & diff_FlowConstant.y == 1 & passed_ss_FlowConstant.y == 1 & lfc <= 0) ~ "Less",
                        TRUE ~ "No"),
         Flow = case_when((Effect == "Intermittent.c" | Effect == "Intermittent.w") ~ "Intermittent",
                          (Effect == "Stochastic.c" | Effect == "Stochastic.w") ~ "Stochastic",
                          (Effect == "Constant.c" | Effect == "Constant.w") ~ "Constant",
                          TRUE ~ "NA"),
         Temp = paste0("", str_sub(Effect, start = -2)),
         Temperature = case_when((Temp == ".c") ~ "control",
                                 (Temp == ".w") ~ "warm",
                                 TRUE ~ "NA"),
         Treatment = case_when((Flow == "Intermittent" & Temperature == "control") ~ "Ic",
                               (Flow == "Intermittent" & Temperature == "warm") ~ "Iw",
                               (Flow == "Stochastic" & Temperature == "control") ~ "Sc",
                               (Flow == "Stochastic" & Temperature == "warm") ~ "Sw",
                               (Flow == "Constant" & Temperature == "control") ~ "Cc",
                               (Flow == "Constant" & Temperature == "warm") ~ "Cw",
                               TRUE ~ "NA")) %>%
  dplyr::rename("Genus" = "taxon") %>%
  select(Genus, Effect, lfc, p_value, DA, Flow, Temperature, Treatment)
df_volcano_16S_flow$Treatment <- factor(df_volcano_16S_flow$Treatment, levels = c("Ic", "Iw", "Sc", "Sw", "Cc", "Cw"))
df_volcano_16S_flow$DA <- factor(df_volcano_16S_flow$DA, levels = c("More", "Less", "No"))

#Count how many Genera and families DA in flow
df_volcano_16S_flow %>%
  left_join(taxonomy_16S) %>% 
  filter(DA != "No") %>%
  summarise(DA_genera = n_distinct(Genus), DA_phyla = n_distinct(Phylum))

DA_16S <- df_volcano_16S_flow %>%
  left_join(taxonomy_16S) %>% 
  filter(DA != "No") %>%
  group_by(Flow) %>% 
  filter(Flow == "Intermittent") %>% 
  summarise(DA_genera = n_distinct(Genus), DA_phyla = n_distinct(Phylum))

DA_intermittent <- df_volcano_16S_flow %>%
  left_join(taxonomy_16S) %>% 
  filter(DA != "No") %>%
  group_by(Flow) %>% 
  filter(Flow == "Intermittent") %>%
  distinct(Genus, .keep_all = TRUE) %>% 
  group_by(DA) %>%
  add_count(DA, sort = FALSE, name= "DA_frequency")

df_volcano_16S_flow %>%
  left_join(taxonomy_16S) %>% 
  filter(DA != "No") %>%
  group_by(Treatment) %>% 
  summarise(DA_genera = n_distinct(Genus), DA_phyla = n_distinct(Phylum))

volcano_DA_bacteria_flow <- df_volcano_16S_flow %>%
  filter(p_value != 1) %>% 
  ggplot(aes(x=lfc, y=-log10(p_value), col=Treatment, alpha= DA)) + 
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype="dashed", color = "firebrick", linewidth = 0.8, alpha = 0.5) +
  geom_label_repel(aes(label = Genus), size = 5, fontface = "italic", show.legend = FALSE, box.padding = 0.3, point.padding = 0.3, max.overlaps= 30, max.time = 10, data = (df_volcano_16S_flow %>% filter(DA != "No" & abs(lfc) > 0.80))) +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 14, by = 2)) +
  scale_color_manual(values=c("#b03a2e","#ec7063","#1e8449","#52be80", "#2874a6", "#5dade2"),
                     labels= c(bquote("N"[flow]*"C"[temp]*"  vs "*"I"[flow]*"C"[temp]), 
                               bquote("N"[flow]*"W"[temp]*" vs "*"I"[flow]*"W"[temp]), 
                               bquote("N"[flow]*"C"[temp]*"  vs "*"S"[flow]*"C"[temp]), 
                               bquote("N"[flow]*"W"[temp]*" vs "*"S"[flow]*"W"[temp]), 
                               bquote("N"[flow]*"C"[temp]*"  vs "*"C"[flow]*"C"[temp]), 
                               bquote("N"[flow]*"W"[temp]*" vs "*"C"[flow]*"W"[temp]))) +
  scale_alpha_manual(values=c(1, 1, 0.1), guide = "none") +
  labs(x = "log Fold-Change", y = "-log p-value", color= "Flow comparision") +
  theme(legend.title = element_text(size = 14, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text =  element_text(size = 14),
        legend.position = "right",
        legend.justification = "left") +
  xlim(-3,3)
volcano_DA_bacteria_flow

df_volcano_16S_temperature <- left_join(res_dunn_16S_tempinnatural_genus, res_dunn_16S_tempinintermittent_genus, by = "taxon") %>% 
  left_join(res_dunn_16S_tempinstochastic_genus,  by = "taxon") %>% 
  left_join(res_dunn_16S_tempinconstant_genus,  by = "taxon") %>% 
  mutate(Natural = round(lfc_Temperaturewarm.x, 2),
         Intermittent = round(lfc_Temperaturewarm.y, 2),
         Stochastic = round(lfc_Temperaturewarm.x.x, 2),
         Constant = round(lfc_Temperaturewarm.y.y, 2)) %>% 
  pivot_longer(cols = "Natural":"Constant", names_to = "Effect", values_to = "lfc") %>% 
  mutate(p_value = case_when(Effect == "Natural" ~ q_Temperaturewarm.x,
                             Effect == "Intermittent" ~ q_Temperaturewarm.y,
                             Effect == "Stochastic" ~ q_Temperaturewarm.x.x,
                             Effect == "Constant" ~ q_Temperaturewarm.y.y,
                             TRUE ~ 0),
         DA = case_when((Effect == "Natural" & diff_Temperaturewarm.x == 1 & passed_ss_Temperaturewarm.x == 1 & lfc >= 0) ~ "More",
                        (Effect == "Intermittent" & diff_Temperaturewarm.y == 1 & passed_ss_Temperaturewarm.y == 1 & lfc >= 0) ~ "More",
                        (Effect == "Stochastic" & diff_Temperaturewarm.x.x == 1 & passed_ss_Temperaturewarm.x.x == 1 & lfc >= 0) ~ "More",
                        (Effect == "Constant" & diff_Temperaturewarm.y.y == 1 & passed_ss_Temperaturewarm.y.y == 1 & lfc >= 0) ~ "More",
                        (Effect == "Natural" & diff_Temperaturewarm.x == 1 & passed_ss_Temperaturewarm.x == 1 & lfc <= 0) ~ "Less",
                        (Effect == "Intermittent" & diff_Temperaturewarm.y == 1 & passed_ss_Temperaturewarm.y == 1 & lfc <= 0) ~ "Less",
                        (Effect == "Stochastic" & diff_Temperaturewarm.x.x == 1 & passed_ss_Temperaturewarm.x.x == 1 & lfc <= 0) ~ "Less",
                        (Effect == "Constant" & diff_Temperaturewarm.y.y == 1 & passed_ss_Temperaturewarm.y.y == 1 & lfc <= 0) ~ "Less",
                        TRUE ~ "No")) %>% 
  dplyr::rename("Genus" = "taxon") %>% 
  select(Genus, Effect, lfc, p_value, DA)
df_volcano_16S_temperature$Effect <- factor(df_volcano_16S_temperature$Effect, levels = c("Natural", "Intermittent", "Stochastic", "Constant"))
df_volcano_16S_temperature$DA <- factor(df_volcano_16S_temperature$DA, levels = c("More", "Less", "No"))

#Count how many Genera and families DA in flow
df_volcano_16S_temperature %>%
  left_join(taxonomy_16S) %>% 
  filter(DA != "No") %>%
  summarise(DA_genera = n_distinct(Genus), DA_phyla = n_distinct(Phylum))

df_volcano_16S_temperature %>%
  left_join(taxonomy_16S) %>% 
  filter(DA != "No") %>%
  group_by(Effect) %>% 
  summarise(DA_genera = n_distinct(Genus), DA_phyla = n_distinct(Phylum))

DA_temperature <- df_volcano_16S_temperature %>%
  left_join(taxonomy_16S) %>% 
  filter(DA != "No") %>%
  distinct(Genus, .keep_all = TRUE) %>% 
  group_by(DA) %>%
  add_count(DA, sort = FALSE, name= "DA_frequency")

volcano_DA_bacteria_temperature <- df_volcano_16S_temperature %>% 
  filter(p_value != 1) %>% 
  ggplot(aes(x=lfc, y=-log10(p_value), col=Effect, alpha= DA)) + 
  geom_point(size = 3) +
  scale_y_continuous(breaks = seq(0, 3, by = 1)) +
  geom_vline(xintercept = 0, linetype="dashed", color = "firebrick", linewidth = 0.8, alpha = 0.5) +
  geom_label_repel(aes(label = Genus), size = 5, fontface = "italic", show.legend = FALSE,  box.padding = 0.8, point.padding = 0.8, data = df_volcano_16S_temperature[df_volcano_16S_temperature$DA %in% c("More", "Less"),]) +
  theme_classic() +
  scale_color_manual(values=c("#85929e", "#ec7063", "#52be80","#5dade2"),
                     labels= c(bquote("N"[flow]*"C"[temp]*" vs "*"N"[flow]*"W"[temp]), 
                               bquote("I"[flow]*"C"[temp]*"  vs  "*"I"[flow]*"W"[temp]), 
                               bquote("S"[flow]*"C"[temp]*" vs "*"S"[flow]*"W"[temp]), 
                               bquote("C"[flow]*"C"[temp]*" vs "*"C"[flow]*"W"[temp]))) +
  scale_alpha_manual(values=c(1, 1, 0.1), guide = "none") +
  labs(x = "log Fold-Change", y = "-log p-value", color= "Temperature comparison") +
  xlim(-1.5,1.5) +
  theme(legend.title = element_text(size = 14, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text =  element_text(size = 14),
        legend.position = "right")
volcano_DA_bacteria_temperature

Figure4 <- ggarrange(volcano_DA_bacteria_flow, volcano_DA_bacteria_temperature, ncol = 1, nrow = 2, labels = c("A", "B"), font.label = list(size = 18), align = "v", heights = c(1.3, 1))
Figure4
ggsave("Figure4.pdf", Figure4, width = 35,
       height = 45,
       units = "cm")

##18S

#Glom at the class level
Sentinel18S_2022_eukaryote_filt_class <- tax_glom(Sentinel18S_2022_eukaryote_filt, "Class")
Sentinel18S_2022_eukaryote_filt_class <- Sentinel18S_2022_eukaryote_filt_class %>% 
  ps_mutate(Succession = as.numeric(Succession), FlumeID = as.numeric(FlumeID))
Sentinel18S_2022_eukaryote_filt_class

#Flow DA
#Filter per temperature
Sentinel18S_2022_eukaryote_filt_class_flowincontrol <- ps_filter(Sentinel18S_2022_eukaryote_filt_class, Temperature == "control")
Sentinel18S_2022_eukaryote_filt_class_flowinwarm <- ps_filter(Sentinel18S_2022_eukaryote_filt_class, Temperature == "warm")

set.seed(19950930)
TSE_eukaryote_2022_class_flowincontrol <- mia::makeTreeSummarizedExperimentFromPhyloseq(Sentinel18S_2022_eukaryote_filt_class_flowincontrol)
TSE_eukaryote_2022_class_flowincontrol$Flow <- factor(TSE_eukaryote_2022_class_flowincontrol$Flow, levels = c("Natural", "Intermittent", "Stochastic", "Constant"))
TSE_eukaryote_2022_class_flowincontrol$FlumeID <- factor(TSE_eukaryote_2022_class_flowincontrol$FlumeID, levels = c("1", "2", "3", "4", "5", "6","7", "8", "9", "10", "11", "12"))

TSE_eukaryote_2022_class_flowinwarm <- mia::makeTreeSummarizedExperimentFromPhyloseq(Sentinel18S_2022_eukaryote_filt_class_flowinwarm)
TSE_eukaryote_2022_class_flowinwarm$Flow <- factor(TSE_eukaryote_2022_class_flowinwarm$Flow, levels = c("Natural", "Intermittent", "Stochastic", "Constant"))
TSE_eukaryote_2022_class_flowinwarm$FlumeID <- factor(TSE_eukaryote_2022_class_flowinwarm$FlumeID, levels = c("1", "2", "3", "4", "5", "6","7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))

ancombc_18S_flowincontrol_class <- ancombc2(data = TSE_eukaryote_2022_class_flowincontrol, assay_name = "counts", tax_level = "Class",
                                            fix_formula = "Flow", rand_formula = "((1|Succession) + (1|FlumeID))",
                                            p_adj_method = "BH",
                                            prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                            group = "Flow", struc_zero = TRUE, neg_lb = TRUE,
                                            alpha = 0.05, n_cl = 2, verbose = TRUE,
                                            global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                                            iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                                            em_control = list(tol = 1e-5, max_iter = 100),
                                            lme_control = lme4::lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                            mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
                                            trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100))
res_dunn_18S_flowincontrol_class <- ancombc_18S_flowincontrol_class$res_dunn

ancombc_18S_flowinwarm_class <- ancombc2(data = TSE_eukaryote_2022_class_flowinwarm, assay_name = "counts", tax_level = "Class",
                                         fix_formula = "Flow", rand_formula = "((1|Succession) + (1|FlumeID))",
                                         p_adj_method = "BH",
                                         prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                         group = "Flow", struc_zero = TRUE, neg_lb = TRUE,
                                         alpha = 0.05, n_cl = 2, verbose = TRUE,
                                         global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                                         iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                                         em_control = list(tol = 1e-5, max_iter = 100),
                                         lme_control = lme4::lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                         mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
                                         trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100))
res_dunn_18S_flowinwarm_class <- ancombc_18S_flowinwarm_class$res_dunn


#Temperature DA

#Filter per Flow
Sentinel18S_2022_eukaryote_filt_class_tempinnatural <- ps_filter(Sentinel18S_2022_eukaryote_filt_class, Flow == "Natural")
Sentinel18S_2022_eukaryote_filt_class_tempinintermittent <- ps_filter(Sentinel18S_2022_eukaryote_filt_class, Flow == "Intermittent")
Sentinel18S_2022_eukaryote_filt_class_tempinstochastic <- ps_filter(Sentinel18S_2022_eukaryote_filt_class, Flow == "Stochastic")
Sentinel18S_2022_eukaryote_filt_class_tempinconstant <- ps_filter(Sentinel18S_2022_eukaryote_filt_class, Flow == "Constant")

set.seed(19950930)
TSE_eukaryote_2022_class_tempinnatural <- mia::makeTreeSummarizedExperimentFromPhyloseq(Sentinel18S_2022_eukaryote_filt_class_tempinnatural)
TSE_eukaryote_2022_class_tempinnatural$Temperature <- factor(TSE_eukaryote_2022_class_tempinnatural$Temperature, levels = c("control", "warm"))
TSE_eukaryote_2022_class_tempinnatural$FlumeID <- factor(TSE_eukaryote_2022_class_tempinnatural$FlumeID, levels = c("1", "2", "3", "4", "5", "6","7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))

TSE_eukaryote_2022_class_tempinintermittent <- mia::makeTreeSummarizedExperimentFromPhyloseq(Sentinel18S_2022_eukaryote_filt_class_tempinintermittent)
TSE_eukaryote_2022_class_tempinintermittent$Temperature <- factor(TSE_eukaryote_2022_class_tempinintermittent$Temperature, levels = c("control", "warm"))
TSE_eukaryote_2022_class_tempinintermittent$FlumeID <- factor(TSE_eukaryote_2022_class_tempinintermittent$FlumeID, levels = c("1", "2", "3", "4", "5", "6","7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))

TSE_eukaryote_2022_class_tempinstochastic <- mia::makeTreeSummarizedExperimentFromPhyloseq(Sentinel18S_2022_eukaryote_filt_class_tempinstochastic)
TSE_eukaryote_2022_class_tempinstochastic$Temperature <- factor(TSE_eukaryote_2022_class_tempinstochastic$Temperature, levels = c("control", "warm"))
TSE_eukaryote_2022_class_tempinstochastic$FlumeID <- factor(TSE_eukaryote_2022_class_tempinstochastic$FlumeID, levels = c("1", "2", "3", "4", "5", "6","7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))

TSE_eukaryote_2022_class_tempinconstant <- mia::makeTreeSummarizedExperimentFromPhyloseq(Sentinel18S_2022_eukaryote_filt_class_tempinconstant)
TSE_eukaryote_2022_class_tempinconstant$Temperature <- factor(TSE_eukaryote_2022_class_tempinconstant$Temperature, levels = c("control", "warm"))
TSE_eukaryote_2022_class_tempinconstant$FlumeID <- factor(TSE_eukaryote_2022_class_tempinconstant$FlumeID, levels = c("1", "2", "3", "4", "5", "6","7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))

ancombc_18S_tempinnatural_class <- ancombc2(data = TSE_eukaryote_2022_class_tempinnatural, assay_name = "counts", tax_level = "Class",
                                            fix_formula = "Temperature", rand_formula = "((1|Succession) + (1|FlumeID))",
                                            p_adj_method = "BH",
                                            prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                            group = "Temperature", struc_zero = TRUE, neg_lb = TRUE,
                                            alpha = 0.05, n_cl = 2, verbose = TRUE,
                                            global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                                            iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                                            em_control = list(tol = 1e-5, max_iter = 100),
                                            lme_control = lme4::lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                            mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
                                            trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100))
res_dunn_18S_tempinnatural_class <- ancombc_18S_tempinnatural_class$res

ancombc_18S_tempinintermittent_class <- ancombc2(data = TSE_eukaryote_2022_class_tempinintermittent, assay_name = "counts", tax_level = "Class",
                                                 fix_formula = "Temperature", rand_formula = "((1|Succession) + (1|FlumeID))",
                                                 p_adj_method = "BH",
                                                 prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                                 group = "Temperature", struc_zero = TRUE, neg_lb = TRUE,
                                                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                                                 global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                                                 iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                                                 em_control = list(tol = 1e-5, max_iter = 100),
                                                 lme_control = lme4::lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                                 mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
                                                 trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100))
res_dunn_18S_tempinintermittent_class <- ancombc_18S_tempinintermittent_class$res

ancombc_18S_tempinstochastic_class <- ancombc2(data = TSE_eukaryote_2022_class_tempinstochastic, assay_name = "counts", tax_level = "Class",
                                               fix_formula = "Temperature", rand_formula = "((1|Succession) + (1|FlumeID))",
                                               p_adj_method = "BH",
                                               prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                               group = "Temperature", struc_zero = TRUE, neg_lb = TRUE,
                                               alpha = 0.05, n_cl = 2, verbose = TRUE,
                                               global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                                               iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                                               em_control = list(tol = 1e-5, max_iter = 100),
                                               lme_control = lme4::lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                               mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
                                               trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100))
res_dunn_18S_tempinstochastic_class <- ancombc_18S_tempinstochastic_class$res

ancombc_18S_tempinconstant_class <- ancombc2(data = TSE_eukaryote_2022_class_tempinconstant, assay_name = "counts", tax_level = "Class",
                                             fix_formula = "Temperature", rand_formula = "((1|Succession) + (1|FlumeID))",
                                             p_adj_method = "BH",
                                             prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                                             group = "Temperature", struc_zero = TRUE, neg_lb = TRUE,
                                             alpha = 0.05, n_cl = 2, verbose = TRUE,
                                             global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                                             iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                                             em_control = list(tol = 1e-5, max_iter = 100),
                                             lme_control = lme4::lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                             mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
                                             trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100))
res_dunn_18S_tempinconstant_class <- ancombc_18S_tempinconstant_class$res

df_volcano_18S_flow <- left_join(res_dunn_18S_flowincontrol_class, res_dunn_18S_flowinwarm_class, by = "taxon") %>% 
  mutate(Intermittent.c = round(lfc_FlowIntermittent.x, 2),
         Stochastic.c = round(lfc_FlowStochastic.x, 2),
         Constant.c = round(lfc_FlowConstant.x, 2), 
         Intermittent.w = round(lfc_FlowIntermittent.y, 2),
         Stochastic.w = round(lfc_FlowStochastic.y, 2),
         Constant.w = round(lfc_FlowConstant.y, 2)) %>% 
  pivot_longer(cols = "Intermittent.c":"Constant.w", names_to = "Effect", values_to = "lfc") %>% 
  mutate(p_value = case_when(Effect == "Intermittent.c" ~ q_FlowIntermittent.x,
                             Effect == "Intermittent.w" ~ q_FlowIntermittent.y,
                             Effect == "Stochastic.c" ~ q_FlowStochastic.x,
                             Effect == "Stochastic.w" ~ q_FlowStochastic.y,
                             Effect == "Constant.c" ~ q_FlowConstant.x,
                             Effect == "Constant.w" ~ q_FlowConstant.y,
                             TRUE ~ 0),
         DA = case_when((Effect == "Intermittent.c" & diff_FlowIntermittent.x == 1 & passed_ss_FlowIntermittent.x == 1 & lfc >= 0) ~ "More",
                        (Effect == "Intermittent.w" & diff_FlowIntermittent.y == 1 & passed_ss_FlowIntermittent.y == 1 & lfc >= 0) ~ "More",
                        (Effect == "Stochastic.c" & diff_FlowStochastic.x == 1 & passed_ss_FlowStochastic.x == 1 & lfc >= 0) ~ "More",
                        (Effect == "Stochastic.w" & diff_FlowStochastic.y == 1 & passed_ss_FlowStochastic.y == 1 & lfc >= 0) ~ "More",
                        (Effect == "Constant.c" & diff_FlowConstant.x == 1 & passed_ss_FlowConstant.x == 1 & lfc >= 0) ~ "More",
                        (Effect == "Constant.w" & diff_FlowConstant.y == 1 & passed_ss_FlowConstant.y == 1 & lfc >= 0) ~ "More",
                        (Effect == "Intermittent.c" & diff_FlowIntermittent.x == 1 & passed_ss_FlowIntermittent.x == 1 & lfc <= 0) ~ "Less",
                        (Effect == "Intermittent.w" & diff_FlowIntermittent.y == 1 & passed_ss_FlowIntermittent.y == 1 & lfc <= 0) ~ "Less",
                        (Effect == "Stochastic.c" & diff_FlowStochastic.x == 1 & passed_ss_FlowStochastic.x == 1 & lfc <= 0) ~ "Less",
                        (Effect == "Stochastic.w" & diff_FlowStochastic.y == 1 & passed_ss_FlowStochastic.y == 1 & lfc <= 0) ~ "Less",
                        (Effect == "Constant.c" & diff_FlowConstant.x == 1 & passed_ss_FlowConstant.x == 1 & lfc <= 0) ~ "Less",
                        (Effect == "Constant.w" & diff_FlowConstant.y == 1 & passed_ss_FlowConstant.y == 1 & lfc <= 0) ~ "Less",
                        TRUE ~ "No"),
         Flow = case_when((Effect == "Intermittent.c" | Effect == "Intermittent.w") ~ "Intermittent",
                          (Effect == "Stochastic.c" | Effect == "Stochastic.w") ~ "Stochastic",
                          (Effect == "Constant.c" | Effect == "Constant.w") ~ "Constant",
                          TRUE ~ "NA"),
         Temp = paste0("", str_sub(Effect, start = -2)),
         Temperature = case_when((Temp == ".c") ~ "control",
                                 (Temp == ".w") ~ "warm",
                                 TRUE ~ "NA"),
         Treatment = case_when((Flow == "Intermittent" & Temperature == "control") ~ "Ic",
                               (Flow == "Intermittent" & Temperature == "warm") ~ "Iw",
                               (Flow == "Stochastic" & Temperature == "control") ~ "Sc",
                               (Flow == "Stochastic" & Temperature == "warm") ~ "Sw",
                               (Flow == "Constant" & Temperature == "control") ~ "Cc",
                               (Flow == "Constant" & Temperature == "warm") ~ "Cw",
                               TRUE ~ "NA")) %>%
  dplyr::rename("Class" = "taxon") %>%
  select(Class, Effect, lfc, p_value, DA, Flow, Temperature, Treatment)
df_volcano_18S_flow$Treatment <- factor(df_volcano_18S_flow$Treatment, levels = c("Ic", "Iw", "Sc", "Sw", "Cc", "Cw"))
df_volcano_18S_flow$DA <- factor(df_volcano_18S_flow$DA, levels = c("More", "Less", "No"))

volcano_DA_eukaryote_flow <- df_volcano_18S_flow %>%
  ggplot(aes(x=lfc, y=-log10(p_value), col=Treatment, alpha= DA)) + 
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype="dashed", color = "firebrick", linewidth = 0.8, alpha = 0.5) +
  geom_label_repel(aes(label = Class), size = 5, show.legend = FALSE, box.padding = 0.3, point.padding = 0.3, data = df_volcano_18S_flow[df_volcano_18S_flow$DA %in% c("More", "Less"),]) +
  theme_classic() +
  scale_color_manual(values=c("#b03a2e","#ec7063","#1e8449","#52be80", "#2874a6", "#5dade2"),
                     labels= c(bquote("N"[flow]*"C"[temp]*"  vs "*"I"[flow]*"C"[temp]), 
                               bquote("N"[flow]*"W"[temp]*" vs "*"I"[flow]*"W"[temp]), 
                               bquote("N"[flow]*"C"[temp]*"  vs "*"S"[flow]*"C"[temp]), 
                               bquote("N"[flow]*"W"[temp]*" vs "*"S"[flow]*"W"[temp]), 
                               bquote("N"[flow]*"C"[temp]*"  vs "*"C"[flow]*"C"[temp]), 
                               bquote("N"[flow]*"W"[temp]*" vs "*"C"[flow]*"W"[temp]))) +
  scale_alpha_manual(values=c(1, 1, 0.1), guide = "none") +
  labs(x = "log Fold-Change", y = "-log p-value", color= "Flow effect in") +
  theme(legend.title = element_text(size = 14, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text =  element_text(size = 14),
        legend.position = "right",
        legend.justification = "left") +
  xlim(-3,3) +
  ylim(0,5)
volcano_DA_eukaryote_flow

df_volcano_18S_temperature <- left_join(res_dunn_18S_tempinnatural_class, res_dunn_18S_tempinintermittent_class, by = "taxon") %>% 
  left_join(res_dunn_18S_tempinstochastic_class,  by = "taxon") %>% 
  left_join(res_dunn_18S_tempinconstant_class,  by = "taxon") %>% 
  mutate(Natural = round(lfc_Temperaturewarm.x, 2),
         Intermittent = round(lfc_Temperaturewarm.y, 2),
         Stochastic = round(lfc_Temperaturewarm.x.x, 2),
         Constant = round(lfc_Temperaturewarm.y.y, 2)) %>% 
  pivot_longer(cols = "Natural":"Constant", names_to = "Effect", values_to = "lfc") %>% 
  mutate(p_value = case_when(Effect == "Natural" ~ q_Temperaturewarm.x,
                             Effect == "Intermittent" ~ q_Temperaturewarm.y,
                             Effect == "Stochastic" ~ q_Temperaturewarm.x.x,
                             Effect == "Constant" ~ q_Temperaturewarm.y.y,
                             TRUE ~ 0),
         DA = case_when((Effect == "Natural" & diff_Temperaturewarm.x == 1 & passed_ss_Temperaturewarm.x == 1 & lfc >= 0) ~ "More",
                        (Effect == "Intermittent" & diff_Temperaturewarm.y == 1 & passed_ss_Temperaturewarm.y == 1 & lfc >= 0) ~ "More",
                        (Effect == "Stochastic" & diff_Temperaturewarm.x.x == 1 & passed_ss_Temperaturewarm.x.x == 1 & lfc >= 0) ~ "More",
                        (Effect == "Constant" & diff_Temperaturewarm.y.y == 1 & passed_ss_Temperaturewarm.y.y == 1 & lfc >= 0) ~ "More",
                        (Effect == "Natural" & diff_Temperaturewarm.x == 1 & passed_ss_Temperaturewarm.x == 1 & lfc <= 0) ~ "Less",
                        (Effect == "Intermittent" & diff_Temperaturewarm.y == 1 & passed_ss_Temperaturewarm.y == 1 & lfc <= 0) ~ "Less",
                        (Effect == "Stochastic" & diff_Temperaturewarm.x.x == 1 & passed_ss_Temperaturewarm.x.x == 1 & lfc <= 0) ~ "Less",
                        (Effect == "Constant" & diff_Temperaturewarm.y.y == 1 & passed_ss_Temperaturewarm.y.y == 1 & lfc <= 0) ~ "Less",
                        TRUE ~ "No")) %>% 
  dplyr::rename("Class" = "taxon") %>% 
  select(Class, Effect, lfc, p_value, DA)
df_volcano_18S_temperature$Effect <- factor(df_volcano_18S_temperature$Effect, levels = c("Natural", "Intermittent", "Stochastic", "Constant"))
df_volcano_18S_temperature$DA <- factor(df_volcano_18S_temperature$DA, levels = c("More", "Less", "No"))

#Count how many Genera and families DA in flow
taxonomy_18S_class <- taxonomy_18S %>% select(-Order, -Family, -Genus) %>% distinct()

test <- df_volcano_18S_flow %>%
  left_join(taxonomy_18S_class) %>% 
  filter(DA != "No") %>%
  summarise(DA_classes = n_distinct(Class), DA_phyla = n_distinct(Phylum))

df_volcano_18S_flow %>%
  left_join(taxonomy_18S) %>% 
  filter(DA != "No") %>%
  group_by(Flow) %>% 
  summarise(DA_genera = n_distinct(Class), DA_phyla = n_distinct(Phylum))

df_volcano_18S_flow %>%
  left_join(taxonomy_18S) %>% 
  filter(DA != "No") %>%
  group_by(Treatment) %>% 
  summarise(DA_genera = n_distinct(Class), DA_phyla = n_distinct(Phylum))

df_volcano_18S_temperature %>%
  left_join(taxonomy_18S_class) %>% 
  filter(DA != "No") %>%
  group_by(Treatment) %>% 
  summarise(DA_genera = n_distinct(Class), DA_phyla = n_distinct(Phylum))

volcano_DA_eukaryote_temperature <- df_volcano_18S_temperature %>%
  ggplot(aes(x=lfc, y=-log10(p_value), colour=Effect, alpha= DA)) + 
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype="dashed", color = "firebrick", linewidth = 0.8, alpha = 0.5) +
  geom_label_repel(aes(label = Class), size = 5, show.legend = FALSE, box.padding = 0.3, point.padding = 0.3, data = df_volcano_18S_temperature[df_volcano_18S_temperature$DA %in% c("More", "Less"),]) +
  theme_classic() +
  scale_color_manual(values=c("#85929e", "#ec7063", "#52be80","#5dade2"),
                     labels= c(bquote("N"[flow]*"C"[temp]*" vs "*"N"[flow]*"W"[temp]), 
                               bquote("I"[flow]*"C"[temp]*"  vs  "*"I"[flow]*"W"[temp]), 
                               bquote("S"[flow]*"C"[temp]*" vs "*"S"[flow]*"W"[temp]), 
                               bquote("C"[flow]*"C"[temp]*" vs "*"C"[flow]*"W"[temp]))) +
  scale_alpha_manual(values=c(1, 1, 0.1), guide = "none") +
  labs(x = "log Fold-Change", y = "-log p-value", color= "Temperature comparison") +
  xlim(-3,3) +
  ylim(0,5) +
  theme(legend.title = element_text(size = 14, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text =  element_text(size = 12),
        legend.position = "right",
        legend.justification = "left")
volcano_DA_eukaryote_temperature

Figure5 <- ggarrange(volcano_DA_eukaryote_flow, volcano_DA_eukaryote_temperature, ncol = 1, nrow = 2, labels = c("A", "B"), font.label = list(size = 18), align = "v", heights = c(1.3, 1))
Figure5
ggsave("Figure5.pdf", Figure5, width = 35,
       height = 45,
       units = "cm")


#Figure S15
drought_tolerant <- c("Deinococcus", "Nocardioides", "Massilia", "Peredibacter",
                      "Hymenobacter", "Cavicella", "Thermomonas", "Sphingomonas",
                      "Parablastomonas")

Genus_16S_drought <- Sentinel16S_2022_bacteria_filt %>% 
  ps_filter(Flow == "Intermittent" | Flow == "Natural") %>% 
  tax_glom(taxrank = "Genus") %>% 
  transform_sample_counts(function(x) x / sum(x))

Genus_16S_drought2 <- subset_taxa(Genus_16S_drought, Genus %in% drought_tolerant) %>% ps_melt()
Genus_16S_drought2$Time <- as.factor(Genus_16S_drought2$Time) 
Genus_16S_drought3 <- Genus_16S_drought2 %>% 
  select(Sample, sample_Sample, Time, Succession, Abundance, Flow, Genus) %>% 
  group_by(Flow, Time, Genus) %>%
  mutate(Ab_mean = mean(Abundance)+0.00001) %>% 
  ungroup() %>% 
  group_by(Flow, Genus) %>% 
  mutate(Ab_norm = Ab_mean / Ab_mean[Time == "1"],
         Ab_norm = replace_na(Ab_norm, 0),
         Ab_log = log(Ab_norm))

Plottestdrought <- ggplot(Genus_16S_drought3, aes(x = Succession, y = Ab_norm)) +
  geom_line(aes(color = Genus),linewidth = 1.25, se = FALSE, span = 0.5, lineend = "round") +
  scale_color_viridis(discrete=TRUE) +
  theme_classic() +
  labs(x = "Days of growth", y = "Relative fold increase", color = "Differentially more\nabundant genera") +
  facet_grid(cols = vars(factor(Flow, levels = c("Natural", "Intermittent")))) +
  geom_vline(xintercept = c(24, 38, 87), linetype="dashed", 
             color = "#b03a2e", size=1) +
  theme(legend.title = element_text(size = 16, face="bold"), 
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 14), 
        axis.text.y =  element_text(size = 14, angle=90, hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(size=16, face="bold"),
        legend.position = "right")
Plottestdrought

ggsave("FigureS15.pdf", Plottestdrought, width = 30,
       height = 15,
       units = "cm")

### Community assembly ------------------------------------------------------

BMNTD <- read.csv("BMNTD.csv",sep=",")

BMNTD_long <- BMNTD %>% 
  select(-Flume, -Temperature, -Flow) %>% 
  melt(variable.name = "Process", value.name = "Contribution") %>% 
  mutate(Contribution = Contribution*100,
         Process = factor(Process, levels=c("DL", "DR", "HD", "HeS", "HoS")),
         Treatment = factor(Treatment, levels=c("Nc","Nw", "Ic", "Iw", "Sc", "Sw", "Cc", "Cw")))

# Calculate average and sd
mean(BMNTD$DR,na.rm = TRUE)
sd(BMNTD$DR,na.rm = TRUE)
mean(BMNTD$HoS,na.rm = TRUE)
sd(BMNTD$HoS,na.rm = TRUE)
mean(BMNTD$DL,na.rm = TRUE)
sd(BMNTD$DL,na.rm = TRUE)
mean(BMNTD$HD,na.rm = TRUE)
sd(BMNTD$HD,na.rm = TRUE)
mean(BMNTD$HeS,na.rm = TRUE)
sd(BMNTD$HeS,na.rm = TRUE)

plot_BMNTD_overall <- BMNTD_long %>% ggplot(aes(x = Process, y = Contribution)) +
  geom_boxplot(fill = "lightgray") +
  ylab("Contribution (%)") +
  xlab("Assembly processes") +
  labs(color = "Flow regimes") +
  scale_x_discrete(labels=c("DL" = "Dispersal\nlimitation", "DR" = "Drift", "HD" = "Homogenizing\ndispersal","HeS" = "heterogeneous\nselection", "HoS" = "Homogeneous\nselection")) +
  theme_classic() +
  theme(legend.title = element_text(size = 16, face="bold"),
        legend.text = element_text(size = 16), 
        axis.title.y = element_text(size=16, face="bold"), 
        axis.title.x = element_text(size=16, face="bold"),
        axis.text = element_text(size = 14),
        strip.background = element_blank(),
        axis.line=element_line(),
        legend.position = "right",
        legend.justification = "left")+
  stat_summary(fun.y=median, geom="point", shape=5, size=3, color="black", fill="black") 
plot_BMNTD_overall

#Differences between treatments?
library(multcompView)

DR <- BMNTD_long %>% filter(Process == "DR")
DR.aov <- aov(Contribution ~ Treatment, data = DR)
summary(DR.aov) #There is a significant difference
shapiro.test(residuals(DR.aov)) #Normality can be assume
tukey_DR <- TukeyHSD(DR.aov)
letters_DR <- multcompLetters4(DR.aov, TukeyHSD(DR.aov), threshold = 0.05)
letters_DR_df <- data.frame(Treatment = names(letters_DR$Treatment$Letters), Signif = letters_DR$Treatment$Letters, row.names = NULL) %>% mutate(Process = "DR")

DL <- BMNTD_long %>% filter(Process == "DL")
DL.aov <- aov(Contribution ~ Treatment, data = DL)
summary(DL.aov) #There is a significant difference
shapiro.test(residuals(DL.aov)) #Normality can be assume
tukey_DL <- TukeyHSD(DL.aov)
letters_DL <- multcompLetters4(DL.aov, TukeyHSD(DL.aov), threshold = 0.05, reversed = TRUE)
letters_DL_df <- data.frame(Treatment = names(letters_DL$Treatment$Letters), Signif = letters_HD$Treatment$Letters, row.names = NULL) %>% mutate(Process = "DL")

HD <- BMNTD_long %>% filter(Process == "HD")
HD.aov <- aov(Contribution ~ Treatment, data = HD)
summary(HD.aov) #There is a significant difference
shapiro.test(residuals(HD.aov)) #Normality can be assume
tukey_HD <- TukeyHSD(HD.aov)
letters_HD <- multcompLetters4(HD.aov, TukeyHSD(HD.aov), threshold = 0.05, reversed = TRUE)
letters_HD_df <- data.frame(Treatment = names(letters_HD$Treatment$Letters), Signif = letters_HD$Treatment$Letters, row.names = NULL) %>% mutate(Process = "HD")

HoS <- BMNTD_long %>% filter(Process == "HoS")
HoS.aov <- aov(Contribution ~ Treatment, data = HoS)
summary(HoS.aov) #There is a significant difference
shapiro.test(residuals(HoS.aov)) #Normality can be assume
tukey_HoS <- TukeyHSD(HoS.aov)
letters_HoS <- multcompLetters4(HoS.aov, TukeyHSD(HoS.aov), threshold = 0.05, reversed = TRUE)
letters_HoS_df <- data.frame(Treatment = names(letters_HoS$Treatment$Letters), Signif = letters_HoS$Treatment$Letters, row.names = NULL) %>% mutate(Process = "HoS")

Treatment_sign <- rbind(letters_DL_df, letters_DR_df, letters_HD_df, letters_HoS_df)

BMNTD_long <- BMNTD_long %>%
  left_join(Treatment_sign %>% select(Treatment, Process, Signif), 
            by = c("Treatment", "Process"))

summary_sign <- BMNTD_long %>%
  filter(Process != "HeS") %>%
  group_by(Treatment, Process) %>%
  summarize(Max_Contribution = max(Contribution, na.rm = TRUE),
            Signif = (Signif),
            .groups = 'drop') %>% 
  distinct() %>% 
  mutate(Treatment = factor(Treatment, levels=c("Nc","Nw", "Ic", "Iw", "Sc", "Sw", "Cc", "Cw")))

Process_labels = as_labeller(c("DL" = "Dispersal limitation",
                               "DR" = "Drift", 
                               "HD" = "Homogenizing dispersal",
                               "HoS" = "Homogeneous selection"))

plot_BMNTD_treatment <- BMNTD_long %>% 
  filter(Process != "HeS") %>% 
  ggplot(aes(x = Treatment, y = Contribution, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#283747","#85929e","#b03a2e","#ec7063","#1e8449","#52be80","#2874a6","#5dade2"),
                    labels= c(bquote("N"[flow]*"C"[temp]), bquote("N"[flow]*"W"[temp]), bquote("I"[flow]*"C"[temp]), bquote("I"[flow]*"W"[temp]), bquote("S"[flow]*"C"[temp]), bquote("S"[flow]*"W"[temp]), bquote("C"[flow]*"C"[temp]), bquote("C"[flow]*"W"[temp]))) +
  ylab("Contribution (%)") +
  labs(color = "Flow regimes") +
  theme_classic() +
  theme(legend.title = element_text(size = 16, face="bold"),
        legend.text = element_text(size = 16), 
        axis.title.y = element_text(size=16, face="bold"), 
        axis.title.x = element_text(size=16, face="bold", vjust = -5),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        axis.line=element_line(),
        legend.position = "right",
        legend.justification = "left",
        strip.text = element_text(size = 14, face = "bold"),
        plot.margin = margin(0, 0, 1, 0, "cm")) +
  facet_wrap(~Process, scales = "free", labeller = Process_labels) +
  stat_summary(fun = median, geom = "point", shape = 5, size = 3, color = "black") +
  geom_text(data = summary_sign, aes(x = Treatment, y = Max_Contribution, label = Signif), 
            vjust = -0.2, size = 5, color = "black")
plot_BMNTD_treatment

meanDR_Nc <- BMNTD_long %>% filter(Treatment == "Nc" & Process == "DR")
mean(meanDR_Nc$Contribution)
sd(meanDR_Nc$Contribution)
meanDR_others <- BMNTD_long %>% filter(Treatment != "Nc" & Process == "DR")
mean(meanDR_others$Contribution)
sd(meanDR_others$Contribution)

Figure6 <- ggarrange(plot_BMNTD_overall, plot_BMNTD_treatment, ncol = 2, nrow = 1, labels = c("A", "B"), font.label = list(size = 18), align = "h", common.legend = F, legend ="right")
Figure6 <- plot_grid(plot_BMNTD_overall, plot_BMNTD_treatment, ncol = 2, nrow = 1, align = "h", axis = "bt", labels = c("A", "B"), font.label = list(size = 18))
Figure6

ggsave("Figure6.pdf", Figure6, width = 45,
       height = 25,
       units = "cm")


### Functional metabolic diversity ------------------------------------------

##Correlation biolog vs treatment
Biolog <- sample_data(Sentinel16S_2022_bacteria_filt) %>% data.frame() %>% drop_na(AWCD, R, H, E) %>%
  select(Temperature, Flow, Treatment, Time, Succession, AWCD, R, H, E)
Biolog$FlowTime <- str_c(Biolog$Flow,"_", Biolog$Time)
Biolog$TAD <- Biolog$TAD %>% replace(is.na(.), 0)

#R
shapiro.test(Biolog$R) #Not normal
kruskal.test(R ~ Succession, data = Biolog) #Yes
kruskal.test(R ~ Flow, data = Biolog) #No
kruskal.test(R ~ Temperature, data = Biolog) #Yes

range(Biolog$R, na.rm = TRUE) #20 to 31
mean(Biolog$R) #26.3535
sd(Biolog$R)

R.warm <- Biolog %>% 
  filter(Temperature == "warm")
a <- mean(R.warm$R)  

R.control <- Biolog %>% 
  filter(Temperature == "control")
b <- mean(R.control$R)  
c <- (a - b)

#H
shapiro.test(Biolog$H) #Not normal
kruskal.test(H ~ Succession, data = Biolog) #Yes
kruskal.test(H ~ Flow, data = Biolog) #No
kruskal.test(H ~ Temperature, data = Biolog) #Yes

range(Biolog$H, na.rm = TRUE) #2.969 to 3.306
mean(Biolog$H) #3.176
sd(Biolog$H)

e <- mean(R.warm$H)  
f <- mean(R.control$H)  
g <- (e - f)

#E
shapiro.test(Biolog$E) #Not normal
kruskal.test(E ~ Succession, data = Biolog) #Yes
kruskal.test(E ~ Flow, data = Biolog) #No
kruskal.test(E ~ Temperature, data = Biolog) #Yes

range(Biolog$E, na.rm = TRUE) #0.106 to 0.152
mean(Biolog$E) #0.121
sd(Biolog$E)

h <- mean(R.warm$E)  
i <- mean(R.control$E)  
j <- (h - i)

test <- Biolog %>%
  ggplot(aes(x=Succession,y=R,group=Succession, colour = Temperature))+
           geom_boxplot()



Figure7a <- Biolog %>% mutate(Treatment = factor(Treatment, levels=c("Nc","Nw", "Ic", "Iw", "Sc", "Sw", "Cc", "Cw"))) %>% 
  ggplot(aes(x=Temperature,y=R,fill=Treatment)) + 
  scale_fill_manual(values = c("#283747","#85929e","#b03a2e","#ec7063","#1e8449","#52be80","#2874a6","#5dade2"),
                    labels= c(bquote("N"[flow]*"C"[temp]), bquote("N"[flow]*"W"[temp]), bquote("I"[flow]*"C"[temp]), bquote("I"[flow]*"W"[temp]), bquote("S"[flow]*"C"[temp]), bquote("S"[flow]*"W"[temp]), bquote("C"[flow]*"C"[temp]), bquote("C"[flow]*"W"[temp]))) +
  scale_y_continuous(breaks=seq(20,30,2)) +   
  theme_classic() +
  geom_boxplot() +
  ylab("Substrate Richness") + 
  theme(strip.text = element_text(size = 16),
        legend.title = element_text(size = 14, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16, face="bold"),
        axis.text.y =  element_text(size = 14),
        axis.text.x = element_text(size = 14),
        plot.margin = margin(1, 0, 0, 0, "cm")) +
  stat_summary(fun.y=median, geom="point", shape=5, size=3, color="black", fill="black") +
  geom_signif(comparisons = list(c("control", "warm")), map_signif_level=TRUE, textsize=7, col = 1)
Figure7a
  
Figure7b <- Biolog %>% mutate(Treatment = factor(Treatment, levels=c("Nc","Nw", "Ic", "Iw", "Sc", "Sw", "Cc", "Cw"))) %>% 
  ggplot(aes(x=Temperature,y=H,fill=Treatment)) + 
  scale_fill_manual(values = c("#283747","#85929e","#b03a2e","#ec7063","#1e8449","#52be80","#2874a6","#5dade2"),
                    labels= c(bquote("N"[flow]*"C"[temp]), bquote("N"[flow]*"W"[temp]), bquote("I"[flow]*"C"[temp]), bquote("I"[flow]*"W"[temp]), bquote("S"[flow]*"C"[temp]), bquote("S"[flow]*"W"[temp]), bquote("C"[flow]*"C"[temp]), bquote("C"[flow]*"W"[temp]))) +
  theme_classic() +
  geom_boxplot() +
  ylab("Shannon diversity index") + 
  theme(strip.text = element_text(size = 16),
        legend.title = element_text(size = 14, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16, face="bold"),
        axis.text.y =  element_text(size = 14),
        axis.text.x = element_text(size = 14),
        plot.margin = margin(1, 0, 0, 0, "cm")) +
  stat_summary(fun.y=median, geom="point", shape=5, size=3, color="black", fill="black") +
  geom_signif(test="wilcox.test", comparisons = list(c("control", "warm")), map_signif_level=TRUE, textsize=7, col = 1)
Figure7b
  
Figure7 <- ggarrange(Figure7a, Figure7b, ncol = 2, nrow = 1, labels = c("A", "B"), font.label = list(size = 18), common.legend = T, legend = "right")
Figure7
ggsave("Figure7.pdf", Figure7, width = 30,
         height = 20,
         units = "cm")


## Environmental parameters ------------------------------------------------
metadata_env <- metadata16SFull %>%
  column_to_rownames("NAME") %>% 
  filter(Year == 2022) %>% 
  filter(Treatment == "Nw" | Treatment == "Nc" | Treatment == "River") %>% 
  select(-Primers, -Run, -Type, -TAD, -PO4, -NO3, -NH4, -NO2, -ChlA, -BA, -Time, -Sample, -FlumeID, -Flow, -Year, -Velocity) %>% 
  distinct() %>% 
  mutate(Ammonium = Ammonium/1000,
         Bromide = Bromide/1000,
         Calcium = Calcium/1000,
         Chloride = Chloride/1000,
         DOC = DOC/1000,
         Fluoride = Fluoride/1000,
         Lithium = Lithium/1000,
         Magnesium = Magnesium/1000,
         Nitrate = Nitrate/1000,
         Nitrite = Nitrite/1000, 
         Potassium = Potassium/1000, 
         Sodium = Sodium/1000, 
         Strontium = Strontium/1000,
         Sulfate = Sulfate/1000)

parameters <- c("pH", "C25", "O2", "Turbidity", "Fluoride", "Chloride", "Nitrite", 
                "Bromide", "Nitrate", "Sulfate", "Lithium", "Sodium", "Ammonium", "Magnesium", 
                "Potassium", "Calcium", "Strontium", "DOC")

metadata_env_long <- metadata_env %>%
  pivot_longer(cols = all_of(parameters), names_to = "Parameter", values_to = "Value") 

labels_param <- c("Ammonium (ppm)", "Bromide (ppm)", "Conductivity (µS/cm)", "Calcium (ppm)", "Chloride (ppm)", 
                  "DOC (ppm)", "Fluoride (ppm)", "Lithium (ppm)", "Magnesium (ppm)", "Nitrate (ppm)", "Nitrite (ppm)", 
                  "O2 (mg/L)", "pH", "Potassium (ppm)", "Sodium (ppm)", "Strontium (ppm)", "Sulphate (ppm)", "Turbidity (NTU)")
names(labels_param) <- c("Ammonium", "Bromide", "C25", "Calcium", "Chloride", 
                         "DOC", "Fluoride", "Lithium", "Magnesium", "Nitrate", "Nitrite", 
                         "O2", "pH", "Potassium", "Sodium", "Strontium", "Sulfate", "Turbidity")

FigureS2 <- ggplot(metadata_env_long, aes(x = Succession, y = Value, color = Treatment)) +
  geom_point() + 
  geom_line() +
  labs(x = "Succession (days)",
       y = "Value",
       color = "Measurement site") +
  theme_classic() +
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16, face="bold"),
        axis.text.y =  element_text(size = 12),
        axis.text.x = element_text(size = 14),
        legend.position = "bottom") +
  scale_color_manual(values = c("royalblue", "firebrick", "darkblue" ), labels = c("Control temp. header tank", "Warm temp. header tank", "Streamwater")) +
  facet_wrap(~ Parameter, labeller = as_labeller(labels_param), scales = "free_y", ncol = 3)
FigureS2 

ggsave("FigureS2.pdf", FigureS2, width = 25,
       height = 25,
       units = "cm")

##Supplementary Table 1
#Data
Suppl_table_1_values <- metadata_env %>% data_summary(varname="pH", groupnames=c("Treatment"))
#Stats (tested for each parameters)
shapiro.test(metadata_env$Temp) #Not normal if p < 0.05
#Not normal
kruskal.test(Turbidity ~ Treatment, data = metadata_env) # Not different if p > 0.05
#Which one differ
wilcox <- pairwise.wilcox.test(metadata_env$pH, metadata_env$Treatment, p.adjust.method = "BH")

#Normal
annova <- aov(Temp ~ Treatment, data = metadata_env)
#Which one differ
pairwise <- pairwise.t.test(metadata_env$Temp, metadata_env$Treatment, p.adj = "BH")
