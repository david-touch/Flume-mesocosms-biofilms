
#Flume-mesocosms-biofilms

#16S 2022 data preparation

#Load libraries 
library(phyloseq)
library(ggplot2)
library(gridExtra)
library(vegan)
library(grid)
library(tidyverse)
library(dbplyr)
library(ggpubr)
library(agricolae)
library(microViz)
library(ape)
library(stringr)
library(base)
library(gridExtra)
library(grid)
library(picante)
library(ranacapa)

#Import files and format changes
ASV16SFull = read.csv("16SFullMerged_table_tsv.csv",sep=",")
ASV_16SFull <- subset(ASV16SFull, select = -c(taxonomy))
ASV_16SFull <- ASV_16SFull %>% 
  column_to_rownames(var="OTUID")
ASV16S <- ASV_16SFull %>% dplyr::select(1:789)
ASV16S <- subset(ASV16S, select = -c(X63A_63A_NSE724.NSE510)) #Negative control
ASV16S <- ASV16S[(rowSums(ASV16S) > 0),]
ASV16S_water <- ASV_16SFull %>% dplyr::select(790:822)
ASV16S_water <- ASV16S_water[(rowSums(ASV16S_water) > 0),]
tree16SFull <- read.tree("16SFullMergedtree.nwk")

metadata16SFull = read.csv("metadata16SFull.csv", sep=",")
biolog_sampledata = read.csv("biolog_sampledata.csv", sep = ",")
row.names(metadata16SFull) <- NULL
metadata16S <- metadata16SFull[1:789,]
metadata16S <- metadata16S %>% mutate(FlumeID_Time_Year = paste(FlumeID, Time, Year, sep = '_')) %>% left_join(biolog_sampledata) %>% select(-c(FlumeID_Time_Year)) %>% 
  column_to_rownames(var="NAME")
tree16S <- phy_tree(tree16SFull)
metadata16S_water <- metadata16SFull[790:822,] 
row.names(metadata16S_water) <- NULL
metadata16S_water <- subset(metadata16S_water, select = -c(Flow,Temperature, TAD, ChlA, BA, Run)) %>% 
  column_to_rownames(var="NAME")
metadata16S_water$Time <- sprintf("%02d", as.numeric(metadata16S_water$Time))

#Parse taxonomy file
taxonomy16SFull <- subset(ASV16SFull, select = c(OTUID,taxonomy)) %>%
  separate(taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "; ") %>%
  column_to_rownames(var="OTUID")
taxonomy16SFull <- as.matrix(taxonomy16SFull)
taxonomy16SFull <- gsub(".__","",taxonomy16SFull)

##16S Full Phyloseq (without water)
#Forcing data frame
ASV_table_16S <- data.frame(ASV16S)
sample_data_16S <- data.frame(metadata16S)

#Change to phyloseq object
ASV_table_16S  <- otu_table(ASV_table_16S, taxa_are_rows = TRUE)
sample_data_16S <- sample_data(sample_data_16S)
tax_table_16SFull <- tax_table(taxonomy16SFull)

#Merging dataframe
Sentinel16S <- merge_phyloseq(ASV_table_16S, sample_data_16S, tax_table_16SFull, tree16S)
Sentinel16S 

#Rarefaction curve
Sentinel16S_2022 <- Sentinel16S %>% ps_filter(Year == 2022)

factor(sample_data(Sentinel16S_2022)$Treatment, levels = c("Nc", "Nw", "Ic", "Iw", "Sc", "Sw", "Cc", "Cw"))
rarecurve16S <- ggrare(Sentinel16S_2022, step = 1000, color = "Treatment",
                       plot = TRUE, parallel = FALSE, se = F)
rarecurve_16S <- rarecurve16S +
  theme_classic() +
  facet_wrap(~factor(Treatment, c("Nc","Nw","Ic", "Iw", "Sc", "Sw", "Cc", "Cw")), ncol = 2, scales='free') +
  scale_color_manual(values = c("#283747","#85929e","#b03a2e","#ec7063","#1e8449","#52be80", "#2874a6","#5dade2"), breaks=c("Nc","Nw","Ic", "Iw", "Sc", "Sw", "Cc", "Cw"),
                     labels= c(bquote("N"[flow]*"C"[temp]), bquote("N"[flow]*"W"[temp]), bquote("I"[flow]*"C"[temp]), bquote("I"[flow]*"W"[temp]), bquote("S"[flow]*"C"[temp]), bquote("S"[flow]*"W"[temp]), bquote("C"[flow]*"C"[temp]), bquote("C"[flow]*"W"[temp]))) +
  ylab("Species Richness (bacterial ASVs)") +
  theme(legend.title = element_text(size = 14, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 12), 
        axis.text.y =  element_text(size = 12, angle=90, hjust = 0.5),
        strip.text.x = element_blank()) +
  guides(colour = guide_legend(override.aes = list(linewidth= 3, linetype=2, lineend = "round"))) +
  xlim(0, 150000) +
  scale_y_continuous(breaks = c(0, 1000, 2000), limits = c(0, 2500))
rarecurve_16S

FigureS3 <- ggarrange(rarecurve_16S, rarecurve_18S, ncol = 1, nrow = 2, labels = c("A", "B"), font.label = list(size = 18), common.legend = TRUE, legend = "right")
FigureS3
ggsave("FigureS3.pdf", FigureS3, width = 30,
       height = 35,
       units = "cm")


#Validation with microViz
phyloseq_validate(Sentinel16S, remove_undetected = TRUE)
Sentinel16S
Sentinel16S <- tax_fix(Sentinel16S)
Sentinel16S
samdat_tbl(Sentinel16S)

#Filter unassigned Kingdom, chloroplast, mitochondria and Archaea
Sentinel16S.1 <- subset_taxa(Sentinel16S, !is.na(Kingdom) & !Kingdom %in% c("", "Unassigned") & !Genus %in% c("Bacteria", "Bacteria Kingdom"))
Sentinel16S.1
Sentinel16S.c <- subset_taxa(Sentinel16S.1, (Order!="Chloroplast") | is.na(Order))
Sentinel16S.c
Sentinel16S.cm <- subset_taxa(Sentinel16S.c, (Family!="Mitochondria") | is.na(Family))
Sentinel16S.cm
Sentinel16S.cma <- subset_taxa(Sentinel16S.cm, (Kingdom!="Archaea") | is.na(Kingdom))
Sentinel16S.cma 
Sentinel16S_bacteria <- subset_taxa(Sentinel16S.cma, (Kingdom!="Eukaryota") | is.na(Kingdom))
Sentinel16S_bacteria

##Fixing taxonomy
Sentinel16S_bacteria <- Sentinel16S_bacteria %>% tax_fix(
  min_length = 4,
  unknowns = c("Gammaproteobacteria_Incertae_Sedis", "Incertae_Sedis", "Synechococcales_Incertae_Sedis", "Cyanobacteriales_Incertae_Sedis", "Rhizobiales_Incertae_Sedis", "Incertae_Sedis Genus", "Oxyphotobacteria_Incertae_Sedis", "uncultured", "Unknown_Family"),
  sep = " ", anon_unique = TRUE,
  suffix_rank = "classified")

#Create the "Drought" variable
Sentinel16S_bacteria <- Sentinel16S_bacteria %>% 
  ps_mutate(Drought = if_else(Flow == "Intermittent" & Time >= "05", "Yes", "No"))

saveRDS(Sentinel16S_bacteria, file = "Sentinel16S_bacteria.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

#Select 2022 and filter for more than 1 read in at least 2 samples
Sentinel16S_2022_bacteria <- Sentinel16S_bacteria %>% 
  ps_filter(Year == 2022)
Sentinel16S_2022_bacteria
Sentinel16S_2022_bacteria_filt <- Sentinel16S_2022_bacteria %>% 
  filter_taxa(function(x) sum(x > 1) > (0.002*length(x)), TRUE)
Sentinel16S_2022_bacteria_filt

saveRDS(Sentinel16S_2022_bacteria, file = "Sentinel16S_2022_bacteria.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
saveRDS(Sentinel16S_2022_bacteria_filt, file = "Sentinel16S_2022_bacteria_filt.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

##Merging by Samples

#Create the "Drought" variable
Sentinel16S_2022_bacteria_filt <- Sentinel16S_2022_bacteria_filt %>% 
  ps_mutate(Drought = if_else(Flow == "Intermittent" & Time >= "05", "Yes", "No"))

sample_data(Sentinel16S_2022_bacteria_filt)$Sample <- as.factor(sample_data(Sentinel16S_2022_bacteria_filt)$Sample)
sample_data(Sentinel16S_2022_bacteria_filt)$Treatment <- as.factor(sample_data(Sentinel16S_2022_bacteria_filt)$Treatment)
sample_data(Sentinel16S_2022_bacteria_filt)$Flow <- as.factor(sample_data(Sentinel16S_2022_bacteria_filt)$Flow)
sample_data(Sentinel16S_2022_bacteria_filt)$Temperature <- as.factor(sample_data(Sentinel16S_2022_bacteria_filt)$Temperature)
sample_data(Sentinel16S_2022_bacteria_filt)$Drought <- as.factor(sample_data(Sentinel16S_2022_bacteria_filt)$Drought)

Sentinel16S_2022_bacteria_filt_merged <- merge_samples(Sentinel16S_2022_bacteria_filt, "Sample", fun = "mean") %>% 
  ps_select(-c("Sample", "Run", "FlumeID")) %>% 
  ps_mutate(Flow = str_replace(Flow, "1", "Constant"), Flow = str_replace(Flow, "2", "Intermittent"), 
            Flow = str_replace(Flow, "3", "Natural"), Flow = str_replace(Flow, "4", "Stochastic"),
            Temperature = str_replace(Temperature, "1", "control"), Temperature = str_replace(Temperature, "2", "warm"),
            Treatment = str_replace(Treatment, "1", "Cc"), Treatment = str_replace(Treatment, "2", "Cw"), 
            Treatment = str_replace(Treatment, "3", "Ic"), Treatment = str_replace(Treatment, "4", "Iw"), 
            Treatment = str_replace(Treatment, "5", "Nc"), Treatment = str_replace(Treatment, "6", "Nw"), 
            Treatment = str_replace(Treatment, "7", "Sc"), Treatment = str_replace(Treatment, "8", "Sw"),
            Drought = str_replace(Drought, "1", "No"), Drought = str_replace(Drought, "2", "Yes")) %>% 
  ps_mutate(Primers = "16S")
otu_table(Sentinel16S_2022_bacteria_filt_merged) <- t(otu_table(Sentinel16S_2022_bacteria_filt_merged))

saveRDS(Sentinel16S_2022_bacteria_filt_merged, file = "Sentinel16S_2022_bacteria_filt_merged.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

saveRDS(Sentinel16S_2022_bacteria_filt_merged, file = "Sentinel16S_2022_bacteria_filt_merged.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)


## Phyloseq 16S water

#Create dataframe
ASV_table_16S_water <- data.frame(ASV16S_water)
sample_data_16S_water <- data.frame(metadata16S_water)

#Change to phyloseq object
ASV_table_16S_water  <- otu_table(ASV_table_16S_water, taxa_are_rows = TRUE)
sample_data_16S_water <- sample_data(sample_data_16S_water)

#Merging dataframe
Sentinel16S_water <- merge_phyloseq(ASV_table_16S_water, sample_data_16S_water, tax_table_16SFull, tree16S)
Sentinel16S_water 

#Validation with microViz
phyloseq_validate(Sentinel16S_water, remove_undetected = TRUE)
Sentinel16S_water
Sentinel16S_water <- tax_fix(Sentinel16S_water)
Sentinel16S_water
samdat_tbl(Sentinel16S_water)

#Filter unassigned Kingdom, chloroplast, mitochondria and Archaea
Sentinel16S_water.1 <- subset_taxa(Sentinel16S_water, !is.na(Kingdom) & !Kingdom %in% c("", "Unassigned"))
Sentinel16S_water.1
Sentinel16S_water.c <- subset_taxa(Sentinel16S_water.1, (Order!="Chloroplast") | is.na(Order))
Sentinel16S_water.c
Sentinel16S_water.cm <- subset_taxa(Sentinel16S_water.c, (Family!="Mitochondria") | is.na(Family))
Sentinel16S_water.cm
Sentinel16S_water.cma <- subset_taxa(Sentinel16S_water.cm, (Kingdom!="Archaea") | is.na(Kingdom))
Sentinel16S_water.cma 
Sentinel16S_water_bacteria <- subset_taxa(Sentinel16S_water.cma, (Kingdom!="Eukaryota") | is.na(Kingdom))
Sentinel16S_water_bacteria <- filter_taxa(Sentinel16S_water_bacteria, function(x) sum(x) > 1, TRUE)
Sentinel16S_water_bacteria

Sentinel16S_water_bacteria <- Sentinel16S_water_bacteria %>% tax_fix(
  min_length = 4,
  unknowns = c("Gammaproteobacteria_Incertae_Sedis", "Incertae_Sedis", "Synechococcales_Incertae_Sedis", "Cyanobacteriales_Incertae_Sedis", "Rhizobiales_Incertae_Sedis", "Incertae_Sedis Genus", "Oxyphotobacteria_Incertae_Sedis", "uncultured", "Unknown_Family"),
  sep = " ", anon_unique = TRUE,
  suffix_rank = "classified")

saveRDS(Sentinel16S_water_bacteria, file = "Sentinel16S_water_bacteria.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

Sentinel16S_2022_water_bacteria <- Sentinel16S_water_bacteria %>% ps_filter(Year == 2022)
Sentinel16S_2022_water_bacteria

## Venn Diagram of Water samples compared to biofilm samples
library(MicEco)

phyloseq_2022_water_venn_tax_16S <- tax_table(Sentinel16S_2022_water_bacteria)
phyloseq_2022_water_venn_asv_16S <- otu_table(Sentinel16S_2022_water_bacteria)
phyloseq_2022_water_venn_meta_16S <- sample_data(Sentinel16S_2022_water_bacteria)

phyloseq_2022_biofilm_venn_tax_16S <- tax_table(Sentinel16S_2022_bacteria_filt)
phyloseq_2022_biofilm_venn_asv_16S <- otu_table(Sentinel16S_2022_bacteria_filt)
phyloseq_2022_biofilm_venn_meta_16S <- sample_data(Sentinel16S_2022_bacteria_filt)

taxvenn_16S <- merge_phyloseq(phyloseq_2022_water_venn_tax_16S, phyloseq_2022_biofilm_venn_tax_16S)
asvvenn_16S <- merge_phyloseq(phyloseq_2022_water_venn_asv_16S, phyloseq_2022_biofilm_venn_asv_16S)
metavenn_16S <- merge_phyloseq(phyloseq_2022_water_venn_meta_16S, phyloseq_2022_biofilm_venn_meta_16S)

Venn_16S_2022 <- merge_phyloseq(taxvenn_16S, asvvenn_16S, metavenn_16S)

MicEco_16S <- ps_venn(Venn_16S_2022, group = "Type")
MicEco_16S

MicEco_16S_abundance <- ps_venn(Venn_16S_2022, group = "Type", weight = TRUE)
MicEco_16S_abundance

abundance.round <- MicEco_16S_abundance$data$original.values
abundance.round <- round(as.numeric(abundance.round), 4)
MicEco_16S_abundance$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$label <- 10.05
MicEco_16S_abundance$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$label <- 10.50
MicEco_16S_abundance$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$label <- 79.45

FigureS6 <- ggarrange(MicEco_16S, MicEco_16S_abundance, ncol = 1, nrow = 2, labels = c("A", "B"), font.label = list(size = 18))
FigureS6
ggsave("FigureS6.pdf", FigureS6, width = 20,
       height = 20,
       units = "cm")

##Count number of reads
Sentinel16S_2022_bacteria_filt
reads_16S <- sum(unlist(otu_table(Sentinel16S_2022_bacteria_filt)), na.rm = TRUE)
Sentinel16S_2022_water_bacteria
reads_16S_water <- sum(unlist(otu_table(Sentinel16S_2022_water_bacteria)), na.rm = TRUE)

#Plot Graph-Bar for water samples
Graphbar_16S_phylum_water <- Sentinel16S_2022_water_bacteria %>%
  merge_samples(group = "Succession") %>% 
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 14,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    bar_width = 0.9) + 
  labs(x = "Days of growth") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(size=16, face="bold"),
        axis.title.x = element_text(size=16, face="bold"),
        strip.text = element_text(size = 16), 
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  guides(fill = guide_legend(ncol = 3))
Graphbar_16S_phylum_water

Graphbar_16S_family_water <- Sentinel16S_2022_water_bacteria %>%
  merge_samples(group = "Succession") %>% 
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Family", n_taxa = 14,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    bar_width = 0.9) + 
  labs(x = "Days of growth") +
  xlim("5", "10", "14", "19", "24", "26","28", "33", "38", "40", "42", "47", "54", "61", "67", "75", "82", "87", "89", "91", "96", "103") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16, face = "bold", angle = 90, hjust = 0.5, vjust = 0.5),
        axis.title.y = element_text(size=16, face="bold"),
        axis.title.x = element_text(size=16, face="bold"),
        strip.text = element_text(size = 16), 
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  guides(fill = guide_legend(ncol = 3))
Graphbar_16S_family_water

#Figure S6
left_titleBact <- text_grob("Bacteria", size = 18, rot = 90, hjust = 0.5, vjust = 0.5, face ="bold")
left_titlePhoto <- text_grob("Phototrophic eukaryotes", size = 18, rot = 90, hjust = 0.5, vjust = 0.5, face ="bold")
left_titleOthers <- text_grob("Non-phototrophic eukaryotes", size = 18, rot = 90, hjust = 0.5, vjust = 0.5, face ="bold")

FigureS7A <- ggarrange(left_titleBact, Graphbar_16S_phylum_water, Graphbar_16S_family_water, ncol = 3, nrow = 1, labels = c(" ", "A", "B"), font.label = list(size = 18), widths = c(0.05, 1, 1))
FigureS7B <- ggarrange(left_titlePhoto, Graphbar_18S_photo_phylum_water, Graphbar_18S_photo_class_water, ncol = 3, nrow = 1, labels = c(" ", "C", "D"), font.label = list(size = 18), widths = c(0.05, 1, 1))
FigureS7C <- ggarrange(left_titleOthers, Graphbar_18S_others_phylum_water, Graphbar_18S_others_class_water, ncol = 3, nrow = 1, labels = c(" ", "E", "F"), font.label = list(size = 18), widths = c(0.05, 1, 1))
FigureS7 <- ggarrange(FigureS7A, FigureS7B, FigureS7C, ncol = 1, nrow = 3)
FigureS7  
ggsave("FigureS7.pdf", FigureS7, width = 50,
       height = 70,
       units = "cm")


### NMDS biofilm and water --------------------------------------------------

##Hellinger transformation
Sentinel16S_2022_bacteria_filt_hell <- microbiome::transform(Sentinel16S_2022_bacteria_filt, 'hell')
Sentinel16S_2022_bacteria_filt_merged_hell <- microbiome::transform(Sentinel16S_2022_bacteria_filt_merged, 'hell')
Sentinel16S_2022_water_bacteria_hell <- microbiome::transform(Sentinel16S_2022_water_bacteria, 'hell')


phyloseq_2022_water_16S_tax <- tax_table(Sentinel16S_2022_water_bacteria_hell)
phyloseq_2022_water_16S_asv <- otu_table(Sentinel16S_2022_water_bacteria_hell)
phyloseq_2022_water_16S_meta <- sample_data(Sentinel16S_2022_water_bacteria_hell)

phyloseq_2022_biofilm_16S_tax <- tax_table(Sentinel16S_2022_bacteria_filt_merged_hell)
phyloseq_2022_biofilm_16S_asv <- otu_table(Sentinel16S_2022_bacteria_filt_merged_hell)
phyloseq_2022_biofilm_16S_meta <- sample_data(Sentinel16S_2022_bacteria_filt_merged_hell) %>% data.frame() %>% mutate(Type = as.character(Type), Type = if_else(is.na(Type), "Biofilm", Type)) %>% sample_data()

tax_16S_all <- merge_phyloseq(phyloseq_2022_water_16S_tax, phyloseq_2022_biofilm_16S_tax)
asv_16S_all <- merge_phyloseq(phyloseq_2022_water_16S_asv, phyloseq_2022_biofilm_16S_asv)
meta_16S_all <- merge_phyloseq(phyloseq_2022_water_16S_meta, phyloseq_2022_biofilm_16S_meta)

NMDS_all_16S_2022 <- merge_phyloseq(tax_16S_all, asv_16S_all, meta_16S_all)

matrix_dist_16S_2022_all <- as.matrix(t(data.frame(otu_table(NMDS_all_16S_2022))))
metadata_dist_16S_2022_all <- data.frame(sample_data(NMDS_all_16S_2022)) %>% rownames_to_column("Samples")
set.seed(19950930)
dist_16S_2022_bray_all <- vegdist(matrix_dist_16S_2022_all, method ="bray") %>% 
  as.matrix() %>% 
  data.frame() %>% 
  rownames_to_column("Samples")
dist_16S_2022_all <- dist_16S_2022_bray_all %>% 
  dplyr::select(all_of(.[["Samples"]])) %>% 
  as.dist()
nmds_16S_2022_all <- metaMDS(dist_16S_2022_all, trymax=100)
stress_16S_all <- nmds_16S_2022_all$stress
stress_16S_all

scores_nmds_16S_2022_all <- scores(nmds_16S_2022_all) %>% 
  as_tibble(rownames = "Samples") %>% 
  inner_join(., metadata_dist_16S_2022_all, by="Samples")

scores_nmds_16S_2022_all  <- scores_nmds_16S_2022_all[order(scores_nmds_16S_2022_all$Treatment, scores_nmds_16S_2022_all$Time),]

#Differences?
adonis_bc_16S_all <- adonis2(matrix_dist_16S_2022_all ~ Type, data=metadata_dist_16S_2022_all, permutations=999, method="bray")
adonis_bc_16S_all #Type are significant 

NMDS_16S_2022_all <- ggplot(scores_nmds_16S_2022_all, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = Succession, shape = Type), size = 4, alpha = 0.8) +
  scale_color_viridis(discrete = FALSE, option = "viridis", breaks = c(0, 20, 40, 60, 80, 100)) +
  theme_classic() +
  guides(colour = guide_colourbar(barheight = 10)) +
  labs(color = "Days of growth", caption = "Stress 0.087") +
  scale_x_continuous(limits = c(-0.6, 0.5)) +
  theme(legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12),
        strip.background = element_blank(),
        plot.caption = element_text(size = 16, hjust = 1, vjust = 20),
        strip.text = element_text(size = 16, face = "bold"),
        axis.line = element_line())
NMDS_16S_2022_all

FigureS15 <- ggarrange(NMDS_16S_2022_all, NMDS_18S_photo_2022_all, NMDS_18S_others_2022_all, ncol = 1, nrow = 3, labels = c("A", "B", "C"), font.label = list(size = 18), align = "v")
FigureS15

ggsave("Figure15.pdf", FigureS15, width = 30,
       height = 40,
       units = "cm")



