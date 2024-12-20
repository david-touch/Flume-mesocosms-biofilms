
#Flume-mesocosms-biofilms

#18S 2022 data preparation

#Load libraries 
library(phyloseq)
library(ggplot2)
library(gridExtra)
library(vegan)
library(grid)
library(tidyverse)
library(dbplyr)
library(ggpubr)
#library(agricolae)
library(microViz)
library(ape)
library(stringr)
library(base)
library(gridExtra)
library(grid)
library(picante)


#Import files and format changes
ASV18SFull = read.csv("18SFullMerged_table_tsv.csv",sep=",")
ASV_18SFull <- subset(ASV18SFull, select = -c(taxonomy))
ASV_18SFull <- ASV_18SFull %>% 
  column_to_rownames(var="OTUID")
ASV18S <- ASV_18SFull %>% dplyr::select(1:658)
ASV18S <- subset(ASV18S, select = -c(X63A_63A_NSE724.NSE521)) #Negative control
ASV18S <- ASV18S[(rowSums(ASV18S) > 0),]
ASV18S_water <- ASV_18SFull %>% dplyr::select(659:691)
ASV18S_water <- ASV18S_water[(rowSums(ASV18S_water) > 0),]

metadata18SFull = read.csv("metadata18SFull.csv", sep=",")
row.names(metadata18SFull) <- NULL
metadata18S <- metadata18SFull[1:658,] %>% 
  column_to_rownames(var="NAME")
metadata18S$Time <- sprintf("%02d", as.numeric(metadata18S$Time))
metadata18S_water <- metadata18SFull[659:691,] 
row.names(metadata18S_water) <- NULL
metadata18S_water <- subset(metadata18S_water, select = -c(Flow,Temperature, Treatment, FlumeID, TAD, ChlA, BA, Run)) %>% 
  column_to_rownames(var="NAME")
metadata18S_water$Time <- sprintf("%02d", as.numeric(metadata18S_water$Time))


#Parse taxonomy file
taxonomy18SFull <- subset(ASV18SFull, select = c(OTUID,taxonomy)) %>%
  separate(taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "; ") %>%
  column_to_rownames(var="OTUID")
taxonomy18SFull <- as.matrix(taxonomy18SFull)
taxonomy18SFull <- gsub(".__","", taxonomy18SFull)

##18S Full Phyloseq
#Forcing data frame
ASV_table_18S <- data.frame(ASV18S)
sample_data_18S <- data.frame(metadata18S)

#Change to phyloseq object
ASV_table_18S  <- otu_table(ASV_table_18S, taxa_are_rows = TRUE)
sample_data_18S <- sample_data(sample_data_18S)
tax_table_18SFull <- tax_table(taxonomy18SFull)

#Merging dataframe
Sentinel18S <- merge_phyloseq(ASV_table_18S, sample_data_18S, tax_table_18SFull)
Sentinel18S 

#Rarefaction curve
Sentinel18S_2022 <- Sentinel18S %>% ps_filter(Year == 2022)

factor(sample_data(Sentinel18S_2022)$Treatment, levels = c("Nc", "Nw", "Ic", "Iw", "Sc", "Sw", "Cc", "Cw"))
rarecurve18S <- ggrare(Sentinel18S_2022, step = 1000, color = "Treatment",
                       plot = TRUE, parallel = FALSE, se = F)
rarecurve_18S <- rarecurve18S +
  theme_classic() +
  facet_wrap(~factor(Treatment, c("Nc","Nw","Ic", "Iw", "Sc", "Sw", "Cc", "Cw")), ncol = 2, scales='free') +
  scale_color_manual(values = c("#283747","#85929e","#b03a2e","#ec7063","#1e8449","#52be80", "#2874a6","#5dade2"), breaks=c("Nc","Nw","Ic", "Iw", "Sc", "Sw", "Cc", "Cw"),
                     labels= c(bquote("N"[flow]*"C"[temp]), bquote("N"[flow]*"W"[temp]), bquote("I"[flow]*"C"[temp]), bquote("I"[flow]*"W"[temp]), bquote("S"[flow]*"C"[temp]), bquote("S"[flow]*"W"[temp]), bquote("C"[flow]*"C"[temp]), bquote("C"[flow]*"W"[temp]))) +
  ylab("Species Richness (eukaryotic ASVs)") +
  theme(legend.title = element_text(size = 14, face="bold"), 
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 16, face="bold"),
        axis.text.x =  element_text(size = 12), 
        axis.text.y =  element_text(size = 12, angle=90, hjust = 0.5),
        strip.text.x = element_blank()) +
  guides(colour = guide_legend(override.aes = list(linewidth= 3, linetype=2, lineend = "round"))) +
  xlim(0, 150000) +
  scale_y_continuous(breaks = c(0, 500, 800), limits = c(0, 1000))
rarecurve_18S

#Validation with microViz
phyloseq_validate(Sentinel18S, remove_undetected = TRUE)
Sentinel18S
Sentinel18S <- tax_fix(Sentinel18S, min_length = 3)
Sentinel18S
samdat_tbl(Sentinel18S)

#Filter bacteria, Unassigned Kingdom, Arthropoda and not resolved Genus
Sentinel18S.1 <- subset_taxa(Sentinel18S, !is.na(Kingdom) & !Kingdom %in% c("", "Unassigned") & !Genus %in% c("Eukaryota Kingdom", "Eukaryota"))
Sentinel18S.1
Sentinel18S.a <- subset_taxa(Sentinel18S.1, !Phylum %in% c("Arthropoda"))
Sentinel18S.a
Sentinel18S_eukaryote <- subset_taxa(Sentinel18S.a, (Kingdom!="Bacteria") | is.na(Kingdom))
Sentinel18S_eukaryote

Sentinel18S_eukaryote <- Sentinel18S_eukaryote %>% tax_fix(
  min_length = 4,
  unknowns = c("Incertae_Sedis", "uncultured", "LKM11"),
  sep = " ", anon_unique = TRUE,
  suffix_rank = "classified")
tax_table(Sentinel18S_eukaryote) <- tax_table(Sentinel18S_eukaryote) %>% 
  data.frame() %>% 
  mutate(Order = replace(Order, Order == "Bacillariophyceae", "Bacillariophyceae Class"),
         Family = replace(Family, Family == "Bacillariophyceae", "Bacillariophyceae Class"),
         Genus = replace(Genus, Genus == "Bacillariophyceae", "Bacillariophyceae Class")) %>% 
  mutate(Family = replace(Family, Family == "Monhysterida", "Monhysterida Order"),
         Genus = replace(Genus, Genus == "Monhysterida", "Monhysterida Order")) %>% 
  mutate(Genus = replace(Genus, Genus == "Vampyrellidae", "Vampyrellidae Family")) %>%
  mutate(Class = replace(Class, Class == "Labyrinthulomycetes", "Labyrinthulomycetes Phylum"),
         Order = replace(Order, Order == "Labyrinthulomycetes", "Labyrinthulomycetes Phylum"),
         Family = replace(Family, Family == "Labyrinthulomycetes", "Labyrinthulomycetes Phylum"),
         Genus = replace(Genus, Genus == "Labyrinthulomycetes", "Labyrinthulomycetes Phylum")) %>% 
  mutate(Genus = replace(Genus, Genus == "Bacillariophyceae Family", "Bacillariophyceae Class"),
         Species = replace(Genus, Species == "Bacillariophyceae Family", "Bacillariophyceae Class"),
         Family = replace(Family, Family == "Labyrinthulomycetes Order", "Labyrinthulomycetes Phylum")) %>% 
  mutate(Genus = replace(Genus, Genus == "Labyrinthulomycetes Order", "Labyrinthulomycetes Phylum"),
         Species = replace(Genus, Species == "Labyrinthulomycetes Order", "Labyrinthulomycetes Phylum"),
         Phylum = replace(Phylum, Phylum == "Labyrinthulomycetes", "Eukaryota Kingdom")) %>% 
  mutate(Class = replace(Class, Class == "Cercozoa", "Cercozoa Phylum"),
         Order = replace(Order, Order == "Cercozoa", "Cercozoa Phylum"),
         Family = replace(Family, Family == "Cercozoa", "Cercozoa Phylum"),
         Genus = replace(Genus, Genus == "Cercozoa", "Cercozoa Phylum")) %>% 
  mutate(Class = replace(Class, Class == "Aphelidea", "Aphelidea Phylum"),
         Order = replace(Order, Order == "Aphelidea", "Aphelidea Phylum"),
         Family = replace(Family, Family == "Aphelidea", "Aphelidea Phylum"),
         Genus = replace(Genus, Genus == "Aphelidea", "Aphelidea Phylum")) %>% 
  mutate(Genus = replace(Genus, Genus == "Aphelidea Familly", "Aphelidea Phylum")) %>% 
  as.matrix() %>% 
  tax_table()
Sentinel18S_eukaryote

#Create the "Drought" variable
Sentinel18S_eukaryote <- Sentinel18S_eukaryote %>% 
  ps_mutate(Drought = if_else(Flow == "Intermittent" & Time >= "05", "Yes", "No"))

Sentinel18S_2022_eukaryote_filt <- Sentinel18S_eukaryote %>% 
  ps_filter(Year == 2022) %>% 
  filter_taxa(function(x) sum(x > 1) > (0.002*length(x)), TRUE)

#Divide phototrophs
Sentinel18S_photo <- subset_taxa(Sentinel18S_eukaryote, Phylum %in% c("Dinoflagellata", "Diatomea", "Ochrophyta", "Chlorophyta",
                                                                      "Prymnesiophyceae", "Cryptophyceae", "Porphyridiophyceae",
                                                                      "Rhodellophyceae", "Charophyta", "Chlorokybophyceae", "Haptophyta",
                                                                      "Klebsormidiophyceae", "Pavlovophyceae", "Florideophycidae")) 
Sentinel18S_photo

Sentinel18S_others <- subset_taxa(Sentinel18S_eukaryote, !Phylum %in% c("Dinoflagellata", "Diatomea", "Ochrophyta", "Chlorophyta",
                                                                        "Prymnesiophyceae", "Cryptophyceae", "Porphyridiophyceae",
                                                                        "Rhodellophyceae", "Charophyta", "Chlorokybophyceae", "Haptophyta",
                                                                        "Klebsormidiophyceae", "Pavlovophyceae", "Florideophycidae"))
Sentinel18S_others

saveRDS(Sentinel18S_photo, file = "Sentinel18S_photo.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

saveRDS(Sentinel18S_others, file = "Sentinel18S_others.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

#Select 2022 and filter for more than 1 read in at least 2 samples
Sentinel18S_2022_photo <- Sentinel18S_photo %>% 
  ps_filter(Year == 2022)
Sentinel18S_2022_photo
Sentinel18S_2022_photo_filt <- Sentinel18S_2022_photo %>% 
  filter_taxa(function(x) sum(x > 1) > (0.002*length(x)), TRUE)
Sentinel18S_2022_photo_filt

saveRDS(Sentinel18S_2022_photo, file = "Sentinel18S_2022_photo.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
saveRDS(Sentinel18S_2022_photo_filt, file = "Sentinel18S_2022_photo_filt.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

Sentinel18S_2022_others <- Sentinel18S_others %>% 
  ps_filter(Year == 2022)
Sentinel18S_2022_others
Sentinel18S_2022_others_filt <- Sentinel18S_2022_others %>% 
  filter_taxa(function(x) sum(x > 1) > (0.002*length(x)), TRUE)
Sentinel18S_2022_others_filt

saveRDS(Sentinel18S_2022_others, file = "Sentinel18S_2022_others.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
saveRDS(Sentinel18S_2022_others_filt, file = "Sentinel18S_2022_others_filt.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)


##Merging by Samples
sample_data(Sentinel18S_2022_photo_filt)$Sample <- as.factor(sample_data(Sentinel18S_2022_photo_filt)$Sample)
sample_data(Sentinel18S_2022_photo_filt)$Treatment <- as.factor(sample_data(Sentinel18S_2022_photo_filt)$Treatment)
sample_data(Sentinel18S_2022_photo_filt)$Flow <- as.factor(sample_data(Sentinel18S_2022_photo_filt)$Flow)
sample_data(Sentinel18S_2022_photo_filt)$Temperature <- as.factor(sample_data(Sentinel18S_2022_photo_filt)$Temperature)
sample_data(Sentinel18S_2022_photo_filt)$Drought <- as.factor(sample_data(Sentinel18S_2022_photo_filt)$Drought)

Sentinel18S_2022_photo_filt_merged <- merge_samples(Sentinel18S_2022_photo_filt, "Sample", fun = "mean") %>% 
  ps_select(-c("Sample", "Run", "FlumeID")) %>% 
  ps_mutate(Flow = str_replace(Flow, "1", "Constant"), Flow = str_replace(Flow, "2", "Intermittent"), 
            Flow = str_replace(Flow, "3", "Natural"), Flow = str_replace(Flow, "4", "Stochastic"),
            Temperature = str_replace(Temperature, "1", "control"), Temperature = str_replace(Temperature, "2", "warm"),
            Treatment = str_replace(Treatment, "1", "Cc"), Treatment = str_replace(Treatment, "2", "Cw"), 
            Treatment = str_replace(Treatment, "3", "Ic"), Treatment = str_replace(Treatment, "4", "Iw"), 
            Treatment = str_replace(Treatment, "5", "Nc"), Treatment = str_replace(Treatment, "6", "Nw"), 
            Treatment = str_replace(Treatment, "7", "Sc"), Treatment = str_replace(Treatment, "8", "Sw"),
            Drought = str_replace(Drought, "1", "No"), Drought = str_replace(Drought, "2", "Yes")) %>% 
  ps_mutate(Primers = "18S")
otu_table(Sentinel18S_2022_photo_filt_merged) <- t(otu_table(Sentinel18S_2022_photo_filt_merged))
Sentinel18S_2022_photo_filt_merged

saveRDS(Sentinel18S_2022_photo_filt_merged, file = "Sentinel18S_2022_photo_filt_merged.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)


sample_data(Sentinel18S_2022_others_filt)$Sample <- as.factor(sample_data(Sentinel18S_2022_others_filt)$Sample)
sample_data(Sentinel18S_2022_others_filt)$Treatment <- as.factor(sample_data(Sentinel18S_2022_others_filt)$Treatment)
sample_data(Sentinel18S_2022_others_filt)$Flow <- as.factor(sample_data(Sentinel18S_2022_others_filt)$Flow)
sample_data(Sentinel18S_2022_others_filt)$Temperature <- as.factor(sample_data(Sentinel18S_2022_others_filt)$Temperature)
sample_data(Sentinel18S_2022_others_filt)$Drought <- as.factor(sample_data(Sentinel18S_2022_others_filt)$Drought)

Sentinel18S_2022_others_filt_merged <- merge_samples(Sentinel18S_2022_others_filt, "Sample", fun = "mean") %>% 
  ps_select(-c("Sample", "Run", "FlumeID")) %>% 
  ps_mutate(Flow = str_replace(Flow, "1", "Constant"), Flow = str_replace(Flow, "2", "Intermittent"), 
            Flow = str_replace(Flow, "3", "Natural"), Flow = str_replace(Flow, "4", "Stochastic"),
            Temperature = str_replace(Temperature, "1", "control"), Temperature = str_replace(Temperature, "2", "warm"),
            Treatment = str_replace(Treatment, "1", "Cc"), Treatment = str_replace(Treatment, "2", "Cw"), 
            Treatment = str_replace(Treatment, "3", "Ic"), Treatment = str_replace(Treatment, "4", "Iw"), 
            Treatment = str_replace(Treatment, "5", "Nc"), Treatment = str_replace(Treatment, "6", "Nw"), 
            Treatment = str_replace(Treatment, "7", "Sc"), Treatment = str_replace(Treatment, "8", "Sw"),
            Drought = str_replace(Drought, "1", "No"), Drought = str_replace(Drought, "2", "Yes")) %>% 
  ps_mutate(Primers = "18S")
otu_table(Sentinel18S_2022_others_filt_merged) <- t(otu_table(Sentinel18S_2022_others_filt_merged))
Sentinel18S_2022_others_filt_merged

saveRDS(Sentinel18S_2022_others_filt_merged, file = "Sentinel18S_2022_others_filt_merged.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

## Phyloseq 18S Water

#Create dataframe
ASV_table_18S_water <- data.frame(ASV18S_water)
sample_data_18S_water <- data.frame(metadata18S_water)

#Change to phyloseq object
ASV_table_18S_water  <- otu_table(ASV_table_18S_water, taxa_are_rows = TRUE)
sample_data_18S_water <- sample_data(sample_data_18S_water)

#Merging dataframe
Sentinel18S_water <- merge_phyloseq(ASV_table_18S_water, sample_data_18S_water, tax_table_18SFull)
Sentinel18S_water 

#Validation with microViz
phyloseq_validate(Sentinel18S_water, remove_undetected = TRUE)
Sentinel18S_water
Sentinel18S_water <- tax_fix(Sentinel18S_water)
Sentinel18S_water
samdat_tbl(Sentinel18S_water)

#Filter bacteria, Unassigned Kingdom, Arthropoda and not resolved Genus
Sentinel18S_water_eukaryote <- subset_taxa(Sentinel18S_water, !is.na(Kingdom) & !Kingdom %in% c("", "Bacteria","Unassigned") & !Phylum %in% c("Arthropoda") & !Genus %in% c("Eukaryota Kingdom", "Eukaryota"))
Sentinel18S_water_eukaryote <- filter_taxa(Sentinel18S_water_eukaryote, function(x) sum(x) > 1, TRUE)
Sentinel18S_water_eukaryote

Sentinel18S_2022_water_eukaryote <- Sentinel18S_water_eukaryote %>% ps_filter(Year == 2022) %>% tax_fix(unknowns = c("Incertae_Sedis", "Incertae_Sedis Family", "uncultured","Labyrinthulomycetes"))
Sentinel18S_2022_water_eukaryote

tax_table(Sentinel18S_2022_water_eukaryote) <- tax_table(Sentinel18S_2022_water_eukaryote) %>% 
  data.frame() %>% 
  mutate(Class = replace(Class, Class == "Eukaryota Kingdom Phylum", "Eukaryota Kingdom")) %>% 
  as.matrix() %>% 
  tax_table()

## Venn Diagram of Water samples compared to biofilm samples
library(MicEco)

phyloseq_2022_water_venn_tax_18S <- tax_table(Sentinel18S_2022_water_eukaryote)
phyloseq_2022_water_venn_asv_18S <- otu_table(Sentinel18S_2022_water_eukaryote)
phyloseq_2022_water_venn_meta_18S <- sample_data(Sentinel18S_2022_water_eukaryote)

Sentinel18S_2022_eukaryote_filt_merged <- merge_phyloseq(Sentinel18S_2022_photo_filt_merged, Sentinel18S_2022_others_filt_merged)
phyloseq_2022_biofilm_venn_tax_18S <- tax_table(Sentinel18S_2022_eukaryote_filt)
phyloseq_2022_biofilm_venn_asv_18S <- otu_table(Sentinel18S_2022_eukaryote_filt)
phyloseq_2022_biofilm_venn_meta_18S <- sample_data(Sentinel18S_2022_eukaryote_filt)

taxvenn_18S <- merge_phyloseq(phyloseq_2022_water_venn_tax_18S, phyloseq_2022_biofilm_venn_tax_18S)
asvvenn_18S <- merge_phyloseq(phyloseq_2022_water_venn_asv_18S, phyloseq_2022_biofilm_venn_asv_18S)
metavenn_18S <- merge_phyloseq(phyloseq_2022_water_venn_meta_18S, phyloseq_2022_biofilm_venn_meta_18S)

Venn_18S_2022 <- merge_phyloseq(taxvenn_18S, asvvenn_18S, metavenn_18S)

MicEco_18S <- ps_venn(Venn_18S_2022, group = "Type")
MicEco_18S

MicEco_18S_abundance <- ps_venn(Venn_18S_2022, group = "Type", weight = TRUE)
MicEco_18S_abundance

abundance.round <- MicEco_18S_abundance$data$original.values
abundance.round <- round(as.numeric(abundance.round), 4)
MicEco_18S_abundance$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.1$children$tag.quantity.1$label <- 3.61
MicEco_18S_abundance$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.2$children$tag.quantity.2$label <- 3.04
MicEco_18S_abundance$children$canvas.grob$children$diagram.grob.1$children$tags$children$tag.number.3$children$tag.quantity.3$label <- 93.35

FigureS10 <- ggarrange(MicEco_18S, MicEco_18S_abundance, ncol = 1, nrow = 2, labels = c("A", "B"), font.label = list(size = 18))
FigureS10
ggsave("FigureS10.pdf", FigureS10, width = 20,
       height = 20,
       units = "cm")

##Count number of reads
Sentinel18S_2022_eukaryote_filt
reads_18S <- sum(unlist(otu_table(Sentinel18S_2022_eukaryote_filt)), na.rm = TRUE)
Sentinel18S_2022_photo_filt
reads_18S_photo <- sum(unlist(otu_table(Sentinel18S_2022_photo_filt)), na.rm = TRUE)
percent_photo <- (reads_18S_photo/reads_18S)*100
Sentinel18S_2022_others_filt
reads_18S_others <- sum(unlist(otu_table(Sentinel18S_2022_others_filt)), na.rm = TRUE)
percent_others <- (reads_18S_others/reads_18S)*100


Sentinel18S_2022_water_eukaryote
reads_18S_water <- sum(unlist(otu_table(Sentinel18S_2022_water_eukaryote)), na.rm = TRUE)

##Divide water phyloseq by phototrophs
Sentinel18S_2022_water_photo <- subset_taxa(Sentinel18S_2022_water_eukaryote, Phylum %in% c("Dinoflagellata", "Diatomea", "Ochrophyta", "Chlorophyta",
                                                                      "Prymnesiophyceae", "Cryptophyceae", "Porphyridiophyceae",
                                                                      "Rhodellophyceae", "Charophyta", "Chlorokybophyceae", "Haptophyta",
                                                                      "Klebsormidiophyceae", "Pavlovophyceae", "Florideophycidae")) 
Sentinel18S_2022_water_photo

Sentinel18S_2022_water_others <- subset_taxa(Sentinel18S_2022_water_eukaryote, !Phylum %in% c("Dinoflagellata", "Diatomea", "Ochrophyta", "Chlorophyta",
                                                                        "Prymnesiophyceae", "Cryptophyceae", "Porphyridiophyceae",
                                                                        "Rhodellophyceae", "Charophyta", "Chlorokybophyceae", "Haptophyta",
                                                                        "Klebsormidiophyceae", "Pavlovophyceae", "Florideophycidae"))
Sentinel18S_2022_water_others




#Plot Graph-Bar for water samples (photo)
Graphbar_18S_photo_phylum_water <- Sentinel18S_2022_water_photo %>%
  merge_samples(group = "Succession") %>%
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 19,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    palette = distinct_palette(19, pal = "kelly"),
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
        plot.margin = margin(1, 1, 1.5, 1, "cm")) +
  guides(fill = guide_legend(ncol = 3))
Graphbar_18S_photo_phylum_water

Graphbar_18S_photo_class_water <- Sentinel18S_2022_water_photo %>%
  merge_samples(group = "Succession") %>%
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
Graphbar_18S_photo_class_water

#Creating new palette
myPal_phylum <- tax_palette(data = Sentinel18S_2022_water_others, rank = "Phylum", n = 20, pal = "greenArmytage", add = c(Other = "lightgray"))
myPal_phylum[1] <- "darkorange"
myPal_class <- tax_palette(data = Sentinel18S_2022_water_others, rank = "Class", n = 20, pal = "greenArmytage", add = c(Other = "lightgray"))
myPal_class[1] <- "darkorange"

#Plot Graph-Bar for water samples (photo)
Graphbar_18S_others_phylum_water <- Sentinel18S_2022_water_others %>%
  merge_samples(group = "Succession") %>%
  ps_arrange(Succession) %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 14,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    palette = myPal_phylum,
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
Graphbar_18S_others_phylum_water

Graphbar_18S_others_class_water <- Sentinel18S_2022_water_others %>%
  merge_samples(group = "Succession") %>%
  ps_arrange(Succession) %>%
  tax_fix(unknowns = c("Incertae_Sedis", "Labyrinthulomycetes", "uncultured")) %>% 
  comp_barplot(
    tax_level = "Class", n_taxa = 14,
    x = "Time",
    bar_outline_colour = "grey5",
    merge_other = FALSE,
    sample_order = "default",
    palette = myPal_class,
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
Graphbar_18S_others_class_water


### NMDS biofilm and water --------------------------------------------------

##Hellinger transformation
Sentinel18S_2022_eukaryote_filt_hell <- microbiome::transform(Sentinel18S_2022_eukaryote_filt, 'hell')
Sentinel18S_2022_eukaryote_filt_merged_hell <- microbiome::transform(Sentinel18S_2022_eukaryote_filt_merged, 'hell')
Sentinel18S_2022_photo_filt_hell <- microbiome::transform(Sentinel18S_2022_photo_filt, 'hell')
Sentinel18S_2022_photo_filt_merged_hell <- microbiome::transform(Sentinel18S_2022_photo_filt_merged, 'hell')
Sentinel18S_2022_others_filt_hell <- microbiome::transform(Sentinel18S_2022_others_filt, 'hell')
Sentinel18S_2022_others_filt_merged_hell <- microbiome::transform(Sentinel18S_2022_others_filt_merged, 'hell')

Sentinel18S_2022_water_photo_hell <- microbiome::transform(Sentinel18S_2022_water_photo, 'hell')
Sentinel18S_2022_water_others_hell <- microbiome::transform(Sentinel18S_2022_water_others, 'hell')

#Phototrophic eukaryotes
phyloseq_2022_water_18S_photo_tax <- tax_table(Sentinel18S_2022_water_photo_hell)
phyloseq_2022_water_18S_photo_asv <- otu_table(Sentinel18S_2022_water_photo_hell)
phyloseq_2022_water_18S_photo_meta <- sample_data(Sentinel18S_2022_water_photo_hell)

phyloseq_2022_biofilm_18S_photo_tax <- tax_table(Sentinel18S_2022_photo_filt_merged_hell)
phyloseq_2022_biofilm_18S_photo_asv <- otu_table(Sentinel18S_2022_photo_filt_merged_hell)
phyloseq_2022_biofilm_18S_photo_meta <- sample_data(Sentinel18S_2022_photo_filt_merged_hell) %>% data.frame() %>% mutate(Type = as.character(Type), Type = if_else(is.na(Type), "Biofilm", Type)) %>% sample_data()

tax_18S_photo_all <- merge_phyloseq(phyloseq_2022_water_18S_photo_tax, phyloseq_2022_biofilm_18S_photo_tax)
asv_18S_photo_all <- merge_phyloseq(phyloseq_2022_water_18S_photo_asv, phyloseq_2022_biofilm_18S_photo_asv)
meta_18S_photo_all <- merge_phyloseq(phyloseq_2022_water_18S_photo_meta, phyloseq_2022_biofilm_18S_photo_meta)

NMDS_all_18S_photo_2022 <- merge_phyloseq(tax_18S_photo_all, asv_18S_photo_all, meta_18S_photo_all)

matrix_dist_18S_photo_2022_all <- as.matrix(t(data.frame(otu_table(NMDS_all_18S_photo_2022))))
metadata_dist_18S_photo_2022_all <- data.frame(sample_data(NMDS_all_18S_photo_2022)) %>% rownames_to_column("Samples")
set.seed(19950930)
dist_18S_photo_2022_bray_all <- vegdist(matrix_dist_18S_photo_2022_all, method ="bray") %>% 
  as.matrix() %>% 
  data.frame() %>% 
  rownames_to_column("Samples")
dist_18S_photo_2022_all <- dist_18S_photo_2022_bray_all %>% 
  dplyr::select(all_of(.[["Samples"]])) %>% 
  as.dist()
nmds_18S_photo_2022_all <- metaMDS(dist_18S_photo_2022_all, trymax=100)
stress_18S_photo_all <- nmds_18S_photo_2022_all$stress
stress_18S_photo_all

scores_nmds_18S_photo_2022_all <- scores(nmds_18S_photo_2022_all) %>% 
  as_tibble(rownames = "Samples") %>% 
  inner_join(., metadata_dist_18S_photo_2022_all, by="Samples")

scores_nmds_18S_photo_2022_all  <- scores_nmds_18S_photo_2022_all[order(scores_nmds_18S_photo_2022_all$Treatment, scores_nmds_18S_photo_2022_all$Time),]

#Differences?
adonis_bc_18S_photo_all <- adonis2(matrix_dist_18S_photo_2022_all ~ Type, data=metadata_dist_18S_photo_2022_all, permutations=999, method="bray")
adonis_bc_18S_photo_all #Type are significant 

NMDS_18S_photo_2022_all <- ggplot(scores_nmds_18S_photo_2022_all, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = Succession, shape = Type), size = 4, alpha = 0.8) +
  scale_color_viridis(discrete = FALSE, option = "magma", breaks = c(0, 20, 40, 60, 80, 100)) +
  theme_classic() +
  guides(colour = guide_colourbar(barheight = 10)) +
  labs(color = "Days of growth", caption = "Stress 0.143") +
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
NMDS_18S_photo_2022_all

#Non-phototrophic eukaryotes
phyloseq_2022_water_18S_others_tax <- tax_table(Sentinel18S_2022_water_others_hell)
phyloseq_2022_water_18S_others_asv <- otu_table(Sentinel18S_2022_water_others_hell)
phyloseq_2022_water_18S_others_meta <- sample_data(Sentinel18S_2022_water_others_hell)

phyloseq_2022_biofilm_18S_others_tax <- tax_table(Sentinel18S_2022_others_filt_merged_hell)
phyloseq_2022_biofilm_18S_others_asv <- otu_table(Sentinel18S_2022_others_filt_merged_hell)
phyloseq_2022_biofilm_18S_others_meta <- sample_data(Sentinel18S_2022_others_filt_merged_hell) %>% data.frame() %>% mutate(Type = as.character(Type), Type = if_else(is.na(Type), "Biofilm", Type)) %>% sample_data()

tax_18S_others_all <- merge_phyloseq(phyloseq_2022_water_18S_others_tax, phyloseq_2022_biofilm_18S_others_tax)
asv_18S_others_all <- merge_phyloseq(phyloseq_2022_water_18S_others_asv, phyloseq_2022_biofilm_18S_others_asv)
meta_18S_others_all <- merge_phyloseq(phyloseq_2022_water_18S_others_meta, phyloseq_2022_biofilm_18S_others_meta)

NMDS_all_18S_others_2022 <- merge_phyloseq(tax_18S_others_all, asv_18S_others_all, meta_18S_others_all)

matrix_dist_18S_others_2022_all <- as.matrix(t(data.frame(otu_table(NMDS_all_18S_others_2022))))
metadata_dist_18S_others_2022_all <- data.frame(sample_data(NMDS_all_18S_others_2022)) %>% rownames_to_column("Samples")
set.seed(19950930)
dist_18S_others_2022_bray_all <- vegdist(matrix_dist_18S_others_2022_all, method ="bray") %>% 
  as.matrix() %>% 
  data.frame() %>% 
  rownames_to_column("Samples")
dist_18S_others_2022_all <- dist_18S_others_2022_bray_all %>% 
  dplyr::select(all_of(.[["Samples"]])) %>% 
  as.dist()
nmds_18S_others_2022_all <- metaMDS(dist_18S_others_2022_all, trymax=100)
stress_18S_others_all <- nmds_18S_others_2022_all$stress
stress_18S_others_all

scores_nmds_18S_others_2022_all <- scores(nmds_18S_others_2022_all) %>% 
  as_tibble(rownames = "Samples") %>% 
  inner_join(., metadata_dist_18S_others_2022_all, by="Samples")

scores_nmds_18S_others_2022_all  <- scores_nmds_18S_others_2022_all[order(scores_nmds_18S_others_2022_all$Treatment, scores_nmds_18S_others_2022_all$Time),]

#Differences?
adonis_bc_18S_others_all <- adonis2(matrix_dist_18S_others_2022_all ~ Type, data=metadata_dist_18S_others_2022_all, permutations=999, method="bray")
adonis_bc_18S_others_all #Type are significant 

NMDS_18S_others_2022_all <- ggplot(scores_nmds_18S_others_2022_all, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour = Succession, shape = Type), size = 4, alpha = 0.8) +
  scale_color_viridis(discrete = FALSE, option = "mako", breaks = c(0, 20, 40, 60, 80, 100)) +
  theme_classic() +
  guides(colour = guide_colourbar(barheight = 10)) +
  labs(color = "Days of growth", caption = "Stress 0.089") +
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
NMDS_18S_others_2022_all

