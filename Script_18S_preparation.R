
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
library(agricolae)
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
metadata18S_water <- subset(metadata18S_water, select = -c(Flow,Temperature, Treatment, Succession, FlumeID, TAD, ChlA, BA, Run)) %>% 
  column_to_rownames(var="NAME")
metadata18S_water$Time <- sprintf("%02d", as.numeric(metadata18S_water$Time))


#Parse taxonomy file
taxonomy18SFull <- subset(ASV18SFull, select = c(OTUID,taxonomy)) %>%
  separate(taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "; ") %>%
  column_to_rownames(var="OTUID")
taxonomy18SFull <- as.matrix(taxonomy18SFull)
taxonomy18SFull <- gsub("d__","",
                        gsub("p__","",
                             gsub("c__","",
                                  gsub("o__","",
                                       gsub("f__","",
                                            gsub("g__","",
                                                 gsub("s__","", taxonomy18SFull)))))))

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

Sentinel18S_eukaryote <- subset_taxa(Sentinel18S, !Phylum %in% c("Arthropoda") & !Genus %in% c("Eukaryota Kingdom", "Eukaryota"))
Sentinel18S_eukaryote

#Create the "Drought" variable
Sentinel18S_eukaryote <- Sentinel18S_eukaryote %>% 
  ps_mutate(Drought = if_else(Flow == "Intermittent" & Time >= "05", "Yes", "No"))

saveRDS(Sentinel18S_eukaryote, file = "Sentinel18S_eukaryote.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

#Select 2022 and filter for more than 1 read in at least 2 samples
Sentinel18S_2022_eukaryote <- Sentinel18S_eukaryote %>% 
  ps_filter(Year == 2022)
Sentinel18S_2022_eukaryote
Sentinel18S_2022_eukaryote_filt <- Sentinel18S_2022_eukaryote %>% 
  filter_taxa(function(x) sum(x > 1) > (0.002*length(x)), TRUE)
Sentinel18S_2022_eukaryote_filt

saveRDS(Sentinel18S_2022_eukaryote, file = "Sentinel18S_2022_eukaryote.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
saveRDS(Sentinel18S_2022_eukaryote_filt, file = "Sentinel18S_2022_eukaryote_filt.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

##Merging by Samples
sample_data(Sentinel18S_2022_eukaryote_filt)$Sample <- as.factor(sample_data(Sentinel18S_2022_eukaryote_filt)$Sample)
sample_data(Sentinel18S_2022_eukaryote_filt)$Treatment <- as.factor(sample_data(Sentinel18S_2022_eukaryote_filt)$Treatment)
sample_data(Sentinel18S_2022_eukaryote_filt)$Flow <- as.factor(sample_data(Sentinel18S_2022_eukaryote_filt)$Flow)
sample_data(Sentinel18S_2022_eukaryote_filt)$Temperature <- as.factor(sample_data(Sentinel18S_2022_eukaryote_filt)$Temperature)
sample_data(Sentinel18S_2022_eukaryote_filt)$Drought <- as.factor(sample_data(Sentinel18S_2022_eukaryote_filt)$Drought)

Sentinel18S_2022_eukaryote_filt_merged <- merge_samples(Sentinel18S_2022_eukaryote_filt, "Sample", fun = "mean") %>% 
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
otu_table(Sentinel18S_2022_eukaryote_filt_merged) <- t(otu_table(Sentinel18S_2022_eukaryote_filt_merged))
Sentinel18S_2022_eukaryote_filt_merged

saveRDS(Sentinel18S_2022_eukaryote_filt_merged, file = "Sentinel18S_2022_eukaryote_filt_merged.rds", ascii = FALSE, version = NULL,
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
Sentinel18S_water_eukaryote

Sentinel18S_2022_water_eukaryote <- Sentinel18S_water_eukaryote %>% ps_filter(Year == 2022)
Sentinel18S_2022_water_eukaryote


## Venn Diagram of Water samples compared to biofilm samples
library(MicEco)

phyloseq_2022_water_venn_tax_18S <- tax_table(Sentinel18S_2022_water_eukaryote)
phyloseq_2022_water_venn_asv_18S <- otu_table(Sentinel18S_2022_water_eukaryote)
phyloseq_2022_water_venn_meta_18S <- sample_data(Sentinel18S_2022_water_eukaryote)

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

FigureS7 <- ggarrange(MicEco_18S, MicEco_18S_abundance, ncol = 1, nrow = 2, labels = c("A", "B"), font.label = list(size = 18))
FigureS7
ggsave("FigureS7.pdf", FigureS7, width = 20,
       height = 20,
       units = "cm")
