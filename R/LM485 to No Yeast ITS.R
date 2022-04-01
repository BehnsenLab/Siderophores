library("ggplot2")
library("ggpubr")
library("ggsignif")
library("ggfortify")
library("RColorBrewer")
library("gt")
library("knitr")
library("stats")
library("dplyr")
library("plyr")
library("dada2")
library("phyloseq")
library("png")
library("microbiome")
library("ape")
library("reshape2")
library("data.table")
library("vegan")
library("ade4")
library("ggh4x")
library("writexl")
library("tidyr")
library("stringr")

################### ITS Phyloseq Object Setup ###################
dataITS <- read.delim("SampleData_ITS2.txt")
colnames(dataITS) <- gsub("ITS.Behnsen.MY.", "MY", colnames(dataITS))
colnames(dataITS) <- gsub("ITS.Negative.Control.sample2", "Negative Control 1", colnames(dataITS))
colnames(dataITS) <- gsub("ITS.Negative.Ctrl.2.Underhill", "Negative Control 2", colnames(dataITS))
colnames(dataITS) <- gsub("ITS.Positive.Ctrl.2.Underhill", "Positive 2", colnames(dataITS))
dataITS1 <- filter(dataITS, Taxonomy != "k__Fungi ; NA ; NA ; NA ; NA ; NA ; NA")

asvdataITS <- dataITS %>%
  select(!c(Taxonomy))
row.names(asvdataITS) <- asvdataITS$ID
asvdataITS <- asvdataITS %>% select(-ID)
asvdataITS_matrix = as.matrix(asvdataITS)

taxtreeITS <- dataITS1 %>%
  select(ID, Taxonomy)
taxtreeITS = separate(data = taxtreeITS, col = Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\;" )

taxtreeITS <- taxtreeITS %>%
  mutate_at("Kingdom", str_replace_all, "k__", "") %>%
  mutate_at("Phylum", str_replace_all, "p__", "") %>%
  mutate_at("Class", str_replace_all, "c__", "") %>%
  mutate_at("Order", str_replace_all, "o__", "") %>%
  mutate_at("Family", str_replace_all, "f__", "") %>%
  mutate_at("Genus", str_replace_all, "g__", "") %>%
  mutate_at("Species", str_replace_all, "s__", "")
row.names(taxtreeITS) <- taxtreeITS$ID
taxtreeITS <- taxtreeITS %>% select(-ID)
taxtreeITS_matrix = as.matrix(taxtreeITS)

data31 <- read.delim("samplelabels.txt")
row.names(data31) <- data31$Samples
data31$Diet[data31$Diet == "chow 1"] <- "chow LM485"
data31$Diet[data31$Diet == "chow 2"] <- "chow 2914"
data31$Diet[data31$Diet == "50ppm iron"] <- "purified 50ppm iron"
data31$Bug[data31$Bug == "Salmonella"] <- "STM"

physeq2 = phyloseq(otu_table(asvdataITS_matrix, taxa_are_rows = TRUE), tax_table(taxtreeITS_matrix), sample_data(data31))
random_tree2 = rtree(ntaxa(physeq2), rooted = TRUE, tip.label = taxa_names(physeq2))
physeq19 = merge_phyloseq(physeq2, random_tree2)

noyeast <- c("MY17", "MY18", "MY19", "MY20", "MY21", "MY35", "MY36", "MY37", "MY38", "MY39")
physeq18_noyeast <- subset_samples(physeq19, (Samples %in% c(LM485_50ppm, noyeast)))

################### Theme Setup ###################
cleanbg <- theme(panel.background = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 text = element_text(family = "Arial", size = 50))

rel.abund.theme <- theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.75),
        panel.background = element_blank(),
        text = element_text(family = "Arial", size = 50),
        legend.text.align = 0)

dietcol <- scale_fill_manual(values = c("#1F78B4", "#33A02C"))

Dietgroup <- c("chow LM485",
               "chow LM485",
               "chow LM485",
               "chow LM485",
               "chow LM485",
               "chow 2914",
               "chow 2914",
               "chow 2914",
               "chow 2914",
               "chow 2914")

box_aes <- theme(panel.background = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 aspect.ratio = 1)

cellhosttheme <- theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.75),
                       text = element_text(family = "Arial", size = 50))

color.pal.ITS <- read.delim("Genus color palette ITS.txt", sep = "")

################### Stacked Relative Abundance ###################
gen.glom.ITS <- tax_glom(physeq18_noyeast, taxrank = rank_names(physeq18_noyeast)[6], NArm = FALSE)
glom_rel_gen_noyeastITS<- psmelt(phyloseq::transform_sample_counts(gen.glom.ITS, function(x){x / sum(x)}))
glom_rel_gen_noyeastITS$Genus[glom_rel_gen_noyeastITS$Abundance<0.02] <- "Other <2%"

gen.legend.ITS.noyeast <- glom_rel_gen_noyeastITS %>%
  select(Genus) %>%
  filter(Genus !=  "Other <2%") %>%
  arrange(Genus)
gen.legend.ITS.noyeast.list <- as.list(unique(gen.legend.ITS.noyeast$Genus))
italic.gen.legend.ITS.noyeast <- mixedFontLabel(gen.legend.ITS.noyeast.list, italic = TRUE)

gen.unique.ny <- unique(glom_rel_gen_noyeastITS$Genus)
gen.unique.df.ny <- as.data.frame(gen.unique.ny)
colnames(gen.unique.df.ny) <- "Genus"
gen.unique.df.ny <- gen.unique.df.ny %>%
  arrange(Genus)

joined.col.ny <- dplyr::inner_join(gen.unique.df.ny, color.pal.ITS, by = "Genus")
col.chow.ny<- joined.col.ny$pal.generator %>%
  as.character(joined.col.ny$pal.generator)
col.chow.ny <- dput(col.chow.ny)

relabund_ITS_noyeast_Gen <- ggplot(glom_rel_gen_noyeastITS, aes(x = Mouse, y = Abundance, fill = Genus)) +
  geom_bar(aes(), color = "black", stat = "identity", position = "stack", width = 0.95) +
  facet_wrap(~factor(Diet, levels = c("chow LM485", "chow 2914")), scales = "free", nrow = 1) +
  labs(x = "", y = "Relative Abundance\n", title = "ITS") +
  rel.abund.theme +
  scale_fill_manual(values = col.chow.ny, labels = c(italic.gen.legend.ITS.noyeast, "Other <2%")) +
  scale_y_continuous(expand = c(0.005,0.005)) +
  theme(legend.position = "right") +
  guides(fill = guide_legend(ncol = 2))

png(filename = "Relative Abundance_ITS_NoYeast_Genera.png", width = 12000, height = 9600, units = "px", res = 300)
plot(relabund_ITS_noyeast_Gen)
dev.off()

relabund.ITS.noyeast.genera.legend.fig <-get_legend(relabund_ITS_noyeast_Gen)
png(filename = "Relative Abundance_ITS_NoYeast_Genera legend.png", width = 12000, height = 9600, units = "px", res = 300)
plot(relabund.ITS.noyeast.genera.legend.fig)
dev.off()

relabund.ITS.noyeast <- glom_rel_gen_noyeastITS %>%
  select(OTU, Abundance, Mouse, Diet, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  filter(Abundance > 0) %>%
  arrange(Mouse, desc(Diet))

write_xlsx(relabund.ITS.noyeast, "Table S5 - Relative Abundance No Yeast ITS.xlsx")

################### Alpha Diversity  ###################
tab_18S_noyeast <- microbiome::alpha(physeq18_noyeast, index = "all")
tab_18S_noyeast$Dietgroup <- Dietgroup
write_xlsx(tab_18S_noyeast, "Alpha Diversity No Yeast ITS.xlsx")

################### Beta Diversity  ###################
sample_data(physeq18_noyeast)$Diet <- factor(sample_data(physeq18_noyeast)$Diet, levels = c("chow LM485", "chow 2914"))
DistBCITSnoyeast = distance(physeq18_noyeast, method = "bray")
ordBCITSnoyeast = ordinate(physeq18_noyeast, method="PCoA", distance = DistBCITSnoyeast)
plot_scree(ordBCITSnoyeast)
PCoAITSnoyeast <- plot_ordination(physeq18_noyeast, ordBCITSnoyeast, color = "Diet", shape = "Diet", label = "Mouse") +
  stat_ellipse() +
  cleanbg +
  scale_color_manual(values = c("#1F78B4", "#33A02C")) +
  ggtitle("ITS")

png(filename = "PCoA ITS chow LM485 No Yeast.png", width = 4800, height = 3600, units = "px", res = 300)
plot(PCoAITSnoyeast)
dev.off()

################### References ###################
citation(package="ggplot2")
citation(package="ggplot2")
citation(package="ggpubr")
citation(package="ggsignif")
citation(package="ggfortify")
citation(package="RColorBrewer")
citation(package="gt")
citation(package="knitr")
citation(package="stats")
citation(package="dplyr")
citation(package="plyr")
citation(package="dada2")
citation(package="phyloseq")
citation(package="png")
citation(package="microbiome")
citation(package="ape")
citation(package="reshape2")
citation(package="data.table")
citation(package="vegan")
citation(package="ade4")
citation(package="ggh4x")
citation(package="writexl")
citation(package="tidyr")
citation(package="stringr")
