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

################### 16S Phyloseq Object Setup ###################
data16S <- read.delim("SampleData_16S2.txt")
colnames(data16S) <- gsub("X16S.Behnsen.MY.", "MY", colnames(data16S))
colnames(data16S) <- gsub("X16S.Negative.Control.sample2", "Negative Control 1", colnames(data16S))
colnames(data16S) <- gsub("X16S.Negative.Ctrl.2.Underhill", "Negative Control 2", colnames(data16S))
colnames(data16S) <- gsub("X16S.Positive.Ctrl.2.Underhill", "Positive Control 2", colnames(data16S))
data16S <- filter(data16S, Taxonomy != "k__Bacteria ; p__ NA ; c__ NA ; o__ NA ; f__ NA ; g__ NA ; s__ NA")
data16S <- filter(data16S, Taxonomy != "k__Eukaryota ; p__ NA ; c__ NA ; o__ NA ; f__ NA ; g__ NA ; s__ NA")
data16S <- filter(data16S, Taxonomy != "k__NA ; p__ NA ; c__ NA ; o__ NA ; f__ NA ; g__ NA ; s__ NA")

asvdata16S <- data16S %>%
  select(!c(Taxonomy))
row.names(asvdata16S) <- asvdata16S$ID
asvdata16S <- asvdata16S %>% select(-ID)
asvdata16S_matrix = as.matrix(asvdata16S)

taxtree16S <- data16S %>%
  select(ID, Taxonomy)
taxtree16S = separate(data = taxtree16S, col = Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\;" )
taxtree16S <- taxtree16S %>%
  mutate_at("Kingdom", str_replace_all, "k__", "") %>%
  mutate_at("Phylum", str_replace_all, "p__", "") %>%
  mutate_at("Class", str_replace_all, "c__", "") %>%
  mutate_at("Order", str_replace_all, "o__", "") %>%
  mutate_at("Family", str_replace_all, "f__", "") %>%
  mutate_at("Genus", str_replace_all, "g__", "") %>%
  mutate_at("Species", str_replace_all, "s__", "")

row.names(taxtree16S) <- taxtree16S$ID
taxtree16S <- taxtree16S %>% select(-ID)
taxtree16S_matrix = as.matrix(taxtree16S)

data32 <- read.delim("samplelabels.txt")
row.names(data32) <- data32$Samples
data32$Diet[data32$Diet == "chow 1"] <- "chow LM485"
data32$Diet[data32$Diet == "50ppm iron"] <- "purified 50ppm iron"
data32$Diet[data32$Diet == "chow 2"] <- "chow 2914"
data32$Bug[data32$Bug == "Salmonella"] <- "STM"

physeq17 = phyloseq(otu_table(asvdata16S_matrix, taxa_are_rows = TRUE), tax_table(taxtree16S_matrix), sample_data(data32))
random_tree16S = rtree(ntaxa(physeq17), rooted=TRUE, tip.label=taxa_names(physeq17))
physeq17 = merge_phyloseq(physeq17, random_tree16S)

noyeast <- c("MY17", "MY18", "MY19", "MY20", "MY21", "MY35", "MY36", "MY37", "MY38", "MY39")
physeq16_noyeast <- subset_samples(physeq17, (Samples %in% c(LM485_50ppm, noyeast)))

################### Theme Setup ###################
cleanbg <- theme(panel.background = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 text = element_text(family = "Arial", size = 50))

rel.abund.theme <- theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.75),
        panel.background = element_blank(),
        text = element_text(family = "Arial", size = 50),
        legend.text.align = 0)

dietcol <- scale_fill_manual(values = c("#1F78B4", "#33A02C"))

Dietgroup <- c("chow LM485", "chow LM485", "chow LM485", "chow LM485", "chow LM485", "chow 2914", "chow 2914", "chow 2914", "chow 2914", "chow 2914")

box_aes <- theme(panel.background = element_blank(),
                 panel.border = element_rect(fill = NA, color = "black"),
                 aspect.ratio = 1)

cellhosttheme <- theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.75),
                       text = element_text(family = "Arial", size = 50))

color.pal.16S <- read.delim("Genus color palette 16S.txt", sep = "")

################### Stacked Relative Abundance ###################
gen.glom.16 <- tax_glom(physeq16_noyeast, taxrank = rank_names(physeq16_noyeast)[6], NArm = FALSE)
glom_rel_gen_noyeast16 <- psmelt(phyloseq::transform_sample_counts(gen.glom.16, function(x){x / sum(x)}))
glom_rel_gen_noyeast16$Genus[glom_rel_gen_noyeast16$Abundance<0.02] <- "Other <2%"

gen.legend.16S.noyeast <- glom_rel_gen_noyeast16 %>%
  select(Genus) %>%
  filter(Genus !=  "Other <2%") %>%
  arrange(Genus)
glom_rel_gen_noyeast16_list <- as.list(unique(gen.legend.16S.noyeast$Genus))
italic.gen.legend.16S.noyeast <- mixedFontLabel(glom_rel_gen_noyeast16_list, italic = TRUE)

gen.unique <- unique(glom_rel_gen_noyeast16$Genus)
gen.unique.df <- as.data.frame(gen.unique)
colnames(gen.unique.df) <- "Genus"
gen.unique.df <- gen.unique.df %>%
  arrange(Genus)

joined.col <- dplyr::inner_join(gen.unique.df, color.pal.16S, by = "Genus")
col.chow <- joined.col$pal.generator %>%
  as.character(joined.col$pal.generator)
col.chow <- dput(col.chow)

relabund_16S_noyeast_Gen <- ggplot(glom_rel_gen_noyeast16, aes(x = Mouse, y = Abundance, fill = Genus)) +
  geom_bar(color = "black", stat = "identity", position = "stack", width = 0.95) +
  facet_wrap(~factor(Diet, levels = c("chow LM485", "chow 2914")), scales = "free", nrow = 1) +
  labs(x = "", y = "Relative Abundance\n", title = "16S") +
  rel.abund.theme +
  scale_fill_manual(values = col.chow, labels = c(italic.gen.legend.16S.noyeast, "Other <2%")) +
  scale_y_continuous(expand = c(0.005,0.005)) +
  theme(legend.position = "right") +
  guides(fill = guide_legend(ncol = 1))

png(filename = "Relative Abundance_16S_NoYeast_Genera.png", width = 12000, height = 9600, units = "px", res = 300)
plot(relabund_16S_noyeast_Gen)
dev.off()

relabund.16S.noyeast.genera.legend.fig <-get_legend(relabund_16S_noyeast_Gen)
png(filename = "Relative Abundance_16S_NoYeast_Genera legend.png", width = 12000, height = 9600, units = "px", res = 300)
plot(relabund.16S.noyeast.genera.legend.fig)
dev.off()

relabund.16S.noyeast <- glom_rel_gen_noyeast16 %>%
  select(OTU, Abundance, Mouse, Diet, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  filter(Abundance > 0) %>%
  arrange(Mouse, desc(Diet))

write_xlsx(relabund.16S.noyeast, "Table S4 - Relative Abundance No Yeast 16S.xlsx")

################### Alpha Diversity  ###################
tab_16S_noyeast <- microbiome::alpha(physeq16_noyeast, index = "all")
tab_16S_noyeast$Dietgroup <- Dietgroup
write_xlsx(tab_16S_noyeast, "Alpha Diversity No Yeast 16S.xlsx")

################### Beta Diversity  ###################
sample_data(physeq16_noyeast)$Diet <- factor(sample_data(physeq16_noyeast)$Diet, levels = c("chow LM485", "chow 2914"))
DistBCnoyeast = distance(physeq16_noyeast, method = "bray")
ordBCnoyeast = ordinate(physeq16_noyeast, method="PCoA", distance = DistBCnoyeast)
plot_scree(ordBCnoyeast)

PCoA16Snoyeast <- plot_ordination(physeq16_noyeast, ordBCnoyeast, color = "Diet", shape = "Diet", label = "Mouse") +
  stat_ellipse() +
  cleanbg +
  scale_color_manual(values = c("#1F78B4", "#33A02C")) +
  ggtitle("16S")
png(filename = "PCoA 16S No Yeast.png", width = 4800, height = 3600, units = "px", res = 300)
plot(PCoA16Snoyeast)
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
