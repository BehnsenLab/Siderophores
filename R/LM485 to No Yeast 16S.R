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

asvdata16S <- data16S %>% select(!c(Taxonomy))
row.names(asvdata16S) <- asvdata16S$ID
asvdata16S <- asvdata16S %>% select(-ID)
asvdata16S_matrix = as.matrix(asvdata16S)

taxtree16S <- data16S %>% select(ID, Taxonomy)
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

data31 <- read.delim("samplelabels.txt")
row.names(data31) <- data31$Samples
data31$Diet[data31$Diet == "chow 1"] <- "chow LM485"
data31$Diet[data31$Diet == "chow 2"] <- "chow 2914"
data31$Bug[data31$Bug == "Salmonella"] <- "STM"

physeq = phyloseq(otu_table(asvdata16S_matrix, taxa_are_rows = TRUE), tax_table(taxtree16S_matrix), sample_data(data31))
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
physeq16 = merge_phyloseq(physeq, random_tree)

JAX <- c("MY12", "MY13", "MY14", "MY15", "MY16", "MY30", "MY31", "MY32", "MY33", "MY34", "MY63", "MY64")
physeq16 <- subset_samples(physeq16, !(Samples %in% JAX))

noyeast <- c("MY17", "MY18", "MY19", "MY20", "MY21", "MY35", "MY36", "MY37", "MY38", "MY39")
physeq16_noyeast <- subset_samples(physeq16, (Samples %in% noyeast))

################### Theme Setup ###################
display.brewer.all(colorblindFriendly = TRUE)
pal <- brewer.pal(n = 12, name = "Paired")
pal1 <- brewer.pal(n = 7, name = "Set2")
pal2 <- brewer.pal(n = 7, name = "Dark2")

pal_dark1.25 <- pal %>% 
  adjust_luminance(-1.25)
pal_dark2 <- pal %>% 
  adjust_luminance(-2)

pal1_dark1.25 <- pal1 %>% 
  adjust_luminance(-1.25)
pal1_dark2 <- pal1 %>% 
  adjust_luminance(-2)

pal2_dark1.25 <- pal2 %>% 
  adjust_luminance(-1.25)
pal2_dark2 <- pal2 %>% 
  adjust_luminance(-2)

pal_med_1.25 <- pal %>% 
  adjust_luminance(1.25)
pal_light2 <- pal %>% 
  adjust_luminance(2)

pal1_med_1.25 <- pal1 %>% 
  adjust_luminance(1.25)
pal1_light2 <- pal1 %>% 
  adjust_luminance(2)

pal2_med_1.25 <- pal2 %>% 
  adjust_luminance(1.25)
pal2_light2 <- pal2 %>% 
  adjust_luminance(2)

palettemix1 <- c(pal_dark1.25[1:2], 
                 pal_med_1.25[1:2], 
                 pal_dark1.25[3:4], 
                 pal_med_1.25[3:4], 
                 pal_dark1.25[5:6], 
                 pal_med_1.25[5:6], 
                 pal_dark1.25[7:8], 
                 pal_med_1.25[7:8], 
                 pal_dark1.25[9:10], 
                 pal_med_1.25[9:10], 
                 pal_dark1.25[11:12], 
                 pal_med_1.25[11:12], 
                 pal_dark2[1:2], 
                 pal_light2[1:2], 
                 pal_dark2[3:4], 
                 pal_light2[3:4], 
                 pal_dark2[5:6], 
                 pal_light2[5:6], 
                 pal_dark2[7:8], 
                 pal_light2[7:8], 
                 pal_dark2[9:10], 
                 pal_light2[9:10], 
                 pal_dark2[11:12], 
                 pal_light2[11:12])

palettemix2 <- c(pal1_dark1.25[1:2], 
                 pal1_med_1.25[1:2], 
                 pal1_dark1.25[3:4], 
                 pal1_med_1.25[3:4], 
                 pal1_dark1.25[5:6], 
                 pal1_med_1.25[5:6], 
                 pal1_dark1.25[7], 
                 pal1_med_1.25[7], 
                 pal1_dark2[1:2], 
                 pal1_light2[1:2], 
                 pal1_dark2[3:4], 
                 pal1_light2[3:4], 
                 pal1_dark2[5:6], 
                 pal1_light2[5:6], 
                 pal1_dark2[7], 
                 pal1_light2[7])

palettemix3 <- c(pal2_dark1.25[1:2], 
                 pal2_med_1.25[1:2], 
                 pal2_dark1.25[3:4], 
                 pal2_med_1.25[3:4], 
                 pal2_dark1.25[5:6], 
                 pal2_med_1.25[5:6], 
                 pal2_dark1.25[7], 
                 pal2_light2[7],
                 pal2_dark2[1:2], 
                 pal2_light2[1:2], 
                 pal2_dark2[3:4], 
                 pal2_light2[3:4], 
                 pal2_dark2[5:6], 
                 pal2_light2[5:6], 
                 pal2_dark2[7], 
                 pal2_light2[7])

custom.pal <- c(palettemix1, palettemix2, palettemix3)

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

################### Stacked Relative Abundance ###################
gen.glom.16 <- tax_glom(physeq16_noyeast, taxrank = rank_names(physeq16_noyeast)[6], NArm = FALSE)
glom_rel_gen_noyeast16 <- psmelt(phyloseq::transform_sample_counts(gen.glom.16, function(x){x / sum(x)}))
glom_rel_gen_noyeast16$Genus[glom_rel_gen_noyeast16$Abundance<0.01] <- "Other <1%"

gen.legend.16S.noyeast <- glom_rel_gen_noyeast16 %>% 
  select(Genus) %>% 
  filter(Genus !=  "Other <1%") %>% 
  arrange(Genus)
glom_rel_gen_noyeast16_list <- as.list(unique(gen.legend.16S.noyeast$Genus))
italic.gen.legend.16S.noyeast <- mixedFontLabel(glom_rel_gen_noyeast16_list, italic = TRUE)

#Color setup 
write_xlsx(as.data.frame(glom_rel_gen_noyeast16_list), "gen_legend_16Snoyeast_palette.xlsx")

col13 <- c( "#1F78B4", 
            "#A6761D", 
            "#B2DF8A", 
            "#E6AB02", 
            "#FF7F00", 
            "#CAB2D6", 
            "#66A61E", 
            "#B15928", 
            "#66C2A5", 
            "#FC8D62", 
            "#8DA0CB", 
            "#FFD92F", 
            "#E5C494", 
            "#D95F02", 
            "#666666")  

relabund_16S_noyeast_Gen <- ggplot(glom_rel_gen_noyeast16, aes(x = Mouse, y = Abundance, fill = Genus)) +
  geom_bar(color = "black", stat = "identity", position = "stack", width = 0.95) +
  facet_wrap(~factor(Diet, levels = c("chow LM485", "chow 2914")), scales = "free", nrow = 1) + 
  labs(x = "", y = "Relative Abundance\n", title = "16S") + 
  rel.abund.theme + 
  scale_fill_manual(values = col13, labels = c(italic.gen.legend.16S.noyeast, "Other <1%")) +
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

write_xlsx(relabund.16S.noyeast, "Table S5 - Relative Abundance No Yeast 16S.xlsx")

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

