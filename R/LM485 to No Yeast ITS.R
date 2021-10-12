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
data31$Bug[data31$Bug == "Salmonella"] <- "STM"

physeq2 = phyloseq(otu_table(asvdataITS_matrix, taxa_are_rows = TRUE), tax_table(taxtreeITS_matrix), sample_data(data31))
random_tree2 = rtree(ntaxa(physeq2), rooted = TRUE, tip.label = taxa_names(physeq2))
physeq18 = merge_phyloseq(physeq2, random_tree2)

JAX <- c("MY12", "MY13", "MY14", "MY15", "MY16", "MY30", "MY31", "MY32", "MY33", "MY34", "MY63", "MY64")
physeq18 <- subset_samples(physeq18, !(Samples %in% JAX))

noyeast <- c("MY17", "MY18", "MY19", "MY20", "MY21", "MY35", "MY36", "MY37", "MY38", "MY39")
physeq18_noyeast <- subset_samples(physeq18, (Samples %in% noyeast))

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

################### Stacked Relative Abundance ###################
gen.glom.ITS <- tax_glom(physeq18_noyeast, taxrank = rank_names(physeq18_noyeast)[6], NArm = FALSE)
glom_rel_gen_noyeastITS<- psmelt(phyloseq::transform_sample_counts(gen.glom.ITS, function(x){x / sum(x)}))
glom_rel_gen_noyeastITS$Genus[glom_rel_gen_noyeastITS$Abundance < 0.01] <- "Other <1%"

gen.legend.ITS.noyeast <- glom_rel_gen_noyeastITS %>% 
  select(Genus) %>% 
  filter(Genus !=  "Other <1%") %>% 
  arrange(Genus)
gen.legend.ITS.noyeast.list <- as.list(unique(gen.legend.ITS.noyeast$Genus))
italic.gen.legend.ITS.noyeast <- mixedFontLabel(gen.legend.ITS.noyeast.list, italic = TRUE)

#Color setup 
write_xlsx(as.data.frame(gen.legend.ITS.noyeast.list), "gen_legend_ITSSnoyeast_palette.xlsx")

col53 <- c("#A43E00FF", 
           "#5C889CFF", 
           "#FFC29BFF", 
           "#3B547FFF", 
           "#983577FF", 
           "#FFC1F9FF", 
           "#003B80FF", 
           "#609100FF", 
           "#C4A000FF", 
           "#CBF2FFFF", 
           "#85C2FFFF", 
           "#FFF15BFF", 
           "#9E7F4AFF", 
           "#719D41FF", 
           "#FFD796FF", 
           "#004E30FF", 
           "#D0FEAAFF", 
           "#771100FF", 
           "#FFD1ABFF", 
           "#750055FF", 
           "#81E47EFF", 
           "#A94C4BFF", 
           "#DAECFFFF", 
           "#FFD2FFFF", 
           "#316200FF", 
           "#937200FF", 
           "#FFF763FF", 
           "#6F5100FF", 
           "#8E0000FF", 
           "#00552CFF", 
           "#FFCAC9FF", 
           "#FF7A7BFF", 
           "#B47B00FF", 
           "#A33000FF", 
           "#FFE296FF", 
           "#9B0048FF", 
           "#FFB86DFF", 
           "#FFDA66FF", 
           "#80678BFF", 
           "#5C3100FF", 
           "#003905FF", 
           "#5B00A6FF", 
           "#93F9D1FF", 
           "#F7DFFFFF", 
           "#FFBE8FFF", 
           "#B38FE4FF", 
           "#F2F28BFF", 
           "#100091FF", 
           "#8A0034FF", 
           "#DBD6FFFF", 
           "#FFA3EDFF", 
           "#003900FF", 
           "#666666")
relabund_ITS_noyeast_Gen <- ggplot(glom_rel_gen_noyeastITS, aes(x = Mouse, y = Abundance, fill = Genus)) + 
  geom_bar(aes(), color = "black", stat = "identity", position = "stack", width = 0.95) + 
  facet_wrap(~factor(Diet, levels = c("chow LM485", "chow 2914")), scales = "free", nrow = 1) + 
  labs(x = "", y = "Relative Abundance\n", title = "ITS") + 
  rel.abund.theme + 
  scale_fill_manual(values = col53, labels = c(italic.gen.legend.ITS.noyeast, "Other <1%")) +
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
