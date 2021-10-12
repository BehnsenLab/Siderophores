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
data32$Bug[data32$Bug == "Salmonella"] <- "STM"

physeqFe = phyloseq(otu_table(asvdata16S_matrix, taxa_are_rows = TRUE), tax_table(taxtree16S_matrix), sample_data(data32))
random_treeFe = rtree(ntaxa(physeqFe), rooted=TRUE, tip.label=taxa_names(physeqFe))
physeq16Fe = merge_phyloseq(physeqFe, random_treeFe)

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

dietcol <- scale_fill_manual(values = c("#1F78B4", "#E31A1C"))

Dietgroup <- c("chow LM485", "chow LM485", "purified 50ppm iron", "purified 50ppm iron", "purified 50ppm iron", "purified 50ppm iron", "purified 50ppm iron", "chow LM485", "chow LM485", "chow LM485")

box_aes <- theme(panel.background = element_blank(), 
                 panel.border = element_rect(fill = NA, color = "black"), 
                 text = element_text(family = "Arial", size = 50), 
                 aspect.ratio = 1)

JAX <- c("MY12", "MY13", "MY14", "MY15", "MY16", "MY30", "MY31", "MY32", "MY33", "MY34", "MY63", "MY64")
physeq16Fe <- subset_samples(physeq16Fe, !(Samples %in% JAX))

LM485_50ppm <- c("MY7", "MY8", "MY9", "MY10", "MY11","MY22", "MY23", "MY24", "MY25", "MY26")
physeq16_Fe <- subset_samples(physeq16Fe, (Samples %in% LM485_50ppm))

################### Stacked Relative Abundance  ################### 
gen.glom.16Fe <- tax_glom(physeq16_Fe, taxrank = rank_names(physeq16_Fe)[6], NArm = FALSE)
glom_rel_gen_16Fe <- psmelt(phyloseq::transform_sample_counts(gen.glom.16Fe, function(x){x / sum(x)}))
glom_rel_gen_16Fe$Genus[glom_rel_gen_16Fe$Abundance<0.01] <- "Other <1%"
gen.legend.16Fe <- glom_rel_gen_16Fe %>% 
  select(Genus) %>% 
  filter(Genus !=  "Other <1%") %>% 
  arrange(Genus)
gen.legend.16Fe.list <- as.list(unique(gen.legend.16Fe$Genus))
italic.gen.legend.16Fe <- mixedFontLabel(gen.legend.16Fe.list, italic = TRUE)

#Color palette
gen.legend.16Fe.list2 <- gen.legend.16Fe.list
names(gen.legend.16Fe.list2) <- as.list(c(pal, pal1, pal2[1:2]))
write_xlsx(as.data.frame(gen.legend.16Fe.list2), "gen.legend.16SFe palette.xlsx")

col22 <- c(pal, pal1, pal2[1:2], "#666666")

relabund_16S_50ppm_Gen <- ggplot(glom_rel_gen_16Fe, aes(x = Mouse, y = Abundance, fill = Genus)) + 
  geom_bar(aes(), color = "black", stat = "identity", position = "stack", width = 0.95) + 
  facet_wrap(~factor(glom_rel_gen_16Fe$Diet, levels = c("chow LM485", "purified 50ppm iron")), scales = "free", nrow = 1) + 
  labs(x = "", y = "Relative Abundance\n", title = "16S") + 
  rel.abund.theme + 
  scale_fill_manual(values = col22, labels = c(italic.gen.legend.16Fe, "Other <1%")) +
  scale_y_continuous(expand = c(0.005,0.005)) +
  theme(legend.position = "right") +
  guides(fill = guide_legend(ncol = 1))

png(filename = "Relative Abundance_16S_50ppm_Genera.png", width = 12000, height = 9600, units = "px", res = 300)
plot(relabund_16S_50ppm_Gen)
dev.off()

relabund.16S.50ppm.genera.legend.fig <-get_legend(relabund_16S_50ppm_Gen)
png(filename = "Relative Abundance_16S_50ppm_Genera legend.png", width = 12000, height = 9600, units = "px", res = 300)
plot(relabund.16S.50ppm.genera.legend.fig)
dev.off()

relabund.16S.purified <- glom_rel_gen_16Fe %>% 
  select(OTU, Abundance, Mouse, Diet, Kingdom, Phylum, Class, Order, Family, Genus) %>% 
  filter(Abundance > 0) %>% 
  arrange(Mouse, desc(Diet))

write_xlsx(relabund.16S.purified, "Table S2 - Relative Abundance 50ppm 16S.xlsx")

################### Alpha Diversity  ################### 
tab_16S_50ppm <- microbiome::alpha(physeq16_Fe, index = "all")
tab_16S_50ppm$Dietgroup <- Dietgroup
write_xlsx(tab_16S_50ppm, "Alpha Diversity 50ppm 16S.xlsx")

################### Beta Diversity  ################### 
sample_data(physeq16_Fe)$Diet <- factor(sample_data(physeq16_Fe)$Diet, levels = c("chow LM485", "purified 50ppm iron"))
DistBC = distance(physeq16_Fe, method = "bray")
ordBC = ordinate(physeq16_Fe, method="PCoA", distance = DistBC)
plot_scree(ordBC)
PCoA16S <- plot_ordination(physeq16_Fe, ordBC, color = "Diet", shape = "Diet", label = "Mouse") + 
  stat_ellipse() + 
  cleanbg + 
  scale_color_manual(values = c("#1F78B4", "#E31A1C")) + 
  ggtitle("16S")
png(filename = "PCoA 16S chow LM485 to 50ppm.png", width = 4800, height = 3600, units = "px", res = 300)
plot(PCoA16S)
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
