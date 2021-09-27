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
pal_darker <- pal %>% adjust_luminance(-2)
pal_lighter <- pal %>% adjust_luminance(2)
pal_med_dark <- pal %>% adjust_luminance(-1.5)

pyramid.theme <- theme(axis.text.y = element_blank(), 
                       axis.ticks.y = element_blank(),
                       text = element_text(family = "Arial", size = 25), 
                       strip.text.y = element_text(angle = 0, face = "italic"), 
                       panel.spacing.x = unit(0, "lines"), 
                       strip.background = element_rect(color = "black"), 
                       panel.background = element_part_rect("tblr", color = "black", fill = "white"))

cleanbg <- theme(panel.background = element_blank(), 
                 panel.border = element_rect(fill = NA, color = "black"), 
                 text = element_text(family = "Arial", size = 25))

dietcol <- scale_fill_manual(values = c("#1F78B4", "#33A02C"))

Dietgroup <- c("chow LM485", "chow LM485", "chow LM485", "chow LM485", "chow LM485", "chow 2914", "chow 2914", "chow 2914", "chow 2914", "chow 2914")

box_aes <- theme(panel.background = element_blank(), 
                 panel.border = element_rect(fill = NA, color = "black"), 
                 aspect.ratio = 1)

cellhosttheme <- theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.75), 
                       text = element_text(family = "Arial", size = 25))

################### Pyramid Plots ###################
dummy16 <- data.frame(Mouse = c("UIC127", "UIC127"), AbundInv = c(-1,1), Diet = c("chow LM485", "chow 2914"))
dummy16$Diet <- factor(dummy16$Diet, levels = c("chow LM485", "chow 2914"))

gen.glom.ITS <- tax_glom(physeq18_noyeast, taxrank = rank_names(physeq18_noyeast)[6], NArm = FALSE)
glom_rel <- psmelt(phyloseq::transform_sample_counts(gen.glom.ITS, function(x){x / sum(x)}))
glom_rel$Genus[glom_rel$Abundance < 0.05] <- "Other <5%"
glom_rel2 <- glom_rel %>% mutate(AbundInv = ifelse(Diet == "chow LM485", Abundance*-1, Abundance))
glom_rel2$Diet <- factor(glom_rel2$Diet, levels = c("chow LM485", "chow 2914"))

ITSpyramid <- ggplot(glom_rel2, aes(x = Mouse, y = AbundInv, fill = Mouse)) + 
  ylab("Relative Abundance") +
  geom_bar(stat = "identity", position = "stack", width = 1) + 
  facet_grid(Genus~Diet, scales = "free_x",  space = "free_x") + 
  coord_flip() + 
  pyramid.theme +
  geom_hline(yintercept = 0) + 
  scale_y_continuous(labels = abs, expand = c(0,0)) +
  geom_blank(data = dummy16)

png(filename = "ITS abund pyramid chart chow 2914.png", width = 3600, height = 2400, units = "px", res = 300)
plot(ITSpyramid)
dev.off()

###Table
Abund.table.noyeast.ITS <- glom_rel2 %>% 
  select(OTU, Abundance, Mouse, Diet, Kingdom, Phylum, Class, Order, Family, Genus) %>% 
  filter(Abundance > 0) %>% 
  arrange(Mouse, desc(Diet))
write_xlsx(Abund.table.noyeast.ITS, "Relative Abundance No Yeast ITS Pyramid.xlsx")

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

png(filename = "PCoA ITS chow LM485 No Yeast.png", width = 2400, height = 2400, units = "px", res = 300)
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
