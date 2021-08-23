library("ggplot2")
library("ggpubr")
library("RColorBrewer")
library("gt")
library("knitr")
library("stats")
library("dplyr")
library("plyr")
library("phyloseq")
library("microbiome")
library("ape")
library("ggh4x")
library("writexl")

################### Setup ###################
otudat1 <- read.delim("OTU_ITS_Reads.txt")
row.names(otudat1) <- otudat1$OTU
otudat1 <- otudat1 %>% select (-OTU)
otu18s_matrix1 = as.matrix(otudat1)

otuda2 <- read.delim("OTU_ITS_Tree.txt")
row.names(otuda2) <- otuda2$OTU
otuda2$Genus[otuda2$Genus == " "] <- "Unclassified"
otuda2 <- otuda2 %>% select (-OTU)
otuda2 <- select(otuda2, -Species)
otu18s_matrix2 = as.matrix(otuda2)

otudata3 <- read.delim("samplelabels.txt")
row.names(otudata3) <- otudata3$Samples
otudata3$Diet[otudata3$Diet == "chow 1"] <- "chow LM485"
otudata3$Diet[otudata3$Diet == "chow 2"] <- "chow 2914"
otudata3$Bug[otudata3$Bug == "Salmonella"] <- "STM"

physeq2 = phyloseq(otu_table(otu18s_matrix1, taxa_are_rows = TRUE), tax_table(otu18s_matrix2), sample_data(otudata3))
random_tree2 = rtree(ntaxa(physeq2), rooted = TRUE, tip.label = taxa_names(physeq2))
physeq18 = merge_phyloseq(physeq2, random_tree2)

JAX <- c("MY12", "MY13", "MY14", "MY15", "MY16", "MY30", "MY31", "MY32", "MY33", "MY34", "MY63", "MY64")
physeq18 <- subset_samples(physeq18, !(Samples %in% JAX))

noyeast <- c("MY17", "MY18", "MY19", "MY20", "MY21", "MY35", "MY36", "MY37", "MY38", "MY39")
physeq18_noyeast <- subset_samples(physeq18, (Samples %in% noyeast))

display.brewer.all(colorblindFriendly = TRUE)
pal <- brewer.pal(n = 12, name = "Paired")
pal_darker <- pal %>% adjust_luminance(-2)
pal_lighter <- pal %>% adjust_luminance(2)


################### Pyramid Plots ###################
pyramid.theme <- theme(axis.text.y = element_blank(), 
                       axis.ticks.y = element_blank(),
                       text = element_text(family = "Arial", size = 25), 
                       strip.text.y = element_text(angle = 0, face = "italic"), 
                       panel.spacing.x = unit(0, "lines"), 
                       strip.background = element_rect(color = "black"), 
                       panel.background = element_part_rect("tblr", color = "black", fill = "white"))

glom_rel <- psmelt(phyloseq::transform_sample_counts(gen.glom.ITS, function(x){x / sum(x)}))
glom_rel$Genus[glom_rel$Abundance < 0.05] <- "Other <5%"
glom_rel2 <- glom_rel %>% mutate(AbundInv = ifelse(Diet == "chow LM485", Abundance*-1, Abundance))
glom_rel2$Diet <- factor(glom_rel2$Diet, levels = c("chow LM485", "chow 2914"))

dummy16 <- data.frame(Mouse = c("UIC127", "UIC127"), AbundInv = c(-1,1), Diet = c("chow LM485", "chow 2914"))
dummy16$Diet <- factor(dummy16$Diet, levels = c("chow LM485", "chow 2914"))


ITSpyramid <- ggplot(glom_rel2, aes(x = Mouse, y = AbundInv, fill = Mouse)) + 
  ylab("Relative Abundance") +
  geom_bar(stat = "identity", position = "dodge", width = 1) + 
  facet_grid(Genus~Diet, scales = "free_x",  space = "free_x") + 
  coord_flip() + 
  pyramid.theme +
  geom_hline(yintercept = 0) + 
  scale_y_continuous(labels = abs, expand = c(0,0)) +
  geom_blank(data = dummy16)

png(filename = "ITS abund pyramid chart chow 2914.png", width = 3600, height = 2400, units = "px", res = 300)
plot(ITSpyramid)
dev.off()

Abund.table.noyeast.ITS <- glom_rel2 %>% select(OTU, Abundance, Mouse, Diet, Kingdom, Phylum, Class, Order, Family, Genus) %>% filter(Abundance > 0) %>% arrange(Mouse, desc(Diet))

write_xlsx(Abund.table.noyeast.ITS, "Relative Abundance No Yeast ITS Pyramid.xlsx")
################### Alpha Diversity ###################
tab_18S_noyeast <- microbiome::alpha(physeq18_noyeast, index = "all")
write.table(tab_18S_noyeast, file = "Alpha Diversity No Yeast ITS.csv", sep = ",")

################### Beta Diversity ###################
cleanbg <- theme(panel.background = element_blank(), 
                 panel.border = element_rect(fill = NA, color = "black"), 
                 text = element_text(family = "Arial", size = 25))

Dietgroup <- c("chow LM485", "chow LM485", "chow LM485", "chow LM485", "chow LM485", "chow 2914", "chow 2914", "chow 2914", "chow 2914", "chow 2914")

sample_data(physeq18_noyeast)$Diet <- Dietgroup
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
citation(package="ggpubr")
citation(package="RColorBrewer")
citation(package="gt")
citation(package="knitr")
citation(package="stats")
citation(package="dplyr")
citation(package="plyr")
citation(package="phyloseq")
citation(package="microbiome")
citation(package="ape")
citation(package="ggh4x")
citation(package="writexl")
