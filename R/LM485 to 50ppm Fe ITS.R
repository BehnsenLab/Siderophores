
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
otudata3$Diet[otudata3$Diet == "50ppm iron"] <- "purified 50ppm iron"
otudata3$Bug[otudata3$Bug == "Salmonella"] <- "STM"


physeq2 = phyloseq(otu_table(otu18s_matrix1, taxa_are_rows = TRUE), tax_table(otu18s_matrix2), sample_data(otudata3))
random_tree2 = rtree(ntaxa(physeq2), rooted=TRUE, tip.label=taxa_names(physeq2))
physeq18 = merge_phyloseq(physeq2, random_tree2)

display.brewer.all(colorblindFriendly = TRUE)
pal <- brewer.pal(n = 12, name = "Paired")
pal_darker <- pal %>% adjust_luminance(-2)
pal_lighter <- pal %>% adjust_luminance(2)

JAX <- c("MY12", "MY13", "MY14", "MY15", "MY16", "MY30", "MY31", "MY32", "MY33", "MY34", "MY63", "MY64")
physeq18 <- subset_samples(physeq18, !(Samples %in% JAX))

LM485_50ppm <- c("MY7", "MY8", "MY9", "MY10", "MY11","MY22", "MY23", "MY24", "MY25", "MY26")
physeq18_Fe <- subset_samples(physeq18, (Samples %in% LM485_50ppm))

################### Pyramid Plots ###################
pyramid.theme <- theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        text = element_text(family = "Arial", size = 25), 
        strip.text.y = element_text(angle = 0, face = "italic"), 
        panel.spacing.x = unit(0, "lines"), 
        strip.background = element_rect(color = "black"), 
        panel.background = element_part_rect("tblr", color = "black", fill = "white"))

glom_rel4 <- psmelt(phyloseq::transform_sample_counts(gen.glom.ITSFe, function(x){x / sum(x)}))
glom_rel4$Genus[glom_rel4$Abundance < 0.05] <- "Other <5%"

glom_rel5 <- glom_rel4 %>% mutate(AbundInv = ifelse(Diet == "chow LM485", Abundance*-1, Abundance))
glom_rel5$Diet <- factor(glom_rel5$Diet, levels = c("chow LM485", "purified 50ppm iron"))

dummy <- data.frame(Mouse = c("UIC107", "UIC107"), AbundInv = c(-1,1), Diet = c("chow LM485", "purified 50ppm iron"))
dummy$Diet <- factor(dummy$Diet, levels = c("chow LM485", "purified 50ppm iron"))

ITSpyramid2 <- ggplot(glom_rel5, aes(x = Mouse, y = AbundInv, fill = Mouse)) + 
  ylab("Relative Abundance") +
  geom_bar(stat = "identity", position = "dodge", width = 1) + 
  facet_grid(Genus~Diet, scales = "free_x",  space = "free_x")+ 
  coord_flip() + 
  pyramid.theme +
  geom_hline(yintercept = 0) + 
  scale_y_continuous(labels = abs, expand = c(0,0)) +
  geom_blank(data = dummy) 

png(filename = "ITS abund pyramid chart purified iron.png", width = 3600, height = 2400, units = "px", res = 300)
plot(ITSpyramid2)
dev.off()

Abund.table.purified.ITS <- glom_rel5 %>% select(OTU, Abundance, Mouse, Diet, Kingdom, Phylum, Class, Order, Family, Genus) %>% filter(Abundance > 0) %>% arrange(Mouse, desc(Diet))

write_xlsx(Abund.table.purified.ITS, "Relative Abundance 50ppm ITS Pyramid.xlsx")

################### Alpha Diversity ###################
tab_18S_50ppm <- microbiome::alpha(physeq18_Fe, index = "all")
write.table(tab_18S_50ppm, file = "Alpha Diversity 50ppm ITS.csv", sep = ",")

################### Beta Diversity ###################
cleanbg <- theme(panel.background = element_blank(), 
                 panel.border = element_rect(fill = NA, color = "black"), 
                 text = element_text(family = "Arial", size = 25))
Dietgroup <- c("chow LM485", "chow LM485", "chow LM485", "chow LM485", "chow LM485", "purified 50ppm iron", "purified 50ppm iron", "purified 50ppm iron", "purified 50ppm iron", "purified 50ppm iron")

sample_data(physeq18_Fe)$Diet <- Dietgroup
sample_data(physeq18_Fe)$Diet <- factor(sample_data(physeq18_Fe)$Diet, levels = c("chow LM485", "purified 50ppm iron"))
DistBCITS = distance(physeq18_Fe, method = "bray")
ordBCITS = ordinate(physeq18_Fe, method="PCoA", distance = DistBCITS)
plot_scree(ordBCITS)
PCoAITS <- plot_ordination(physeq18_Fe, ordBCITS, color = "Diet", shape = "Diet", label = "Mouse") + 
  stat_ellipse() + 
  cleanbg + 
  scale_color_manual(values = c("#1F78B4", "#E31A1C")) + 
  ggtitle("ITS")
png(filename = "PCoA ITS chow LM485 to 50ppm.png", width = 2400, height = 2400, units = "px", res = 300)
plot(PCoAITS)
dev.off()

################### References ###################
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
citation(package="phyloseq")
citation(package="microbiome")
citation(package="ape")
citation(package="ggh4x")
citation(package="writexl")

