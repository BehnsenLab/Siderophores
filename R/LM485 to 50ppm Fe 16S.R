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
otudata <- read.delim("Behnsen_OTU 16S_v2 samples.txt")
row.names(otudata) <- otudata$OTU
otudata <- otudata %>% select (-OTU) 
otu16s_matrix1 = as.matrix(otudata)

otudata2 <- read.delim("Behnsen_OTU 16S_v2 tree.txt")
row.names(otudata2) <- otudata2$OTU
otudata2$Genus[otudata2$Genus == " "]
otudata2 <- otudata2 %>% select (-OTU)
otudata2 <- select(otudata2, -Species)
otu16s_matrix2 = as.matrix(otudata2)

otudata3 <- read.delim("samplelabels.txt")
row.names(otudata3) <- otudata3$Samples
otudata3$Diet[otudata3$Diet == "chow 1"] <- "chow LM485"
otudata3$Diet[otudata3$Diet == "50ppm iron"] <- "purified 50ppm iron"
otudata3$Bug[otudata3$Bug == "Salmonella"] <- "STM"

physeq = phyloseq(otu_table(otu16s_matrix1, taxa_are_rows = TRUE), tax_table(otu16s_matrix2), sample_data(otudata3))
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
physeq16 = merge_phyloseq(physeq, random_tree)

display.brewer.all(colorblindFriendly = TRUE)
pal <- brewer.pal(n = 12, name = "Paired")
pal_darker <- pal %>% adjust_luminance(-2)
pal_lighter <- pal %>% adjust_luminance(2)

JAX <- c("MY12", "MY13", "MY14", "MY15", "MY16", "MY30", "MY31", "MY32", "MY33", "MY34", "MY63", "MY64")
physeq16 <- subset_samples(physeq16, !(Samples %in% JAX))

LM485_50ppm <- c("MY7", "MY8", "MY9", "MY10", "MY11","MY22", "MY23", "MY24", "MY25", "MY26")
physeq16_Fe <- subset_samples(physeq16, (Samples %in% LM485_50ppm))

################### Pyramid Plots ###################
pyramid.theme <- theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        text = element_text(family = "Arial", size = 25), 
        strip.text.y = element_text(angle = 0, face = "italic"), 
        panel.spacing.x = unit(0, "lines"), 
        strip.background = element_rect(color = "black"), 
        panel.background = element_part_rect("tblr", color = "black", fill = "white"))

glom_rel_16 <- psmelt(phyloseq::transform_sample_counts(gen.glom.16Fe, function(x){x / sum(x)}))
glom_rel_16$Genus[glom_rel_16$Abundance < 0.05] <- "Other <5%"

glom_rel16 <- glom_rel_16 %>% mutate(AbundInv = ifelse(Diet == "chow LM485", Abundance*-1, Abundance))
glom_rel16$Diet <- factor(glom_rel_16$Diet, levels = c("chow LM485", "purified 50ppm iron"))

dummy <- data.frame(Mouse = c("UIC107", "UIC107"), AbundInv = c(-1,1), Diet = c("chow LM485", "purified 50ppm iron"))
dummy$Diet <- factor(dummy$Diet, levels = c("chow LM485", "purified 50ppm iron"))

pyramid_16S <- ggplot(glom_rel16, aes(x = Mouse, y = AbundInv, fill = Mouse)) + 
  ylab("Relative Abundance") +
  geom_bar(stat = "identity", position = "dodge", width = 1) + 
  facet_grid(Genus~Diet, scales = "free_x",  space = "free_x")+ 
  coord_flip() + 
  pyramid.theme +
  geom_hline(yintercept = 0) + 
  scale_y_continuous(labels = abs, expand = c(0,0)) +
  geom_blank(data = dummy) 

png(filename = "16S abund pyramid chart purified iron.png", width = 3600, height = 2400, units = "px", res = 300)
plot(pyramid_16S)
dev.off()

Abund.table.purified.16S <- glom_rel16 %>% select(OTU, Abundance, Mouse, Diet, Kingdom, Phylum, Class, Order, Family, Genus) %>% filter(Abundance > 0) %>% arrange(Mouse, desc(Diet))
write_xlsx(Abund.table.purified.16S, "Relative Abundance 50ppm 16S Pyramid.xlsx")

################### Alpha Diversity ###################
tab_16S_50ppm <- microbiome::alpha(physeq16_Fe, index = "all")
write.table(tab_16S_50ppm, file = "Alpha Diversity 50ppm 16S.csv", sep = ",")

################### Beta Diversity ###################
cleanbg <- theme(panel.background = element_blank(), 
                 panel.border = element_rect(fill = NA, color = "black"), 
                 text = element_text(family = "Arial", size = 25))
Dietgroup <- c("chow LM485", "chow LM485", "chow LM485", "chow LM485", "chow LM485", "purified 50ppm iron", "purified 50ppm iron", "purified 50ppm iron", "purified 50ppm iron", "purified 50ppm iron")

sample_data(physeq16_Fe)$Diet <- Dietgroup
sample_data(physeq16_Fe)$Diet <- factor(sample_data(physeq16_Fe)$Diet, levels = c("chow LM485", "purified 50ppm iron"))
DistBC = distance(physeq16_Fe, method = "bray")
ordBC = ordinate(physeq16_Fe, method="PCoA", distance = DistBC)
plot_scree(ordBC)
PCoA16S <- plot_ordination(physeq16_Fe, ordBC, color = "Diet", shape = "Diet", label = "Mouse") + 
  stat_ellipse() + 
  cleanbg + 
  scale_color_manual(values = c("#1F78B4", "#E31A1C")) + 
  ggtitle("16S")
png(filename = "PCoA 16S chow LM485 to 50ppm.png", width = 2400, height = 2400, units = "px", res = 300)
plot(PCoA16S)
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

