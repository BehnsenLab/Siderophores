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

################### Setup ###################
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


JAX <- c("MY12", "MY13", "MY14", "MY15", "MY16", "MY30", "MY31", "MY32", "MY33", "MY34", "MY63", "MY64")
physeq16 <- subset_samples(physeq16, !(Samples %in% JAX))

LM485_50ppm <- c("MY7", "MY8", "MY9", "MY10", "MY11","MY22", "MY23", "MY24", "MY25", "MY26")
physeq16_Fe <- subset_samples(physeq16, (Samples %in% LM485_50ppm))

################### Theme Setup ###################
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
dietcol <- scale_fill_manual(values = c("#1F78B4", "#E31A1C"))
Dietgroup <- c("chow LM485", "chow LM485", "chow LM485", "chow LM485", "chow LM485", "purified 50ppm iron", "purified 50ppm iron", "purified 50ppm iron", "purified 50ppm iron", "purified 50ppm iron")
box_aes <- theme(panel.background = element_blank(), 
                 panel.border = element_rect(fill = NA, color = "black"), 
                 text = element_text(family = "Arial", size = 25), 
                 aspect.ratio = 1)

################### Pyramid Plots ###################
glom_rel_16 <- psmelt(phyloseq::transform_sample_counts(gen.glom.16Fe, function(x){x / sum(x)}))
glom_rel_16$Genus[glom_rel_16$Abundance < 0.05] <- "Other <5%"

glom_rel16 <- glom_rel_16 %>% mutate(AbundInv = ifelse(Diet == "chow LM485", Abundance*-1, Abundance))
glom_rel16$Diet <- factor(glom_rel_16$Diet, levels = c("chow LM485", "purified 50ppm iron"))

dummy <- data.frame(Mouse = c("UIC107", "UIC107"), AbundInv = c(-1,1), Diet = c("chow LM485", "purified 50ppm iron"))
dummy$Diet <- factor(dummy$Diet, levels = c("chow LM485", "purified 50ppm iron"))

pyramid_16S <- ggplot(glom_rel16, aes(x = Mouse, y = AbundInv, fill = Mouse)) + 
  ylab("Relative Abundance") +
  geom_bar(stat = "identity", position = "stack", width = 1) + 
  facet_grid(Genus~Diet, scales = "free_x",  space = "free_x")+ 
  coord_flip() + 
  pyramid.theme +
  geom_hline(yintercept = 0) + 
  scale_y_continuous(labels = abs, expand = c(0,0)) +
  geom_blank(data = dummy) 

png(filename = "16S abund pyramid chart purified iron.png", width = 3600, height = 2400, units = "px", res = 300)
plot(pyramid_16S)
dev.off()

Abund.table.purified.16S <- glom_rel16 %>% 
  select(OTU, Abundance, Mouse, Diet, Kingdom, Phylum, Class, Order, Family, Genus) %>% 
  filter(Abundance > 0) %>% 
  arrange(Mouse, desc(Diet))

write_xlsx(Abund.table.purified.16S, "Relative Abundance 50ppm 16S Pyramid.xlsx")

################### Alpha Diversity ###################
tab_16S_50ppm <- microbiome::alpha(physeq16_Fe, index = "all")
write_xlsx(tab_16S_50ppm, "Alpha Diversity 50ppm 16S.xlsx")

################### Beta Diversity ###################
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
