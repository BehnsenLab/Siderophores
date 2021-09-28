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

asvdataITS <- dataITS1 %>% 
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

data32 <- read.delim("samplelabels.txt")
row.names(data32) <- data32$Samples
data32$Diet[data32$Diet == "chow 1"] <- "chow LM485"
data32$Diet[data32$Diet == "50ppm iron"] <- "purified 50ppm iron"
data32$Bug[data32$Bug == "Salmonella"] <- "STM"

physeq2Fe = phyloseq(otu_table(asvdataITS_matrix, taxa_are_rows = TRUE), tax_table(taxtreeITS_matrix), sample_data(data32))
random_tree2Fe = rtree(ntaxa(physeq2Fe), rooted = TRUE, tip.label = taxa_names(physeq2Fe))
physeq18Fe = merge_phyloseq(physeq2Fe, random_tree2Fe)

JAX <- c("MY12", "MY13", "MY14", "MY15", "MY16", "MY30", "MY31", "MY32", "MY33", "MY34", "MY63", "MY64")
physeq18Fe <- subset_samples(physeq18Fe, !(Samples %in% JAX))

LM485_50ppm <- c("MY7", "MY8", "MY9", "MY10", "MY11","MY22", "MY23", "MY24", "MY25", "MY26")
physeq18_Fe <- subset_samples(physeq18Fe, (Samples %in% LM485_50ppm))

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
dietcol <- scale_fill_manual(values = c("#1F78B4", "#E31A1C"))
Dietgroup <- c("chow LM485", "chow LM485", "purified 50ppm iron", "purified 50ppm iron", "purified 50ppm iron", "purified 50ppm iron", "purified 50ppm iron", "chow LM485", "chow LM485", "chow LM485")
box_aes <- theme(panel.background = element_blank(), 
                 panel.border = element_rect(fill = NA, color = "black"), 
                 text = element_text(family = "Arial", size = 25), 
                 aspect.ratio = 1)



################### Pyramid Plots  ################### 
dummy <- data.frame(Mouse = c("UIC107", "UIC107"), AbundInv = c(-1,1), Diet = c("chow LM485", "purified 50ppm iron"))
dummy$Diet <- factor(dummy$Diet, levels = c("chow LM485", "purified 50ppm iron"))

gen.glom.ITSFe <- tax_glom(physeq18_Fe, taxrank = rank_names(physeq18_Fe)[6], NArm = FALSE)
glom_rel4 <- psmelt(phyloseq::transform_sample_counts(gen.glom.ITSFe, function(x){x / sum(x)}))
glom_rel4$Genus[glom_rel4$Abundance < 0.05] <- "Other <5%"
glom_rel5 <- glom_rel4 %>% mutate(AbundInv = ifelse(Diet == "chow LM485", Abundance*-1, Abundance))
glom_rel5$Diet <- factor(glom_rel5$Diet, levels = c("chow LM485", "purified 50ppm iron"))

ITSpyramid2 <- ggplot(glom_rel5, aes(x = Mouse, y = AbundInv, fill = Mouse)) + 
  ylab("Relative Abundance") +
  geom_bar(stat = "identity", position = "stack", width = 1) + 
  facet_grid(Genus~Diet, scales = "free_x",  space = "free_x")+ 
  coord_flip() + 
  pyramid.theme +
  geom_hline(yintercept = 0) + 
  scale_y_continuous(labels = abs, expand = c(0,0)) +
  geom_blank(data = dummy) 

png(filename = "ITS abund pyramid chart purified iron.png", width = 4800, height = 3600, units = "px", res = 300)
plot(ITSpyramid2)
dev.off()

###Table
Abund.table.purified.ITS <- glom_rel5 %>% 
  select(OTU, Abundance, Mouse, Diet, Kingdom, Phylum, Class, Order, Family, Genus) %>% 
  filter(Abundance > 0) %>% 
  arrange(Mouse, desc(Diet))

write_xlsx(Abund.table.purified.ITS, "Relative Abundance 50ppm ITS Pyramid.xlsx")

################### Alpha Diversity  ################### 
tab_18S_50ppm <- microbiome::alpha(physeq18_Fe, index = "all")
tab_18S_50ppm$Dietgroup <- Dietgroup
write_xlsx(tab_18S_50ppm, "Alpha Diversity 50ppm ITS.xlsx")

################### Beta Diversity  ################### 
sample_data(physeq18_Fe)$Diet <- factor(sample_data(physeq18_Fe)$Diet, levels = c("chow LM485", "purified 50ppm iron"))
DistBCITS = distance(physeq18_Fe, method = "bray")
ordBCITS = ordinate(physeq18_Fe, method="PCoA", distance = DistBCITS)
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
