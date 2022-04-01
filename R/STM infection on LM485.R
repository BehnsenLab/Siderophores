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
library("data.table")
library("ggh4x")
library("tidyr")
library("stringr")
library("cowplot")
library("tidyverse")
library("writexl")

data <- read.delim("SampleData_ITS.txt")
colnames(data) <- gsub("ITS.Behnsen.MY.", "MY", colnames(data))
colnames(data) <- gsub("ITS.Negative.Control.sample2", "Negative Control 1", colnames(data))
colnames(data) <- gsub("ITS.Negative.Ctrl.2.Underhill", "Negative Control 2", colnames(data))
colnames(data) <- gsub("ITS.Positive.Ctrl.2.Underhill", "Positive Control 2", colnames(data))
data <- filter(data, Taxonomy != "k__Fungi ; NA ; NA ; NA ; NA ; NA ; NA")

asvdata <- data %>%
  select(!c(Taxonomy))

row.names(asvdata) <- asvdata$ID
asvdata <- asvdata %>% select(-ID)
asvdata_matrix = as.matrix(asvdata)

taxtree <- data %>%
  select(ID, Taxonomy)
taxtree = separate(data = taxtree, col = Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\;" )

taxtree <- taxtree %>%
  mutate_at("Kingdom", str_replace_all, "k__", "") %>%
  mutate_at("Phylum", str_replace_all, "p__", "") %>%
  mutate_at("Class", str_replace_all, "c__", "") %>%
  mutate_at("Order", str_replace_all, "o__", "") %>%
  mutate_at("Family", str_replace_all, "f__", "") %>%
  mutate_at("Genus", str_replace_all, "g__", "") %>%
  mutate_at("Species", str_replace_all, "s__", "")
row.names(taxtree) <- taxtree$ID
taxtree <- taxtree %>% select(-ID)
taxtree_matrix = as.matrix(taxtree)

data3 <- read.delim("samplelabels.txt")
row.names(data3) <- data3$Samples
data3$Diet[data3$Diet == "chow 1"] <- "chow LM485"
data3$Diet[data3$Diet == "chow 2"] <- "chow 2914"
data3$Diet[data3$Diet == "50ppm iron"] <- "purified 50ppm iron"
data3$Bug[data3$Bug == "Salmonella"] <- "STM"

physeqITS = phyloseq(otu_table(asvdata_matrix, taxa_are_rows = TRUE), tax_table(taxtree_matrix), sample_data(data3))
random_tree = rtree(ntaxa(physeqITS), rooted=TRUE, tip.label=taxa_names(physeqITS))
physeq2ITS = merge_phyloseq(physeqITS, random_tree)

data3 <- read.delim("samplelabels.txt")
row.names(data3) <- data3$Samples
data3$Diet[data3$Diet == "chow 1"] <- "chow LM485"
data3$Diet[data3$Diet == "chow 2"] <- "chow 2914"
data3$Diet[data3$Diet == "50ppm iron"] <- "purified 50ppm iron"
data3$Bug[data3$Bug == "Salmonella"] <- "STM"

physeqITS = phyloseq(otu_table(asvdata_matrix, taxa_are_rows = TRUE), tax_table(taxtree_matrix), sample_data(data3))
random_tree = rtree(ntaxa(physeqITS), rooted=TRUE, tip.label=taxa_names(physeqITS))
physeq2ITS = merge_phyloseq(physeqITS, random_tree)

STMinf_ITS <- c("MY1", "MY2", "MY3", "MY4", "MY5", "MY6", "MY48", "MY49", "MY50", "MY51", "MY52", "MY53")
physeq2STMinf_ITS <- subset_samples(physeq2ITS, (Samples %in% STMinf_ITS))

rel.abund.theme <- theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.75),
        panel.background = element_blank(),
        text = element_text(family = "Arial", size = 50),
        legend.text.align = 0)

STM20 <- prune_taxa(taxa_sums(physeq2STMinf_ITS) >= 4750, physeq2STMinf_ITS)

glomSTM_ITS <- tax_glom(STM20, taxrank = "Genus")
glomSTM_relITS <- phyloseq::transform_sample_counts(glomSTM_ITS, function(x){x / sum(x)})
glomSTM_relITS<-psmelt(glomSTM_relITS)
glomSTM_relITS$Genus <- as.character(glomSTM_relITS$Genus)
glomSTM_relITS$Treatment[glomSTM_relITS$Treatment == "Untreated"] <- "Pre-Infection"
glomSTM_relITS$Treatment[glomSTM_relITS$Treatment == "Strep+STM"] <- "Post-Infection"

STMPre <- glomSTM_relITS %>%
  filter(Treatment == "Pre-Infection")
STMPre$Genus[STMPre$Abundance<0.02] <- "Other <2%"

color.pal.ITS <- read.delim("Genus color palette ITS.txt", sep = "")

gen.legend.STMPre <- STMPre %>%
select(Genus) %>%
filter(Genus !=  "Other <2%") %>%
arrange(Genus)
glom_rel_gen_STMPre_list <- as.list(unique(gen.legend.STMPre$Genus))
italic.gen.legend.STMPre <- mixedFontLabel(glom_rel_gen_STMPre_list, italic = TRUE)

gen.legend.ITS <- STMPre %>%
  select(Genus) %>%
  filter(Genus !=  "Other <2%") %>%
  arrange(Genus)
gen.legend.ITS.list <- as.list(unique(gen.legend.ITS$Genus))
italic.gen.legend.ITS <- mixedFontLabel(gen.legend.ITS.list, italic = TRUE)

gen.unique <- unique(STMPre$Genus)
gen.unique.df <- as.data.frame(gen.unique)
colnames(gen.unique.df) <- "Genus"
gen.unique.df <- gen.unique.df %>%
  arrange(Genus)

  joined.col <- dplyr::inner_join(gen.unique.df, color.pal.ITS, by = "Genus")
  col.chow <- joined.col$pal.generator %>%
    as.character(joined.col$pal.generator)
  col.chow <- dput(col.chow)

  relabund_ITS_STM_Gen <- ggplot(STMPre, aes(x = Mouse, y = Abundance, fill = Genus)) +
    geom_bar(aes(), color = "black", stat = "identity", position = "stack", width = 0.95) +
    facet_wrap(~factor(Diet, levels = c("chow LM485", "purified 50ppm iron", "chow 2914")), scales = "free", nrow = 1) +
    labs(x = "", y = "Relative Abundance\n", title = "ITS") +
    rel.abund.theme +
    scale_fill_manual(values = col.chow, labels = c(italic.gen.legend.ITS, "Other <2%")) +
    scale_y_continuous(expand = c(0.005,0.005)) +
    theme(legend.position = "right") +
    guides(fill = guide_legend(ncol = 3))

png(filename = "RelativeAbundanceITS_STM.png", width = 12000, height = 9600, units = "px", res = 300)
plot(relabundSTM)
dev.off()

relabund.STMPre.legend.fig <-get_legend(relabundSTM)
png(filename = "Relative Abundance_ITS_STMPre_Genera_legend.png", width = 12000, height = 9600, units = "px", res = 300)
plot(relabund.STMPre.legend.fig)
dev.off()

relabund.ITS.STMPre <- STMPre %>%
  select(OTU, Abundance, Mouse, Diet, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  filter(Abundance > 0) %>%
  arrange(Mouse, desc(Diet))

write_xlsx(relabund.ITS.STMPre, "Table S1 - Relative Abundance ITS uninfected SPF mice.xlsx")
