
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


STMinf_ITS <- c("MY1", "MY2", "MY3", "MY4", "MY5", "MY6", "MY48", "MY49", "MY50", "MY51", "MY52", "MY53")
physeq2STMinf_ITS <- subset_samples(physeq2ITS, (Samples %in% STMinf_ITS))

rel.abund.theme <- theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.75), 
        panel.background = element_blank(), 
        text = element_text(family = "Arial", size = 50), 
        legend.text.align = 0)

STM20 <- prune_taxa(taxa_sums(physeq2STMinf_ITS) >= 6650, physeq2STMinf_ITS)

glomSTM_ITS <- tax_glom(STM20, taxrank = "Genus")
glomSTM_relITS <- phyloseq::transform_sample_counts(glomSTM_ITS, function(x){x / sum(x)})
glomSTM_relITS<-psmelt(glomSTM_relITS)
glomSTM_relITS$Genus <- as.character(glomSTM_relITS$Genus)
glomSTM_relITS$Treatment[glomSTM_relITS$Treatment == "Untreated"] <- "Pre-Infection"
glomSTM_relITS$Treatment[glomSTM_relITS$Treatment == "Strep+STM"] <- "Post-Infection"

STMPre <- glomSTM_relITS %>% 
  filter(Treatment == "Pre-Infection")
STMPre$Genus[STMPre$Abundance<0.01] <- "Other <1%"

gen.legend.STMPre <- STMPre %>% 
select(Genus) %>% 
filter(Genus !=  "Other <1%") %>% 
arrange(Genus)
glom_rel_gen_STMPre_list <- as.list(unique(gen.legend.STMPre$Genus))
italic.gen.legend.STMPre <- mixedFontLabel(glom_rel_gen_STMPre_list, italic = TRUE)

write_xlsx(as.data.frame(glom_rel_gen_STMPre_list), "gen.legend.STMPre palette.xlsx")

col20 <- c("#A43E00FF", 
           "#FFC29BFF", 
           "#3B547FFF", 
           "#C8DAFFFF", 
           "#609100FF", 
           "#C4A000FF", 
           "#9E7F4AFF", 
           "#FFD1ABFF", 
           "#FFF763FF", 
           "#B47B00FF", 
           "#9B0048FF", 
           "#FFE876FF", 
           "#5C3100FF", 
           "#003905FF", 
           "#93F9D1FF", 
           "#FFBE8FFF", 
           "#100091FF", 
           "#DBD6FFFF", 
           "#FFA3EDFF", 
           "#BAFB8DFF", 
           "#666666")

relabundSTM <- ggplot(STMPre, aes(x = Mouse, y = Abundance, fill = Genus)) +
  geom_bar(aes(), color = "black", stat = "identity", position = "stack") + 
  labs(x = "", y = "Relative Abundance\n", title = "ITS") + 
  rel.abund.theme + 
  scale_fill_manual(values = col20, labels = c(italic.gen.legend.STMPre, "Other <1%")) + 
  scale_y_continuous(expand = c(0.005,0.005))

png(filename = "RelativeAbundanceITS_STM.png", width = 12000, height = 9600, units = "px", res = 300)
plot(relabundSTM)
dev.off()

relabund.STMPre.legend.fig <-get_legend(relabundSTM)
png(filename = "Relative Abundance_ITS_STMPre_Genera_legend.png", width = 12000, height = 9600, units = "px", res = 300)
plot(relabund.STMPre.legend.fig)
dev.off()
