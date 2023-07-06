#install DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")

#packages used
library(tidyverse)
#install.packages("matrixStats")
library(matrixStats)
library(DESeq2)
library(phyloseq)


#phyloseq object
phyloseq_object_all <- readRDS("/home/lgschaer/old/Chesapeake/March_2023_Analysis/phyloseq_out/ps_unrarefied_03212023.rds")
head(sample_data(phyloseq_object_all))

#filter out eukaryotes and mitochondria
#will also remove inocula samples and subset to only include the final transfer
justbacteria <- phyloseq_object_all %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "Mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  ) %>%
  subset_samples(location_sampled == "openwater") %>%
  subset_samples(location_description == "port")
justbacteria

#saveRDS file
saveRDS(justbacteria, file = "/home/lgschaer/old/Chesapeake/March_2023_Analysis/deseq_out/openwater_port_samples_unrarefied_ps.rds")
head(sample_data(justbacteria))

#add a pseudo count of one to the OTU table
justbacteria@otu_table <- as.matrix(justbacteria@otu_table)+1

#subset phyloseq object to make sample comparison categories, Media_Carbon = DCPET_BH or TPA_BH
sample_info <- as.data.frame(sample_data(justbacteria))
unique(sample_info$location)


## make sure there are no zero counts
any(taxa_sums(justbacteria)==0)

# Convert phyloseq object to deseq2 format
ds.ps <- phyloseq_to_deseq2(justbacteria, ~ location)

ds.ps$location<-relevel(ds.ps$location, "Baltimore")

# Run DESeq2 analysis (all taxa at once!)
dds_ps <- DESeq(ds.ps)

# Investigate results
resultsNames(dds_ps)

# Put DESeq output into data frames
res <- as.data.frame(results(dds_ps, contrast=c("location","Norfolk","Baltimore"))) %>% mutate(Comparison="Baltimore vs. Norfolk") %>% rownames_to_column(var = "taxon")
head(res)

# Add taxonomy to DESeq output
res_tax <- as.data.frame(tax_table(justbacteria)) %>% rownames_to_column(var = "taxon") %>% full_join(res)
head(res_tax)

##check dimensions
dim(res)
dim(res_tax)

# Join everything together
enriched <- res_tax %>% 
  filter(!is.na(padj)) %>%
  mutate(
    threshold = ifelse(padj <= 0.001 & abs(log2FoldChange) >= 2, "Enriched", "Not_Enriched"),
    Enriched_Genus = ifelse(threshold == "Enriched", as.character(Genus), "Not Enriched"),
    Enriched_Genus = ifelse(is.na(Enriched_Genus), "Unclassified Genus", Enriched_Genus)
  ) 
head(enriched)
#View(enriched)

# Are any ASVs enriched?
any(enriched$threshold == "Enriched")
sum(enriched$threshold == "Enriched")

# Save a csv of the results
write.csv(enriched,"/home/lgschaer/old/Chesapeake/March_2023_Analysis/deseq_out/enriched_w_taxonomy_03222023.csv")

enriched_counts <- enriched %>%
  mutate(Category = ifelse(log2FoldChange > 2, "Right", "not_sig"),
         Category = ifelse(log2FoldChange < -2 & Category == "not_sig", "Left", Category)) %>%
  group_by(Comparison, Category) %>%
  summarise(Count = n()) %>%
  filter(Category != "not_sig")
enriched_counts

#any(is.na(enriched_w_tax))
#View(enriched_counts)

# Re-Order data to organize legend
sort(unique(enriched$Enriched_Genus))
length(unique(enriched$Enriched_Genus))

#enriched_w_tax$Enriched_Genus <- factor(enriched_w_tax$Enriched_Genus, 
#                                       levels = c("Achromobacter",      "Aminobacter",        "Ancylobacter",       "Aquamicrobium",      "Bauldia",            "Bosea",
#                                                 "Brevundimonas",      "Bryobacter",         "Candidimonas",       "Chelatococcus",     
#                                                "Chitinophaga",       "Cryobacterium",      "Devosia",            "Hydrogenophaga",     "Hyphomonas",         "Legionella",         
#                                               "Luteimonas",         "Mesorhizobium",      "Microbacterium",     "Millisia",           "Orrella",            "Paramesorhizobium",
#                                              "Parapedobacter",     "Parapusillimonas",   "Parvibaculum",       "Pedomicrobium",      "Pelagibacterium",   
#                                             "Persicitalea",       "Planktosalinus",     "Pseudaminobacter",   "Pseudolabrys",       "Pseudomonas",        "Pseudoxanthomonas", 
#                                            "Pusillimonas",       "Rhodobacter",        "Rhodococcus",        "Shinella",           "SN8",                "Sphingobacterium",  
#                                           "Tepidimonas",        "Thermomonas",        "Tianweitania",       "Variovorax",         "Verticiella",        "Youhaiella",   
#                                          "Unclassified Genus", "Not Enriched"))
# Save color palette

cl <- colors(distinct = TRUE)
set.seed(15887) # to set random generator seed
colors11 <- sample(cl, length(unique(enriched$Enriched_Genus)))

#colors11 <- c(
# "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
#"magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
#  "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
# "magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","firebrick","lightblue",
#"magenta","lightpink","yellow","orange","cyan","blue","mediumseagreen","gray63","white"
#)

shapes <- c(
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24, 
  21, 22, 23, 24, 
  21, 22, 23, 24,
  21, 22, 23, 24,
  21, 22, 23, 24,
  22, 21, 21
)
length(shapes)

#colors11 <- c(
# "orangered",      "purple",        "green",           "cyan",          "orange",        "khaki4",             "mediumslateblue",
#"mediumpurple1",  "darkmagenta",   "darkgreen",       "wheat2",        "yellow",        "lawngreen",          "plum",  
#  "royalblue",      "magenta",       "mediumseagreen",  "palegoldenrod", "grey47",        "chocolate4",         "darkorange3",        
# "lightblue",      "firebrick",     "yellowgreen",     "turquoise3",    "purple4",       "blue",               "red",            
#"lightcyan",       "coral1",       "cyan",            "goldenrod",     "yellowgreen",   "turquoise3",    "purple4",       "blue",               "red",            
#  "lightcyan",       "coral1",       "cyan",            "goldenrod",     "black",         "white"   
#) 
#length(colors11)

#colors11 <- c(
# "palegoldenrod","palegoldenrod","palegoldenrod","palegoldenrod",
#"orange",  "orange",  "orange",  "orange",  
#  "firebrick","firebrick","firebrick","firebrick",
# "pink","lightpink","lightpink","lightpink",
#"purple","purple","purple","purple",
#  "darkblue","darkblue","darkblue","darkblue",
# "darkcyan", "darkcyan", "darkcyan", "darkcyan", 
#"lightblue","lightblue","lightblue","lightblue",
#  "olivedrab2", "olivedrab2", "olivedrab2", "olivedrab2", 
# "darkgreen", "darkgreen", "darkgreen", "darkgreen", 
#"grey77","grey77","white"
#)

#volcano plot:
ggplot(data=enriched, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(fill=Enriched_Genus), shape = 21, color = "black", size=6) +
  facet_grid(cols = vars(Comparison))+
  #scale_fill_manual(values=colors11) +
  #scale_shape_manual(values=shapes) +
  labs(x = "log2 fold change", 
       y = "-log10 p-value") +
  theme_classic(base_size = 14)+
  geom_hline(yintercept = -log10(0.001), colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -2, colour="#990000", linetype="dashed")+ 
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 25, face = "bold", angle = 0),
        legend.position = "bottom",
        title = element_text(size = 18))#+
#guides(fill = guide_legend(override.aes = list(shape = shapes)))

# another way to plot
library(ggh4x)

unique(enriched$Comparison)

enriched_filt <- enriched %>%
  filter(threshold == "Enriched") %>%
  #mutate(Genus2 = ifelse(is.na(Genus), "", Genus),
   #      Species2 = ifelse(is.na(Species), "", Species)) %>%
  #unite(Taxonomy, Genus2, Species2, sep = " ") %>%
  mutate(#Taxonomy = ifelse(Taxonomy == " ", "Unclassified Organism", Taxonomy),
         Where_Enriched = case_when(
           (Comparison == "Baltimore vs. Norfolk" & log2FoldChange > 2) ~ "Norfolk",
           (Comparison == "Baltimore vs. Norfolk" & log2FoldChange < -2) ~ "Baltimore"
         )) %>%
  group_by(Comparison, Enriched_Genus, Where_Enriched) #%>%
# summarise(Count = n())
head(enriched_filt)
#View(enriched_filt)

colors <- c("olivedrab", "lightgoldenrod")

max <- 0.001
min <- round(min(enriched_filt$padj), digits = 3)
#min <- 0.5

ggplot(data=enriched_filt, aes(y=Enriched_Genus, x=log2FoldChange)) +
  geom_point(aes(fill=padj, shape = Where_Enriched), color = "black", size = 6) +
  facet_nested(cols = vars(Comparison), rows = vars(Class), space = "free_y", scales = "free_y")+
  scale_fill_gradientn(colors = colors, breaks = seq(min, max, by = 0.0001),#)+#,
                       limits = c(min, max), labels = as.character(seq(min, max, by = 0.0001))) +
  scale_shape_manual(values=c(21, 22, 23)) +
  xlab("log2 fold change")+ 
  ylab("Significantly Enriched Genera") +
  theme_linedraw(base_size = 14)+
  theme(axis.text.x = element_text(size = 18, angle = 0, vjust = 1, hjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        legend.text = element_text(size = 18),
        #legend.title = element_blank(),
        strip.text.x = element_text(size = 25, face = "bold", angle = 0),
        strip.text.y = element_text(size = 18, face = "bold", angle = 0),
        legend.position = "bottom",
        title = element_text(size = 18)) +
  guides(fill = guide_colourbar(barwidth = 40, barheight = 1, title = "Adjusted P-value"))


ggplot(data=enriched_filt, aes(y=Enriched_Genus, x=log2FoldChange)) +
  facet_nested(cols = vars(Comparison), rows = vars(Class), space = "free_y", scales = "free_y")+
  geom_col(aes(fill = Where_Enriched, y=Enriched_Genus, x=log2FoldChange), color = "black", show.legend = FALSE)+
  scale_fill_manual(values = c("darkcyan", "chocolate1"))+
  xlab("log2 fold change")+ 
  ylab("Significantly Enriched Taxa") +
  theme_linedraw(base_size = 14)+
  theme(axis.text.x = element_text(size = 16, angle = 0, vjust = 1, hjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 1, vjust = 0.5),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        #legend.text = element_text(size = 28),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 20, face = "bold", angle = 0),
        strip.text.y = element_text(size = 20, face = "bold", angle = 0),
        legend.position = "bottom",
        title = element_text(size = 18))

# Getting a list of differentially abundant taxa

head(enriched_filt)

diffAbundant <- enriched_filt %>%
  ungroup() %>%
  filter(threshold == "Enriched") %>%
  select(taxon, Genus, Where_Enriched) %>%
  unique()
head(diffAbundant)

write.csv(diffAbundant,"/home/lgschaer/old/Chesapeake/March_2023_Analysis/deseq_out/list_of_diff_abundant_tax.csv")
