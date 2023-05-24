
###PHYLOSEQ ANALYSIS

library(phyloseq)
library(tidyverse)
library(csv)


##### I USED THE FOLLOWING CODE TO GENERATE A CORRECT METADATA TABLE
##### Skip the code with hash marks in front of it and skip to loading the corrected file

#sample_names <- sample_names %>%
 # mutate(old_fastq_name = sample.names)
#head(sample_names)
#dim(sample_names)

#sample_key <- as.csv("/home/lgschaer/old/Chesapeake/fastq_naming_information.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE) %>%
#  rownames_to_column(var = "full_old_fq") %>%
#  separate(full_old_fq, into = c("old_fastq_name", "file_info"), sep = "_S") %>%
#  separate(file_info, into = c(NA, NA, "read", NA), sep = "_") %>%
#  filter(read != "R2") %>%
#  mutate(is.duplicated = ifelse(duplicated(new_fastq_name), "duplicate", "not_duplicate")) %>%
#  filter(is.duplicated != "duplicate")
#head(sample_key)
#dim(sample_key)

#duplicated(sample_key$new_fastq_name)

#write_csv(sample_names, "/home/lgschaer/old/Chesapeake/March_2023_Analysis/phyloseq_out/sample_names_03212023.csv")
#sample_names <- as.csv("/home/lgschaer/old/Chesapeake/March_2023_Analysis/phyloseq_out/sample_names_03212023.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
  

#geo_summary <- as.csv("/home/lgschaer/old/Chesapeake/March_2023_Analysis/geo_summary_LGS.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE) %>%
 # mutate(
  #  latitude = as.character(latitude),
   # longitude = as.character(longitude)
  #)
#head(geo_summary)
#dim(geo_summary)

#sdata <- as.csv("/home/lgschaer/old/Chesapeake/cb_metadata_rg_LGS.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE) %>%
 # select(-geographic_order) %>%
  #mutate(location_sampled_w_boat_name = ifelse(location_sampled == "openwater", location_sampled, boat_name),
   #      latitude = as.character(latitude),
    #     longitude = as.character(longitude),
     #    rownames = SampleID) %>%
#  left_join(geo_summary, by = c("latitude", "longitude")) %>%
 # column_to_rownames(var = "SampleID")
#head(sdata)
#dim(sdata)
#View(sdata)

#write_csv(sdata, "/home/lgschaer/old/Chesapeake/March_2023_Analysis/phyloseq_out/full_cleaned_metadata.csv")

#geo_summary <- select(sdata, latitude, longitude, geographic_order) %>%
 # filter(longitude != "blank") %>%
  #unique()
#head(geo_summary)
#dim(geo_summary)

#write_csv(geo_summary, "/home/lgschaer/old/Chesapeake/March_2023_Analysis/geo_summary.csv")


#sdataFilt <- sdata %>%
 # filter(is.na(new_fastq_name))
#head(sdataFilt)
#View(sdataFilt)

#write_csv(sdataFilt, "/home/lgschaer/old/Chesapeake/March_2023_Analysis/phyloseq_out/extra_metadata_rows.csv")

#dim(sample_key)
#dim(sdata)

#sdata2 <- sdata %>%
 # column_to_rownames(var = "sample.names")
#head(sdata2)

##### START HERE #####

#load metadata
sdata <- as.csv("/home/lgschaer/old/Chesapeake/March_2023_Analysis/phyloseq_out/full_cleaned_metadata_edited04112023.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(sdata)

#load sequence table
seqtab <- readRDS("/home/lgschaer/old/Chesapeake/March_2023_Analysis/dada2_out/seqtab.rds")
colnames(seqtab) <- NULL
dim(seqtab)

#load taxa table
taxa <- readRDS("/home/lgschaer/old/Chesapeake/March_2023_Analysis/dada2_out/taxtab.rds")   
rownames(taxa) <- NULL
taxa[2:5,2:5]


### Section 2: making a phyloseq object

#make the first phyloseq object
samdata = sample_data(sdata)
seqtab = otu_table(seqtab, taxa_are_rows = FALSE)
taxtab = tax_table(taxa)

#sample_names(seqtab)
#sample_names(samdata)

#taxa_names(seqtab)
#taxa_names(taxtab)

#combine all components into a phyloseq object
ps = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata))
ps


write_rds(ps, "/home/lgschaer/old/Chesapeake/March_2023_Analysis/phyloseq_out/ps_unrarefied_03212023.rds")


#Filter out eukaryotes and mitochondria
ps_filt <- ps %>%
  subset_samples(sample_description != "blank") %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "Mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  )
ps_filt

### Section 3: Normalize the data

#hist(sample_sums(ps_filt))

samplesover1000_all <- subset_samples(ps_filt, sample_sums(ps_filt) > 2000)
any(taxa_sums(samplesover1000_all) == 0)
sum(taxa_sums(samplesover1000_all) == 0)
prune_samplesover1000_all <- prune_taxa(taxa_sums(samplesover1000_all) > 0, samplesover1000_all)
any(taxa_sums(prune_samplesover1000_all) == 0)

#for reproducible data
set.seed(81)

#hist(sample_sums(ps_filt))

rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))
pps<- rarefy_samplesover1000_all 
pps
min(sample_sums(prune_samplesover1000_all))

#mdata <- sample_data(pps) %>%
 # group_by(sample_type, transit, location, boat, filter_type, swab_type) %>%
  #summarise(
   # n = n()
#  )
#head(mdata)

#write_csv(mdata, "/home/lgschaer/old/Chesapeake/March_2023_Analysis/phyloseq_out/n_by_category.csv")

###Section 4: plotting

head(sample_data(pps))

colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(sdata$sample_detail_description)))
sample_colors <- c(color_list)

#violin plot
pps %>%                                                              #phyloseq object
  plot_richness(
    x = "sample_detail_description",                                                 #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_boxplot(aes(x = sample_detail_description, fill = sample_detail_description), show.legend = FALSE)+           #make violin plot, set fill aes to sampletype
  theme_linedraw()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, vjust = 0.5, hjust = 1, angle = 90),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = sample_colors)+ 
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position


# Is there a difference in alpha diversity between the ports?

portsOnly <- subset_samples(pps, location_sampled == "openwater" & location_description == "port")

portMeta <- as.data.frame(sample_data(portsOnly))
head(portMeta)

colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(portMeta$location)))
sample_colors <- c(color_list)


sample_colors <- c("darkcyan", "chocolate2")


#violin plot
portsOnly %>%                                                              #phyloseq object
  plot_richness(
    x = "location",                                                 #compare diversity of datatype
    measures = c("Observed", "Shannon")) +                           #choose diversity measures
  geom_violin(aes(x = location, fill = location), show.legend = FALSE)+           #make violin plot, set fill aes to sampletype
  geom_jitter(aes(x = location, shape = sample_description), fill = "white", size = 4, width = 0.15, height = 0) +
  theme_linedraw()+                                                   #change theme to classic
  xlab(NULL)+                                                        #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                  #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_blank())+                     #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+        #adjust headings
  scale_fill_manual(values = sample_colors)+ 
  scale_shape_manual(values = c(22, 21))+
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5)) #change title size, face and position

# Is there a statistical difference in alpha diversity between ports?

#sdata2 <- sdata %>% 
 # rownames_to_column(var = "SampleID")
#head(sdata2)

alphadiv <- estimate_richness(portsOnly, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sdata, by = "SampleID") %>%
  filter(location_description == "port" & (location == "Norfolk" | location == "Baltimore") & location_sampled == "openwater")
head(alphadiv)


alphadiv2 <- alphadiv %>%
  group_by(location) %>%
  summarise(
    maxO = max(Observed), 
    meanO = mean(Observed),
    minO = min(Observed),
    maxS = max(Shannon), 
    meanS = mean(Shannon),
    minS = min(Shannon)
  )
View(alphadiv2)


#Kruskal-Wallis Test
library(FSA)
set.seed(81)

##NORFOLK VS BALTIMORE, openwater only
#Observed
kruskal.test(Observed ~ location, data = alphadiv) 

#Shannon
kruskal.test(Shannon ~ location, data = alphadiv)


alphadiv <- write_csv(alphadiv, "/home/lgschaer/old/Chesapeake/March_2023_Analysis/phyloseq_out/alpha_diversity.csv")


#Making a PCoA plot
### UNIFRAC
library(ape)

#adding a phylogenetic tree to phyloseq object using ape library
random_tree = rtree(ntaxa(pps), rooted=TRUE, tip.label=taxa_names(pps))
plot(random_tree)

justbacteria3 = merge_phyloseq(pps, samdata, random_tree)
justbacteria3

#ordination
uni_distance <- ordinate(
  physeq = justbacteria3, 
  method = "PCoA", 
  distance = "unifrac"
)

#plot

colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(sdata$sample_detail_description)))
sample_colors <- c(color_list)
sample_colors <- c("#C74718", "#FFC957", "#A9BA6B", "#1E8B6C", "#94CCD8", "#00008B")
#sample_colors <- c("#1E8B6C", "#A9BA6B", "#FFC957", "#C74718", "#F4A0AA", "#551A8B")

plot_ordination(
  physeq = justbacteria3,                                                          #phyloseq object
  ordination = uni_distance)+                                                #ordination
  geom_point(aes(fill = sample_detail_description, shape = boat_name), size = 6) +                         #sets fill color to sampletype
  scale_fill_manual(values = sample_colors) +
  scale_shape_manual(values = c(22,21, 23))+
  theme_linedraw() +                                                      #changes theme, removes grey background
  theme(                             
    legend.title = element_blank(),                                      #removes legend title
    legend.position = "right",
    legend.text = element_text(size = 20, face = "bold"),                                 
    axis.text.y.left = element_text(size = 10),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 20),
    strip.text = element_text(face = "bold", size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = c(21))))          #fills legend points based on the fill command


# Is there a difference in beta diversity between the boats and the open water?
#View(sdata)
head(sdata)

boatPort <- subset_samples(pps, !is.na(location_sampled) & location_sampled != "blank" & !is.na(boat_name))
boatPort

#adding a phylogenetic tree to phyloseq object using ape library
random_tree = rtree(ntaxa(boatPort), rooted=TRUE, tip.label=taxa_names(portsOnly))
plot(random_tree)

boatPort2 = merge_phyloseq(boatPort, samdata, random_tree)
boatPort2

#ordination
uni_distance2 <- ordinate(
  physeq = boatPort2, 
  method = "PCoA", 
  distance = "unifrac"
)

#plot
head(sample_data(boatPort2))

sample_colors <- c("sienna3","darkseagreen", "aliceblue")

plot_ordination(
  physeq = boatPort2,                                                          #phyloseq object
  ordination = uni_distance2)+                                                 #ordination
  geom_point(aes(fill = location_sampled_w_boat_name, shape = sample_type), size = 6) +     #sets fill color to sampletype
  scale_fill_manual(values = sample_colors) +
  scale_shape_manual(values = c(23, 21, 22))+
  theme_linedraw() +                                                           #changes theme, removes grey background
  theme(                             
    legend.title = element_blank(),                                            #removes legend title
    legend.position = "bottom",
    legend.text = element_text(size = 20, face = "bold"),                                 
    axis.text.y.left = element_text(size = 10),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 20),
    strip.text = element_text(face = "bold", size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = c(21))))              #fills legend points based on the fill command


#PERMANOVA

library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(phyloseq)

# Adonis test

#Comparing boat samples to the open water samples
head(sdata)

simple_sdata <- as_tibble(sample_data(boatPort2)) %>% 
  #rownames_to_column(var = "SampleID") %>%
  select(rownames, location_sampled_w_boat_name, sample_type) %>%
  unite(sample_info_by_boat, location_sampled_w_boat_name, sample_type, sep = "_")
head(simple_sdata)
colnames(simple_sdata)

adonisData <- otu_table(boatPort2) %>%
  as.data.frame() %>%
  rownames_to_column(var = "rownames") %>%                                          #change rownames to a column so there is a common variable to join by
  left_join(simple_sdata, by = "rownames") %>%                                            #join sample data to the sequence table
  #filter(!is.na(sample_info_by_boat)) %>%
  select(-c("rownames")) %>%   #remove all metadata columns except the one to be used to compare
  select(sample_info_by_boat, everything())
dim(adonisData)
head(adonisData[,1:10])

any(is.na(adonisData$sample_info_by_boat))
is.numeric(adonisData[,2:27601])

#similarity euclidean from vegdist and holm correction
#pairwise.adonis(x=adonisData[,2:3722],factors=adonisData$Substrate,sim.function='vegdist',
#               sim.method='bray',p.adjust.m='holm')

pairwise.adonis2(adonisData[,2:27601]~sample_info_by_boat,method="bray",data=adonisData, na.rm = TRUE)


# Is there a difference in beta diversity between the two ports?

#adding a phylogenetic tree to phyloseq object using ape library
random_tree = rtree(ntaxa(portsOnly), rooted=TRUE, tip.label=taxa_names(portsOnly))
plot(random_tree)

portsOnly3 = merge_phyloseq(portsOnly, samdata, random_tree)
portsOnly3

#ordination
uni_distance2 <- ordinate(
  physeq = portsOnly3, 
  method = "PCoA", 
  distance = "unifrac"
)

#plot

sample_colors <- c("darkcyan", "chocolate2")

plot_ordination(
  physeq = portsOnly3,                                                          #phyloseq object
  ordination = uni_distance2)+                                                #ordination
  geom_point(aes(fill = location, shape = sample_description), size = 6) +                         #sets fill color to sampletype
  scale_fill_manual(values = sample_colors) +
  scale_shape_manual(values = c(22,21, 23))+
  theme_linedraw() +                                                      #changes theme, removes grey background
  theme(                             
    legend.title = element_blank(),                                      #removes legend title
    legend.position = "right",
    legend.text = element_text(size = 20, face = "bold"),                                 
    axis.text.y.left = element_text(size = 10),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 20),
    strip.text = element_text(face = "bold", size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = c(21))))          #fills legend points based on the fill command


#head(sdata2)
#unifrac_distance <- as.data.frame(uni_distance$vectors) %>%
 # rownames_to_column(var = "SampleID") %>%
  #left_join(sdata2, by = "SampleID") %>%
  #select(SampleID, Microbial_Consortia, Time, Media, row_number, everything())
#head(unifrac_distance)

#write_csv(unifrac_distance, "/home/lgschaer/old/KcBaruah/phyloseq_out/unifrac_distance.csv")

#PERMANOVA

library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(phyloseq)

# Adonis test

#Comparing openwater port samples (Baltimore vs. Norfolk)
head(sdata)

simple_sdata <- as_tibble(sample_data(portsOnly)) %>% 
  dplyr::select(rownames, location)
head(simple_sdata)
#class(simple_sdata)

adonisData <- as.data.frame(otu_table(portsOnly)) %>% 
  rownames_to_column(var = "rownames") %>%
  left_join(simple_sdata, by = "rownames") %>%
  #filter(!is.na(location)) %>%
  select(-c("rownames")) %>%   
  select(location, everything()) 
dim(adonisData)
head(adonisData[,1:10])
class(adonisData)

any(is.na(adonisData))

#similarity euclidean from vegdist and holm correction
#pairwise.adonis(x=adonisData[,2:3722],factors=adonisData$Substrate,sim.function='vegdist',
#               sim.method='bray',p.adjust.m='holm')

pairwise.adonis2(adonisData[,2:27600]~location,method="bray",data=adonisData)


#EXPLORING TAXA


#Summarize abundance of each class
genusabundance <- pps %>%
  tax_glom(taxrank = "Genus") %>%                      # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus) 
head(genusabundance)

#save color palatte
colors10 <- c(
  "black",   "darkcyan",     "orchid1",   "green",       "blue",   "olivedrab",
  "grey47",  "cyan",    "coral1",     "yellow",    "darkgreen",   "palegoldenrod",    
  "grey77",  "darkblue",     "orange",    "red",         "mediumpurple1", "tan4",   "purple4",
  "dodgerblue",    "white", "firebrick", "yellowgreen", "magenta", "blue", "green", "red", "orchid", "lightblue"
) 

#Select and summarize necessary variables
head(genusabundance)

all <- genusabundance %>%
  select(Phylum, Class, Family, Genus, Sample, location_sampled, sample_description, sample_type, sample_detail_description, boat_name,
         location_description, location, true_geographic_order, km_to_baltimore, Abundance) %>%
  filter(Abundance > 0) %>%
  mutate(
    Phylum = as.character(Phylum),
    Class = as.character(Class),
    Family = as.character(Family),
    Genus = as.character(Genus))
head(all)
#View(all)

head(all)





# Is the microbial community of the water distinct between the two ports?

head(all)

genus_ports <- all %>%
  filter(Abundance > 0) %>%
  filter(location_description == "port" & location_sampled == "openwater") %>%
  group_by(sample_detail_description, sample_description, location) %>%
  mutate(totalSum = sum(Abundance)) %>%
  ungroup() %>%
  dplyr::group_by(sample_detail_description, sample_description, location, Phylum, Class, Family, Genus, totalSum) %>%
  #unite(Genus_Species, Genus, Species, sep = " ", remove = FALSE) %>%
  summarise(
    Abundance = sum(Abundance),
    Genus = ifelse(Abundance < 0.12, "< 12 %", Genus)
    #Genus_Species = ifelse(Abundance < 0.03, "< 3 %", Genus_Species)
  ) %>%               #change Genus label to group low abundance taxa together
  group_by(sample_detail_description, sample_description, location, Genus, totalSum) %>%  #now group and summarize again to group newly labeled low abundance taxa together
  summarise(
    Abundance = sum(Abundance),
    RelAb = Abundance/totalSum) %>%
  unique()
head(genus_ports)
#View(genus)

colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(genus_ports$Genus)))
colors <- c("black", color_list)
length(colors)


colors10 <- c(
  "black",  "lightblue",  "firebrick", "darkcyan", "chocolate1",    "orchid1",   "green",       "blue",   "olivedrab",
  "grey47",  "cyan",    "coral1",     "yellow", "maroon3",   "darkgreen",   "palegoldenrod",    
  "grey77",  "darkblue",  "#A765A5",   "orange",    "red",         "mediumpurple1", "tan4",   "purple4",
  "dodgerblue",    "white", "firebrick", "yellowgreen", "magenta", "blue", "#00008B", "green", "red", "orchid", "lightblue"
) 


ggplot(genus_ports)+
  geom_col(mapping = aes(x = location, y = RelAb, fill = Genus), color = "black", position = "stack", show.legend = TRUE)+
  facet_grid(cols = vars(sample_description))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors10)+                        
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "right",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 18, face = "bold", angle = 0),
        strip.text.y = element_text(size = 18, face = "bold", angle = 270),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=1,byrow=FALSE))


# Is the microbial community composition of the boat different than the open water?

genus <- all %>%
  filter(Abundance > 0) %>%
  #filter(location_sampled == "openwater") %>%
  filter(!is.na(sample_description) & sample_description != "blank") %>%
  group_by(location_sampled, sample_description, sample_type, location, true_geographic_order, km_to_baltimore) %>%
  mutate(totalSum = sum(Abundance)) %>%
  ungroup() %>%
  dplyr::group_by(location_sampled, sample_description, sample_type, location, true_geographic_order, km_to_baltimore, Phylum, Class, Family, Genus, totalSum) %>%
  #unite(Genus_Species, Genus, Species, sep = " ", remove = FALSE) %>%
  summarise(
    Abundance = sum(Abundance),
    Phylum = ifelse(Abundance < 0.01, "< 1 %", Phylum)) %>% 
  group_by(location_sampled, sample_description, sample_type, location, true_geographic_order, km_to_baltimore, Phylum, totalSum) %>%  #now group and summarize again to group newly labeled low abundance taxa together
  summarise(
    Abundance = sum(Abundance),
    RelAb = Abundance/totalSum) %>%
  unique()
head(genus)
#View(genus)

colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(genus$Phylum)))
colors <- c("black", color_list)
length(colors)

#unique(genus$location_sampled)
#unique(genus$sample_description)
#unique(genus$sample_type)

ggplot(genus)+
  geom_col(mapping = aes(x = true_geographic_order, y = RelAb, fill = Phylum), color = "black", position = "stack", show.legend = TRUE)+
  facet_nested(rows = vars(location_sampled, sample_type, sample_description), space = "free", scales = "free")+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors)+                        
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "right",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 18, face = "bold", angle = 0),
        strip.text.y = element_text(size = 18, face = "bold", angle = 270),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=1,byrow=TRUE))



#write_csv(genus, "/home/lgschaer/old/KBaruah/phyloseq_out/gecnus_relative_abundances_w_species.csv")




# Is the microbial community composition of the boat different than the open water?

genus <- all %>%
  filter(Abundance > 0) %>%
  filter(location_sampled != "blank") %>%
  filter(!is.na(sample_description)) %>%
  group_by(location_sampled, sample_description, sample_type) %>%
  mutate(totalSum = sum(Abundance)) %>%
  ungroup() %>%
  dplyr::group_by(location_sampled, sample_description, sample_type, Phylum, Class, Family, Genus, totalSum) %>%
  #unite(Genus_Species, Genus, Species, sep = " ", remove = FALSE) %>%
  summarise(
    Abundance = sum(Abundance),
    Class = ifelse(Abundance < 0.10, "< 10 %", Class)) %>% 
  group_by(location_sampled, sample_description, sample_type, Class, totalSum) %>%  #now group and summarize again to group newly labeled low abundance taxa together
  summarise(
    Abundance = sum(Abundance),
    RelAb = Abundance/totalSum) %>%
  unique()
head(genus)
#View(genus)

colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(genus$Class)))
colors <- c("black", color_list)
length(colors)

unique(genus$location_sampled)
unique(genus$sample_description)
unique(genus$sample_type)

ggplot(genus)+
  geom_col(mapping = aes(x = sample_description, y = RelAb, fill = Class), color = "black", position = "stack", show.legend = TRUE)+
  facet_nested(cols = vars(location_sampled, sample_type), space = "free", scales = "free")+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors)+                        
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 18, face = "bold", angle = 0),
        strip.text.y = element_text(size = 18, face = "bold", angle = 270),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(nrow=5,byrow=TRUE))



#write_csv(genus, "/home/lgschaer/old/KBaruah/phyloseq_out/gecnus_relative_abundances_w_species.csv")


# How does the relative abundance of differentially abundant taxa change during the voyage?


#Summarize abundance of each class
asvabundance <- pps %>%
  #tax_glom(taxrank = "Genus") %>%                      # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() 
head(asvabundance)

#save color palatte
colors10 <- c(
  "black",   "darkcyan",     "orchid1",   "green",       "blue",   "olivedrab",
  "grey47",  "cyan",    "coral1",     "yellow",    "darkgreen",   "palegoldenrod",    
  "grey77",  "darkblue",     "orange",    "red",         "mediumpurple1", "tan4",   "purple4",
  "dodgerblue",    "white", "firebrick", "yellowgreen", "magenta", "blue", "green", "red", "orchid", "lightblue"
) 

#Select and summarize necessary variables
head(asvabundance)



diffAbundant <- as.csv("/home/lgschaer/old/Chesapeake/March_2023_Analysis/deseq_out/list_of_diff_abundant_tax.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(diffAbundant)


all <- asvabundance %>%
  filter(Abundance > 0) %>%
  filter(sample_description != "blank") %>%
  mutate(
    taxon = OTU,
    Phylum = as.character(Phylum),
    Class = as.character(Class),
    Family = as.character(Family),
    Genus = as.character(Genus)) %>%
  select(Phylum, Class, Family, Genus, taxon, Sample, location_sampled, sample_description, sample_type, sample_detail_description, boat_name,
         location_description, location, true_geographic_order, km_to_baltimore, Abundance) %>%
  left_join(diffAbundant, by = c("Genus", "taxon")) %>%
  mutate(Where_Enriched = ifelse(is.na(Where_Enriched), "Not_Enriched", Where_Enriched))
head(all)
#View(all)

dim(asvabundance)
dim(diffAbundant)
dim(all)


tax_list <- unique(diffAbundant$Genus)

genusDiff <- all %>%
  filter(Abundance > 0) %>%
  filter(location_sampled != "blank") %>%
  filter(!is.na(sample_description)) %>%
  group_by(location_sampled, sample_description, sample_type, true_geographic_order) %>%
  mutate(totalSum = sum(Abundance)) %>%
  ungroup() %>%
  dplyr::group_by(location_sampled, sample_description, sample_type, true_geographic_order, Where_Enriched,
                  #Phylum, Class, Family, Genus, taxon, 
                  totalSum) %>%
  #unite(Genus_Species, Genus, Species, sep = " ", remove = FALSE) %>%
  summarise(
    Abundance = sum(Abundance),
    #Phylum = ifelse(Abundance < 0.01, "< 1 %", Phylum)
    #Genus = ifelse(Genus %in% tax_list, Genus, "(Not_Enriched)"),#
    ) %>% 
  group_by(location_sampled, sample_description, sample_type, true_geographic_order, Where_Enriched, totalSum) %>%  #now group and summarize again to group newly labeled low abundance taxa together
  summarise(
    Abundance = sum(Abundance),
    RelAb = Abundance/totalSum) %>%
  unique()
head(genusDiff)
#View(genusDiff)

#length(unique(genusDiff$Genus))
#length(tax_list)

#colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
#color_list <- colFunc(length(unique(genus$Genus)))
#colors <- c("black", color_list)
#length(colors)

#unique(genus$location_sampled)
#unique(genus$sample_description)
#unique(genus$sample_type)

#unique(genusDiff$Where_Enriched)

colors <- c("darkcyan", "chocolate1", "black")

ggplot(genusDiff)+
  geom_col(mapping = aes(x = true_geographic_order, y = RelAb, fill = Where_Enriched), color = "black", position = "stack", show.legend = TRUE)+
  facet_nested(rows = vars(location_sampled, sample_type, sample_description), space = "free", scales = "free")+
  ylab("Proportion of Community") +
  xlab("Baltimore --> Norfolk") +
  scale_fill_manual(values = colors)+                        
  #xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.title.x = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 18, face = "bold", angle = 0),
        strip.text.y = element_text(size = 18, face = "bold", angle = 270),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))



