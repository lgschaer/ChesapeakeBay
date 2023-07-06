library("biomformat")
library(vctrs)
library(tidyverse)
library(csv)
library(ggh4x)

# Preparing the input files for Sourcetracker

ps <- readRDS("/home/lgschaer/old/Chesapeake/March_2023_Analysis/phyloseq_out/ps_unrarefied_03212023.rds")
ps

# making map.txt

ps2 <- ps %>%
  subset_samples(sample_description != "blank") %>% 
  subset_samples(sample_sums(ps_filt) > 2000)

map <- as_tibble(sample_data(ps2)) %>%
  select(SampleID, location_sampled, sample_type, sample_description, sample_detail_description, boat_name, location_description, location, true_geographic_order) %>%
  unite(sample_detail_description_boat_name, sample_detail_description, boat_name, sep = "_", remove = FALSE) %>%
  mutate(
    sample_detail_description_boat_name = ifelse(location_sampled == "openwater", sample_detail_description, sample_detail_description_boat_name),
    SourceSink = ifelse(location_sampled == "openwater" & (location == "Baltimore" | location == "Norfolk"), "source", "sink"),
    Env = ifelse(SourceSink == "source", paste0(location, "_", sample_description), sample_detail_description_boat_name)
  ) %>%
  select(SampleID, SourceSink, Env) %>%
  filter(!is.na(Env)) 
head(map)
dim(map)
#colnames(map) <- c("#SampleID", "SourceSink", "Env")
rownames(map) <- NULL
#View(map)
head(map)

write.table(map, "/home/lgschaer/old/Chesapeake/March_2023_Analysis/sourcetracker/inputs/map.txt", sep="\t", row.names=FALSE, quote=FALSE)


# making otu.biom
count_table <- (otu_table(ps2))
count_table[1:5,1:5]

#transpose
t_count_table <- t(count_table)
t_count_table[,5:1][1:5,]
class(t_count_table)

#write table
write.table(t_count_table, "/home/lgschaer/old/Chesapeake/March_2023_Analysis/sourcetracker/inputs/otu.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

##next step: use putty to convert to .biom table

####### Processing Sourcetracker Output

sdataOG <- as.csv("/home/lgschaer/old/Chesapeake/March_2023_Analysis/phyloseq_out/full_cleaned_metadata_edited04112023.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE) 
colnames(sdataOG)

sdata <- sdataOG%>%
  mutate(geographic_order = true_geographic_order) %>%
  select(SampleID, location_sampled, sample_type, sample_description, sample_detail_description, latitude, longitude, geographic_order, station, boat_name, transit, location_description, location)
head(sdata)

mxProp <- read.table("/home/lgschaer/old/Chesapeake/March_2023_Analysis/sourcetracker/out_sourcePortsByFilter_sinkBoats_04112023MD//mixing_proportions.txt")
head(mxProp)
dim(mxProp)

mxSD <- read.table("/home/lgschaer/old/Chesapeake/March_2023_Analysis/sourcetracker/out_sourcePortsByFilter_sinkBoats_04112023MD/mixing_proportions_stds.txt") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = c("Baltimore_PES_0.2_um", "Baltimore_glass_fiber", "Norfolk_PES_0.2_um", "Norfolk_glass_fiber", "Unknown"), names_to = "Source", values_to = "StdDev")
head(mxSD)
dim(mxSD)

#colnames(mxSD)


## Is the boat microbial community sourced from the water?
mxProp2 <- mxProp %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sdata, by = "SampleID") %>%
  pivot_longer(cols = c("Baltimore_PES_0.2_um", "Baltimore_glass_fiber", "Norfolk_PES_0.2_um", "Norfolk_glass_fiber", "Unknown"), names_to = "Source", values_to = "Proportion") %>%
  left_join(mxSD, by = c("SampleID", "Source")) %>%
  filter(location_sampled == "boat") %>%
  group_by(sample_description, Source, sample_type, boat_name) %>%
  summarise(
    meanProp = mean(Proportion),
    maxSD = max(StdDev)
  )
head(mxProp2)
#View(mxProp2)

mxProp2$meanProp

colors <- c("darkcyan", "cyan", "chocolate1", "orange", "grey")

ggplot(mxProp2, aes(x=boat_name, y = meanProp))+
  geom_col(mapping = aes(fill = Source), color = "black", position = "fill", show.legend = TRUE)+
  facet_nested(cols = vars(sample_type, sample_description), space = "free", scales = "free")+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors) +
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 18, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))


mxProp2B <- mxProp %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sdata, by = "SampleID") %>%
  pivot_longer(cols = c("Baltimore_PES_0.2_um", "Baltimore_glass_fiber", "Norfolk_PES_0.2_um", "Norfolk_glass_fiber", "Unknown"), names_to = "Source", values_to = "Proportion") %>%
  left_join(mxSD, by = c("SampleID", "Source")) %>%
  filter(location_sampled == "boat") %>%
  mutate(Source2 = ifelse(Source == "Unknown", Source, "Port"))%>%
  group_by(sample_description, Source2, sample_type, boat_name) %>%
  summarise(
    meanProp = mean(Proportion),
    maxSD = max(StdDev)
  )
head(mxProp2B)
#View(mxProp2B)


mxProp3 <- mxProp %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sdata, by = "SampleID") %>%
  pivot_longer(cols = c("Baltimore_PES_0.2_um", "Baltimore_glass_fiber", "Norfolk_PES_0.2_um", "Norfolk_glass_fiber", "Unknown"), names_to = "Source", values_to = "Proportion") %>%
  mutate(
    Source = case_when(
    Source == "Baltimore_PES_0.2_um" ~ "Baltimore 0.2 um PES filter",
    Source == "Norfolk_PES_0.2_um" ~ "Norfolk 0.2 um PES filter",
    Source == "Baltimore_glass_fiber" ~ "Baltimore glass fiber filter",
    Source == "Norfolk_glass_fiber" ~ "Norfolk glass fiber filter",
    Source == "Unknown" ~ "Unknown",
  )) %>%
  left_join(mxSD, by = c("SampleID", "Source")) %>%
  filter(location_sampled == "boat") %>%
  group_by(sample_description, Source, sample_type) %>%
  summarise(
    meanProp = mean(Proportion),
    maxSD = max(StdDev)
  )
head(mxProp3)
#View(mxProp3)


colors <- c("darkcyan", "cyan", "chocolate1", "orange", "grey")

ggplot(mxProp3, aes(x=sample_description, y = meanProp))+
  geom_col(mapping = aes(fill = Source), color = "black", position = "fill", show.legend = TRUE)+
  facet_nested(cols = vars(sample_type), space = "free", scales = "free")+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors) +
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 18, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(nrow=3,byrow=TRUE))


## Do boat microbial signatures persist from one location to another?




mxProp3a <- mxProp %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sdata, by = "SampleID") %>%
  pivot_longer(cols = c("Baltimore_PES_0.2_um", "Baltimore_glass_fiber", "Norfolk_PES_0.2_um", "Norfolk_glass_fiber", "Unknown"), names_to = "Source", values_to = "Proportion") %>%
  filter(!is.na(station)) %>%
  mutate(Source = factor(Source, levels = c("Baltimore_glass_fiber", "Baltimore_PES_0.2_um", "Unknown", "Norfolk_PES_0.2_um", "Norfolk_glass_fiber")),
         station = as.numeric(station),
         direction = ifelse(transit == "t1" | transit == "t3", "Departing", "Returning"),
         direction_abbr = ifelse(direction == "Departing", "D", "R")) %>%
  unite(direction_station, direction_abbr, station, sep = "", remove = FALSE) %>%
  mutate(direction_station = factor(direction_station, levels = c("D1", "D2", "D3", "D4", "D5", "D6", "D8", "D10", "D11", "D12", "D13", "D14", "D15", "D16", "D17",
                                                                  "D18", "D19", "D21", "D22", "D23", "D24", "D25", "D26", "D27", "D28", "R28", "R26", "R24",  "R22",
                                                                  "R20", "R19", "R17", "R15", "R13", "R12", "R9", "R7")),
         sample_type = factor(sample_type, levels = c("filter", "bilge", "surface")),
         sample = case_when(sample_detail_description == "boat_bilge_biofilm" ~ "Bilge Biofilm",
                            sample_detail_description == "boat_bilge_water" ~ "Bilge Water",
                            sample_detail_description == "boat_surface_hull" ~ "Hull",
                            sample_detail_description == "boat_surface_transom" ~ "Transom",
                            sample_detail_description == "openwater_filter_glass_fiber" ~ "Open Water (GF)",
                            sample_detail_description == "openwater_filter_PES_0.2_um" ~ "Open Water (PES)"),
         sample = factor(sample, levels= c("Open Water (PES)", "Open Water (GF)", "Hull", "Transom", "Bilge Water", "Bilge Biofilm"))) %>%
  select(SampleID, location_sampled, sample, sample_type, sample_description, sample_detail_description, geographic_order, direction_station, boat_name, location_description, location, Source, Proportion) %>%
  separate(Source, into = c("Source", NA), sep = "_") %>%
  filter(Source != "Unknown")
head(mxProp3a)
#View(mxProp3)


#colors <- c("darkcyan", "cyan", "grey", "orange", "chocolate1")
colors <- c("skyblue1", "skyblue3", "maroon4", "maroon2", "olivedrab1", "olivedrab4")


ggplot(mxProp3a, aes(x=direction_station, y = Proportion))+
  geom_boxplot(mapping = aes(fill = sample), color = "black", linewidth = 0.1, position = "dodge", show.legend = TRUE)+
  facet_nested(cols = vars(boat_name), rows = vars(Source), space = "free", scales = "free")+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors) +
  xlab("Baltimore to Norfolk")+
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
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=3,byrow=TRUE))



head(mxProp)

mxProp3b <- mxProp %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sdata, by = "SampleID") %>%
  pivot_longer(cols = c("Baltimore_PES_0.2_um", "Baltimore_glass_fiber", "Norfolk_PES_0.2_um", "Norfolk_glass_fiber", "Unknown"), names_to = "Source", values_to = "Proportion") %>%
  filter(!is.na(station)) %>%
  mutate(Source = factor(Source, levels = c("Baltimore_glass_fiber", "Baltimore_PES_0.2_um", "Unknown", "Norfolk_PES_0.2_um", "Norfolk_glass_fiber")),
    station = as.numeric(station),
    direction = ifelse(transit == "t1" | transit == "t3", "Departing", "Returning"),
    direction_abbr = ifelse(direction == "Departing", "D", "R")) %>%
  unite(direction_station, direction_abbr, station, sep = "", remove = FALSE) %>%
  mutate(direction_station = factor(direction_station, levels = c("D1", "D2", "D3", "D4", "D5", "D6", "D8", "D10", "D11", "D12", "D13", "D14", "D15", "D16", "D17",
                                                                  "D18", "D19", "D21", "D22", "D23", "D24", "D25", "D26", "D27", "D28", "R28", "R26", "R24",  "R22",
                                                                  "R20", "R19", "R17", "R15", "R13", "R12", "R9", "R7"))) %>%
  #filter(location_sampled == "boat") %>%
  group_by(sample_type, sample_description, boat_name, direction_station, Source) %>%
  summarise(
    meanProp = mean(Proportion)
  ) 
head(mxProp3b)
#View(mxProp3b)

#write_csv(mxProp3b, "/home/lgschaer/old/Chesapeake/sourcetracker_output_06012023.csv")

colors <- c("darkcyan", "cyan", "grey", "orange", "chocolate1")

ggplot(mxProp3b, aes(x=direction_station, y = meanProp))+
  geom_col(mapping = aes(fill = Source), color = "black", linewidth = 0.1, position = "fill", show.legend = TRUE)+
  facet_nested(cols = vars(boat_name), rows = vars(sample_type, sample_description), space = "free", scales = "free")+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors) +
  xlab("Baltimore to Norfolk")+
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
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=3,byrow=TRUE))



#library(ggridges)

#data(iris)

#mxProp3 <- mxProp3 %>% mutate(station = as.factor(station),
 #                             meanProp = as.factor(meanProp))

#ggplot(mxProp3, aes(x = meanProp, y = Source)) +
 # geom_density_ridges(aes(fill = Source)) +
  #scale_fill_manual(values = c("darkcyan", "orange", "grey"))

head(mxProp)

mxProp4 <- mxProp %>% 
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sdata, by = "SampleID") %>%
  pivot_longer(cols = c("Baltimore_PES_0.2_um", "Baltimore_glass_fiber", "Norfolk_PES_0.2_um", "Norfolk_glass_fiber", "Unknown"), names_to = "Source", values_to = "Proportion") %>%
  filter(!is.na(station)) %>%
  mutate(station = as.numeric(station),
         direction = ifelse(transit == "t1" | transit == "t3", "Departing", "Returning")) %>%
  #filter(location_sampled == "boat") %>%
  group_by(sample_type, sample_description, station, boat_name, transit, geographic_order, direction, Source) %>%
  summarise(
    meanProp = mean(Proportion)
  ) %>%
  #filter(station != 13 | station != 14 | station != 10) %>%
  mutate(direction_abbr = ifelse(direction == "Departing", "D", "R")) %>%
  unite(direction_station, direction_abbr, station, sep = "", remove = FALSE) %>%
  mutate(
    Source = factor(Source, levels = c("Baltimore_glass_fiber", "Baltimore_PES_0.2_um", "Unknown", "Norfolk_PES_0.2_um", "Norfolk_glass_fiber")),
    direction_station = factor(direction_station, levels = c("D1", "D2", "D3", "D4", "D5", "D6", "D8", "D10", "D11", "D12", "D13", "D14", "D15", "D16", "D17",
                                                                  "D18", "D19", "D21", "D22", "D23", "D24", "D25", "D26", "D27", "D28", "R28", "R26", "R24",  "R22",
                                                                  "R20", "R19", "R17", "R15", "R13", "R12", "R9", "R7")))
head(mxProp4)
#View(mxProp4)
length(unique(mxProp4$direction_station))


mxProp4B <- mxProp4 %>%
  filter(sample_type != "filter")


colors <- c("darkcyan", "cyan", "grey", "orange", "chocolate1")

ggplot(mxProp4B, aes(x=direction_station, y = meanProp))+
  geom_col(mapping = aes(fill = Source), linewidth = 0.1, color = "black")+
  ylim(0,1)+
  facet_nested(cols = vars(boat_name, direction), rows = vars(sample_type, sample_description), space = "free_x", scales = "free_x")+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors) +
  xlab("Baltimore to Norfolk")+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 16, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        panel.background = element_rect(fill = "grey"),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

mxProp4B_stats <- mxProp4B %>%
  filter(Source != "Unknown") %>%
  pivot_wider(names_from = Source, values_from = meanProp) %>%
  mutate(
    PortWater = Baltimore_glass_fiber + Baltimore_PES_0.2_um + Norfolk_glass_fiber + Norfolk_PES_0.2_um
  ) %>%
  pivot_longer(cols = c("PortWater"), names_to = "Source", values_to = "meanProp") %>%
  group_by(sample_type, sample_description, Source, boat_name) %>%
  summarise(
    minProp = min(meanProp),
    maxProp = max(meanProp),
    sdProp = sd(meanProp),
    meanProp = mean(meanProp)
  )
head(mxProp4B_stats)
View(mxProp4B_stats)


colnames(mxProp4)

mxProp5 <- mxProp4 %>%
  filter(sample_type == "filter") %>%
  filter(boat_name != "Callinectes" | direction_station != "D10") %>%  #station D_10 removed for the departing voyage due to suspected incorrect labeling
  group_by(direction_station, Source, direction, sample_type, sample_description) %>%
  summarise(meanProp = mean(meanProp))


ggplot(mxProp5, aes(x=direction_station, y = meanProp))+
  geom_col(mapping = aes(fill = Source), linewidth = 0.2, color = "black")+
  ylim(0,1)+
  facet_nested(cols = vars(direction), rows = vars(sample_type, sample_description), space = "free_x", scales = "free_x")+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors) +
  xlab("Baltimore to Norfolk")+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 18, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.spacing.x = unit(0.4, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
       # panel.background = element_rect(fill = "grey"),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

head(mxProp5)

mxProp5_stats <- mxProp5 %>%
  filter(Source != "Unknown") %>%
  filter(direction_station == "D5" |
           direction_station == "D24" |
           direction_station == "R24" |
           direction_station == "R7") %>%
  mutate(Source = as.character(Source),
         #Source = ifelse(str_detect(Source, "Baltimore."), "Baltimore", Source),
         #Source = ifelse(str_detect(Source, "Norfolk."), "Norfolk", Source)
         ) %>%
  pivot_wider(names_from = Source, values_from = meanProp) %>%
  mutate(
    Baltimore = Baltimore_glass_fiber + Baltimore_PES_0.2_um,
    Norfolk = Norfolk_glass_fiber + Norfolk_PES_0.2_um
  ) %>%
  pivot_longer(cols = c("Baltimore", "Norfolk"), names_to = "Source", values_to = "meanProp") %>%
  group_by(direction_station, Source) %>%
  summarise(
    sdProp = sd(meanProp),
    meanProp = mean(meanProp))
head(mxProp5_stats)
View(mxProp5_stats)

## Looking at the feature tables from Sourcetracker

fttable <- read.table("/home/lgschaer/old/Chesapeake/March_2023_Analysis/sourcetracker/out_sourcePorts_sinkBoats_fixedmetadata/feature_tables/c11t3bw.feature_table.txt")
head(fttable[3:10])
dim(fttable)

fttable2 <- read.table("/home/lgschaer/old/Chesapeake/March_2023_Analysis/sourcetracker/out_sourcePorts_sinkBoats_fixedmetadata/feature_tables/c12t3bw.feature_table.txt")
head(fttable2[3:10])
dim(fttable2)

setwd("/home/lgschaer/old/Chesapeake/March_2023_Analysis/sourcetracker/out_sourcePorts_sinkBoats_fixedmetadata/feature_tables/")
getwd()

file_names <- dir("./") #where you have your files
file_names

#file_names2 <- file_names[1:4]

feature_tables <- do.call(rbind,lapply(file_names,read_tsv))
head(feature_tables[1:5, 1:10])
dim(feature_tables)


feature_tables2 <- feature_tables %>%
  rownames_to_column(var = "row_num") %>%
  unite(SampleID, row_num, ...1, sep = "_") %>%
  separate(SampleID, into = c("row_num", "Condition"), sep = "_")
head(feature_tables2[1:5, 1:10])

write_csv(feature_tables2, "/home/lgschaer/old/Chesapeake/March_2023_Analysis/sourcetracker/out_sourcePorts_sinkBoats_fixedmetadata/feature_tables_combined.csv")
dim(feature_tables2)

#adding sampleIDs

#fruits <- c("banana", "apple", "orange", "cherry", "pineapple")
#length(fruits)

#fruits2 <- vec_rep_each(fruits, 3)
#length(fruits2)
#fruits2

file_names
length(file_names)

file_names2 <- vec_rep_each(file_names, 3)
length(file_names2)
head(file_names2)

#### Adding metadata

#load metadata
sdata <- as.csv("/home/lgschaer/old/Chesapeake/March_2023_Analysis/phyloseq_out/full_cleaned_metadata_edited04062023.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE) %>%
  select(-c("npoc_um", "tn_um", "phosphate_umoll", "silicate_umoll", "nn_umoll", "ammonia_umoll", "tds_mgl", "temp_f", "temp_c", "odo_mgl",    "ph",   "orp"))
head(sdata)
dim(sdata)

#### Adding Taxa Names
ASV <- colnames(feature_tables2[,3:55147])
length(ASV)
head(ASV)

#load taxa table
taxa <- readRDS("/home/lgschaer/old/Chesapeake/March_2023_Analysis/dada2_out/taxtab.rds")   
rownames(taxa) <- NULL
taxa[2:5,2:5]
dim(taxa)

taxa2 <- cbind(ASV, taxa) %>% as.data.frame()
taxa2[1:5,1:5]

feature_tables3 <- cbind(file_names2, feature_tables2)

feature_table_summary <- feature_tables3 %>%
  column_to_rownames(var = "row_num") %>%
  separate(file_names2, into = c("SampleID", NA), sep = ".feature") %>%
  pivot_longer(cols = c(everything(), -SampleID, -Condition), names_to = "ASV", values_to = "Count") %>%
  left_join(sdata, by = "SampleID") %>%
  left_join(taxa2, by = "ASV") %>%
  filter(Count != 0) %>%
  group_by(SampleID, sample_detail_description, Condition, boat_name) %>%
  mutate(totalSum = sum(Count),
         Genus = ifelse(is.na(Genus), paste0("Unclassified Genus from ", Class), Genus)) %>%
  group_by(SampleID, sample_detail_description, Condition, Class, Genus, boat_name) %>%
  summarise(RelAb = Count/totalSum) %>%
  group_by(SampleID, sample_detail_description, Condition, Class, Genus, boat_name) %>%
  summarise(RelAb = sum(RelAb),
            Genus = ifelse(RelAb < 0.04, "< 4 %", Genus)) %>%
  group_by(SampleID, sample_detail_description, Condition, Class, Genus, boat_name) %>%
  summarise(RelAb = sum(RelAb)) %>%
  select(SampleID, sample_detail_description, boat_name, Class, Genus, RelAb, everything())
head(feature_table_summary)
dim(feature_table_summary)
colnames(feature_table_summary)


write_csv(feature_table_summary, "/home/lgschaer/old/Chesapeake/March_2023_Analysis/sourcetracker/out_sourcePorts_sinkBoats_fixedmetadata/feature_tables_combined.csv")


colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(feature_table_summary$Genus)))
colors <- c("black", color_list)
length(colors)

ggplot(feature_table_summary, aes(x = SampleID, y = RelAb, fill = Genus))+
  facet_nested(cols = vars(sample_detail_description, boat_name), rows = vars(Condition), space = "free", scales = "free")+
  geom_col(color = "black", position = "stack", show.legend = FALSE)+
  scale_fill_manual(values = colors) +
  #scale_y_reverse()+
  ylab(NULL)+
  xlab(NULL)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 15, color = "black"),
        #axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1, color = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 20),
        #legend.text = element_text(size = 15, angle = 0, vjust = 0.5, hjust = 0),
        #legend.position = "right",
        #legend.title = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 1, face = "bold"),
        strip.background = element_rect(fill = "white"),
        strip.text.y = element_text(size = 17, face = "bold", angle = 270, hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text(size = 16, face = "bold", angle = 0),
        title = element_text(size = 18))+
  guides(fill=guide_legend(ncol=1,byrow=FALSE))


feature_table_summary2 <- feature_tables3 %>%
  column_to_rownames(var = "row_num") %>%
  separate(file_names2, into = c("SampleID", NA), sep = ".feature") %>%
  pivot_longer(cols = c(everything(), -SampleID, -Condition), names_to = "ASV", values_to = "Count") %>%
  left_join(sdata, by = "SampleID") %>%
  left_join(taxa2, by = "ASV") %>%
  filter(Count != 0) %>%
  group_by(sample_detail_description, Condition) %>%
  mutate(totalSum = sum(Count),
         Genus = ifelse(is.na(Genus), paste0("Unclassified Genus from ", Class), Genus)) %>%
  group_by(sample_detail_description, Condition, Class, Genus) %>%
  summarise(RelAb = Count/totalSum) %>%
  group_by(sample_detail_description, Condition, Class, Genus) %>%
  summarise(RelAb = sum(RelAb),
            Genus = ifelse(RelAb < 0.04, "< 4 %", Genus)) %>%
  group_by(sample_detail_description, Condition, Class, Genus) %>%
  summarise(RelAb = sum(RelAb)) %>%
  select(sample_detail_description, Class, Genus, RelAb, everything())
head(feature_table_summary2)



write_csv(feature_table_summary2, "/home/lgschaer/old/Chesapeake/March_2023_Analysis/sourcetracker/out_sourcePorts_sinkBoats_fixedmetadata/feature_tables_combined_sample_detail_description.csv")


colFunc <- colorRampPalette(c("purple4", "lightpink", "firebrick", "orange", "lemonchiffon", "olivedrab4", "darkcyan", "lightblue", "darkblue"))
color_list <- colFunc(length(unique(feature_table_summary2$Genus)))
colors <- c("black", color_list)
length(colors)

ggplot(feature_table_summary2, aes(x = sample_detail_description, y = RelAb, fill = Genus))+
  facet_nested(cols = vars(Condition), space = "free", scales = "free")+
  geom_col(color = "black", position = "stack")+
  scale_fill_manual(values = colors) +
  #scale_y_reverse()+
  ylab(NULL)+
  xlab(NULL)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1, color = "black"),
        #axis.text.x = element_blank(),
        axis.title.x = element_text(size = 20),
        legend.text = element_text(size = 15, angle = 0, vjust = 0.5, hjust = 0),
        legend.position = "right",
        #legend.title = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 1, face = "bold"),
        strip.background = element_rect(fill = "white"),
        strip.text.y = element_text(size = 17, face = "bold", angle = 270, hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text(size = 16, face = "bold", angle = 0),
        title = element_text(size = 18))+
  guides(fill=guide_legend(ncol=1,byrow=FALSE))

