library("biomformat")
library(vctrs)
library(tidyverse)
library(csv)
library(ggh4x)
#install.packages("pamr", dependencies=TRUE)
library(pamr)
library(caret)
#install.packages("ggpmisc", dependencies = TRUE)
#library(ggpmisc)

# Preparing for machine learning - prepare training set

#ps_global <- readRDS("/home/lgschaer/old/Chesapeake/bothasvs.rds")
#ps_global

ps_ches <- readRDS("/home/lgschaer/old/Chesapeake/March_2023_Analysis/phyloseq_out/ps_rarefied_nochloro_justbacteria_06132023.rds")
ps_ches

ps_training <- subset_samples(ps_ches, (location == "Baltimore" | location == "Norfolk") & sample_type == "filter") %>%
  subset_samples(sample_sums(ps_ches) > 0) %>%
  subset_taxa(taxa_sums(ps_ches) > 0)
ps_training
head(sample_data(ps_training))

training_sd <- as_tibble(sample_data(ps_training)) %>%
  select(SampleID, location)
head(training_sd)

tr_class <- training_sd$location

training_otu <- as.data.frame(otu_table(ps_training)) %>%
  rownames_to_column(var = "SampleID")
training_otu[1:5,1:5]

training <- training_sd %>%
  select(-location) %>%
  left_join(training_otu) %>%
  column_to_rownames(var = "SampleID")
training[1:3,1:5]

dim(training)
length(class)



#### prepare testing set

ps_testing <- subset_samples(ps_ches, (location != "Baltimore" | location != "Norfolk")) %>%
  subset_samples(sample_sums(ps_ches) > 0) %>%
  subset_taxa(taxa_sums(ps_ches) > 0)
ps_testing
ps_ches
ps_training
head(sample_data(ps_testing))

testing_sd <- as_tibble(sample_data(ps_testing)) %>%
  select(SampleID, location)
head(testing_sd)

ts_class <- testing_sd$location
ts_sid <- testing_sd$SampleID

testing_otu <- as.data.frame(otu_table(ps_testing)) %>%
  rownames_to_column(var = "SampleID")
testing_otu[1:5,1:5]

testing <- testing_sd %>%
  select(-location) %>%
  left_join(testing_otu) %>%
  column_to_rownames(var = "SampleID")
testing[1:3,1:5]

dim(testing)
length(class)




#find any predictors that are zero variance or near zero variance
nzv_tr <- nearZeroVar(training)
nzv_ts <- nearZeroVar(testing)
nzv <- unique(c(nzv_tr, nzv_ts))

length(nzv_tr)
length(nzv_ts)
length(nzv)

dim(training)
dim(testing)

training_filt <- training[-nzv]
nearZeroVar(training_filt)

testing_filt <- testing[-nzv]
nearZeroVar(testing_filt)

dim(training_filt)
dim(testing_filt)

#Nearest Shrunken Centroids Model

nscGrid <- data.frame(.threshold = 0:30)
set.seed(476)
ctrl <- trainControl(classProbs = TRUE, method = "repeatedcv", number = 3, repeats = 25)
nscTuned <- train(x = training_filt,
                  y = tr_class,
                  method = "pam",
                  preProcess = c("center", "scale"),
                  tuneGrid = nscGrid,
                  metric = "Accuracy",
                  trControl = ctrl)
nscTuned

varImp(nscTuned)

nscTop20 <- c("sp89", "sp24", "sp39", "sp1", "sp22", "sp3", "sp37", "sp38", "sp11", "sp260", "sp150", "sp12", "sp17", "sp65", "sp25", 
              "sp9", "sp2", "sp6", "sp223", "sp4")


sdata <- as_tibble(sample_data(ps_ches)) %>% select(c("SampleID", "location"))
head(sdata)


nscPredict <- predict(nscTuned, newdata = testing_filt, type = "prob") %>%
  rownames_to_column(var = "SampleID") %>%
  #left_join(sdata, by = "SampleID") %>%
  pivot_longer(cols = c("Baltimore", "Norfolk"), names_to = "nsc_Port", values_to = "nsc_Prob")# %>%
  #filter(location == "Baltimore")
head(nscPredict)
#View(nscPredict)


#Random Forest
dim(training)

#rfGrid <- data.frame(.mtry = c(2, 500, 5000, 10000, 50000))
#set.seed(476)
#ctrl <- trainControl(classProbs = TRUE, method = "repeatedcv", number = 3, repeats = 25)
#rfTuned <- train(x = training,
 #                 y = tr_class,
  #                method = "rf",
   #               tuneGrid = rfGrid,
    #              metric = "Accuracy",
     #             trControl = ctrl)
#rfTuned

#varImp(rfTuned)

#rfPredict <- predict(rfTuned, newdata = testing, type = "prob") %>%
 # rownames_to_column(var = "SampleID2") %>%
  #pivot_longer(cols = c("Baltimore", "Norfolk"), names_to = "rf_Port", values_to = "rf_Prob")
#summary(rfPredict)
#rfPredict

#dim(training)

rfGrid <- data.frame(.mtry = c(2, 10, 25, 50, 75, 100))
set.seed(476)
ctrl <- trainControl(classProbs = TRUE, method = "repeatedcv", number = 3, repeats = 25)
rfTuned_nonzv <- train(x = training_filt,
                 y = tr_class,
                 method = "rf",
                 tuneGrid = rfGrid,
                 metric = "Accuracy",
                 trControl = ctrl)
rfTuned_nonzv

varImp(rfTuned_nonzv)

rfTop20 <- c("sp3", "sp89", "sp39", "sp1", "sp24", "sp22", "sp11", "sp38", "sp37", "sp4", "sp2", "sp150", "sp25",
             "sp12", "sp260", "sp65", "sp9", "sp769", "sp17", "sp6")

rfPredict_nonzv <- predict(rfTuned_nonzv, testing_filt, type = "prob") %>%
  rownames_to_column(var = "SampleID3") %>%
  pivot_longer(cols = c("Baltimore", "Norfolk"), names_to = "rfnzv_Port", values_to = "rfnzv_Prob")
rfPredict_nonzv
print(rfPredict_nonzv)
#confusionMatrix(ts_class, rfPredict_nonzv)


#combining all the results!
dim(nscPredict)
#dim(rfPredict)
dim(rfPredict_nonzv)

sdata <- as.data.frame(sample_data(ps_ches))
colnames(sdata)

results_temp <- nscPredict %>%
  cbind(rfPredict_nonzv) %>%
  select(-c("SampleID3")) %>%
  left_join(sdata) %>%
  filter(!(is.na(true_geographic_order))) %>%
  mutate(
    direction = ifelse(transit == "t1" | transit == "t3", "Departing", "Returning"),
    direction_abbr = ifelse(direction == "Departing", "D", "R")) %>%
  unite(direction_station, direction_abbr, station, sep = "", remove = FALSE) %>%
  mutate(direction_station = factor(direction_station, levels = c("D1", "D2", "D3", "D4", "D5", "D6", "D8", "D10", "D11", "D12", "D13", "D14", "D15", "D16", "D17",
                                                                  "D18", "D19", "D21", "D22", "D23", "D24", "D25", "D26", "D27", "D28", "R28", "R26", "R24",  "R22",
                                                                  "R20", "R19", "R17", "R15", "R13", "R12", "R9", "R7"))) %>%
  filter(sample_description != "blank" & !is.na(sample_description)) %>%
  mutate(
    rfnzv = "rfnzv",
    nsc = "nsc") %>%
  unite(nsc_Port, nsc_Port, nsc, sep = "_", remove = TRUE) %>%
  unite(rfnzv_Port, rfnzv_Port, rfnzv, sep = "_", remove = TRUE) %>%
  pivot_longer(cols = c("nsc_Port", "rfnzv_Port"), names_to = "Model", values_to = "Predicted_Port") %>%
  pivot_longer(cols = c("nsc_Prob", "rfnzv_Prob"), names_to = "Probability_Type", values_to = "Probability") %>%
  separate(Predicted_Port, into = c("Predicted_Port", "Predicted_Port_Model"), sep = "_") %>%
  separate(Probability_Type, into = c("Probability_Type", NA), sep = "_") %>%
  filter(Predicted_Port_Model == Probability_Type)

results <- results_temp %>%
  group_by(direction_station, Predicted_Port, boat_name, sample_type, Model) %>%
  summarise(
    Probability = mean(Probability)
  )
head(results)
dim(results)
colnames(results)
#View(results)

#sum(is.na(results$true_geographic_order))

colors <- c("darkcyan", "chocolate1")

ggplot(results, aes(x = direction_station, y = Probability, fill = Predicted_Port))+
  facet_nested(rows = vars(boat_name, Model), cols = vars(sample_type), scales = "free", space = "free")+
  geom_col(color = "black")+
  scale_fill_manual(values = colors) +
  #xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        legend.text = element_text(size = 15),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=1,byrow=TRUE))


# show the open water signature and the boat samples at the same time
#both models plus sourcetracker

#load sourcetracker output
sourcetracker <- as.csv("/home/lgschaer/old/Chesapeake/sourcetracker_output_06012023.csv") %>%
  mutate(
    Model = "SourceTracker2",
    Source = ifelse(str_detect(Source, "Baltimore."), "Baltimore", Source),
    Source = ifelse(str_detect(Source, "Norfolk."), "Norfolk", Source)) %>%
  select(SampleID, sample_type, sample_description, sample_detail_description, direction_station, boat_name, location, Model, Source, Proportion) %>%
  pivot_wider(names_from = Source, values_from = Proportion, values_fn = sum) %>%
  pivot_longer(cols = c("Baltimore", "Norfolk", "Unknown"), names_to = "Predicted_Port", values_to = "Probability")
head(sourcetracker)
View(sourcetracker)


head(results_temp)
head(results2)

colnames(results_temp)
colnames(sourcetracker)

dim(sourcetracker)
dim(results_temp)
dim(results2)

unique(results2$Model)

results2 <- results_temp %>%
  full_join(sourcetracker, by = c("SampleID", "sample_type", "sample_description", "sample_detail_description", "direction_station", "boat_name", "location", "Model", "Predicted_Port", "Probability")) %>%
  mutate(sample_desc = ifelse(sample_type == "filter", "filter", sample_detail_description)) %>%
  group_by(direction_station, Predicted_Port, boat_name, sample_type, sample_desc, Model) %>%
  summarise(
    Probability = mean(Probability)
  ) %>%
  pivot_wider(names_from = sample_desc, values_from = Probability) %>%
  pivot_longer(cols = c("boat_bilge_biofilm", "boat_bilge_water", "boat_surface_hull", "boat_surface_transom"), names_to = "sample_desc", values_to = "Probability") %>%
  mutate(
    direction_station = factor(direction_station, levels = c("D1", "D2", "D3", "D4", "D5", "D6", "D8", "D10", "D11", "D12", "D13", "D14", "D15", "D16", "D17",
                                                                    "D18", "D19", "D21", "D22", "D23", "D24", "D25", "D26", "D27", "D28", "R28", "R26", "R24",  "R22",
                                                                    "R20", "R19", "R17", "R15", "R13", "R12", "R9", "R7")),
    direction_station_int = as.integer(direction_station),
    Model = case_when(
      Model == "nsc_Port" ~ "Nearest Shrunken Centroids", 
      Model == "rfnzv_Port" ~ "Random Forest",
      Model == "SourceTracker2" ~ "SourceTracker2"),
    sample_desc = case_when(
      sample_desc == "boat_bilge_biofilm" ~ "Bilge Biofilm",
      sample_desc == "boat_bilge_water" ~ "Bilge Water",
      sample_desc == "boat_surface_hull" ~ "Hull",
      sample_desc == "boat_surface_transom" ~ "Transom"
    ),
    Signature = "Open Water Signature",
    Boat_Signature = "Boat Sample Signature"
  ) %>%
  filter(Predicted_Port == "Baltimore")
head(results2)
colnames(results2)
#View(results2)




#colors <- c("chocolate1", "darkcyan")
colors <- c("lightblue")

ggplot(results2, aes(x = direction_station_int))+
  facet_nested(rows = vars(Model), 
               cols = vars(sample_desc, boat_name), 
               scales = "free", space = "free")+
  geom_vline(aes(xintercept = c(22)), linetype = "dotted", color = "chocolate1", size = 1)+
  geom_vline(aes(xintercept = c(27)), linetype = "dotted", color = "chocolate1", size = 1)+
  geom_area(aes(y = filter, fill = Signature))+
  geom_point(aes(y = Probability, group = direction_station, shape = Boat_Signature), fill = "darkcyan", size = 3, color = "black")+
  scale_fill_manual(values = colors)+
  scale_shape_manual(values = c(21))+
  ylab("Probability of being sourced from Baltimore")+
  xlab("Baltimore --> Norfolk --> Baltimore")+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.title.x = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 15, angle = 0, vjust = 1, hjust = 0.5, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0))+
  guides(fill = guide_legend(ncol=1,byrow=TRUE))


head(results2)
colnames(results2)

results_stats <- results2 %>%
  filter(sample_desc == "Bilge Biofilm") %>%
  filter(Predicted_Port == "Baltimore") %>%
  group_by(Model) %>%
  summarise(mean_Prob = mean(Probability, na.rm = TRUE))
head(results_stats)
View(results_stats)


#Quantifying the difference between the water signature and the boat signature

head(results)

qSig <- results %>%
  pivot_wider(names_from = sample_type, values_from = Probability) %>%
  mutate(
    surface_sig = surface - filter,
    bilge_sig = bilge - filter
  ) %>%
  select(-c("bilge", "filter", "surface")) %>%
  pivot_longer(cols = c("surface_sig", "bilge_sig"), names_to = "Signature", values_to = "Difference") %>%
  separate(Signature, into = c("Signature", NA), sep = "_")
head(qSig)

colors <- c("chocolate1", "darkcyan")

ggplot(qSig, aes(x = direction_station, y = Difference, fill = Predicted_Port))+
  facet_nested(rows = vars(boat_name, Model), cols = vars(Signature, Predicted_Port))+
  geom_col(color = "black")+
  scale_fill_manual(values = colors) +
  #xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        legend.text = element_text(size = 15),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=1,byrow=TRUE))

#install.packages("ggTimeSeries")
#library(ggTimeSeries)


#ggplot_waterfall(dtData = qSig, 'direction_station', 'Difference')+
 # facet_nested(rows = vars(boat_name, Model), cols = vars(Predicted_Port, Signature), scales = "free", space = "free")+
  #scale_fill_manual(values = colors) +
  #xlab(NULL)+
#  theme_linedraw()+
 # theme(axis.text.y = element_text(size = 15, color = "black"),
  #      axis.title.y = element_text(size = 15, color = "black"),
   #     axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 0.5, color = "black"),
    #    legend.text = element_text(size = 15),
     #   legend.position = "bottom",
      #  legend.spacing.x = unit(0.1, 'mm'),
       # legend.spacing.y = unit(0.05, 'mm'),
        #plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
#        strip.text = element_text(size = 18, face = "bold", angle = 0),
 #       legend.title = element_text(face="bold", size = 22))+
  #guides(fill=guide_legend(ncol=1,byrow=TRUE))

#ggplot(qSig, aes(x = direction_station, y = Difference, group = Predicted_Port, fill = Predicted_Port)) +
 # stat_steamgraph()


#head(qSig)


#results <- nscPredict %>%
 # cbind(rfPredict_nonzv) %>%
  #select(-c("SampleID3")) %>%
#  left_join(sdata) %>%
 # filter(!(is.na(true_geographic_order))) %>%
  #mutate(
   # direction = ifelse(transit == "t1" | transit == "t3", "Departing", "Returning"),
    #direction_abbr = ifelse(direction == "Departing", "D", "R")) %>%
#  unite(direction_station, direction_abbr, station, sep = "", remove = FALSE) %>%
 # mutate(direction_station = factor(direction_station, levels = c("D1", "D2", "D3", "D4", "D5", "D6", "D8", "D10", "D11", "D12", "D13", "D14", "D15", "D16", "D17",
  #                                                                "D18", "D19", "D21", "D22", "D23", "D24", "D25", "D26", "D27", "D28", "R28", "R26", "R24",  "R22",
   #                                                               "R20", "R19", "R17", "R15", "R13", "R12", "R9", "R7"))) %>%
#  filter(sample_description != "blank" & !is.na(sample_description)) %>%
 # mutate(
  #  rfnzv = "rfnzv",
   # nsc = "nsc") %>%
#  unite(nsc_Port, nsc_Port, nsc, sep = "_", remove = TRUE) %>%
 # unite(rfnzv_Port, rfnzv_Port, rfnzv, sep = "_", remove = TRUE) %>%
  #pivot_longer(cols = c("nsc_Port", "rfnzv_Port"), names_to = "Model", values_to = "Predicted_Port") %>%
#  pivot_longer(cols = c("nsc_Prob", "rfnzv_Prob"), names_to = "Probability_Type", values_to = "Probability") %>%
 # separate(Predicted_Port, into = c("Predicted_Port", "Predicted_Port_Model"), sep = "_") %>%
  #separate(Probability_Type, into = c("Probability_Type", NA), sep = "_") %>%
#  filter(Predicted_Port_Model == Probability_Type) %>%
 # group_by(SampleID, direction_station, Predicted_Port, boat_name, sample_type, Model) %>%
  #summarise(
   # Probability = mean(Probability)
#  )
#head(results)


  
qSig2 <- results %>%
  pivot_wider(names_from = sample_type, values_from = Probability) %>%
  mutate(
    surface_sig = surface - filter,
    bilge_sig = bilge - filter
  ) %>%
  select(-c("bilge", "filter", "surface")) %>%
  pivot_longer(cols = c("surface_sig", "bilge_sig"), names_to = "Signature", values_to = "Difference") %>%
  separate(Signature, into = c("Signature", NA), sep = "_") %>%
  filter(Predicted_Port == "Baltimore") %>%
  mutate(Prediction = ifelse(Difference >= 0, "Higher than open water", "Lower than open water"),
         Model_Name = ifelse(Model == "nsc_Port", "Nearest Shrunken Centroids", "Random Forest"),
         Signature = ifelse(Signature == "bilge", "Bilge", "Surface"))
head(qSig2)
#colnames(qSig2)

colors <- c("olivedrab", "firebrick")

ggplot(qSig2, aes(y = direction_station, x = Difference, fill = Prediction))+
  facet_nested(cols = vars(Model_Name, boat_name, Signature))+
  geom_col(color = "black")+
  scale_fill_manual(values = colors) +
  ylab("Sampling Station\n(D = Departing Voyage, R = Returning Voyage")+
  xlab("Difference in Baltimore signature between boat samples and open water")+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 20, face = "bold", color = "black"),
        axis.text.x = element_text(size = 15, angle = 0, hjust = 1, vjust = 0.5, color = "black"),
        axis.title.x = element_text(size = 20, face = "bold", color = "black"),
        legend.text = element_text(size = 15),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", color = "black", angle = 0),
        strip.background = element_rect(fill = "white"),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=2,byrow=TRUE))


## What is the signature?

# from deseq2
deseq_out <- read_csv("/home/lgschaer/old/Chesapeake/March_2023_Analysis/deseq_out/enriched_w_taxonomy_03222023.csv")
head(deseq_out)
colnames(deseq_out)

# from ML models
nscTop20
rfTop20

signature_list <- as.data.frame(nscTop20) %>%
  mutate(ta1 = ifelse(nscTop20 %in% rfTop20, nscTop20, NA),
         taxon = ta1) %>%
  filter(!is.na(ta1)) %>%
  left_join(deseq_out, by = "taxon") %>%
  mutate(Where_Enriched = case_when(
    (log2FoldChange > 2) ~ "Norfolk",
    (log2FoldChange < -2) ~ "Baltimore")) %>%
  filter(threshold == "Enriched") %>%
  select(ta1, Phylum, Class, Genus, Species, Where_Enriched, log2FoldChange, padj, threshold)
signature_list
dim(signature_list)

#load phyloseq
pps_rare <- readRDS("/home/lgschaer/old/Chesapeake/March_2023_Analysis/phyloseq_out/ps_rarefied_nochloro_justbacteria_06132023.rds")
pps_rare

new_taxtab <- tax_table(pps_rare) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ta1") %>%
  mutate(SampleID = ta1) %>%
  column_to_rownames(var = "SampleID") %>%
  select(ta1)
head(new_taxtab)
dim(new_taxtab)

pps_newtax <- merge_phyloseq(pps_rare, tax_table(new_taxtab))
head(tax_table(pps_newtax))

#subset_pps_newtax <- subset_taxa(pps_newtax, ta1 %in% signature_list$taxon)
#subset_pps_newtax
#tax_table(subset_pps_newtax)


# Make abundance plots of the port microbial signatures
signature_abundance <- pps_newtax %>%
  #subset_taxa(ta1 %in% signature_list$taxon) %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(ta1 %in% signature_list$ta1) %>%
  select(-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
head(signature_abundance)
dim(signature_abundance)

# should be 13 ASVs
length(unique(signature_abundance$ta1))


head(signature_abundance)
View(signature_list)

unique(signature_list$ta1)
unique(signature_abundance$ta1)

filt_sig_abundance <- signature_abundance %>%
  select(ta1, Sample, location_sampled, transit, station, sample_description, sample_type, sample_detail_description, boat_name,
         location_description, location, true_geographic_order, km_to_baltimore, Abundance) %>%
  left_join(signature_list, by = c("ta1")) %>%
  filter(Abundance > 0) %>%
  mutate(
    Where_Enriched = factor(Where_Enriched, levels = c("Baltimore", "Norfolk")),
    as.character(true_geographic_order),
    direction = ifelse(transit == "t1" | transit == "t3", "Departing", "Returning"),
    direction_abbr = ifelse(direction == "Departing", "D", "R")) %>%
  unite(direction_station, direction_abbr, station, sep = "", remove = FALSE) %>%
  mutate(direction_station = factor(direction_station, levels = c("D1", "D2", "D3", "D4", "D5", "D6", "D8", "D10", "D11", "D12", "D13", "D14", "D15", "D16", "D17",
                                                                  "D18", "D19", "D21", "D22", "D23", "D24", "D25", "D26", "D27", "D28", "R28", "R26", "R24",  "R22",
                                                                  "R20", "R19", "R17", "R15", "R13", "R12", "R9", "R7")),
         direction_station_int = as.integer(direction_station)) %>%
  filter(!is.na(sample_description) & sample_description != "blank") %>%
  filter(!is.na(direction_station)) %>%
  group_by(Sample, Where_Enriched, Genus, Species, location_sampled, sample_description, sample_type, location, direction_station) %>%
  mutate(totalSum = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Genus = ifelse(is.na(Genus), "", Genus)) %>%
  unite(Tax_Label, ta1, Class, Genus, sep = " ", remove = FALSE) %>%
  mutate(Tax_Label = factor(Tax_Label, levels = c("sp11 Acidimicrobiia CL500-29 marine group", "sp17 Acidimicrobiia CL500-29 marine group",
                                                  "sp22 Cyanobacteriia Cyanobium PCC-6307",    "sp38 Actinobacteria hgcI clade",           
                                                  "sp39 Cyanobacteriia Cyanobium PCC-6307",    "sp89 Cyanobacteriia Cyanobium PCC-6307",   
                                                  "sp1 Cyanobacteriia Cyanobium PCC-6307",     "sp12 Acidimicrobiia ",                     
                                                  "sp150 Cyanobacteriia Synechococcus CC9902", "sp24 Planctomycetes ",                     
                                                  "sp3 Cyanobacteriia Cyanobium PCC-6307",     "sp37 Planctomycetes ",                     
                                                  "sp65 Planctomycetes "))) %>%
  dplyr::group_by(Where_Enriched, Tax_Label, location_sampled, sample_description, sample_type, direction_station, ta1) %>%
  summarise(Abundance = sum(Abundance)) %>%
  arrange(Where_Enriched)
head(filt_sig_abundance)
#View(filt_sig_abundance)
colnames(filt_sig_abundance)

unique(filt_sig_abundance$Tax_Label)

colFunc1 <- colorRampPalette(c("lightblue", "dodgerblue", "cyan", "darkcyan"))
color_list1 <- colFunc1(6)
colFunc2 <- colorRampPalette(c("lightgoldenrod", "orange", "chocolate1", "firebrick"))
color_list2 <- colFunc2(7) #length(unique(filt_sig_abundance$ta1))
colors <- c(color_list1, color_list2)
#length(colors)


ggplot(filt_sig_abundance)+
  geom_col(mapping = aes(x = Where_Enriched, y = Abundance, fill = Tax_Label), color = "black", position = "stack", show.legend = TRUE)+
  facet_nested(rows = vars(location_sampled, sample_type, sample_description),
               cols = vars(direction_station))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors)+   
  scale_x_discrete(limits = c("Baltimore", "Norfolk"), labels = c("Baltimore" = "B", "Norfolk" = "N"))+
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 18, color = "black", face = "bold", angle = 0),
        strip.text.y = element_text(size = 18, color = "black", face = "bold", angle = 270),
        strip.background.x = element_rect(fill = "white"),
        strip.background.y = element_rect(fill = "white"),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=3,byrow=FALSE))

filt_sig_abundance2 <- filt_sig_abundance %>%
  group_by(Where_Enriched, direction_station, location_sampled, sample_description, sample_type) %>%
  mutate(
    sAbundance = sum(Abundance),
    direction_station_int = as.integer(direction_station)
  )
colnames(filt_sig_abundance2)

area_colors <- c("aliceblue", "lemonchiffon")
point_colors <- c("lightblue", "firebrick", "olivedrab4", "maroon3", "cyan", 
                  "darkcyan", "orange", "purple", "darkgreen", "tan",
                  "darkblue", "tomato", "darkgoldenrod", "pink", "green")

ggplot(filt_sig_abundance2, aes(x = direction_station_int))+
  geom_area(aes(x = direction_station_int, y = sAbundance, fill = Where_Enriched), color = "black", show.legend = FALSE)+
  geom_point(aes(x = direction_station_int, y = Abundance, fill = Tax_Label), color = "black", shape = 21, size = 4, show.legend = TRUE)+
  facet_nested(rows = vars(location_sampled, sample_type, sample_description), cols = vars(Where_Enriched))+
  ylab("Proportion of Community") +
  scale_fill_manual(values = c(area_colors, point_colors))+   
  xlab(NULL)+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 25, color = "black"),
        axis.text.x = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.x = element_text(size = 25, color = "black", face = "bold", angle = 0),
        strip.text.y = element_text(size = 25, color = "black", face = "bold", angle = 270),
        strip.background.x = element_rect(fill = "white"),
        strip.background.y = element_rect(fill = "white"),
        legend.title = element_text(face="bold", size = 22))+
  guides(fill=guide_legend(ncol=3,byrow=FALSE, title = ""))


### Trying to predict continuous variable: distance from Baltimore THIS DID NOT GIVE US ACCURATE MODELS


# format data for machine learning

ps_ches <- readRDS("/home/lgschaer/old/Chesapeake/March_2023_Analysis/phyloseq_out/ps_unrarefied_03212023.rds")
ps_ches

ps_subs <- subset_samples(ps_ches, (sample_type != "blank" & (sample_type == "surface" | sample_type == "filter")) & !is.na(sample_description)) %>%
  subset_samples(sample_sums(ps_ches) > 0) %>%
  subset_taxa(taxa_sums(ps_ches) > 0)
ps_subs
head(sample_data(ps_subs))

ches_sd <- as_tibble(sample_data(ps_subs)) %>%
  select(SampleID, km_to_baltimore)
head(ches_sd)

ches_otu <- as.data.frame(otu_table(ps_subs)) %>%
  rownames_to_column(var = "SampleID")
ches_otu[1:5,1:5]

ml_data <- ches_sd %>%
  #select(-location) %>%
  filter(!is.na(km_to_baltimore)) %>%
  left_join(ches_otu) %>%
  column_to_rownames(var = "SampleID")
ml_data[1:3,1:5]
#dim(ml_data)


# subset into training and testing set

split <- 0.80
trainIndex <- createDataPartition(ml_data$km_to_baltimore, p=split, list=FALSE)
data_train <- ml_data[ trainIndex,]
y_train <- data_train$km_to_baltimore
data_test <- ml_data[-trainIndex,]
y_test <- data_test$km_to_baltimore

train <- data_train %>% select(-km_to_baltimore)
test <- data_test %>% select(-km_to_baltimore)

length(y_train)
length(y_test)

## Remove near zero variance variables

nzv <- nearZeroVar(train)
train_nzv <- train[-nzv]
test_nzv <- test[-nzv]
dim(train)
dim(train_nzv)


## Remove highly correlated variables *******worked after I removed nzv predictors
tooHigh <- findCorrelation(cor(train_nzv), cutoff = .75)
trainXnnet <- train_nzv[, -tooHigh]
testXnnet <- test_nzv[, -tooHigh]


# Neural Network

nnetGrid <- expand.grid(.decay = c(.1, 1, 10),
                             .size = c(1:7),
                             .bag = FALSE) #see next chapter for more info about bagging
ctrl <- trainControl(classProbs = TRUE, method = "repeatedcv", number = 3, repeats = 25)
set.seed(100)
nnetTune <- train(train_nzv, y_train,
                    method = "avNNet",
                    tuneGrid = nnetGrid,
                    trControl = ctrl,
                    preProc = c("center", "scale"),
                    linout = TRUE,
                    trace = FALSE,
                    MaxNWts = 10 * (ncol(trainXnnet) + 1) + 10 + 1,
                    maxit = 500)
nnetTune

nnetPredict <- predict(nnetTune, test_nzv) %>%
  cbind(y_test) %>%
  as_tibble()
colnames(nnetPredict) <- c("predicted", "real")
#rownames_to_column(var = "SampleID2") %>%
#pivot_longer(cols = c("Baltimore", "Norfolk"), names_to = "pls_Port", values_to = "pls_Prob")
head(nnetPredict)


ggplot(nnetPredict, aes(x = real, y = predicted))+
  geom_point(shape = 21, fill = "darkcyan", size = 5)+
  geom_line(aes(x = real, y = real), color = "darkred", linetype = "dotted", size = 2)+
  theme_linedraw()

summary(lm(data = nnetPredict, x = nnetPredict$real, y = nnetPredict$predicted))

# Partial Least Squares

set.seed(100)
plsTune <- train(train_nzv, y_train,
                   method = "pls",
                   tuneLength = 20,    #tuning parameter
                   trControl = ctrl,
                   preProc = c("center", "scale"))
plsTune

plsPredict <- predict(plsTune, test_nzv) %>%
  cbind(y_test) %>%
  as_tibble()
colnames(plsPredict) <- c("predicted", "real")
  #rownames_to_column(var = "SampleID2") %>%
  #pivot_longer(cols = c("Baltimore", "Norfolk"), names_to = "pls_Port", values_to = "pls_Prob")
head(plsPredict)


ggplot(plsPredict, aes(x = real, y = predicted))+
  geom_point(shape = 21, fill = "darkcyan", size = 5)+
  geom_line(aes(x = real, y = real), color = "darkred", linetype = "dotted", size = 2)+
  theme_linedraw()

summary(lm(data = plsPredict, x = plsPredict$real, y = plsPredict$predicted))


# KNN

set.seed(100)
knnTune <- train(train_nzv, y_train,
                   method = "knn",
                   preProc = c("center", "scale"),
                   tuneGrid = data.frame(.k = 1:20),
                   trControl = ctrl)
knnTune

knnPredict <- predict(knnTune, test_nzv) %>%
  cbind(y_test) %>%
  as_tibble()
colnames(knnPredict) <- c("predicted", "real")
head(knnPredict)


ggplot(knnPredict, aes(x = real, y = predicted))+
  geom_point(shape = 21, fill = "darkcyan", size = 5)+
  geom_line(aes(x = real, y = real), color = "darkred", linetype = "dotted", size = 2)+
  theme_linedraw()

summary(lm(data = knnPredict, x = knnPredict$real, y = knnPredict$predicted))



# Random Forest Regression

set.seed(100)
rfGrid <- data.frame(.mtry = c(50, 75, 100, 125, 150, 200, 250, 300, 350))
rfRegTune <- train(train_nzv, y_train,
                 method = "rf",
                 preProc = c("center", "scale"),
                 tuneGrid = rfGrid,
                 trControl = ctrl)
rfRegTune

rfRegPredict <- predict(rfRegTune, test_nzv) %>%
  cbind(y_test) %>%
  as_tibble()
colnames(rfRegPredict) <- c("predicted", "real")
head(rfRegPredict)


ggplot(rfRegPredict, aes(x = real, y = predicted))+
  geom_point(shape = 21, fill = "darkcyan", size = 5)+
  geom_line(aes(x = predicted, y = predicted), color = "darkred", linetype = "dotted", size = 2)+
  theme_linedraw()

summary(lm(data = rfRegPredict, x = rfRegPredict$real, y = rfRegPredict$predicted))

