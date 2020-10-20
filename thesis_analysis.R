# Libraries
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(broom)
source("R/readMappingData.R")
library(ggplot2)
library(gridExtra)
library(formattable)
library(openxlsx)
library(stringr)
set.seed(1)

# Data
arg_phages <- readRDS("data/arg_phages.RDS")
arg_plasmids <- readRDS("data/arg_plasmids.RDS")
arg_TEs <- readRDS("data/arg_TEs.RDS")
all_args <- read.csv("data/nonsubsampled_merged_bedtools.csv", stringsAsFactors = FALSE)
all_args$ID <- gsub(".bam.bedtools.coverage.txt", "", gsub("db/MAPPING_DATA/NONSUBSAMPLED/", "", all_args$filename))
all_args_abud <- readMappingData("data/nonsubsampled_merged_bedtools.csv", without_US_duplicates = FALSE)
antibiotic_use <- read.csv("db/resistanceMap_use_190221_mod_200904.csv", stringsAsFactors = FALSE)
card_aro_categories_index <- read.csv("db/CARD_DB/card-data/aro_categories_index.csv", stringsAsFactors = FALSE, sep = "\t")
card_index <- read.csv("db/CARD_DB/card-data/aro_index.csv", stringsAsFactors = FALSE, sep = "\t")
card_metadata <- full_join(card_index, card_aro_categories_index, by = "Protein.Accession")
card_metadata <- card_metadata[!duplicated(card_metadata$ARO.Accession),]
card_metadata_class <- card_metadata %>% select(ARO.Name, Drug.Class)

# Metadata
metadata <- read.csv("db/metadata.csv", stringsAsFactors = FALSE)
metadata_phages <- readRDS("db/metadata_phages.RDS")
metadata_plasmids <- readRDS("db/metadata_plasmids.RDS")
metadata_TEs <- readRDS("db/metadata_TEs.RDS")

# Combine
arg_phages <- arg_phages %>%
  select(-qseqid) %>%
  rename(mge_name = vcontact_cluster_name) %>%
  mutate(mge = "phage") %>%
  as.data.frame()
arg_plasmids <- arg_plasmids %>%
  select(-qseqid) %>%
  rename(mge_name = plasmid_name_cluster) %>%
  mutate(mge = "plasmid") %>%
  as.data.frame()
arg_TEs <- arg_TEs %>%
  select(-qseqid) %>%
  rename(mge_name = itr_cluster) %>%
  mutate(mge_name = as.character(mge_name)) %>%
  mutate(mge = "IS") %>%
  as.data.frame()
arg_mges <- rbind(arg_phages, arg_plasmids, arg_TEs)

# Add ARG class and mechanism
card_aro_categories_index <- read.csv("db/CARD_DB/card-data/aro_categories_index.csv", stringsAsFactors = FALSE, sep = "\t")
card_index <- read.csv("db/CARD_DB/card-data/aro_index.csv", stringsAsFactors = FALSE, sep = "\t")
card_metadata <- full_join(card_index, card_aro_categories_index, by = "Protein.Accession")
card_metadata <- card_metadata[!duplicated(card_metadata$ARO.Accession),]
card_metadata$Drug.Class.alt <- card_metadata$Drug.Class
card_metadata$Drug.Class.alt[card_metadata$Drug.Class.alt == "macrolide antibiotic;lincosamide antibiotic;streptogramin antibiotic"] <- "MLS antibiotic"
card_metadata$Drug.Class.alt[sapply(card_metadata$Drug.Class.alt, function(x) str_count(x, ";")) > 2] <- "multidrug"
card_metadata$Drug.Class.alt[card_metadata$Resistance.Mechanism == "antibiotic efflux"] <- paste(card_metadata$Drug.Class.alt[card_metadata$Resistance.Mechanism == "antibiotic efflux"], "efflux")

# Redo the AMR metadata
arg_mges <- arg_mges %>% 
  select(-c(AMR.Gene.Family, Drug.Class.alt, Resistance.Mechanism)) %>%
  mutate(ARO.Name = strsplit(ARO.Name, ",")) %>%
  unnest() %>%
  left_join(card_metadata)

# Get common IDs
ids <- Reduce(intersect, list(unique(metadata_phages$ID), unique(metadata_plasmids$ID), unique(metadata_TEs$ID)))

# Join metadata and add timepoints
metadata <- metadata %>% group_by(Location, sample_type, Sample.name) %>%
  mutate(timepoint = rank(as.numeric(Visit_Number))) 
metadata_summary <- metadata %>%
  group_by(Location, sample_type, timepoint) %>%
  summarise(n_total = n_distinct(ID))
metadata_all_mges <- metadata[metadata$ID %in% ids,]
arg_all_mges <- inner_join(arg_mges, metadata_all_mges, by = "ID")
arg_mges <- inner_join(arg_mges, metadata, by = "ID")

# Prevalence of ARGs
arg_mges_sing_summary <- arg_all_mges %>%
  filter(timepoint == 1) %>%
  select(Location, sample_type, ID, ARO.Name, Drug.Class.alt, Resistance.Mechanism, AMR.Gene.Family, mge) %>%
  group_by(Location, sample_type, ID, ARO.Name, Drug.Class.alt, Resistance.Mechanism, AMR.Gene.Family) %>%
  unique() %>%
  summarise_all(paste, collapse = " and ") %>%
  ungroup() %>%
  group_by(Location, sample_type, ARO.Name, Drug.Class.alt, Resistance.Mechanism, AMR.Gene.Family, mge) %>%
  summarise(n = n_distinct(ID)) %>%
  inner_join(metadata_summary[metadata_summary$timepoint == 1,]) %>%
  mutate(perc = n/n_total*100) %>%
  mutate(labels = paste0(Drug.Class.alt, " - ", ARO.Name))
arg_mges_sing_summary$labels <- factor(arg_mges_sing_summary$labels, levels = rev(levels(factor(arg_mges_sing_summary$labels))))

# Plot ARG prevalence
tiff("figures/arg_mge_prevalence.tiff", width = 4000, height = 5750, res = 250)
ggplot(arg_mges_sing_summary, aes(labels, perc, fill = mge)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ sample_type + Location, scale = "free", space = "free", switch = "both") +
  theme_bw() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        strip.text.y = element_text(angle = 180, size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position="bottom") +
  ylab("% samples") + xlab("ARG Class - ARG") +
  scale_fill_manual("MGE", values = brewer.pal(length(unique(arg_mges_sing_summary$mge)), "Set2")) +
  scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100))
dev.off()  

# Prevalence of ARG classes with MGEs and antibiotic use
arg_mges_use <- arg_mges %>%
  filter(timepoint == 1) %>%
  select(Location, sample_type, ID, ARO.Name, Drug.Class.alt, Resistance.Mechanism, AMR.Gene.Family, mge) %>%
  mutate(ARO.Name = c(strsplit(ARO.Name, ","))) %>%
  unnest(ARO.Name) %>%
  inner_join(card_metadata_class) %>%
  ungroup() %>%
  mutate(Drug.Class = c(strsplit(Drug.Class, ";"))) %>%
  unnest(Drug.Class) %>%
  group_by(Location, sample_type, Drug.Class, mge) %>%
  summarise(n = n_distinct(ID)) %>%
  inner_join(metadata_summary[metadata_summary$timepoint == 1,]) %>%
  mutate(perc = n/n_total*100) %>%
  inner_join(antibiotic_use, c("Location", "Drug.Class"="CARD.Class")) %>%
  mutate(Location_sampletype = paste(sample_type, "-", Location))

arg_mges_use_summary <- arg_mges_use %>%
  ungroup() %>%
  select(mge, Location, sample_type, n_total) %>%
  unique()

# Add linear model
arg_mges_use_lm <- arg_mges_use %>%
  group_by(Location, sample_type, mge) %>%
  do(mod = lm(perc ~ DDD.Per.1000.Pop, data = .))

class_cols <- c("black", "grey", brewer.pal(length(unique(arg_mges_use$Drug.Class))-2, "Paired"))
names(class_cols) <- unique(arg_mges_use$Drug.Class)

# IS
is_list <- list()
unique_location_sampletype <- unique(arg_mges_use$Location_sampletype[arg_mges_use$mge == "IS"])
for (i in 1:length(unique_location_sampletype)) {
  is_list[[i]] <- ggplot(arg_mges_use[arg_mges_use$mge == "IS" & arg_mges_use$Location_sampletype == unique_location_sampletype[i],], aes(DDD.Per.1000.Pop, perc, colour = Drug.Class)) +
    geom_point(shape = 4, size = 2, stroke = 1.5) +
    theme_bw() +
    ggtitle(unique_location_sampletype[i]) +
    scale_colour_manual("Antibiotic Class", values = class_cols) +
    theme(title = element_text(size = 14), legend.position = "none") +
    ylab("% samples") + xlab("DDD Per 1000")
}

tiff("figures/antibiotic_use_is_prevalence.tiff", width = 3000, height = 2250, res = 200)
grid.arrange(grobs = is_list, ncol = 4)
dev.off()

# Plasmids
plasmid_list <- list()
unique_location_sampletype <- unique(arg_mges_use$Location_sampletype[arg_mges_use$mge == "plasmid"])
for (i in 1:length(unique_location_sampletype)) {
  plasmid_list[[i]] <- ggplot(arg_mges_use[arg_mges_use$mge == "plasmid" & arg_mges_use$Location_sampletype == unique_location_sampletype[i],], aes(DDD.Per.1000.Pop, perc, colour = Drug.Class)) +
    geom_point(shape = 4, size = 2, stroke = 1.5) +
    theme_bw() +
    ggtitle(unique_location_sampletype[i]) +
    scale_colour_manual("Antibiotic Class", values = class_cols) +
    theme(title = element_text(size = 14), legend.position = "none") +
    ylab("% samples") + xlab("DDD Per 1000")
}

tiff("figures/antibiotic_use_plasmid_prevalence.tiff", width = 3000, height = 1500, res = 200)
grid.arrange(grobs = plasmid_list, nrow = 2, ncol = 4)
dev.off()

# Phages
phage_list <- list()
unique_location_sampletype <- unique(arg_mges_use$Location_sampletype[arg_mges_use$mge == "phage"])
for (i in 1:length(unique_location_sampletype)) {
  phage_list[[i]] <- ggplot(arg_mges_use[arg_mges_use$mge == "plasmid" & arg_mges_use$Location_sampletype == unique_location_sampletype[i],], aes(DDD.Per.1000.Pop, perc, colour = Drug.Class)) +
    geom_point(shape = 4, size = 2, stroke = 1.5) +
    theme_bw() +
    ggtitle(unique_location_sampletype[i]) +
    scale_colour_manual("Antibiotic Class", values = class_cols) +
    theme(title = element_text(size = 14), legend.position = "none") +
    ylab("% samples") + xlab("DDD Per 1000")
}

tiff("figures/antibiotic_use_phage_prevalence.tiff", width = 3000, height = 1500, res = 200)
grid.arrange(grobs = phage_list, nrow = 2, ncol = 4)
dev.off()

# Legend for antibiotic classes
tiff("figures/antibiotic_class_legend.tiff", width = 2000, height = 750, res = 150)
ggplot(arg_mges_use, aes(DDD.Per.1000.Pop, perc, colour = Drug.Class)) +
  geom_point(shape = 4, size = 2, stroke = 1.5) +
  theme_bw() +
  scale_colour_manual("Antibiotic Class", values = class_cols) +
  theme(title = element_text(size = 14), 
        legend.text = element_text(size = 12), 
        legend.position = "bottom") +
  ylab("% samples") + xlab("DDD Per 1000")
dev.off()

# Prevalnce of ARGs vs prevalence of MGEs
arg_prev <- all_args_abud %>%
  filter(ID %in% metadata$ID[metadata$timepoint == 1]) %>%
  group_by(Location, sample_type, ARO.Name, Drug.Class) %>%
  summarise(n_samples_args = n_distinct(ID)) 
  
arg_mge_prev <- arg_all_mges %>%
  filter(timepoint == 1) %>%
  select(Location, sample_type, ID, ARO.Name, Drug.Class.alt, Resistance.Mechanism, AMR.Gene.Family, mge) %>%
  group_by(Location, sample_type, ARO.Name, Drug.Class.alt, Resistance.Mechanism, AMR.Gene.Family, mge) %>%
  summarise(n = n_distinct(ID)) %>%
  inner_join(arg_prev) %>%
  group_by(Location, sample_type, ARO.Name) %>%
  pivot_wider(names_from = mge, values_from = n, values_fill = list(n = 0))

lm_prev <- arg_mge_prev %>%
  group_by(Location, sample_type) %>% 
  do(mod = lm(n_samples_args ~ IS + plasmid + phage, data = .)) 

# Risk score
arg_mge_risk <- arg_all_mges %>%
  filter(timepoint == 1) %>%
  select(Location, sample_type, ID, ARO.Name, mge) %>%
  group_by(Location, sample_type, ARO.Name, mge) %>%
  summarise(n = n_distinct(ID)) %>%
  right_join(arg_prev) %>%
  inner_join(metadata_summary[metadata_summary$timepoint == 1,]) %>%
  select(Location, sample_type, ARO.Name, Drug.Class, mge, n_total, n_samples_args, n) %>%
  rename(n_samples_args_mges = n) %>%
  mutate(mge = replace(mge, is.na(mge), "IS,plasmid,phage")) %>%
  mutate(mge = c(strsplit(mge, ","))) %>%
  unnest(mge) %>%
  mutate(n_samples_args_mges = replace(n_samples_args_mges, is.na(n_samples_args_mges), 0)) %>%
  mutate(perc_samples_args = n_samples_args/n_total*100, perc_samples_args_mges = n_samples_args_mges/n_total*100) %>%
  mutate(Drug.Class = c(strsplit(Drug.Class, ";"))) %>%
  unnest(Drug.Class) %>%
  pivot_wider(names_from = mge, values_from = perc_samples_args_mges, values_fill = list(mge_incidence = 0)) %>%
  mutate(plasmid = replace(plasmid, is.na(plasmid), 0), IS = replace(IS, is.na(IS), 0), phage = replace(phage, is.na(phage), 0)) %>%
  mutate(sum_perc_mge = plasmid + IS + phage) %>%
  group_by(Location, sample_type, Drug.Class) %>%
  mutate(mge_rank = rank(-sum_perc_mge, ties.method = "min")) %>%
  group_by(Location, sample_type, Drug.Class, sum_perc_mge) %>%
  mutate(arg_rank = rank(-perc_samples_args, ties.method = "min")) %>%
  mutate(arg_rank = replace(arg_rank, sum_perc_mge != 0, 1)) %>%
  mutate(rank = mge_rank + (arg_rank - 1)) %>%
  ungroup() %>%
  mutate(perc_samples_args = signif(perc_samples_args, 3), IS = signif(IS, 3), plasmid = signif(plasmid, 3), phage = signif(phage, 3))

# Top three in each category
collapse_rows_df <- function(df, variable){
  
  group_var <- enquo(variable)
  
  df %>%
    group_by(!! group_var) %>%
    mutate(groupRow = 1:n()) %>%
    ungroup() %>%
    mutate(!!quo_name(group_var) := ifelse(groupRow == 1, as.character(!! group_var), "")) %>%
    select(-c(groupRow))
}

arg_mge_risk_top <- arg_mge_risk %>%
  filter(rank == 1) %>%
  filter(sum_perc_mge != 0) %>%
  select(Drug.Class, Location, sample_type, ARO.Name, perc_samples_args, IS, plasmid, phage) %>%
  arrange(Drug.Class) %>%
  collapse_rows_df(Drug.Class) 
names(arg_mge_risk_top) <- c("Antibiotic Class", "Country", "GIT Site", "ARG", "Prevalence of ARGs (%)", "Prevalence of ARG-carrying ISs (%)", "Prevalence of ARG-carrying plasmids (%)", "Prevalence of ARG-carrying phages (%)")
write.xlsx(arg_mge_risk_top, file = "data/arg_mge_risk_top.xlsx")

# Make supplementary tables
arg_mge_risk_sup <- arg_mge_risk %>%
  select(Drug.Class, Location, sample_type, ARO.Name, rank, perc_samples_args, IS, plasmid, phage) %>%
  group_by(Drug.Class, Location, sample_type) %>%
  arrange(rank, .by_group = TRUE)
names(arg_mge_risk_sup) <- c("Antibiotic Class", "Country", "GIT Site", "ARG", "Rank", "Prevalence of ARGs (%)", "Prevalence of ARG-carrying ISs (%)", "Prevalence of ARG-carrying plasmids (%)", "Prevalence of ARG-carrying phages (%)")
  
sheet_names <- sort(unique(arg_mge_risk_sup$"Antibiotic Class"))

wb <- createWorkbook()
for (i in 1:length(sheet_names)){
  addWorksheet(wb, sheet_names[i])
  arg_mge_risk_sup_tmp <- arg_mge_risk_sup[arg_mge_risk_sup$"Antibiotic Class" == sheet_names[i],] %>%
    select(-"Antibiotic Class") 
  writeData(wb, sheet = sheet_names[i], arg_mge_risk_sup_tmp)
}
saveWorkbook(wb, "data/arg_mge_risk_supplementary.xlsx", overwrite = TRUE)

# Abundance of ARGs vs incidence of MGE
arg_mge_inc <- arg_mges %>%
  filter(timepoint == 1) %>%
  select(Location, sample_type, ID, ARO.Name, Drug.Class.alt, Resistance.Mechanism, AMR.Gene.Family, mge_name) %>%
  group_by(ID, Location, sample_type, ARO.Name, Drug.Class.alt, Resistance.Mechanism, AMR.Gene.Family, mge_name) %>%
  summarise(mge_incidence = n())
  
arg_abund <- all_args_abud %>%
  filter(ID %in% metadata$ID[metadata$timepoint == 1]) %>%
  group_by(Location, sample_type, ARO.Name) %>%
  mutate(n_samples_args = n_distinct(ID))  %>%
  filter(n_samples_args >= 10) %>%
  left_join(arg_mge_inc) %>%
  mutate(mge_incidence=replace(mge_incidence, !is.na(mge_incidence), 1)) %>%
  mutate(mge_incidence=replace(mge_incidence, is.na(mge_incidence), 0))

arg_no_mge_abund <- arg_abund[is.na(arg_abund$mge_name),]
arg_mge_abund <- arg_abund[!is.na(arg_abund$mge_name),]

arg_mge_num <- arg_mge_abund %>%
  group_by(Location, sample_type, ARO.Name, mge_name) %>%
  summarise(n = n_distinct(ID)) %>%
  filter(n > 10) %>%
  select(-n)

arg_no_mge_num <- arg_no_mge_abund %>%
  group_by(Location, sample_type, ARO.Name, mge_name) %>%
  summarise(n = n_distinct(ID)) %>%
  filter(n > 10) %>%
  select(-c(n, mge_name))

arg_no_mge_abund$mge_name <- paste(unique(arg_mge_abund$mge_name), collapse = ";")
arg_no_mge_abund <- arg_no_mge_abund %>%
  mutate(mge_name = strsplit(mge_name, ";")) %>%
  unnest() %>%
  inner_join(arg_mge_num) %>%
  inner_join(arg_no_mge_num) 

arg_mge_abund <- arg_abund %>%
  inner_join(arg_mge_num) %>%
  inner_join(arg_no_mge_num) %>%
  rbind(arg_no_mge_abund)

# Chi square test
arg_mge_chi <- arg_mge_abund %>%
  group_by(Location, sample_type, ARO.Name, mge_name) %>%
  do(pval_chi = as.numeric(anova(glm(mge_incidence~1, data = ., family = 'binomial'), glm(mge_incidence ~ rpkm, data = ., family = "binomial"), test = "Chisq")$'Pr(>Chi)'[2])) %>%
  mutate(pval_chi = unlist(pval_chi))

arg_mge_glm <- arg_mge_abund %>%
  group_by(Location, sample_type, ARO.Name) %>% 
  do(mod = glm(mge_incidence ~ rpkm, data = .)) %>%
  inner_join(arg_mge_chi) 

# Set coef to 0 if no significant
arg_mge_glm <- arg_mge_glm[arg_mge_glm$pval_chi < 0.05,]
arg_mge_glm$pval <- unlist(lapply(arg_mge_glm$mod, function(x) summary(x)$coefficient[8]))
arg_mge_glm$se_mge <- unlist(lapply(arg_mge_glm$mod, function(x) summary(x)$coefficient[4]))
arg_mge_glm <- arg_mge_glm[arg_mge_glm$pval < 0.05,]

arg_mge_glm_filter <- arg_mge_glm %>%
  select(-mod) %>%
  inner_join(arg_mge_abund) %>%
  mutate(mge = ifelse(grepl("NODE", mge_name) | grepl("PC", mge_name), "Plasmid", "IS")) %>%
  mutate(arg_mge = paste0(ARO.Name, "\n", mge, ": ", mge_name)) %>%
  mutate(Location_sampletype = paste(sample_type, "-", Location))

# Cohort colours
cohort_cols <- c("grey", brewer.pal(9, "Blues")[c(5,7)], "seagreen3", brewer.pal(9, "YlOrRd")[c(3,5,7,9)], brewer.pal(9, "RdPu")[c(3,5,7,9)])
names(cohort_cols) <- sort(unique(paste(metadata$sample_type, "-", metadata$Location)))

# Plot
unique_arg_mges <- unique(arg_mge_glm_filter$arg_mge)
glm_mge <- list()
for (i in 1:length(unique_arg_mges)) {
  glm_mge[[i]] <- ggplot(arg_mge_glm_filter[arg_mge_glm_filter$arg_mge == unique_arg_mges[i],], aes(rpkm, mge_incidence, colour = Location_sampletype)) +
    geom_point(shape = 1) +
    geom_smooth(method = "glm", se = FALSE, method.args = list(family = "binomial")) +
    theme_bw() +
    scale_colour_manual("GIT Site - Country", values = cohort_cols[names(cohort_cols) %in% arg_mge_glm_filter$Location_sampletype[arg_mge_glm_filter$arg_mge == unique_arg_mges[i]]]) +
    ggtitle(unique_arg_mges[i]) +
    xlab("ARG abundance (RPKM)") +
    ylab("Incidence of MGE") +
    scale_y_continuous(breaks = c(0,1)) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), title = element_text(size = 10)) +
    theme(plot.title = element_text(), legend.position = "none")
}

tiff("figures/glm_mges.tiff", width = 2500, height = 3000, res = 250)
grid.arrange(grobs = glm_mge, ncol = 2)
dev.off()

# Plot legend
tiff("figures/glm_legend.tiff", width = 1000, height = 1000, res = 180)
ggplot(arg_mge_glm_filter, aes(rpkm, mge_incidence, colour = Location_sampletype)) +
  geom_point(shape = 1) +
  geom_smooth(method = "glm", se = FALSE, method.args = list(family = "binomial")) +
  theme_bw() +
  scale_colour_manual("GIT Site - Country", values = cohort_cols) +
  xlab("") +
  ylab("") +
  scale_y_continuous(breaks = c(0,1)) +
  theme(plot.title = element_text())
dev.off()

# # Cohort colours
# cohort_cols <- c("grey", brewer.pal(9, "Blues")[c(5,7)], "seagreen3", brewer.pal(9, "YlOrRd")[c(3,5,7,9)], brewer.pal(9, "RdPu")[c(3,5,7,9)])
# names(cohort_cols) <- sort(unique(arg_mge_abund$Location_sampletype)) 
# 
# # ARG abundance vs ISs
# arg_IS_chi <- arg_mge_abund %>%
#   group_by(Location_sampletype, ARO.Name) %>%
#   filter(sum(IS) >= 10 & sum(IS == 0) >= 10) %>%
#   do(pval_chi = as.numeric(anova(glm(IS~1, data = ., family = 'binomial'), glm(IS ~ rpkm, data = ., family = "binomial"), test = "Chisq")$'Pr(>Chi)'[2])) %>%
#   mutate(pval_chi = unlist(pval_chi))
# arg_IS_glm <- arg_mge_abund %>%
#   group_by(Location_sampletype, ARO.Name) %>% 
#   do(mod = glm(IS ~ rpkm, data = .)) %>%
#   inner_join(arg_IS_chi) 
# arg_IS_glm$pval <- unlist(lapply(arg_IS_glm$mod, function(x) summary(x)$coefficient[8]))
# arg_IS_glm$se_IS <- unlist(lapply(arg_IS_glm$mod, function(x) summary(x)$coefficient[4]))
# 
# # Set coef to 0 if no significant
# arg_IS_glm <- arg_IS_glm[arg_IS_glm$pval < 0.05 & arg_IS_glm$pval_chi < 0.05,]
# 
# # Plot
# arg_IS_glm_filter <- arg_IS_glm %>%
#   inner_join(arg_mge_abund)
# unique_args <- unique(arg_IS_glm_filter$ARO.Name)
# unique_args_titles <- unique_args
# unique_args_titles[unique_args_titles == "Escherichia coli ampC beta-lactamase"] <- "Escherichia coli\nampC beta-lactamase"
# unique_args_titles[unique_args_titles == "Campylobacter coli chloramphenicol acetyltransferase"] <- "Campylobacter coli\nchloramphenicol acetyltransferase"
# glm_IS <- list()
# for (i in 1:length(unique_args)) {
# glm_IS[[i]] <- ggplot(arg_IS_glm_filter[arg_IS_glm_filter$ARO.Name == unique_args[i],], aes(rpkm, IS, colour = Location_sampletype)) +
#   geom_point(shape = 1) +
#   geom_smooth(method = "glm", se = FALSE, method.args = list(family = "binomial")) +
#   theme_bw() +
#   scale_colour_manual("GIT Site - Country", values = cohort_cols[names(cohort_cols) %in% arg_IS_glm_filter$Location_sampletype[arg_IS_glm_filter$ARO.Name == unique_args[i]]]) +
#   ggtitle(paste0("IS - ", unique_args_titles[i])) +
#   xlab("ARG abundance (RPKM)") +
#   ylab("Incidence of IS") +
#   scale_y_continuous(breaks = c(0,1)) +
#   theme(plot.title = element_text(), legend.position = "none")
# }
# 
# # tiff("figures/glm_ISs.tiff", width = 1000, height = 1000, res = 220)
# # grid.arrange(grobs = glm_IS, nrow = 1)
# # dev.off()
# 
# # Plot legend
# tiff("figures/glm_legend.tiff", width = 1000, height = 1000, res = 180)
# ggplot(arg_IS_glm_filter, aes(rpkm, IS, colour = Location_sampletype)) +
#   geom_point(shape = 1) +
#   geom_smooth(method = "glm", se = FALSE, method.args = list(family = "binomial")) +
#   theme_bw() +
#   scale_colour_manual("GIT Site - Country", values = cohort_cols) +
#   xlab("") +
#   ylab("") +
#   scale_y_continuous(breaks = c(0,1)) +
#   theme(plot.title = element_text())
# dev.off()
# 
# # ARG abundance vs plasmids
# arg_plasmid_chi <- arg_mge_abund %>%
#   group_by(Location_sampletype, ARO.Name) %>%
#   filter(sum(plasmid) >= 10 & sum(plasmid == 0) >= 10) %>%
#   do(pval_chi = as.numeric(anova(glm(plasmid~1, data = ., family = 'binomial'), glm(plasmid ~ rpkm, data = ., family = "binomial"), test = "Chisq")$'Pr(>Chi)'[2])) %>%
#   mutate(pval_chi = unlist(pval_chi))
# arg_plasmid_glm <- arg_mge_abund %>%
#   group_by(Location_sampletype, ARO.Name) %>% 
#   do(mod = glm(plasmid ~ rpkm, data = .)) %>%
#   inner_join(arg_plasmid_chi) 
# arg_plasmid_glm$pval <- unlist(lapply(arg_plasmid_glm$mod, function(x) summary(x)$coefficient[8]))
# arg_plasmid_glm$se_plasmid <- unlist(lapply(arg_plasmid_glm$mod, function(x) summary(x)$coefficient[4]))
# 
# # Set coef to 0 if no significant
# arg_plasmid_glm <- arg_plasmid_glm[arg_plasmid_glm$pval < 0.05 & arg_plasmid_glm$pval_chi < 0.05,]
# 
# # Plot
# arg_plasmid_glm_filter <- arg_plasmid_glm %>%
#   inner_join(arg_mge_abund)
# unique_args <- unique(arg_plasmid_glm_filter$ARO.Name)
# glm_plasmid <- list()
# for (i in 1:length(unique_args)) {
#   glm_plasmid[[i]] <- ggplot(arg_plasmid_glm_filter[arg_plasmid_glm_filter$ARO.Name == unique_args[i],], aes(rpkm, plasmid, colour = Location_sampletype)) +
#     geom_point(shape = 1) +
#     geom_smooth(method = "glm", se = FALSE, method.args = list(family = "binomial")) +
#     theme_bw() +
#     scale_colour_manual("GIT Site - Country", values = cohort_cols[names(cohort_cols) %in% arg_plasmid_glm_filter$Location_sampletype[arg_plasmid_glm_filter$ARO.Name == unique_args[i]]]) +
#     ggtitle(paste0("Plasmid - ", unique_args[i])) +
#     scale_y_continuous(breaks = c(0,1)) +
#     xlab("ARG abundance (RPKM)") +
#     ylab("Incidence of plasmid") +
#     theme(plot.title = element_text(), legend.position = "none")
# }
# 
# tiff("figures/glm_mges.tiff", width = 4700, height = 3000, res = 200)
# grid.arrange(grobs = c(glm_IS, glm_plasmid), nrow = 3)
# dev.off()

# # ARG abundance vs phages - No phage-ARGs in 5 or greater samples
# arg_phage_chi <- arg_mge_abund %>%
#   group_by(Location_sampletype, ARO.Name) %>%
#   filter(sum(phage) >= 5 & sum(phage == 0) >= 5) %>%
#   do(pval_chi = as.numeric(anova(glm(phage~1, data = ., family = 'binomial'), glm(phage ~ rpkm, data = ., family = "binomial"), test = "Chisq")$'Pr(>Chi)'[2])) %>%
#   mutate(pval_chi = unlist(pval_chi))
# arg_phage_glm <- arg_mge_abund %>%
#   group_by(Location_sampletype, ARO.Name) %>%
#   do(mod = glm(phage ~ rpkm, data = .)) %>%
#   inner_join(arg_phage_chi)
# arg_phage_glm$pval <- unlist(lapply(arg_phage_glm$mod, function(x) summary(x)$coefficient[8]))
# arg_phage_glm$se_phage <- unlist(lapply(arg_phage_glm$mod, function(x) summary(x)$coefficient[4]))
# 
# # Set coef to 0 if no significant
# arg_phage_glm <- arg_phage_glm[arg_phage_glm$pval < 0.05 & arg_phage_glm$pval_chi < 0.05,]
# 
# # Plot
# arg_phage_glm_filter <- arg_phage_glm %>%
#   inner_join(arg_mge_abund)
# unique_args <- unique(arg_phage_glm_filter$ARO.Name)
# glm_phage <- list()
# for (i in 1:length(unique_args)) {
#   glm_phage[[i]] <- ggplot(arg_phage_glm_filter[arg_phage_glm_filter$ARO.Name == unique_args[i],], aes(rpkm, phage, colour = Location_sampletype)) +
#     geom_point(shape = 1) +
#     geom_smooth(method = "glm", se = FALSE, method.args = list(family = "binomial")) +
#     theme_bw() +
#     scale_colour_manual("GIT Site - Country", values = cohort_cols[names(cohort_cols) %in% arg_phage_glm_filter$Location_sampletype[arg_phage_glm_filter$ARO.Name == unique_args[i]]]) +
#     ggtitle(unique_args[i]) +
#     scale_y_continuous(breaks = c(0,1)) +
#     xlab("") +
#     ylab("") +
#     theme(plot.title = element_text(), legend.position = "none")
# }
# 
# tiff("figures/glm_phages.tiff", width = 1000, height = 500, res = 150)
# grid.arrange(grobs = glm_phage, layout_matrix = matrix(seq(1, length(glm_phage)), nrow = 1))
# dev.off()
# 
