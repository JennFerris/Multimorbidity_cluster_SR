# ****************************************************************************************************************************
#   Accompanying code for Ferris, et al. 2025 "A systematic review and meta-analysis of disease clustering in multimorbidity  
#   AUTHOR: Jennifer Ferris, Postdoctoral fellow (jennifer_ferris@sfu.ca)
#   LAST UPDATED: November 26th, 2025
# ****************************************************************************************************************************
#   Code to generate manuscript tables and figures 
#   Prerequisites : MM_SR_meta_analysis.R script run and results files saved in /data directory
# ****************************************************************************************************************************

#--- load packages and set directories -----
pkg <- c("tidyverse",
         "igraph",
         "ggraph",
         "linkcomm", 
         "munsell",
         "viridis",
         "arrow",
         "cowplot",
         "flextable",
         "scatterpie",
         "scales", 
         "clustAnalytics",
         "fpc")

sapply(pkg, require, character.only = T)


dir_data <- "../data/" #- directory where raw data files are stored
dir_output <- "../output/"  #- directory to output results


#- Note: Plot colour schemes from: https://personal.sron.nl/~pault/

#- load functions 
source("00_functions.R") 


#--- table and plot formatting -----

#- flextable settings 
set_flextable_defaults(digits = 3,
                       padding = 0.5,
                       font.family = 'Garamond')

border_small = officer::fp_border(color="grey", width=1)
border_big = officer::fp_border(color="black", width=1.3)
border_null = officer::fp_border(color="grey95", width=0)

format_flextable <- function(mydf, cat_label){
  mydf %>%
    as_grouped_data(groups = "study_characteristic") %>%
    flextable::as_flextable(hide_grouplabel = TRUE) %>%
    set_header_labels(category = cat_label, output = "Number of studies\n n (%)") %>%
    autofit() %>%
    hline(i=c(1), part="header", border=border_small) %>%
    hline_top(border=border_null, part="header") %>%
    hline_bottom(border=border_big, part="body") %>%
    hline_bottom(border=border_big, part="header") %>%
    padding(j = 2, padding.left = 15, part = "all") %>%
    align(i = ~is.na(study_characteristic), j = 1, align = 'right', part = "body") %>%
    bold(i = ~!is.na(study_characteristic),  bold = T, part = "body") %>%
    height(i = ~!is.na(study_characteristic),  height = 0.25, part = "body", unit = "in") %>%
    valign(i = ~!is.na(study_characteristic),  valign = "bottom", part = "body") %>%
    hrule(rule = "atleast", part = "body") 
}


#- set plot colour palettes 

#- scico palettes
cols_rob <- c("Low risk of bias" = "#034A85", 
              "Some concerns" = "#4A8FB2",
              "High risk of bias" = "#C98157", 
              "Very high risk of bias" = "#9D3709")

cols_sex = c("F" = "#52741B", "M" = "#C5AE32", "Both" = "#A1A1A1")

cols_method = c("Cluster-based" = "#C37469", "Factor-based" = "#5AA3DA")


cols_ages = c("All ages" = "#A1A1A1",
              "0-39" = "#D2A4C4",
              "40-59" = "#BA7AA9",
              "40+" = "#9B5A8C",
              "60+" = "#753D68")

#- from paul tol's colourblind friendly scheme: https://sronpersonalpages.nl/~pault/
cols_stable_clusters = c("1" = "#4477AA",
                         "2" =  "#EE6677",
                         "3" = "#228833",
                         "4" = "#CCBB44",
                         "5" = "#66CCEE",
                         "6" = "#AA3377"
)


#- disease label names
disease_labels <- data.frame(
  label = c("Asthma", "Cancer", "Chronic kidney disease", "COPD", "Dementia", 
            "Diabetes", "Epilepsy", "Gout", "Heart failure", "HIV/AIDS", "Hypertension",
            "Ischemic heart disease", "Mood & anxiety disorders", "Multiple sclerosis",
            "Osteoarthritis", "Osteoporosis", "Parkinson's disease", "Rheumatoid arthritis", 
            "Schizophrenia", "Stroke"),
  link_name = c("asthma", "cancer", "CKD", "COPD", "dementia", "DM", "epilepsy", "gout",
                "HF", "HIV", "HTN", "IHD", "MAD", "MS", "OA", "OP", "PD", "RA", "SCZ", "stroke")
)


#- manual network plot layout based on the 'mds' algorithm
layout <- read.csv(paste0(dir_data, "network_layout.csv")) %>%
  arrange(name, .locale = "en")

#--- read in data -----

#- read in study summary data
SR_summary <- read_csv(paste0(dir_data, "data_SR_summary.csv"))

SR_summary_clustermethods <- read_csv(paste0(dir_data, "data_SR_summary_byclustermethod.csv"))

#-- load data with age stratification of sample n's (fig 3)
SR_summary_agestrat <- read_csv(paste0(dir_data, "data_SR_summary_agestrat.csv"))

#- meta-analytic data 
SR_incl_dis <- read_csv(paste0(dir_data, "data_disease_harmonization.csv")) %>%
  mutate(across(everything(), ~ifelse(. == "n/a", NA, .)))

meta_dis <- read_csv(paste0(dir_data, "data_cluster_meta_analytics.csv"))  %>%
  mutate(across(everything(), ~ifelse(. == "n/a", NA, .))) %>%
  mutate(age_bin_2 = factor(age_bin_2, ordered = T, levels = c("0-39", "40-59", "60+", 
                                                               "All ages\n(0 - 100+)", 
                                                               "Adults\n(16+)", 
                                                               "Middle-aged\nand older\n(40+)")))
#- risk of bias data
data_ROB <- read.csv(paste0(dir_data, "data_RoB.csv")) %>%
  rename_with(str_to_lower) 

#- calculate disease co-occurrence across clusters 

diseases = colnames(meta_dis)[12:31]

combo_list_pairs <- 
  combn(diseases, 2) %>% 
  as.data.frame.table() %>%
  pivot_wider(names_from = Var1, values_from = Freq) %>%
  select(-Var2) %>%
  rename(Disease1 = A, Disease2 = B) %>%
  mutate(Disease1 = as.character(Disease1), Disease2 = as.character(Disease2))

#-- load network data
links <- read_csv(paste0(dir_data, "results_network_link_data.csv"))

nodes <- layout %>% 
  select(name, label) %>%
  rename(dis = name) %>%
  mutate(pct_dis = 100) 

net <- graph_from_data_frame(links, directed = F, vertices = nodes)

seq_link <- c(1,3,5,7,9) #- figure legend sequence


#--- overall study summaries -----
#----  Calculate summary data on study characteristics for tables 1 & 2

#- country/continent of publication
total_n = 82 #- for this table only
sum_pub_country <- 
  SR_summary %>%
  select(country_of_data) %>%
  separate_wider_delim(country_of_data, delim = ", ", names_sep = "_", too_few = "align_start") %>%
  pivot_longer(cols = everything(), values_to = "country_of_data") %>%
  mutate(country_of_data = ifelse(country_of_data %in% c("Finland", "France", "Germany", "Japan", "Norway", "Singapore", "Taiwan", "The netherlands", "the Netherlands"),
                                  "Other", country_of_data)) %>%
  select(-name) %>%
  filter(!is.na(country_of_data)) %>%
  count_vals(., country_of_data) %>% 
  arrange(-n, country_of_data) %>%
  mutate(continent_of_data = case_when(
    country_of_data %in% c("Spain", "Austria", "Denmark", "England", "Finland", "France", "Germany", 
                           "Italy", "Norway", "Sweden", "Switzerland",
                           "The netherlands") ~ "Europe",
    country_of_data %in% c("United states") ~ "North America",
    country_of_data %in% c("Australia") ~ "Oceania",
    country_of_data %in% c("China", "Japan", "Singapore", "South korea", "Taiwan") ~ "Asia",
    T ~ country_of_data
  ),
  country_of_data = str_to_title(country_of_data)) %>%
  relocate(continent_of_data, .before = country_of_data) 


sum_pub_continent <- 
  SR_summary %>%
  select(studyid:year, country_of_data) %>%
  separate_wider_delim(country_of_data, delim = ", ", names_sep = "_", too_few = "align_start") %>%
  pivot_longer(cols = starts_with("country"), values_to = "country_of_data") %>%
  select(-name) %>%
  filter(!is.na(country_of_data)) %>%
  mutate(continent_of_data = case_when(
    country_of_data %in% c("Spain", "Austria", "Denmark", "England", "Finland", "France", "Germany", 
                           "Italy", "Norway", "Sweden", "Switzerland", "United Kingdom",
                           "the Netherlands") ~ "Europe",
    country_of_data %in% c("United States") ~ "North America",
    country_of_data %in% c("Australia") ~ "Oceania",
    country_of_data %in% c("China", "Japan", "Singapore", "South Korea", "Taiwan") ~ "Asia",
    T ~ country_of_data
  )) %>%
  count_vals(continent_of_data, caps = F)  %>%
  rename(continent_output = output)


sum_pub_location <- sum_pub_country %>% 
  left_join(sum_pub_continent, by = "continent_of_data") %>%
  arrange(-n.y, -n.x) %>%
  select(-starts_with(c("n", "perc"))) %>%
  mutate(study_characteristic = "Country") %>%
  relocate(study_characteristic, .before = continent_of_data) %>%
  rename(category = country_of_data)

sum_pub_country <- sum_pub_country %>%
  select(-continent_of_data, -n, -perc_total) %>%
  mutate(study_characteristic = "Country") %>%
  relocate(study_characteristic, .before = country_of_data) %>%
  rename(category = country_of_data) %>%
  mutate(category = factor(category, ordered = T, 
                           levels = c("Spain",
                                      "United Kingdom",
                                      "United States",
                                      "Italy",
                                      "Denmark",
                                      "Sweden",
                                      "China",
                                      "Austria",     
                                      "Switzerland",
                                      "Australia",
                                      "South Korea",
                                      "The Netherlands",
                                      "Other"))) %>%
  arrange(category) %>%
  mutate(category = as.character(category))


#- type of study cohort 
total_n = 79
sum_study_design <- 
  SR_summary %>%
  count_vals(., study_design) %>%
  mutate(study_characteristic = "Study design") %>%
  mutate(study_design = case_when(study_design == "Cohort: retrospective" ~ "Retrospective cohort",
                                  study_design == "Cohort: prospective" ~ "Prospective cohort",
                                  T ~ study_design)) %>%
  relocate(study_characteristic, .before = study_design) %>%
  select(-n, -perc_total) %>%
  rename(category = study_design)


#- participant information 
sum_ages <- 
  SR_summary %>%
  mutate(age_bins = case_when(
    ages_included_low < 16 & (ages_included_high >= 65 | is.na(ages_included_high)) ~ "All ages",
    ages_included_low >= 16 & ages_included_low < 40 & (ages_included_high >= 100 | is.na(ages_included_high))  ~ "All adults (16+)",
    ages_included_low >= 60 ~ "Older adults (60+)",
    ages_included_low >= 40 & (ages_included_high >= 100 | is.na(ages_included_high))  ~ "Middle-aged & older adults (40+)",
    ages_included_low >= 16 & ages_included_low < 40 & ages_included_high <= 65 ~ "Young & middle-aged adults (16 - 65)",
    ages_included_low >= 40 & ages_included_high <= 65 ~ "Middle-aged adults (40 - 65)",
    T ~ "All ages"
  )) %>% 
  count_vals(age_bins) %>%
  mutate(age_bins = factor(age_bins, ordered = T, 
                           levels = c("All ages", "All adults (16+)", 
                                      "Young & middle-aged adults (16 - 65)",
                                      "Middle-aged adults (40 - 65)",
                                      "Middle-aged & older adults (40+)",
                                      "Older adults (60+)")
  )) %>%
  arrange(age_bins)  %>%
  mutate(study_characteristic = "Study population") %>%
  relocate(study_characteristic, .before = age_bins) %>%
  select(-n, -perc_total) %>%
  rename(category = age_bins)

sum_sample_size <- SR_summary %>%
  mutate(n_bin = case_when(
    sample_size <10000 ~ "Under 10,000",
    sample_size >= 10000 & sample_size < 500000 ~ "10,000 - 499,999",
    sample_size >= 500000 & sample_size < 1999999 ~ "500,000 - 1,999,999",
    sample_size >= 2000000 ~ "Over 2M",
    is.na(sample_size) ~ "Not stated"
  )) %>% 
  count_vals(n_bin) %>%
  mutate(n_bin = factor(n_bin, ordered = T, 
                        levels = c("Under 10,000",
                                   "10,000 - 499,999",
                                   "500,000 - 1,999,999",
                                   "Over 2m",
                                   "Not stated")
  )) %>%
  arrange(n_bin) %>%
  mutate(study_characteristic = "Sample size") %>%
  relocate(study_characteristic, .before = n_bin) %>%
  select(-n, -perc_total) %>%
  rename(category = n_bin)


#- MM definitions and disease information
sum_MM_def <- 
  SR_summary %>%
  mutate(MM_definition = replace_na(MM_definition, "not stated")) %>%
  count_vals(MM_definition) %>%
  mutate(MM_definition = str_to_sentence(MM_definition),
         MM_definition = factor(MM_definition, ordered = T,
                                levels = c("2+ chronic conditions",
                                           "2+ conditions", 
                                           "3+ chronic conditions",
                                           "Multiple chronic conditions",
                                           "Not stated"))) %>%
  arrange(MM_definition) %>%
  mutate(study_characteristic = "Multimorbidity definition") %>%
  relocate(study_characteristic, .before = MM_definition) %>%
  select(-n, -perc_total) %>%
  rename(category = MM_definition)


sum_HROs <- SR_summary %>%
  count_vals(included_any_HROs) %>%
  arrange(rev(included_any_HROs)) %>%
  mutate(study_characteristic = "Health-related outcomes included?") %>%
  relocate(study_characteristic, .before = included_any_HROs) %>%
  select(-n, -perc_total) %>%
  rename(category = included_any_HROs)


#- risk of bias --
sum_rob <- SR_summary %>%
  count_vals(risk_of_bias_overall) %>%
  mutate(risk_of_bias_overall = factor(risk_of_bias_overall, ordered = T,
                                       levels = c(
                                         "Low risk of bias",
                                         "Some concerns",
                                         "High risk of bias",
                                         "Very high risk of bias"
                                       ))) %>%
  arrange(risk_of_bias_overall) %>%
  mutate(study_characteristic = "Risk of bias") %>%
  relocate(study_characteristic, .before = risk_of_bias_overall) %>%
  select(-n, -perc_total) %>%
  rename(category = risk_of_bias_overall)

#- number of diseases
sum_n_diseases <- 
  SR_summary %>%
  mutate(n_diseases_bin = case_when(
    n_diseases_included < 50 ~ "Under 50",
    n_diseases_included >= 50 & n_diseases_included < 100 ~ "50 - 99",
    n_diseases_included >= 100 & n_diseases_included < 250 ~ "100 - 249",
    n_diseases_included >= 250 & n_diseases_included < 1000 ~ "250 - 1,000",
    n_diseases_included >= 1000 ~ "Over 1,000",
    is.na(n_diseases_included) ~ "Not stated"
  )) %>%
  count_vals(n_diseases_bin) %>%
  mutate(n_diseases_bin = factor(n_diseases_bin, ordered = T,
                                 levels = c("Under 50",
                                            "50 - 99",
                                            "100 - 249",
                                            "250 - 1,000",
                                            "Over 1,000",
                                            "Not stated"))) %>%
  arrange(n_diseases_bin) %>%
  mutate(study_characteristic = "Number of diseases included") %>%
  relocate(study_characteristic, .before = n_diseases_bin) %>%
  select(-n, -perc_total) %>%
  rename(category = n_diseases_bin)


total_n = 102
sum_disease_codes <- 
  SR_summary %>%
  select(disease_codes) %>%
  separate_wider_delim(disease_codes, delim = " and ", names_sep = "_", too_few = "align_start") %>%
  pivot_longer(cols = starts_with("disease_codes"), values_to = "disease_codes") %>%
  select(-name) %>%
  filter(!is.na(disease_codes)) %>%
  mutate(disease_codes = case_when(str_detect(disease_codes, "ICD") ~ "ICD-9 or -10",
                                   str_detect(disease_codes, "hospital") | str_detect(disease_codes, "FAERS") | 
                                     str_detect(disease_codes, "OPSC") | str_detect(disease_codes, "gemsscript") |
                                     str_detect(disease_codes, "Medcodes") ~ "Other",   
                                   disease_codes == "not stated" ~ "Not stated",
                                   T ~ disease_codes)) %>%
  count_vals(disease_codes, caps = F) %>%
  mutate(disease_codes = factor(disease_codes, ordered = T,
                                levels = c("ICD-9 or -10",
                                           "ATC", 
                                           "ICPC",
                                           "READ",
                                           "SNODMED",
                                           "Other",
                                           "Not stated"))) %>%
  arrange(disease_codes) %>%
  mutate(study_characteristic = "Coding system") %>%
  relocate(study_characteristic, .before = disease_codes) %>%
  select(-n, -perc_total) %>%
  rename(category = disease_codes) %>%
  mutate(category = as.character(category))


total_n = 80
sum_disease_groups <- 
  SR_summary %>%
  select(studyid:year, disease_groupings) %>%
  mutate(disease_groupings = case_when(str_detect(disease_groupings, "[345]") ~ "3- to 5-digit ICD codes", 
                                       str_detect(disease_groupings, "Charlson|Elixhauser") ~ "Charlson or Elixhauser index algorithms",
                                       T ~ disease_groupings)) %>%
  separate_wider_delim(disease_groupings, delim = " + ", names_sep = "_", too_few = "align_start") %>%
  pivot_longer(cols = starts_with("disease_groupings"), values_to = "disease_groupings") %>%
  select(-name) %>%
  filter(!is.na(disease_groupings)) %>% 
  mutate(disease_groupings = case_when(
    !str_detect(disease_groupings, "Quality and Outcomes") &
      !str_detect(disease_groupings, "Custom groupings") &
      !str_detect(disease_groupings, "2017") &
      !str_detect(disease_groupings, "ICD") &
      !str_detect(disease_groupings, "Clinical Classification Software") &
      !str_detect(disease_groupings, "Charlson") &
      !str_detect(disease_groupings, "John's Hopkins") ~ "Other",
    T ~ disease_groupings)) %>%
  count_vals(disease_groupings, caps = F) %>%
  mutate(disease_groupings = factor(disease_groupings, ordered = T,
                                    levels = c("Custom groupings",
                                               "3- to 5-digit ICD codes" ,
                                               "Disease groupings defined by Calderón-Larrañaga, 2017",
                                               "Block-level ICD codes" ,
                                               "Agency for Healthcare Research and Quality Clinical Classification Software",
                                               "Charlson or Elixhauser index algorithms",
                                               "John's Hopkins Adjusted Clinical Group Expanded Diagnostic Clusters",
                                               "National Institute for Health and Care Excellence Quality and Outcomes Framework",
                                               "Other"))) %>%
  arrange(disease_groupings) %>%
  mutate(study_characteristic = "Method to group codes into disease entities") %>%
  relocate(study_characteristic, .before = disease_groupings) %>%
  select(-n, -perc_total) %>%
  rename(category = disease_groupings) %>%
  mutate(category = as.character(category))



#- cluster methods
total_n = 84

sum_cluster_methods_count <- SR_summary_clustermethods %>%
  count_vals(., cluster_category, cluster_method, caps = F) %>%
  arrange(-n, cluster_method) %>%
  ungroup() %>%
  mutate(study_characteristic = "Clustering method") %>%
  relocate(study_characteristic, .before = cluster_method) %>%
  select(-cluster_category, -n, -perc_total) %>%
  rename(category = cluster_method)

sum_cluster_temporality <- SR_summary_clustermethods %>%
  count_vals(temporal_analysis) %>%
  mutate(study_characteristic = "Clustering method") %>%
  relocate(study_characteristic, .before = temporal_analysis) %>%
  select(-n, -perc_total) %>%
  #arrange(rev(temporal_analysis)) %>%
  mutate(study_characteristic = "Disease temporality included in clustering?") %>%
  relocate(study_characteristic, .before = temporal_analysis) %>%
  rename(category = temporal_analysis)

sum_cluster_n <-  SR_summary_clustermethods %>%
  mutate(n_cluster_cat = case_when(
    n_clusters < 5 ~ "under 5",
    n_clusters >= 5 & n_clusters < 25 ~ "5 - 25",
    n_clusters >= 25 & n_clusters < 100 ~ "25 - 100",
    n_clusters >= 100 ~ "over 100",
    T ~ "not reported"
  )) %>%
  count_vals(n_cluster_cat) %>%
  mutate(n_cluster_cat = factor(n_cluster_cat, ordered = T,
                                levels = c(
                                  "Under 5",
                                  "5 - 25",
                                  "25 - 100",
                                  "Over 100",
                                  "Not reported"
                                ))) %>%
  arrange(n_cluster_cat) %>%
  mutate(study_characteristic = "Number of disease clusters obtained") %>%
  relocate(study_characteristic, .before = n_cluster_cat) %>%
  select(-n, -perc_total) %>%
  rename(category = n_cluster_cat)


#----- Table 1: summarizing included studies  -----
table_data_sum <- 
  bind_rows(sum_pub_country, 
            sum_study_design, sum_sample_size, 
            sum_ages, sum_MM_def, 
            sum_HROs, sum_rob)

title_table2 <- "Table 1: Summary characteristics of included studies."
footer_table2 <- "Numbers in table may exceed the total number of included studies (n = 79) because some studies used multiple data sources. Percentages may not add up to 100% because of rounding."

table_sum <- format_flextable(table_data_sum, "Study Characteristics") %>%
  add_footer_lines(values = footer_table2) %>%
  add_header_lines(values = title_table2) %>%
  align(align = "left", part = "footer") %>%  
  hline(i = 1, border=border_big, part="footer") %>%
  fontsize(size = 11, part = "footer")  %>%
  fontsize(size = 12, part = "header") %>%
  autofit() %>%
  width(j = 1, width = 2.25)

save_as_docx(table_sum, path = paste0(dir_output, "table1_study_summaries.docx"))


#----- Table 2: clustering information -----

table_data_sum_cluster <-
  bind_rows(sum_disease_codes, sum_disease_groups, 
            sum_n_diseases, sum_cluster_methods_count, 
            sum_cluster_temporality, sum_cluster_n)


title_table3 <- "Table 2: Summary characteristics disease definitions and clustering methodologies in included studies."
footer_table3 <- "Numbers in table may exceed the total number of included studies (n = 79) because some studies used multiple data sources and/or clustering methods. Percentages may not add to 100% because of rounding."

table_sum_clusters <- format_flextable(table_data_sum_cluster, "Disease clustering characteristics") %>%
  add_footer_lines(values = footer_table3) %>%
  add_header_lines(values = title_table3) %>%
  align(align = "left", part = "footer") %>%  
  hline(i = 1, border=border_big, part="footer") %>%
  fontsize(size = 11, part = "footer")  %>%
  fontsize(size = 12, part = "header") %>%
  autofit() %>%
  width(j = 1, width = 2.5) 

save_as_docx(table_sum_clusters, path = paste0(dir_output, "table2_study_summaries_clusters.docx"))


#----- Figure 2: plot study characteristics by method ---------
total_n = 84

plot_data_cluster_methods <- SR_summary_clustermethods %>%
  #-  collapse network studies into a single category to be consistent with results section
  mutate(cluster_method = case_when(str_detect(cluster_method, "Network") ~ "Network studies",
                                    T ~ cluster_method))

plot_data_cluster_methods_sum <- plot_data_cluster_methods %>%
  count_vals(., cluster_method, caps = F) %>%
  arrange(-n, cluster_method) %>%
  select(-output) %>%
  rename(n_method = n) 

plot_data_cluster_methods_names <- plot_data_cluster_methods_sum %>%
  filter(n_method > 5) %>%
  pull(cluster_method) #-- include only the top 5 most common methods

#- methods across time
plot_data_cluster_methods_time <- plot_data_cluster_methods %>%
  filter(cluster_method %in% plot_data_cluster_methods_names) %>%
  mutate(cluster_method = str_replace(cluster_method, " analysis", "\nanalysis"),
         cluster_method = str_replace(cluster_method, " cluster", "\ncluster"))

plot_data_cluster_methods_names_2 <- plot_data_cluster_methods_names %>%
  str_replace(., " analysis", "\nanalysis") %>%
  str_replace(., " cluster", "\ncluster")


plot_data_cluster_methods_time$cluster_method <- ordered(plot_data_cluster_methods_time$cluster_method, 
                                                         levels = plot_data_cluster_methods_names_2)

plot_cluster_methods_time <- 
  ggplot(plot_data_cluster_methods_time, aes(x = year)) +
  geom_histogram(binwidth = 1, linewidth = 0.5) +
  facet_wrap(~cluster_method, nrow = 2) +
  labs(y = "Number of studies", x = "Publication year") +
  theme_minimal() +
  scale_x_continuous(limits = c(2008, 2025), breaks = c(2009, 2014, 2019, 2024)) +
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(linetype = 'dotted', colour = "dark grey", linewidth = 0.15),
        axis.ticks.x = element_blank(), #panel.spacing = unit(2, "mm"),
        plot.margin = margin(t = 10, b = 2, l = 2, r = 5.5, unit = "pt"),
        text = element_text(size = 4.5), 
        axis.title = element_text(size = 4),
        axis.text.x = element_text(size = 3)) 

plot_cluster_methods_time_all <-
  ggplot(SR_summary_clustermethods, aes(x = year)) +
  geom_histogram(binwidth = 1, linewidth = 0.5) +
  labs(y = "Number of studies", x = "Publication year") +
  theme_minimal() +
  scale_x_continuous(limits = c(2008, 2025), breaks = c(2009, 2014, 2019, 2024)) +
  theme(plot.title.position = "plot",
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(linetype = 'dotted', colour = "dark grey", linewidth = 0.15),
        axis.ticks.x = element_blank(), panel.spacing = unit(3, "mm"),
        plot.margin = margin(t = 10, b = 2, l = 5.5, r = 5.5, unit = "pt"),
        text = element_text(size = 5))

plot_cluster_methods_time_2 <- 
  plot_grid(plot_cluster_methods_time_all, plot_cluster_methods_time, nrow = 1, rel_widths = c(1, 1))



#- most common methods by ROB
plot_data_cluster_methods_rob <- plot_data_cluster_methods %>%
  count_vals(., cluster_method, risk_of_bias_overall, caps = F) %>%
  arrange(-n, cluster_method) %>%
  select(-perc_total, -output)


plot_data_cluster_methods_rob <- left_join(plot_data_cluster_methods_rob, plot_data_cluster_methods_sum) %>%
  mutate(perc_method = n/n_method*100) %>%
  filter(n_method > 5) %>%
  mutate(risk_of_bias_overall = 
           factor(str_to_sentence(risk_of_bias_overall), ordered = T,
                  levels = c(
                    "Very high risk of bias",
                    "High risk of bias",
                    "Some concerns",
                    "Low risk of bias"
                  ))) 


plot_cluster_ROB <- 
  ggplot(plot_data_cluster_methods_rob, aes(x = n, y = reorder(cluster_method, perc_total), fill = risk_of_bias_overall)) + 
  geom_col(width = 0.75, alpha = 0.9) +
  labs(y = "", x = "Number of studies", fill = "") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.major.x = element_line(linetype = 'dotted', colour = "dark grey", linewidth = 0.15),
        text = element_text(size = 4.5),
        plot.margin = margin(t = 20,l = 15,r = 15, b = 5.5),
        legend.key.size = unit(3, "mm"),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = cols_rob) +
  guides(fill = guide_legend(reverse = T))



plot_cluster_methods <- plot_grid(plot_cluster_methods_time_2, plot_cluster_ROB, nrow = 2, rel_heights = c(1, 0.8),
                                  labels = c("A: Included publications by year and clustering method",
                                             "B: Risk of bias by clustering method"),
                                  hjust = 0, label_fontface = "plain", label_size = 5, vjust = c(2, 4), label_x = 0.025)


ggsave(paste0(dir_output, "Fig2_clustering_methods.pdf"), 
       plot_cluster_methods, bg = "white", width = 88, height = 70, units = "mm")


#----- Supp Fig 1: Risk of bias data -------

counts_ROB <- count_vals_ROB(overall.risk.of.bias) %>% 
  rename(risk.of.bias = overall.risk.of.bias, overall.risk.of.bias = n) %>%
  left_join(count_vals_ROB(domain.1.rob) %>% rename(risk.of.bias = domain.1.rob, domain.1.rob = n)) %>%
  left_join(count_vals_ROB(domain.2.rob) %>% rename(risk.of.bias = domain.2.rob, domain.2.rob = n)) %>%
  left_join(count_vals_ROB(domain.5.rob) %>% rename(risk.of.bias = domain.5.rob, domain.5.rob = n)) %>%
  left_join(count_vals_ROB(domain.6.rob) %>% rename(risk.of.bias = domain.6.rob, domain.6.rob = n)) %>%
  left_join(count_vals_ROB(domain.7.rob) %>% rename(risk.of.bias = domain.7.rob, domain.7.rob = n)) %>%
  arrange(risk.of.bias) %>%
  mutate(risk.of.bias = str_to_sentence(risk.of.bias)) %>%
  mutate(risk.of.bias = factor(risk.of.bias, ordered = T, 
                               levels = c("Low risk of bias", "Some concerns",
                                          "High risk of bias", "Very high risk of bias"))) 

counts_ROB_pcnt <- counts_ROB %>%
  mutate(across(overall.risk.of.bias:domain.7.rob, ~(./68*100))) 


plot_ROB <- counts_ROB %>%
  pivot_longer(cols = overall.risk.of.bias:domain.7.rob, names_to = "measure", values_to = "n") %>%
  mutate(measure = case_match(measure, 
                              "overall.risk.of.bias" ~ "Overall\nrisk of bias",
                              "domain.1.rob" ~ "Domain 1-\n Confounding",
                              "domain.2.rob" ~ "Domain 2-\n Exposure",
                              "domain.5.rob" ~ "Domain 5-\n Missing data",
                              "domain.6.rob" ~ "Domain 6-\n Clustering",
                              "domain.7.rob" ~ "Domain 7-\n Reporting results"),
         measure = factor(measure, ordered = T, 
                          levels = c("Overall\nrisk of bias",
                                     "Domain 1-\n Confounding",
                                     "Domain 2-\n Exposure",
                                     "Domain 5-\n Missing data",
                                     "Domain 6-\n Clustering",
                                     "Domain 7-\n Reporting results")),) %>%
  ggplot(aes(y = risk.of.bias, x = measure, fill = risk.of.bias, alpha = n, label = n)) +
  geom_tile() +
  geom_text(alpha = 1, 
            aes(fontface = ifelse(measure == "Overall\nRisk of Bias", "bold", "plain"))) +
  scale_alpha(range = c(0.25,1), na.value = 0) +
  scale_fill_manual(values = cols_rob) +
  theme_minimal() +
  theme(panel.grid.major = element_blank()) +
  labs(y = "", x = "", alpha = "number of studies") +
  guides(fill = "none", alpha = "none")


caption_sfig1 <- paste0(
  "Results of Risk of Bias (ROB) assessment of multimorbidity disease clustering studies, adapted from the
ROBINS-E tool.The table cells represent the number of studies in each domain. Columns present the overall
ROB assessment, and each ROB assessment domain. Domain 1- Confounding: assesses confounding that
was not controlled by the design or analysis of the study – and is therefore likely to bias the estimated effect
of exposure on outcome. Domain 2- Exposure: addresses measurement error or misclassification in exposure.
Domain 5- Missing data: whether data were missing and, if so, the potential that missing data led to bias
in analysis. Domain 6- Clustering: custom domain, evaluating the risk of bias in disease clustering results
based on reported methodology and cluster validation.  Domain 7- Reporting results: assessing bias from 
selective reporting of results, which may arise from a desire for findings to be noteworthy.
"
)

caption_sfig1 <- ggdraw() +
  draw_label(caption_sfig1, x = 0.025, size = 10,
             hjust = -0.005, vjust = 0.5)

title_sfig1 <- ggdraw() +
  draw_label("Supplementary Figure 1: Risk of Bias assessment", 
             size = 11,
             hjust = 1, vjust = 1)

plot_ROB_2 <- plot_grid(title_sfig1, plot_ROB, caption_sfig1, ncol = 1, rel_heights = c(0.05, 1, 0.31))

ggsave(paste0(dir_output, "Supp_Fig1_ROB.png"), plot_ROB_2, bg = "white", dpi = 400,
       height = 5.5, width = 7)


#----- Figure 3: Meta-analysis study characteristics -----

#- format plot data -
plot_data_SR_summary_agestrat_sex <- SR_summary_agestrat %>%
  select(studyid, subgoup_name, age_bin, age_bin_2, subgroup_sex, subgroup_n, subgroup_n_F, subgroup_n_M) %>%
  mutate(subgroup_n_F = case_when(is.na(subgroup_n_F) & subgroup_sex == "F" ~ subgroup_n,
                                  T ~ subgroup_n_F),
         subgroup_n_M = case_when(is.na(subgroup_n_M) & subgroup_sex == "M" ~ subgroup_n,
                                  T ~ subgroup_n_M),
         subgroup_n_Both = case_when(is.na(subgroup_n_F) & is.na(subgroup_n_M) ~ subgroup_n,
                                     T ~ NA)) %>%
  select(-subgroup_n, -subgroup_sex) %>%
  pivot_longer(cols = subgroup_n_F:subgroup_n_Both, names_to = "sex", values_to = "sample_size") %>%
  mutate(sex = str_remove(sex, "subgroup_n_"),
         sex = factor(sex, ordered = T, levels = c("F", "M", "Both")),
         sample_size = na_if(sample_size, 0)) %>%
  filter(!is.na(sample_size))


total_n = 73

#- rename some diseases for nicer plotting
plot_data_harmonized_dis <- SR_incl_dis %>% 
  filter(!is.na(harmonized_disease_name), !is.na(study_disease_name)) %>%
  count_vals(., harmonized_disease_name) %>%
  mutate(harmonized_disease_name = case_when(
    str_detect(harmonized_disease_name, "Dementia") ~ "Dementia",
    harmonized_disease_name == "Copd" ~ "COPD",
    str_detect(harmonized_disease_name, "Hiv") ~ "HIV/AIDS",
    str_detect(harmonized_disease_name, "Parkin") ~ "Parkinson's disease",
    harmonized_disease_name == "Diabetes mellitus" ~ "Diabetes",
    harmonized_disease_name == "Mood and anxiety disorders" ~ "Mood & anxiety disorders",
    T ~ harmonized_disease_name
  ))

disease_order <- plot_data_harmonized_dis %>% pull(harmonized_disease_name) %>% rev()

plot_data_harmonized_dis$harmonized_disease_name <- factor(plot_data_harmonized_dis$harmonized_disease_name, levels = disease_order)


#- create plots
plot_study_ages <- 
  ggplot(plot_data_SR_summary_agestrat_sex, aes(y = sample_size, x = age_bin_2, color = sex, fill = sex)) + 
  geom_col() + 
  labs(x = "Age groups", y = "Total study participants", color = "", fill = "")  +
  scale_y_continuous(labels = label_comma()) +
  scale_fill_manual(values = cols_sex, labels = c("Female", "Male", "Sex distribution\nnot reported")) +
  scale_color_manual(values = cols_sex, labels = c("Female", "Male", "Sex distribution\nnot reported")) +
  theme_minimal()  +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.major.y = element_line(linetype = 'dotted', colour = "dark grey", linewidth = 0.15),
        text = element_text(size = 6.5),
        #plot.margin = margin(t = 20,15,15,15),
        legend.key.size = unit(3, "mm"),
        axis.ticks.x = element_blank())


plot_dis_pct <- 
  ggplot(plot_data_harmonized_dis, aes(x = perc_total, y = harmonized_disease_name)) +
  geom_point(size = 0.75) +
  geom_segment(aes(x = 0, xend = perc_total, y = harmonized_disease_name, yend = harmonized_disease_name),
               linewidth = 0.2) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dotted', colour = "dark grey", linewidth = 0.15),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 6.5),
        plot.margin = margin(t = 20, b= 20, l = 0, r = 20, unit = "pt"),
        axis.ticks.x = element_blank()) +
  labs(y = "", x = "Percent of studies (n = 73)" )+
  #title = "Frequency of disease inclusion across studies") +
  scale_x_continuous(limits = c(0, 100), labels = scales::label_percent(scale = 1)) 



plot_meta_sample_char <- 
  plot_grid(plot_study_ages, plot_dis_pct, align = "h", rel_widths = c(1, 0.7),
            labels = c("A: Meta-analysis sample age and sex distribution",
                       "B: Harmonized disease inclusion across studies"),
            hjust = 0, label_fontface = "plain", label_size = 7, vjust = 1.75, label_x = 0.01)


ggsave(paste0(dir_output, "Fig3_meta_analysis_sample.pdf"), 
       plot_meta_sample_char, bg = "white", width = 180, height = 75, units = "mm")


#----- Figure 4: Meta-analytic network and sample characteristics ------

#--- prepare data for node characteristics subplots
#- node degree
net_degree <- as.data.frame(degree(net))

plot_data_degree <- data.frame(row.names(net_degree)) %>% 
  rename(`dis` = `row.names.net_degree.`) %>%
  left_join(disease_labels, by = join_by("dis" == "link_name")) 

plot_data_degree$degree <- net_degree$`degree(net)`
plot_data_degree_order <- plot_data_degree %>% arrange(degree) %>% pull(label)
plot_data_degree$label <- factor(plot_data_degree$label, levels = plot_data_degree_order)

#- link lift
plot_data_links <- links %>%
  pivot_longer(cols = c(from, to), values_to = "dis") %>%
  select(-name) %>%
  group_by(dis) %>%
  summarize(n = n(), avg_weight = mean(weight, na.rm=T)) %>%
  left_join(disease_labels, by = join_by("dis" == "link_name")) %>%
  filter(!is.na(label))

#- order by node degree
plot_data_links$label <- factor(plot_data_links$label, levels = plot_data_degree_order)


#- create disease network plot
plot_disease_network <- plot_network(net) +
  labs(edge_color = "Lift",
       edge_alpha = "Lift", edge_width = "Lift") + 
  theme(plot.title = element_text(size = 12, hjust = 0, face = "plain"),
        plot.caption = element_text(size = 9.5, hjust = 0),
        text = element_text(size  = 9),
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(size = 8)) 


#---- create node subplots
plot_B <- lolipop_plot(plot_data_degree, degree, label) +
  labs(y = "", x = "Node degree") +
  scale_x_continuous(limits = c(0, 20)) + 
  theme(plot.margin = margin(t = 5.5, b= 5.5, l = 20, r = 10, unit = "pt"))

plot_C <- lolipop_plot(plot_data_links, avg_weight, label) +
  labs(y = "", x = "Average link lift") +
  scale_x_continuous(limits = c(0, 4))  + 
  theme(plot.margin = margin(t = 5.5, b= 5.5, l = 10, r = 30, unit = "pt"))


plot_row_top <- plot_grid(plot_disease_network, labels = "A: Meta-analytic disease network", 
                          hjust = 0, label_fontface = "plain", label_size = 10,
                          label_x = 0.025) 

plot_row_bottom <- plot_grid(plot_B, plot_C, nrow = 1, align = "hv", 
                             labels = c("B: Number of connections with other diseases in network", 
                                        "C: Strength of association with other diseases in network"), 
                             hjust = 0, label_fontface = "plain", label_size = 9, vjust = -1,
                             label_x = 0.05) 


plot_fig2 <- plot_grid(plot_row_top, plot_row_bottom, nrow = 2, rel_heights = c(1, 0.5)) 

ggsave(paste0(dir_output, "Fig4_disease_network.pdf"), plot_fig2, bg = "white", dpi = 400,
       height = 180, width = 180, units = "mm")






#----- Figure 5: network with clusters ------

#- Results:  stability estimates -- 
cluster_stability_results <- read.csv(paste0(dir_data, "results_cluster_stability.csv")) %>%
  arrange(-jaccard_stability) %>%
  mutate(stable_cluster = 1:43) %>%
  select(cluster, stable_cluster, jaccard_stability)

clusters_stable <- read.csv(paste0(dir_data, "results_clusters.csv")) %>%
  left_join(cluster_stability_results) %>%
  filter(jaccard_stability > 0.5) %>%
  mutate(cluster_id = as.ordered(cluster),
         stable_cluster_id = as.ordered(stable_cluster))

layout_2 <- read.csv(paste0(dir_data, "network_layout_2.csv")) %>%
  mutate(x = ifelse(name == "stroke", x+0.04, x)) %>%
  mutate(x = ifelse(name == "HF", x-0.01, x)) %>%
  arrange(name, .locale = "en")

links_clust <- links %>%
  left_join(clusters_stable) %>%
  arrange(stable_cluster_id) %>%
  mutate(stable_cluster_id = factor(as.character(stable_cluster_id), ordered = T,
                                    levels = c("1", "2", "3", "4", "5", "6", NA))) %>%
  rename(jaccard = jaccard_stability) %>%
  select(from, to, weight, stable_cluster_id, jaccard)


nodes_clust <- links_clust %>%
  pivot_longer(cols = c(from, to), values_to = "node") %>%
  select(-weight, -name) %>%
  filter(!is.na(stable_cluster_id)) %>%
  group_by(node) %>%
  mutate(total = n()) %>%
  group_by(stable_cluster_id, node, total) %>%
  tally %>%
  mutate(perc_cluster = n/total*100) %>%
  left_join(select(layout_2, name, x, y), by = join_by("node" == "name")) %>%
  rename(x1 = x, y1 = y) %>%
  ungroup() %>%
  arrange(stable_cluster_id) %>%
  select(-n, -total) %>%
  pivot_wider(names_from = stable_cluster_id, values_from = perc_cluster) %>%
  mutate(across(everything(), ~replace_na(., 0))) %>%
  arrange(node, .locale = "en")


net_clust <- graph_from_data_frame(links_clust, directed = F, vertices = nodes)

#- create plots 

legend_data <- data.frame(
  stable_cluster_id = factor(c("1", "2", "3", "4", "5", "6", NA),
                             levels = c("1", "2", "3", "4", "5", "6", NA)),
  x = 0,
  y = 0
)



plot_disease_network_clust <- 
  plot_network_with_clusters(net_clust) + 
  labs(edge_color = "Lift",
       edge_alpha = "Lift", edge_width = "Lift") +
  theme(text = element_text(size = 9),
        plot.margin = margin(b = 0, t = 5.5, l = 5.5, r = 5.5, unit = "mm"),
        legend.key.height = unit(0.75, "lines")) 


nodes_clust_plot <- nodes %>%
  #- change the nudge variable so all disease labels are consistently plotted for cluster output
  mutate(neg_nudge = F) %>%
  left_join(select(layout_2, name, x, y), by = join_by("dis" == "name"))


#- create cluster plots
plot_network_cluster(links_clust, nodes_clust_plot, 1, neg_y_nudge_nodes = c("Diabetes"), x_nudge_nodes = c("Chronic kidney disease"), x_nudge = -0.05)
plot_network_cluster(links_clust, nodes_clust_plot, 2, neg_y_nudge_nodes = c("Diabetes", "Hypertension"), x_nudge_nodes = c("Hypertension"), x_nudge = -0.1, link_strength = 0.05)
plot_network_cluster(links_clust, nodes_clust_plot, 5, neg_y_nudge_nodes = c("Diabetes", "Chronic kidney disease", "COPD", "Osteoarthritis"), y_nudge = 0.15,
                     x_nudge_nodes = c("Chronic kidney disease"), x_nudge = 0.05)

#- custom functions for some plot tweaks
plot_network_cluster_3 <- function(link_df, node_df, clust_num){
  link_strength = 0.1
  seq_links = seq_link 
  round_seq_by = 1
  label_size = 2.25
  col_clusters = T 
  node_size = 2.8
  x_nudge = 0 
  y_nudge = 0.09
  neg_y_nudge_nodes = c("Chronic kidney disease")
  jacc_stability = T
  
  clust_links <- link_df %>%
    filter(stable_cluster_id == clust_num)
  
  clust <- link_df %>%
    pivot_longer(cols = c(from, to), names_to = "node", values_to = "disease") %>% 
    dplyr::select(-node, -weight) 
  
  if(jacc_stability == T){
    clust <- clust %>%
      dplyr::select(-jaccard)
  }
  
  clust <- clust %>%
    distinct() %>%
    filter(stable_cluster_id == clust_num) %>% pull(disease)
  
  clust_nodes <- node_df %>%
    filter(dis %in% clust)
  
  clust_net <- graph_from_data_frame(clust_links, directed = T, vertices = clust_nodes)
  
  clust_col = "#228833"
  
  
  clust_plot <- 
    ggraph(clust_net, layout = "manual", x = clust_nodes$x, y = clust_nodes$y) +
    geom_edge_arc(edge_width = 0.5,
                  strength = link_strength, alpha = 0.75, color = clust_col) + 
    geom_node_point(size = {{node_size}},
                    shape = 21, fill = clust_col, color = "00") +
    geom_node_text(aes(label = str_wrap(label, 20, whitespace_only = F),
                       x = x + ifelse(label == "Chronic kidney disease", -0.45, -0.05),
                       y = y + case_when(label %in% neg_y_nudge_nodes ~ -y_nudge, 
                                         label == "Schizophrenia" ~ 0.12,
                                         T ~ y_nudge), 
                       vjust = ifelse(label %in% neg_y_nudge_nodes, 0.25, 0)),
                   size = label_size,
                   show.legend = F) +
    theme_void() 
  
  
  clust_plot <- clust_plot +
    labs(title = paste0("Cluster ", clust_num)) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, vjust = 0, size = 8, face = "bold"), 
          plot.subtitle = element_text(hjust = 0.5, vjust = 1, size = 7, face = "bold")) +
    scale_x_continuous(expand = ggplot2::expansion(mult = 0.25)) +
    scale_y_continuous(expand = ggplot2::expansion(mult = 0.35))
  
  if(jacc_stability == T){
    jacc = clust_links %>% pull(jaccard) %>% unique() %>% round(., 2)
    
    clust_plot <- clust_plot + 
      labs(subtitle = paste0("Jaccard : ", format(jacc, nsmall = 2)))
  }
  
  assign(paste0("clust_plot_", clust_num), clust_plot, envir = .GlobalEnv)
  
}
plot_network_cluster_3(links_clust, nodes_clust_plot, 3)

plot_network_cluster_4 <- function(link_df, node_df, clust_num){
  link_strength = 0.1
  seq_links = seq_link 
  round_seq_by = 1
  label_size = 2.25
  col_clusters = T 
  node_size = 2.8
  x_nudge = 0 
  y_nudge = 0.09
  neg_y_nudge_nodes = c("Mood & anxiety disorders")
  jacc_stability = T
  
  clust_links <- link_df %>%
    filter(stable_cluster_id == clust_num)
  
  clust <- link_df %>%
    pivot_longer(cols = c(from, to), names_to = "node", values_to = "disease") %>% 
    dplyr::select(-node, -weight) 
  
  if(jacc_stability == T){
    clust <- clust %>%
      dplyr::select(-jaccard)
  }
  
  clust <- clust %>%
    distinct() %>%
    filter(stable_cluster_id == clust_num) %>% pull(disease)
  
  clust_nodes <- node_df %>%
    filter(dis %in% clust)
  
  clust_net <- graph_from_data_frame(clust_links, directed = T, vertices = clust_nodes)
  
  clust_col = "#CCBB44"
  
  
  clust_plot <- 
    ggraph(clust_net, layout = "manual", x = clust_nodes$x, y = clust_nodes$y) +
    geom_edge_arc(edge_width = 0.5,
                  strength = link_strength, alpha = 0.75, color = clust_col) + 
    geom_node_point(size = {{node_size}},
                    shape = 21, fill = clust_col, color = "00") +
    geom_node_text(aes(label = str_wrap(label, 20, whitespace_only = F),
                       x = x + case_when(!label %in% neg_y_nudge_nodes ~ -0.1, 
                                         #label == "COPD" ~ -2,
                                         T ~ 0.075),
                       y = y + ifelse(label %in% neg_y_nudge_nodes, -y_nudge - 0.2, y_nudge), 
                       vjust = case_when(label %in% neg_y_nudge_nodes ~ 0.25, 
                                         label == "COPD" ~ 1.5,
                                         T ~ 0),
                       hjust = ifelse(label == "COPD", 1, 0.5)),
                   size = label_size,
                   show.legend = F) +
    theme_void() 
  
  
  clust_plot <- clust_plot +
    labs(title = paste0("Cluster ", clust_num)) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, vjust = 0, size = 8, face = "bold"), 
          plot.subtitle = element_text(hjust = 0.5, vjust = 1, size = 7, face = "bold")) +
    scale_x_continuous(expand = ggplot2::expansion(mult = 0.25)) +
    scale_y_continuous(expand = ggplot2::expansion(mult = 0.35))
  
  if(jacc_stability == T){
    jacc = clust_links %>% pull(jaccard) %>% unique() %>% round(., 2)
    
    clust_plot <- clust_plot + 
      labs(subtitle = paste0("Jaccard : ", format(jacc, nsmall = 2)))
  }
  
  assign(paste0("clust_plot_", clust_num), clust_plot, envir = .GlobalEnv)
  
}
plot_network_cluster_4(links_clust, nodes_clust_plot, 4)

plot_network_cluster_6 <- function(link_df, node_df, clust_num, x_nudge_nodes = c("Ischemic heart disease", "Diabetes"), x_nudge = -0.1){
  link_strength = 0.1
  auto_nudge = F
  seq_links = seq_link
  round_seq_by = 1
  label_size = 2.25
  col_clusters = T
  node_size = 2.8
  y_nudge = 0.1
  neg_y_nudge_nodes = c("Dementia")
  jacc_stability = T
  
  clust_links <- link_df %>%
    filter(stable_cluster_id == clust_num)
  
  clust <- link_df %>%
    pivot_longer(cols = c(from, to), names_to = "node", values_to = "disease") %>% 
    dplyr::select(-node, -weight) 
  
  if(jacc_stability == T){
    clust <- clust %>%
      dplyr::select(-jaccard)
  }
  
  clust <- clust %>%
    distinct() %>%
    filter(stable_cluster_id == clust_num) %>% pull(disease)
  
  clust_nodes <- node_df %>%
    filter(dis %in% clust)
  
  
  if(auto_nudge == T){
    if(length(unique(clust_nodes$dis)) == 3) {
      y_nudge = 0.015
      #link_strength = 0.05
    }
  }
  
  clust_net <- graph_from_data_frame(clust_links, directed = T, vertices = clust_nodes)
  
  clust_col = "#698896"
  
  if(col_clusters == T){
    clust_col <- if (clust_num == 1){
      "#4477AA"
    }  else if (clust_num == 2){
      "#EE6677"
    }   else if (clust_num == 3){
      "#228833"
    } else if (clust_num == 4){
      "#CCBB44"
    } else if (clust_num == 5){
      "#66CCEE"
    }  else if (clust_num == 6){
      "#AA3377" 
    } 
    
  } 
  
  
  clust_plot <- 
    ggraph(clust_net, layout = "manual", x = clust_nodes$x, y = clust_nodes$y) +
    geom_edge_arc(edge_width = 0.5,
                  strength = link_strength, alpha = 0.75, color = clust_col) + 
    geom_node_point(size = {{node_size}},
                    shape = 21, fill = clust_col, color = "00") +
    geom_node_text(aes(label = str_wrap(label, 20, whitespace_only = F),
                       x = x + case_when(label %in% x_nudge_nodes ~ x_nudge, 
                                         label == "Chronic kidney disease" ~ 0.25,
                                         T ~ 0),
                       y = y + ifelse(label %in% neg_y_nudge_nodes, -y_nudge, y_nudge), 
                       vjust = case_when(label %in% neg_y_nudge_nodes ~ 1, 
                                         label == "Chronic kidney disease" ~ 0.5,
                                         T ~ 0)),
                   size = label_size,
                   show.legend = F) +
    theme_void() 
  
  
  clust_plot <- clust_plot +
    labs(title = paste0("Cluster ", clust_num)) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, vjust = 0, size = 8, face = "bold"), 
          plot.subtitle = element_text(hjust = 0.5, vjust = 1, size = 7, face = "bold")) +
    scale_x_continuous(expand = ggplot2::expansion(mult = 0.25)) +
    scale_y_continuous(expand = ggplot2::expansion(mult = 0.35))
  
  if(jacc_stability == T){
    jacc = clust_links %>% pull(jaccard) %>% unique() %>% round(., 2)
    
    clust_plot <- clust_plot + 
      labs(subtitle = paste0("Jaccard : ", format(jacc, nsmall = 2)))
  }
  
  assign(paste0("clust_plot_", clust_num), clust_plot, envir = .GlobalEnv)
  
}
plot_network_cluster_6(links_clust, nodes_clust_plot, 6)

row_1 <- cowplot::plot_grid(clust_plot_1, clust_plot_2, clust_plot_3,  nrow = 1) 
row_2 <- cowplot::plot_grid(clust_plot_4, clust_plot_5, clust_plot_6, nrow = 1) +
  theme(plot.margin = margin(t = 0.25, b = 0, unit = "mm"))

plot_clusters <- cowplot::plot_grid(row_1, row_2, ncol = 1) +
  theme(plot.margin = margin(t = 0.25, b = 0, l = 12, r =  12, unit = "mm"))


plot_disease_network_clust_all <- cowplot::plot_grid(plot_disease_network_clust, plot_clusters, 
                                                     rel_heights = c(1, 0.9), ncol = 1,
                                                     labels = c("A: Meta-analytic network with stable disease clusters", 
                                                                "B: Stable disease clusters"), 
                                                     hjust = 0, label_fontface = "plain", 
                                                     label_size = 9, vjust = c(1,-1), label_x = 0.05) + 
  theme(plot.margin = margin(2.5,5,2,5, unit = "mm"))


ggsave(paste0(dir_output, "Fig5_disease_network_clusters.pdf"), plot_disease_network_clust_all, bg = "white", dpi = 400,
       height = 200, width = 180, units = "mm")



#----- Supp Table 1: all cluster results ----

#- Results:  stability estimates -- 
cluster_stability_results <- read.csv(paste0(dir_data, "results_cluster_stability.csv")) %>%
  arrange(-jaccard_stability) %>%
  mutate(stable_cluster = 1:43) %>%
  select(cluster, stable_cluster, jaccard_stability)

clusters_stable_all <- read.csv(paste0(dir_data, "results_clusters.csv")) %>%
  left_join(cluster_stability_results) %>%
  mutate(cluster_id = as.ordered(cluster),
         stable_cluster_id = as.ordered(stable_cluster)) %>%
  arrange(stable_cluster_id) 


links_clust_all <- links %>%
  left_join(clusters_stable_all) %>%
  arrange(stable_cluster_id) %>%
  mutate(cluster = stable_cluster_id) %>%
  select(-stable_cluster_id) %>% 
  left_join(disease_labels, by = join_by("from" == "link_name")) %>%
  select(-from) %>%
  rename(from = label) %>%
  left_join(disease_labels, by = join_by("to" == "link_name")) %>%
  select(-to) %>%
  rename(to = label) %>%
  mutate(cluster_links = paste(from, "+", to))


cluster_table_data <- select(links_clust_all, cluster, jaccard_stability, cluster_links) %>%
  mutate(jaccard_stability = as.character(round(jaccard_stability, 2))) %>%
  pivot_longer(jaccard_stability:cluster_links) %>%
  mutate(name = case_when(name == "jaccard_stability" ~ "Jaccard stability", 
                          name == "cluster_links" ~ "Cluster links"),
         cluster = paste0("Cluster ", cluster)) %>%
  distinct() %>%
  mutate(name = case_when(name == "Cluster links" & lag(name) != "Cluster links" ~ name, 
                          name == "Cluster links" & lag(name) == "Cluster links" ~ NA,
                          T ~ name)) %>%
  filter(cluster != "Cluster NA")

#- create table
cluster_table <- cluster_table_data %>%
  as_grouped_data("cluster") %>%
  as_flextable(hide_grouplabel = T) %>%
  set_header_labels(name = "", value = "") %>%
  add_header_lines("Supplementary Table 2: Meta-analytic disease clusters") %>%
  merge_at(i = 1:2, part = "header") %>%
  hline(i = 1, part = "header") %>%
  bold(j = 1, i = ~!is.na(cluster), part = "body") %>%
  fontsize(size = 9, part = "body") %>%
  fontsize(size = 10, part = "header") %>%
  fontsize(size = 10, i = ~!is.na(cluster), part = "body") %>%
  add_footer_lines("Summary information for all disease clusters identified from meta-analytic disease network.") %>%
  fontsize(size = 10, part = "footer") %>%
  autofit() %>%
  width(j = 1, width = 1) 

save_as_docx(cluster_table, path = paste0(dir_output, "supp_table2_all_cluster_stats.docx"))


#----- Supp Fig 2: Sensitivity analyses low ROB network ------

#-- load meta-network for studies with ROB low/some concerns 

links_ROB <- read_csv(paste0(dir_data, "results_links_ROB_sensitivity.csv"))

missing_dis <- anti_join(layout, nodes, by = join_by(name == dis)) %>%
  select(name) %>% rename(dis = name)

net_ROB <- graph_from_data_frame(links_ROB, directed = F, vertices = nodes)

#-- create plots for low ROB network vs full network --- 

plot_disease_network_ROB <- plot_network(net_ROB) +
  labs(edge_color = "Lift",
       edge_alpha = "Lift", edge_width = "Lift") + 
  theme(plot.title = element_text(size = 12, hjust = 0, face = "plain"),
        plot.caption = element_text(size = 9.5, hjust = 0),
        text = element_text(size  = 9),
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(size = 8)) 

title_sfig2_1 <- ggdraw() +
  draw_label("Full network", 
             fontface = "bold", 
             size = 10,
             x = 0.05,
             hjust = 0, vjust = 1)

title_sfig2_2 <- ggdraw() +
  draw_label("Sensitivity analysis: Low ROB network", 
             fontface = "bold", 
             size = 10,
             x = 0.05,
             hjust = 0, vjust = 1)

row_1 <- plot_grid(title_sfig2_1, 
                   plot_disease_network, 
                   ncol = 1, rel_heights = c(0.05, 1))
row_2 <- plot_grid(title_sfig2_2, 
                   plot_disease_network_ROB, 
                   ncol = 1, rel_heights = c(0.05, 1))

plot_nets_ROB <- plot_grid(row_1, row_2, nrow = 2)
plot_nets_ROB <- plot_grid(plot_nets_ROB, #leg, 
                           rel_widths = c(1, 0.1)) +
  theme(plot.margin = margin(l = 1, r = 2, unit = "mm"))

caption_sfig2 <- paste0(
  "Meta-analytic network from the full network dataframe (top; n = 79), compared to network restricted to data to
studies with Risk of Bias ratings of 'Low' or 'Some Concerns' (bottom; n = 53)"
)

caption_sfig2 <- ggdraw() +
  draw_label(caption_sfig2, x = 0.025, size = 11,
             hjust = 0, vjust = 0.5)

title_sfig2 <- ggdraw() +
  draw_label("Supplementary Figure 2: Sensitivity analysis", 
             size = 12,  x = 0.025,
             hjust = 0, vjust = 1)


plot_nets_ROB <- plot_grid(title_sfig2, NULL, plot_nets_ROB, caption_sfig2, ncol = 1, rel_heights = c(0.05, 0.005, 1, 0.05))

ggsave(paste0(dir_output, "Supp_Fig2_disease_network_sensitivity.png"), plot_nets_ROB, bg = "white", dpi = 400,
       height = 10, width = 8)


#----- Supp Fig 3: clusters in sensitivity analyses ------

#- rename old plots and remove --
vars <- paste0("clust_plot_", 1:6)
prefix <- "orig_"

lapply(vars, function(var) {
  new_name <- paste0(prefix, var) 
  assign(new_name, get(var), envir = .GlobalEnv) 
  rm(list = var, envir = .GlobalEnv) 
})


#- filter new links data and rename clusters
new_clust_nums <- c(1, 3, 32, 27, 2, 8)

clusters_ROB <-  read_csv(paste0(dir_data, "results_cluster_ROB_sensitivity.csv")) %>%
  rename(from = node1, to = node2)

sensitivity_cluster_matches <- read_csv(paste0(dir_data, "results_sensitivity_cluster_matches.csv"))

links_ROB_clust <-
  links_ROB %>% left_join(clusters_ROB)

links_ROB_clust <- 
  links_ROB_clust %>% 
  filter(cluster %in% new_clust_nums) %>%
  mutate(new_ROB_cluster_id = case_when(cluster == 1 ~ 1,
                                        cluster == 3 ~ 2,
                                        cluster == 32 ~ 3,
                                        cluster == 27 ~ 4,
                                        cluster == 2 ~ 5,
                                        cluster == 8 ~ 6)) %>%
  select(-cluster) %>%
  rename(stable_cluster_id = new_ROB_cluster_id) %>%
  left_join(sensitivity_cluster_matches, relationship = "many-to-many") %>%
  rename(jaccard = jacc) 


#-- generate new cluster plots -
plot_network_cluster(links_ROB_clust, nodes_clust_plot, 1, neg_y_nudge_nodes = c("Diabetes"), x_nudge_nodes = c("Chronic kidney disease"), x_nudge = -0.05)
plot_network_cluster(links_ROB_clust, nodes_clust_plot, 2, neg_y_nudge_nodes = c("Diabetes", "Hypertension"), x_nudge_nodes = c("Hypertension"), x_nudge = -0.1, link_strength = 0.05)
plot_network_cluster(links_ROB_clust, nodes_clust_plot, 5, neg_y_nudge_nodes = c("Diabetes", "Chronic kidney disease", "COPD", "Osteoarthritis"), y_nudge = 0.1,
                     x_nudge_nodes = c("Chronic kidney disease"), x_nudge = 0.05)

#- custom function for some plot tweaks
plot_network_cluster_3 <- function(link_df, node_df, clust_num){
  link_strength = 0.05
  seq_links = seq_link 
  round_seq_by = 1
  label_size = 2.25
  col_clusters = T 
  node_size = 2.8
  x_nudge = 0 
  y_nudge = 0.09
  neg_y_nudge_nodes = NULL
  jacc_stability = T
  
  clust_links <- link_df %>%
    filter(stable_cluster_id == clust_num)
  
  clust <- link_df %>%
    pivot_longer(cols = c(from, to), names_to = "node", values_to = "disease") %>% 
    dplyr::select(-node, -weight) 
  
  if(jacc_stability == T){
    clust <- clust %>%
      dplyr::select(-jaccard)
  }
  
  clust <- clust %>%
    distinct() %>%
    filter(stable_cluster_id == clust_num) %>% pull(disease)
  
  clust_nodes <- node_df %>%
    filter(dis %in% clust)
  
  clust_net <- graph_from_data_frame(clust_links, directed = T, vertices = clust_nodes)
  
  clust_col = "#228833"
  
  
  clust_plot <- 
    ggraph(clust_net, layout = "manual", x = clust_nodes$x, y = clust_nodes$y) +
    geom_edge_arc(edge_width = 0.5,
                  strength = link_strength, alpha = 0.75, color = clust_col) + 
    geom_node_point(size = {{node_size}},
                    shape = 21, fill = clust_col, color = "00") +
    geom_node_text(aes(label = str_wrap(label, 20, whitespace_only = F),
                       x = x + -0.05,
                       y = y + case_when(label %in% neg_y_nudge_nodes ~ -y_nudge, 
                                         label == "Schizophrenia" ~ 0.12,
                                         T ~ y_nudge), 
                       vjust = ifelse(label %in% neg_y_nudge_nodes, 0.25, 0)),
                   size = label_size,
                   show.legend = F) +
    theme_void() 
  
  
  clust_plot <- clust_plot +
    labs(title = paste0("Cluster ", clust_num)) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, vjust = 0, size = 8, face = "bold"), 
          plot.subtitle = element_text(hjust = 0.5, vjust = 1, size = 7, face = "bold")) +
    scale_x_continuous(expand = ggplot2::expansion(mult = 0.25)) +
    scale_y_continuous(expand = ggplot2::expansion(mult = 0.35))
  
  if(jacc_stability == T){
    jacc = clust_links %>% pull(jaccard) %>% unique() %>% round(., 2)
    
    clust_plot <- clust_plot + 
      labs(subtitle = paste0("Jaccard : ", format(jacc, nsmall = 2)))
  }
  
  assign(paste0("clust_plot_", clust_num), clust_plot, envir = .GlobalEnv)
  
}
plot_network_cluster_3(links_ROB_clust, nodes_clust_plot, 3)

plot_network_cluster_4 <- function(link_df, node_df, clust_num){
  link_strength = 0.1
  seq_links = seq_link 
  round_seq_by = 1
  label_size = 2.25
  col_clusters = T 
  node_size = 2.8
  x_nudge = 0 
  y_nudge = 0.09
  neg_y_nudge_nodes = c("Mood & anxiety disorders")
  jacc_stability = T
  
  clust_links <- link_df %>%
    filter(stable_cluster_id == clust_num)
  
  clust <- link_df %>%
    pivot_longer(cols = c(from, to), names_to = "node", values_to = "disease") %>% 
    dplyr::select(-node, -weight) 
  
  if(jacc_stability == T){
    clust <- clust %>%
      dplyr::select(-jaccard)
  }
  
  clust <- clust %>%
    distinct() %>%
    filter(stable_cluster_id == clust_num) %>% pull(disease)
  
  clust_nodes <- node_df %>%
    filter(dis %in% clust)
  
  clust_net <- graph_from_data_frame(clust_links, directed = T, vertices = clust_nodes)
  
  clust_col = "#CCBB44"
  
  
  clust_plot <- 
    ggraph(clust_net, layout = "manual", x = clust_nodes$x, y = clust_nodes$y) +
    geom_edge_arc(edge_width = 0.5,
                  strength = link_strength, alpha = 0.75, color = clust_col) + 
    geom_node_point(size = {{node_size}},
                    shape = 21, fill = clust_col, color = "00") +
    geom_node_text(aes(label = str_wrap(label, 20, whitespace_only = F),
                       x = x + case_when(!label %in% neg_y_nudge_nodes ~ -0.1, 
                                         #label == "COPD" ~ -2,
                                         T ~ 0.075),
                       y = y + ifelse(label %in% neg_y_nudge_nodes, -y_nudge - 0.2, y_nudge), 
                       vjust = case_when(label %in% neg_y_nudge_nodes ~ 0.25, 
                                         label == "COPD" ~ 1.5,
                                         T ~ 0),
                       hjust = ifelse(label == "COPD", 1, 0.5)),
                   size = label_size,
                   show.legend = F) +
    theme_void() 
  
  
  clust_plot <- clust_plot +
    labs(title = paste0("Cluster ", clust_num)) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, vjust = 0, size = 8, face = "bold"), 
          plot.subtitle = element_text(hjust = 0.5, vjust = 1, size = 7, face = "bold")) +
    scale_x_continuous(expand = ggplot2::expansion(mult = 0.25)) +
    scale_y_continuous(expand = ggplot2::expansion(mult = 0.35))
  
  if(jacc_stability == T){
    jacc = clust_links %>% pull(jaccard) %>% unique() %>% round(., 2)
    
    clust_plot <- clust_plot + 
      labs(subtitle = paste0("Jaccard : ", format(jacc, nsmall = 2)))
  }
  
  assign(paste0("clust_plot_", clust_num), clust_plot, envir = .GlobalEnv)
  
}
plot_network_cluster_4(links_ROB_clust, nodes_clust_plot, 4)

plot_network_cluster_6 <- function(link_df, node_df, clust_num, x_nudge_nodes = c("Ischemic heart disease", "Diabetes"), x_nudge = -0.1){
  link_strength = 0.1
  auto_nudge = F
  seq_links = seq_link
  round_seq_by = 1
  label_size = 2.25
  col_clusters = T
  node_size = 2.8
  y_nudge = 0.1
  neg_y_nudge_nodes = c("Dementia")
  jacc_stability = T
  
  clust_links <- link_df %>%
    filter(stable_cluster_id == clust_num)
  
  clust <- link_df %>%
    pivot_longer(cols = c(from, to), names_to = "node", values_to = "disease") %>% 
    dplyr::select(-node, -weight) 
  
  if(jacc_stability == T){
    clust <- clust %>%
      dplyr::select(-jaccard)
  }
  
  clust <- clust %>%
    distinct() %>%
    filter(stable_cluster_id == clust_num) %>% pull(disease)
  
  clust_nodes <- node_df %>%
    filter(dis %in% clust)
  
  
  if(auto_nudge == T){
    if(length(unique(clust_nodes$dis)) == 3) {
      y_nudge = 0.015
      #link_strength = 0.05
    }
  }
  
  clust_net <- graph_from_data_frame(clust_links, directed = T, vertices = clust_nodes)
  
  clust_col = "#698896"
  
  if(col_clusters == T){
    clust_col <- if (clust_num == 1){
      "#4477AA"
    }  else if (clust_num == 2){
      "#EE6677"
    }   else if (clust_num == 3){
      "#228833"
    } else if (clust_num == 4){
      "#CCBB44"
    } else if (clust_num == 5){
      "#66CCEE"
    }  else if (clust_num == 6){
      "#AA3377" 
    } 
    
  } 
  
  
  clust_plot <- 
    ggraph(clust_net, layout = "manual", x = clust_nodes$x, y = clust_nodes$y) +
    geom_edge_arc(edge_width = 0.5,
                  strength = link_strength, alpha = 0.75, color = clust_col) + 
    geom_node_point(size = {{node_size}},
                    shape = 21, fill = clust_col, color = "00") +
    geom_node_text(aes(label = str_wrap(label, 20, whitespace_only = F),
                       x = x + case_when(label %in% x_nudge_nodes ~ x_nudge, 
                                         label == "Chronic kidney disease" ~ 0.25,
                                         T ~ 0),
                       y = y + ifelse(label %in% neg_y_nudge_nodes, -y_nudge, y_nudge), 
                       vjust = case_when(label %in% neg_y_nudge_nodes ~ 1, 
                                         label == "Chronic kidney disease" ~ 0.5,
                                         T ~ 0)),
                   size = label_size,
                   show.legend = F) +
    theme_void() 
  
  
  clust_plot <- clust_plot +
    labs(title = paste0("Cluster ", clust_num)) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, vjust = 0, size = 8, face = "bold"), 
          plot.subtitle = element_text(hjust = 0.5, vjust = 1, size = 7, face = "bold")) +
    scale_x_continuous(expand = ggplot2::expansion(mult = 0.25)) +
    scale_y_continuous(expand = ggplot2::expansion(mult = 0.35))
  
  if(jacc_stability == T){
    jacc = clust_links %>% pull(jaccard) %>% unique() %>% round(., 2)
    
    clust_plot <- clust_plot + 
      labs(subtitle = paste0("Jaccard : ", format(jacc, nsmall = 2)))
  }
  
  assign(paste0("clust_plot_", clust_num), clust_plot, envir = .GlobalEnv)
  
}
plot_network_cluster_6(links_ROB_clust, nodes_clust_plot, 6)



#- compile figure
row_1 <- cowplot::plot_grid(orig_clust_plot_1 + labs(subtitle = element_blank()), clust_plot_1, NULL, orig_clust_plot_2 + labs(subtitle = element_blank()), clust_plot_2,
                            nrow = 1, 
                            rel_widths = c(1, 1, 0.25, 1, 1)) +
  theme(plot.margin = margin(r = 2, unit = "mm"))

row_2 <- cowplot::plot_grid(orig_clust_plot_3 + labs(subtitle = element_blank()), clust_plot_3, NULL, orig_clust_plot_4 + labs(subtitle = element_blank()), clust_plot_4,
                            nrow = 1,
                            rel_widths = c(1, 1, 0.25, 1, 1)) +
  theme(plot.margin = margin(r = 2, unit = "mm"))

row_3 <- cowplot::plot_grid(orig_clust_plot_5 + labs(subtitle = element_blank()), clust_plot_5, NULL, orig_clust_plot_6 + labs(subtitle = element_blank()), clust_plot_6, 
                            nrow = 1,
                            rel_widths = c(1, 1, 0.25, 1, 1)) +
  theme(plot.margin = margin(r = 2, unit = "mm"))


title_sfig1 <- ggdraw() +
  draw_label("Original analysis", fontface = "bold",
             hjust = 0.5, vjust = 0, size = 10)

title_sfig2 <- ggdraw() +
  draw_label("Sensitivity analysis", fontface = "bold",
             hjust = 0.5, vjust = 0, size = 10)

title_row <- plot_grid(title_sfig1, title_sfig2, NULL, title_sfig1, title_sfig2,
                       nrow = 1,
                       rel_widths = c(1, 1, 0.25, 1, 1))

caption_sfig <- paste0("Meta-analytic disease clusters from the full dataframe ('Original analysis'), compared to the most similar clusters identified 
in a sensitivity analysis restricted to data to studies with Risk of Bias (ROB) ratings of 'Low' or 'Some Concerns' (n = 53). 
The Jaccard index scores the similarity between the original cluster and its' closest match in the low ROB network."
)

caption_sfig <- ggdraw() +
  draw_label(caption_sfig, x = 0.025, size = 10,
             hjust = 0, vjust = 0.4)

title_sfig3 <- ggdraw() +
  draw_label("Supplementary Figure 3: Sensitivity analysis - disease cluster replicability in low ROB studies", 
             size = 11,  x = 0.015,
             hjust = 0, vjust = 1)

plot_clusters_ROB <- cowplot::plot_grid(title_row, row_1, row_2, row_3, ncol = 1, rel_heights = c(0.1, 1, 1, 1)) +
  theme(plot.margin = margin(t = 0.5, unit = "cm"))
plot_clusters_ROB <- plot_grid(title_sfig3, plot_clusters_ROB, caption_sfig, ncol =1, rel_heights = c(0.05, 1, 0.05))

ggsave(paste0(dir_output, "Supp_Fig3_disease_clusters_ROB_sensitivity.png"), plot_clusters_ROB, bg = "white", dpi = 400,
       height = 8, width = 8)


#----- Figure 6: Cluster study characteristics ------
#-- Goal is to describe characteristics of studies where clusters have been observed
#-- Go through each line of meta-analytic dataframe and find matches of clusters that have the anchor disease + n-1 other diseases in the cluster
#-- Then describe the characteristics of studies with matching clusters

#-- create unique identified for meta-analytic cluster data
meta_dis <- meta_dis %>% 
  mutate(clustid = paste0(studyid, "_", subgroup_num, "_", row_number()))

clusters_stable_wide <- clusters_stable %>%
  mutate(member = 1) %>%
  pivot_longer(cols = c(from, to), values_to = "dis") %>%
  select(-name) %>% unique() %>%
  pivot_wider(id_cols = stable_cluster_id, names_from = dis, values_from = member) %>%
  arrange(stable_cluster_id)

meta_dis_cluster_replicability <- map(clusters_stable$stable_cluster, find_cluster_match_jacc) %>% list_rbind() %>%
  mutate(stable_cluster_id = as.ordered(stable_cluster_id)) 


#-- calculate summary info about cluster matches
cluster_replicability_denom <- meta_dis_cluster_replicability %>%
  group_by(stable_cluster_id) %>%
  summarise(n_total = n())

cluster_replicability_mean <- meta_dis_cluster_replicability %>%
  filter(cluster_match == "y") %>%
  group_by(stable_cluster_id) %>%
  summarise(n_matching_clusters = n(),
            jacc_max = max(jacc),
            jacc_mean = mean(jacc),
            jacc_sd = sd(jacc),
            jacc_med = median(jacc),
            jacc_IQR_low = quantile(jacc, 0.25),
            jacc_IQR_high = quantile(jacc, 0.75)) %>%
  left_join(cluster_replicability_denom) %>%
  mutate(perc_total_clusters = n_matching_clusters/n_total*100)



#-- calculate descriptive info for clusters that match SR clusters
meta_dis_cluster_replicability_matches <- meta_dis_cluster_replicability %>% filter(cluster_match == "y") 

cluster_replicability_matches_category <- meta_dis_cluster_replicability_matches %>%
  group_by(stable_cluster_id, category) %>% tally() %>%
  left_join(select(cluster_replicability_mean, stable_cluster_id, n_matching_clusters)) %>%
  mutate(perc_total = n/n_matching_clusters*100,
         feature = "method") %>%
  ungroup()

cluster_replicability_matches_age <- meta_dis_cluster_replicability_matches %>%
  mutate(age_bin_2 = case_when(age_bin_2 %in% c("All ages\n(0 - 100+)", "Adults\n(16+)") ~ "All ages", 
                               age_bin_2 == "Middle-aged\nand older\n(40+)" ~ "40+",
                               T ~ age_bin_2)) %>%
  group_by(stable_cluster_id, age_bin_2) %>% tally() %>%
  left_join(select(cluster_replicability_mean, stable_cluster_id, n_matching_clusters)) %>%
  mutate(perc_total = n/n_matching_clusters*100, 
         feature = "age") %>%
  rename(category = age_bin_2)  %>%
  ungroup()


cluster_replicability_matches_sex <-meta_dis_cluster_replicability_matches %>%
  group_by(stable_cluster_id, subgroup_sex) %>% tally() %>%
  left_join(select(cluster_replicability_mean, stable_cluster_id, n_matching_clusters)) %>%
  mutate(perc_total = n/n_matching_clusters*100, feature = "sex") %>%
  rename(category = subgroup_sex)  %>%
  ungroup()

cluster_replicability_matches <- bind_rows(cluster_replicability_matches_category, cluster_replicability_matches_age, cluster_replicability_matches_sex)

#rm(cluster_replicability_matches_category, cluster_replicability_matches_age, cluster_replicability_matches_sex)

cluster_replicability_matches <- cluster_replicability_matches %>%
  mutate(category = str_to_sentence(category)) %>%
  mutate(category = factor(category, ordered = T, levels = 
                             c("Factor-based", "Cluster-based", 
                               "F", "M", "Both", 
                               "60+", "40+",  "40-59", "0-39", 
                               "All ages")))



clust_1_bars <- stack_clust_and_bar_plots(1) + 
  annotate("text", x = 0.69, y = 0.08, label = "0-39\n6%", size = 2.05, fontface = "bold", color = "black") + 
  annotate("segment", x = 0.69, xend = 0.695, y = 0.13, yend = 0.165, linewidth = 0.3) +
  annotate("text", x = 0.76, y = 0.08, label = "40+\n4%", size = 2.05, fontface = "bold", color = "black") +
  annotate("segment", x = 0.76, xend = 0.755, y = 0.13, yend = 0.165, linewidth = 0.3)

clust_2_bars <- stack_clust_and_bar_plots(2) + 
  annotate("text", x = 0.76, y = 0.08, label = "0-39\n2%", size = 2.05, fontface = "bold", color = "black") + 
  annotate("segment", x = 0.76, xend = 0.773, y = 0.13, yend = 0.165, linewidth = 0.3) +
  annotate("text", x = 0.83, y = 0.08, label = "40+\n4%", size = 2.05, fontface = "bold", color = "black") +
  annotate("segment", x = 0.83, xend = 0.81, y = 0.13, yend = 0.165, linewidth = 0.3)
clust_3_bars <- stack_clust_and_bar_plots(3)
clust_4_bars <- stack_clust_and_bar_plots(4)
clust_5_bars <- stack_clust_and_bar_plots(5)
clust_6_bars <- stack_clust_and_bar_plots(6)

clust_features <- plot_grid(
  plot_grid(clust_1_bars, NULL, clust_2_bars, nrow = 1, rel_widths = c(1, 0.05, 1)),
  plot_grid(clust_3_bars, NULL, clust_4_bars, nrow = 1, rel_widths = c(1, 0.05, 1)),
  plot_grid(clust_5_bars, NULL, clust_6_bars, nrow = 1, rel_widths = c(1, 0.05, 1)),
  ncol = 1
)


ggsave(paste0(dir_output, "Fig6_cluster_features.pdf"), clust_features, bg = "white", dpi = 400,
       height = 140, width = 180, units = "mm")

