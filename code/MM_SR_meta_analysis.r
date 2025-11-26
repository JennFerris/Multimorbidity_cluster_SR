# ****************************************************************************************************************************
#   Accompanying code for Ferris, et al. 2025 "A systematic review and meta-analysis of disease clustering in multimorbidity  
#   AUTHOR: Jennifer Ferris, Postdoctoral fellow (jennifer_ferris@sfu.ca)
#   LAST UPDATED: November 21, 2025
# ****************************************************************************************************************************
#   Code to perform disease-level cluster meta-analysis 
# ****************************************************************************************************************************

#--- load packages and set directories -----
pkg <- c("dplyr",
         "tidyr",
         "readr",
         "tibble",
         "purrr",
         "ggplot2",
         "arrow",
         "stats",
         "clustAnalytics",
         "igraph",
         "linkcomm", 
         "fpc")

sapply(pkg, require, character.only = T)


dir_data <- "../data/" #- directory where raw data files are stored
dir_output <- "../output/"  #- directory to output results

if (!dir.exists(dir_output)) {
  dir.create(dir_output)
} 

#- load functions 
source("00_functions.R") 

#--- Load and wrangle network data ------
meta_dis <- read_csv(paste0(dir_data, "data_cluster_meta_analytics.csv"))  %>%
  mutate(across(everything(), ~ifelse(. == "n/a", NA, .))) %>%
  mutate(age_bin_2 = factor(age_bin_2, ordered = T, levels = c("0-39", "40-59", "60+", 
                                                               "All ages\n(0 - 100+)", 
                                                               "Adults\n(16+)", 
                                                               "Middle-aged\nand older\n(40+)")))
layout <- read_csv(paste0(dir_data, "network_layout.csv")) %>%
  arrange(name, .locale = "en")

#- load node data
nodes <- layout %>% 
  select(name, label) %>%
  rename(dis = name) %>%
  mutate(pct_dis = 100) 

#- Generate network link data
#- calculate disease co-occurrence across clusters 
diseases = colnames(meta_dis)[12:31]

combo_list_pairs <- 
  combn(diseases, 2) %>% 
  as.data.frame.table() %>%
  pivot_wider(names_from = Var1, values_from = Freq) %>%
  select(-Var2) %>%
  rename(Disease1 = A, Disease2 = B) %>%
  mutate(Disease1 = as.character(Disease1), Disease2 = as.character(Disease2))

links.2 <- calc_network_metrics(meta_dis) 

links <- links.2 %>%
  select(dis1, dis2, weight) %>%
  rename(from = dis1, to = dis2) %>%
  filter(weight > 1)

net <- graph_from_data_frame(links, directed = F, vertices = nodes)

#- save link data
write_csv(links, paste0(dir_data, "results_network_link_data.csv"))


#--- Test overall network stability per Arriata (2021) with Louvain algorithm ------

network_assessment <- data.frame(evaluate_significance_r(net, alg_list = list(Louvain = cluster_louvain), n_reps = 1000))

network_assessment$feature <- row.names(network_assessment)

#- output with directional arrows for ease of viewing
filter(network_assessment, feature %in% c("conductance", "norm cut", "max ODF", "average ODF", "flake ODF", "modularity")) %>%
  select(feature, Louvain, Louvain_r, Louvain_percentile) %>%
  mutate(feature = case_when(
    feature == "conductance" ~ paste("Conductance", sprintf('\u2193')),
    feature == "norm cut" ~ paste("Normalized Cut", sprintf('\u2193')),
    feature == "max ODF" ~ paste("Maximum ODF", sprintf('\u2193')),
    feature == "average ODF" ~ paste("Average ODF", sprintf('\u2193')),
    feature == "flake ODF" ~ paste("Flake ODF", sprintf('\u2193')),
    feature == "modularity" ~ paste("Modularity", sprintf('\u2191')),
    feature == "clustering coef" ~ paste("Clustering coefficient", sprintf('\u2191')),
    T ~ feature))



#--- Determine linkcomm dendrogram cut-height ----------

#- generate LC dendrogram
LC <- getLinkCommunities(as.data.frame(links), directed = F,
                         hcmethod = "ward.D2", removetrivial = F)


#- find optimal cutheight
LC_stats <-
  LinkCommunitiesCutStats(as.data.frame(links), 60, directed = F, removetrivial = F, hcmethod = "ward.D2") %>%
  rownames_to_column(var = "measure") %>%
  #filter(measure %in% c("cluster.number", "within.cluster.ss", "avg.silwidth", "ch", "dunn2", "sindex")) %>%
  pivot_longer(cols = starts_with("Test"), names_prefix = "Test.", names_to = "cutheight") %>%
  pivot_wider(names_from = "measure") %>%
  arrange(-ch) %>%
  mutate(pseudo.F = average.between/average.within, cut.height = as.numeric(cutheight))

elbow_plot_2(LC_stats, within.cluster.ss) 
elbow_plot_2(LC_stats, avg.silwidth) 
elbow_plot_2(LC_stats, ch) 
elbow_plot_2(LC_stats, pseudo.F) 

#-- based on plots, there is a range of values where ch peaks
View(LC_stats %>% dplyr::filter(ch >= 0.49) %>% select(cutheight:sindex, pseudo.F))

#-- at these ranges, 0.249 gives smaller clusters and has good score on other metrics
newcutheight <- 0.249
LC_net <- newLinkCommsAt(LC, newcutheight)

#- plots to explore cluster members
plotLinkCommDend(LC_net, labels = T, droptrivial = F)
plotLinkCommMembers(LC_net, maxclusters = 50, nodes = names(LC_net$numclusters))

#- save cluster outputs
clusters <- LC_net$edges %>%
  rename(from = node1, to = node2)

write_csv(clusters, paste0(dir_data, "results_clusters.csv"))


#--- Perform cluster stability testing with bootstrapped networks ------
#- note: this section takes several hours to run

#-- input dataframe is disease data from meta-clustering
meta_dis_clust <- meta_dis %>% select(asthma:stroke)

#- get list of links
links_key <- calc_network_metrics(meta_dis_clust) %>% select(dis1, dis2) %>%
  arrange(dis1, dis2, .locale = "en") %>%
  rename(node1 = dis1, node2 = dis2)


#- note: the output is numbered by cluster; 1st partition == cluster 1
start <- Sys.time()

cluster_stability <- clusterboot_linkcomm(
  meta_dis_clust, 
  B = 100, 
  distances = F, 
  datatomatrix = F, 
  clustermethod = cluster_network
)

Sys.time() - start

cluster_stability_results <- data.frame(
  cluster = seq(1:43), 
  jaccard_stability = cluster_stability$bootmean
  ) %>% 
  mutate(cluster = as.character(cluster))

write_csv(cluster_stability_results, paste0(dir_data, "results_cluster_stability.csv"))



#--- sensitivity analysis in low ROB network --------
#-- filter meta-network for studies with ROB low/some concerns 

meta_dis_ROB <- meta_dis %>% 
  filter(risk_of_bias_overall %in% c("Low Risk of Bias", "Some Concerns"))

links_ROB <- calc_network_metrics(meta_dis_ROB) 

links_ROB <- links_ROB %>%
  select(dis1, dis2, weight) %>%
  rename(from = dis1, to = dis2) 

write_csv(links_ROB, paste0(dir_data, "results_links_ROB_sensitivity.csv"))

#- cluster low ROB network - 
LC_net_ROB <- linkcomm_optimal_cutheight(links_ROB, 60)
clusters_ROB <- LC_net_ROB$edges
write_csv(clusters_ROB, paste0(dir_data, "results_cluster_ROB_sensitivity.csv"))

#-- comparing clusters ---

#- create clusterid for each cluster to test
clusters_ROB <-  clusters_ROB %>%
  rename(from = node1, to = node2)

links_ROB_clust <-
  links_ROB %>% left_join(clusters_ROB)

cluster_stability_results <- read.csv(paste0(dir_data, "results_cluster_stability.csv")) %>%
  arrange(-jaccard_stability) %>%
  mutate(stable_cluster = 1:43) %>%
  select(cluster, stable_cluster, jaccard_stability)

clusters_stable <- read.csv(paste0(dir_data, "results_clusters.csv")) %>%
  left_join(cluster_stability_results) %>%
  filter(jaccard_stability > 0.5) %>%
  mutate(cluster_id = as.ordered(cluster),
         stable_cluster_id = as.ordered(stable_cluster))

#- pivot cluster data to wide format for replicability testing 
clusters_ROB_wide <- 
  clusters_ROB %>%
  mutate(member = 1) %>%
  pivot_longer(cols = c(from, to), values_to = "dis") %>%
  select(-name) %>% unique() %>%
  pivot_wider(id_cols = cluster, names_from = dis, values_from = member) %>%
  rename(lowROB_cluster_id = cluster)

clusters_stable_wide <- clusters_stable %>%
  mutate(member = 1) %>%
  pivot_longer(cols = c(from, to), values_to = "dis") %>%
  select(-name) %>% unique() %>%
  pivot_wider(id_cols = stable_cluster_id, names_from = dis, values_from = member) %>%
  arrange(stable_cluster_id)


#-- jaccard index tested ONLY for diseases within the SR cluster 
find_closest_clusters <- function(cluster_id){
  #- for each stable cluster....
  #-- extract cluster diseases
  clust_diseases <- clusters_stable_wide %>% 
    filter(stable_cluster_id == {{cluster_id}}) %>% 
    select(-stable_cluster_id) %>%
    select_if(~.x == as.numeric(1)) %>% names
  
  #-- score jaccard index for each cluster in low ROB network 
  cluster_results_jaccard_sensitivity <- data.frame(stable_cluster_id = numeric(0), lowROB_cluster_id = numeric(0), jacc = numeric(0))
  
  jaccard_candidate_clusters <- function(test_cluster_id){
    
    #- extract cluster diseases
    test_dis <- clusters_ROB_wide[test_cluster_id,] %>%
      select(-lowROB_cluster_id) %>%
      select_if(~.x == as.numeric(1)) %>% names
    
    #-- calculate jaccard index for each cluster
    add_row(cluster_results_jaccard_sensitivity, 
            stable_cluster_id = as.numeric({{cluster_id}}), 
            lowROB_cluster_id = clusters_ROB_wide[test_cluster_id,] %>% pull(lowROB_cluster_id), 
            jacc = jaccard(clust_diseases, test_dis))
  }
  
  cluster_results_jaccard_sensitivity <- map(1:nrow(clusters_ROB_wide), jaccard_candidate_clusters) %>% list_rbind() 
  
  #-- find highest matching cluster for each cluster
  cluster_results_jaccard_sensitivity %>% group_by(stable_cluster_id) %>% slice_max(order_by = jacc, with_ties = T)
  #cluster_results_jaccard_sensitivity
  
}

sensitivity_cluster_matches <- map(clusters_stable_wide$stable_cluster_id, find_closest_clusters) %>% list_rbind()

write_csv(sensitivity_cluster_matches, paste0(dir_data, "results_sensitivity_cluster_matches.csv"))

