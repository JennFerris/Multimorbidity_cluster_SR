
#- count variables
count_vals <- function(mydf, ..., caps = T){
  
  group_vars <- quos(...)
  
  mydf <- mydf %>%
    dplyr::group_by(...) %>%
    dplyr::tally() %>%
    dplyr::mutate(perc_total = round(n/total_n*100),
                  output = paste0(n, " (", perc_total, ")")) %>%
    dplyr::arrange(-n) 
  
  if(caps == T){
    mydf <- mydf %>% dplyr::mutate(across(c(!!!group_vars), ~str_to_sentence(.))) 
  }
  
  return(mydf)
}

count_vals_ROB <- function(...){
  data_ROB %>%
    dplyr::group_by(...) %>%
    dplyr::tally() 
}

#- calculate meta-analytic network variables
calc_network_metrics <- function(mydf, ...){
  
  grouping_vars <- quos(...)
  
  count_dis_cooccurence <- function(Disease1, Disease2){
    
    Disease1 <- as.name(Disease1)
    Disease2 <- as.name(Disease2)
    
    #- calculate denominator: find number of clusters that included both diseases
    N <- mydf %>% 
      filter(!is.na({{Disease1}}) & !is.na({{Disease2}})) %>% summarise(total_N_dis1dis2 = n())
    
    #-- calculate n for D1 (for studies that included D1 and D2)
    nD1 <- mydf %>% 
      filter({{Disease1}} == 1 & !is.na({{Disease2}})) %>% summarise(n_dis1 = n())
    
    #-- calculate n for D2 (for studies that included D1 and D2)
    nD2 <- mydf %>% 
      filter(!is.na({{Disease1}}) & {{Disease2}} == 1) %>% summarise(n_dis2 = n())
    
    outputdf <- mydf %>% 
      filter({{Disease1}} == 1 & {{Disease2}} == 1) %>% 
      summarise(n_dis1dis2 = n()) %>%
      mutate(
        disease_pair = paste0(Disease1, ", ", Disease2),
        dis1 = paste0(Disease1),
        dis2 = paste0(Disease2)) %>%
      relocate(n_dis1dis2, .after = dis2)
    
    #- create joining variable if no grouping supplied
    if (length(grouping_vars) == 0){ 
      N <- mutate(N, group = "all")
      nD1 <- mutate(nD1, group = "all")
      nD2 <- mutate(nD2, group = "all")
      outputdf <- mutate(outputdf, group = "all")
    }
    #-- calculate n for clusters with Disease1 and Disease 2 and add other info to dataframe
    outputdf %>%
      left_join(N) %>%
      left_join(nD1) %>%
      left_join(nD2) %>%
      select(-any_of("group")) 
    
  }
  
  mydf <- mydf %>% group_by(...)
  
  #- count disease pairs across clusters
  meta_dis_pairs <- pmap(combo_list_pairs, count_dis_cooccurence) %>% list_rbind() 
  
  #-- calculate link Lift
  meta_dis_pairs <- meta_dis_pairs %>%
    mutate(weight = (n_dis1dis2/ (n_dis1 * n_dis2))* total_N_dis1dis2) 
  
  
  return(meta_dis_pairs)
  
  
}


#-- Plots -------

#- set legend sequence
calc_seq <- function(mydf, seq_var, seq_start, seq_length, round_to){
  
  max_val <- max(mydf[seq_var], na.rm = T)
  seq_val <- plyr::round_any(seq(from = seq_start, 
                                 to = plyr::round_any(max_val, round_to, ceiling), 
                                 length.out = seq_length), 
                             round_to, ceiling)
  return(seq_val)
} 


plot_network <- function(mynet, seq_links = seq_link, round_seq_by = 1, node_size = 4.5, 
                         y_nudge = 0.07, label_size = 2.5){
  
  ggraph(mynet, layout = "manual", x = layout$x, y = layout$y) +
    geom_edge_arc(aes(edge_width = weight, edge_alpha = weight, edge_color = weight), 
                  strength = 0.1) +
    scale_edge_color_viridis(
      option = "magma",
      direction = -1,
      end = 0.9,
      begin = 0.1,
      guide = "legend",
      breaks = seq_links, 
      limits = c(min(seq_links), plyr::round_any(max(seq_links), round_seq_by, ceiling))) +
    geom_node_point(size = {{node_size}},
                    shape = 21, fill = "#48b191", color = "00") +
    geom_node_label(aes(label = str_wrap(label, 20, whitespace_only = F), 
                        x = x + case_when(label == "Schizophrenia" ~ -0.05, 
                                          label == "Ischemic heart disease" ~ 0.025, 
                                          T ~ 0),
                        y = y + y_nudge,
                        vjust = 0),
                    alpha = 0.75, fill = "white", label.size = 0,
                    label.padding = unit(0, 'pt'),
                    label.r = unit(0.05, 'lines'),
                    size = label_size, fontface = "bold",
                    show.legend = F) +
    scale_edge_width(range = c(0.2, 0.8), 
                     breaks = seq_links, 
                     limits = c(min(seq_links), plyr::round_any(max(seq_links), round_seq_by, ceiling))) +
    scale_edge_alpha_continuous(range = c(0.5, 0.75), 
                                breaks = seq_links, guide = "none",
                                limits = c(min(seq_links), 
                                           plyr::round_any(max(seq_links), round_seq_by, ceiling))) +
    theme_graph(base_family = "sans") 
  
}

plot_network_with_clusters <- function(mynet, seq_links = seq_link, y_nudge = 0.05){
  
  ggraph(mynet, layout = "manual", x = layout_2$x, y = layout_2$y) +
    geom_edge_arc(
      aes(
        edge_color = stable_cluster_id,
        edge_width = ifelse(is.na(stable_cluster_id), 0, 1),
        edge_alpha = ifelse(is.na(stable_cluster_id), 0, 1)
      ), strength = 0.1) +
    geom_node_point(shape = 21, fill = "#999999", color = "00", size = 1.5) +
    scale_edge_width(range = c(0.1, 0.5), guide = "none") +
    scale_edge_alpha_continuous(range = c(0.4, 0.9), guide = "none") +
    geom_node_label(aes(label = str_wrap(label, 20, whitespace_only = F),
                        x = x + ifelse(label == "Schizophrenia", -0.05, 0),
                        y = y + y_nudge,
                        vjust = 0),
                    size = 2.25, fontface = "bold",
                    show.legend = F, label.size = 0, alpha = 0.75) +
    theme_graph(base_family = "sans") +
    geom_scatterpie(data= nodes_clust, aes(x = x1, y= y1, group = node), pie_scale = 1.25,
                    cols = paste0(1:6), linewidth = 0.05) +
    scale_edge_color_manual(values = cols_stable_clusters,
                            breaks = c("1", "2", "3", "4", "5", "6",  NA)) +
    scale_fill_manual(values = cols_stable_clusters) +
    guides(fill = "none", 
           edge_color = guide_legend(title = "Cluster\nNumber", override.aes = list(edge_width = 2, spacing = unit(0, "mm")))) 
  
  
}

plot_network_with_clusters_weightedlinks <- function(mynet, seq_links = seq_link, round_seq_by = 1, y_nudge = 0.05){
  
  
  ggraph(mynet, layout = "manual", x = layout_2$x, y = layout_2$y) +
    geom_edge_arc(
      aes(
        edge_color = stable_cluster_id,
        edge_width = weight, 
        edge_alpha = ifelse(is.na(stable_cluster_id), 0, 1)
      ), strength = 0.1) +
    geom_node_point(shape = 21, fill = "#999999", color = "00", size = 1.5) +
    scale_edge_width(range = c(0.2, 0.8), 
                     breaks = seq_links, 
                     limits = c(min(seq_links), plyr::round_any(max(seq_links), round_seq_by, ceiling))) +
    scale_edge_alpha_continuous(range = c(0.25, 1), guide = "none") +
    geom_node_label(aes(label = str_wrap(label, 20, whitespace_only = F),
                        x = x + ifelse(label == "Schizophrenia", -0.05, 0),
                        y = y + y_nudge,
                        vjust = 0),
                    size = 2.25, fontface = "bold",
                    show.legend = F, label.size = 0, alpha = 0.75) +
    theme_graph(base_family = "sans") +
    geom_scatterpie(data= nodes_clust, aes(x = x1, y= y1, group = node), pie_scale = 1.25,
                    cols = paste0(1:6), linewidth = 0.05) +
    scale_edge_color_manual(values = cols_stable_clusters,
                            breaks = c("1", "2", "3", "4", "5", "6",  NA)) +
    scale_fill_manual(values = cols_stable_clusters) +
    guides(fill = "none", 
           edge_color = guide_legend(title = "Cluster\nNumber", override.aes = list(edge_width = 2, spacing = unit(0, "mm")))) 
  
  
}

plot_network_cluster <- function(link_df, node_df, clust_num, link_strength = 0.1, auto_nudge = F, 
                                 seq_links = seq_link, round_seq_by = 1,
                                 label_size = 2.25, col_clusters = T, 
                                 node_size = 2.8,
                                 x_nudge = 0, y_nudge = 0.08,
                                 neg_y_nudge_nodes = c("Gout", "Diabetes", "Dementia", "Parkinson's disease", "Asthma"),
                                 x_nudge_nodes = NULL,
                                 jacc_stability = T){
  
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
                       x = x + ifelse(label %in% x_nudge_nodes, x_nudge, 0),
                       y = y + ifelse(label %in% neg_y_nudge_nodes, -y_nudge, y_nudge), 
                       vjust = ifelse(label %in% neg_y_nudge_nodes, 1, 0)),
                   size = label_size,
                   show.legend = F) +
    scale_edge_width(range = c(0.2, 0.8), 
                     breaks = seq_links, 
                     limits = c(min(seq_links), plyr::round_any(max(seq_links), round_seq_by, ceiling))) +
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
  
  
  # ggsave(clust_plot, filename = paste0(dir_plots_clusters,
  #                                      "cluster", clust_num, ".png"),
  #        dpi = 300, bg = "white",
  #        height = 3.5, width = 3.75)
  
  assign(paste0("clust_plot_", clust_num), clust_plot, envir = .GlobalEnv)
  
}

#- Lolipop plot function
lolipop_plot <- function(mydf, xvar, yvar){
  ggplot(mydf, 
         aes(x = {{xvar}}, y = {{yvar}})) + 
    geom_point(size = 1) +
    geom_segment(aes(x = 0, xend = {{xvar}}, y = {{yvar}}, yend = {{yvar}}),
                 linewidth = 0.4) +
    theme_minimal() +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(linetype = 'dotted', colour = "dark grey", linewidth = 0.3),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(size = 8, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 7.5),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(size = 8)) # plot.margin = margin(t = 10, b= 10, l = 10, r = 10)
}



bar_plot_metastudies <- function(feature = c("method", "sex", "age"), cluster){
  
  plot_title = case_when(feature == "sex" ~ "Sex stratification",
                         feature == "age" ~ "Age stratification",
                         feature == "method" ~ "Clustering method")
  
  palette_choice = match.arg(feature)
  
  col_palette = switch(palette_choice,
                       "sex" = cols_sex,
                       "age" = cols_ages,
                       "method" = cols_method)
  
  cluster_replicability_matches <- cluster_replicability_matches %>%
    filter(feature == {{feature}}, stable_cluster_id == {{cluster}}) %>%
    arrange(stable_cluster_id, category) %>%
    group_by(stable_cluster_id) %>%
    mutate(
      cum = cumsum(perc_total),                 # top of each segment
      mid = cum - perc_total/2,                  # midpoint for label
      label_off = if_else({{cluster}} %in% c(1,2) & feature == "age" & category %in% c("0-39", "40+"), 1, 0),
      label_color = if_else(label_off == 1 , "black", "white")
    ) 
  
  off_plot_labels <- cluster_replicability_matches %>% filter(label_off == 1) %>%
    mutate(hjust = if_else(category  == "0-39", 1, 0),
           label_y = if_else(category  == "0-39", 1.4, 0.65),
           dx = if_else(hjust == 1, -5, 4.),   # left/right depending on hjust
           dy = if_else(hjust == 1, -0.5, 0.4))
  
  
  plot <- cluster_replicability_matches %>% 
    ggplot(aes(fill = category, y = stable_cluster_id, x = perc_total)) +
    geom_bar(position = "stack", stat = "identity", color = "black", linewidth = 0.1) +
    geom_text(aes(x = if_else({{cluster}} == 3 & category == "Cluster-based", 100-mid-9, 100-mid), 
                  label = if_else(label_off == 0, paste0(category, "\n", round(perc_total), "%"), ""),
                  color = label_color, 
                  hjust = if_else({{cluster}} == 3 & category == "Cluster-based", 0, 0.5)), 
              fontface = "bold",
              size = 2) +
    scale_x_continuous(labels = scales::percent_format(scale = 1)) + #expand = expansion(mult = c(0.05, 0.05))) +
    theme_minimal(base_size = 14) +
    labs(title = plot_title) +
    theme(plot.margin = margin(0,0,0,0),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 5),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = 0, size = 6, face = "plain", margin = margin(t = 1, b = 0.5))
    ) +
    guides(fill = "none") +
    scale_fill_manual(values = col_palette) +
    scale_color_identity(guide = "none") 
  
  return(plot)
  
}

stack_bar_plots <- function(cluster_num){
  plot1 <- bar_plot_metastudies("method", {{cluster_num}}) + theme(axis.text.x = element_blank())
  plot2 <- bar_plot_metastudies("sex", {{cluster_num}}) + theme(axis.text.x = element_blank())
  plot3 <- bar_plot_metastudies("age", {{cluster_num}}) + theme(axis.text.x = element_blank())
  
  stacked_plot <- plot_grid(plot1, plot2, plot3, ncol = 1, rel_heights = c(1,1,1), align = "hv")  + 
    theme(plot.margin = margin(t = 8, b = 5, l = 0.5, r = 5, unit = "mm"))
  return(stacked_plot)
}

stack_clust_and_bar_plots <- function(cluster_num, plot_subtitle){
  
  n_clusters <- cluster_replicability_matches %>%
    filter(stable_cluster_id == {{cluster_num}}) %>%
    pull(n_matching_clusters) %>% unique()
  
  plot_subtitle <- paste0("matching cluster n: ", n_clusters)
  
  plot_name <- paste0("orig_clust_plot_", {{cluster_num}})
  plot_obj <- get(plot_name, envir = .GlobalEnv)
  
  plot_obj$layers[[3]]$aes_params$size <- 1.9 #- make labels smaller
  
  plot_grid(plot_obj + labs(subtitle = {{plot_subtitle}}) + theme(plot.title = element_text(size = 8), plot.subtitle = element_text(size = 5.5), text = element_text(size = 1)),
            stack_bar_plots({{cluster_num}}), nrow = 1, rel_widths = c(0.8, 1))
}


#--- meta-analyses ------
#- Custom functions for linkcomm dendrogram cut height testing

LinkCommunitiesCutStats <- function (network, ncuts, hcmethod = "average", use.all.edges = FALSE, 
                                     edglim = 10^4, directed = FALSE, dirweight = 0.5, bipartite = FALSE, 
                                     dist = NULL, plot = TRUE, check.duplicates = TRUE, removetrivial = TRUE, 
                                     verbose = TRUE) {
  
  #- Function modified from the LinkComm source code
  #- this function stops at the point of the hierarchical clustering of the links and outputs 
  #- summary statistics at different cutheights instead. 
  #- this is to allow the user to select the summary statistics to determine cut-height (rather than using the linkcomm default)
  #- ncuts = the number of dendrogram cut heights that will be tested (in descending order)
  
  
  if (is.character(network) && !is.matrix(network)) {
    if (file.access(network) == -1) {
      stop(cat("\nfile not found: \"", network, "\"\n", 
               sep = ""))
    }
    else {
      network <- read.table(file = network, header = FALSE)
    }
  }
  x <- network
  rm(network)
  if (ncol(x) == 3) {
    wt <- as.numeric(as.character(x[, 3]))
    if (length(which(is.na(wt) == TRUE)) > 0) {
      stop("\nedge weights must be numerical values\n")
    }
    x <- cbind(as.character(x[, 1]), as.character(x[, 2]))
  }
  else if (ncol(x) == 2) {
    x <- cbind(as.character(x[, 1]), as.character(x[, 2]))
    wt <- NULL
  }
  else {
    stop("\ninput data must be an edge list with 2 or 3 columns\n")
  }
  if (check.duplicates) {
    dret <- edge.duplicates(x, verbose = verbose)
    x <- dret$edges
    if (!is.null(wt)) {
      if (length(dret$inds) > 0) {
        wt <- wt[-dret$inds]
      }
    }
    rm(dret)
  }
  el <- x
  rm(x)
  len <- nrow(el)
  nnodes <- length(unique(c(as.character(el[, 1]), as.character(el[, 
                                                                   2]))))
  intel <- integer.edgelist(el)
  edges <- intel$edges
  node.names <- names(intel$nodes)
  numnodes <- length(node.names)
  if (bipartite) {
    big <- graph.edgelist(as.matrix(el), directed = directed)
    bip.test <- bipartite.mapping(big)
    if (!bip.test$res) {
      stop("\nnetwork is not bi-partite; change bipartite argument to FALSE\n")
    }
    bip <- rep(1, length(bip.test$type))
    bip[which(bip.test$type == FALSE)] <- 0
    names(bip) <- V(big)$name
    bip <- bip[match(node.names, names(bip))]
    rm(big, bip.test)
  }
  else {
    bip <- 0
  }
  rm(intel)
  if (len <= edglim) {
    disk <- FALSE
    if (is.null(dist)) {
      emptyvec <- rep(1, (len * (len - 1))/2)
      if (!is.null(wt)) {
        weighted <- TRUE
      }
      else {
        wt <- 0
        weighted <- FALSE
      }
      if (!use.all.edges) {
        dissvec <- .C("getEdgeSimilarities", as.integer(edges[, 
                                                              1]), as.integer(edges[, 2]), as.integer(len), 
                      rowlen = integer(1), weights = as.double(wt), 
                      as.logical(directed), as.double(dirweight), 
                      as.logical(weighted), as.logical(disk), dissvec = as.double(emptyvec), 
                      as.logical(bipartite), as.logical(verbose))$dissvec
      }
      else {
        dissvec <- .C("getEdgeSimilarities_all", 
                      as.integer(edges[, 1]), as.integer(edges[, 
                                                               2]), as.integer(len), as.integer(numnodes), 
                      rowlen = integer(1), weights = as.double(wt), 
                      as.logical(FALSE), as.double(dirweight), as.logical(weighted), 
                      as.logical(disk), dissvec = as.double(emptyvec), 
                      as.logical(bipartite), as.logical(verbose))$dissvec
      }
      distmatrix <- matrix(1, len, len)
      distmatrix[lower.tri(distmatrix)] <- dissvec
      colnames(distmatrix) <- 1:len
      rownames(distmatrix) <- 1:len
      distobj <- as.dist(distmatrix)
      
      rm(distmatrix)
    }
    else {
      if (!inherits(dist, "dist")) {
        stop("\ndistance matrix must be of class \"dist\" (see ?as.dist)\n")
      }
      else if (attr(dist, which = "Size") != len) {
        stop("\ndistance matrix size must equal the number of edges in the input network\n")
      }
      else if (length(dist) != (len * (len - 1))/2) {
        stop("\ndistance matrix must be the lower triangular matrix of a square matrix\n")
      }
      distobj <- dist
    }
    if (verbose) {
      cat("\n   Hierarchical clustering of edges...")
    }
    hcedges <- hclust(distobj, method = hcmethod)
    hcedges$order <- rev(hcedges$order)
    #rm(distobj)
    if (verbose) {
      cat("\n")
    }
  }
  else {
    disk <- TRUE
    if (!is.null(wt)) {
      weighted <- TRUE
    }
    else {
      wt <- 0
      weighted <- FALSE
    }
    if (!use.all.edges) {
      rowlen <- .C("getEdgeSimilarities", as.integer(edges[,
                                                           1]), as.integer(edges[, 2]), as.integer(len),
                   rowlen = integer(len - 1), weights = as.double(wt),
                   as.logical(directed), as.double(dirweight), as.logical(weighted),
                   as.logical(disk), dissvec = double(1), as.logical(bipartite),
                   as.logical(verbose))$rowlen
    }
    else {
      rowlen <- .C("getEdgeSimilarities_all", as.integer(edges[,
                                                               1]), as.integer(edges[, 2]), as.integer(len),
                   as.integer(numnodes), rowlen = integer(len -
                                                            1), weights = as.double(wt), as.logical(FALSE),
                   as.double(dirweight), as.logical(weighted), as.logical(disk),
                   dissvec = double(1), as.logical(bipartite), as.logical(verbose))$rowlen
    }
    if (verbose) {
      cat("\n")
    }
    hcobj <- .C("hclustLinkComm", as.integer(len),
                as.integer(rowlen), heights = single(len - 1), hca = integer(len -
                                                                               1), hcb = integer(len - 1), as.logical(verbose))
    if (verbose) {
      cat("\n")
    }
    hcedges <- list()
    hcedges$merge <- cbind(hcobj$hca, hcobj$hcb)
    hcedges$height <- hcobj$heights
    hcedges$order <- .C("hclustPlotOrder", as.integer(len),
                        as.integer(hcobj$hca), as.integer(hcobj$hcb), order = integer(len))$order
    hcedges$order <- rev(hcedges$order)
    hcedges$method <- "single"
    class(hcedges) <- "hclust"
    
  }
  hh <- unique(round(hcedges$height, digits = 3))
  
  heights <- hh[order(hh, decreasing = T)] 
  
  #- if you've specified more cuts than there are heights then test one minus the total number of heights
  if(length(heights) <= ncuts){ncuts = length(heights)-2}
  
  stats <- cstats.table.2(distobj, hcedges, ncuts, heights) 
  
  return(stats)
  
}


linkcomm_optimal_cutheight <- function(links_df, ncuts = 50){
  
  #- initialize dendrogram
  LC <- getLinkCommunities(as.data.frame(links_df), directed = F, 
                           hcmethod = "ward.D2", removetrivial = F)
  #- find optimal cutheight
  stats_LC <- 
    LinkCommunitiesCutStats(as.data.frame(links_df), ncuts, directed = F, removetrivial = F, hcmethod = "ward.D2") %>% 
    rownames_to_column(var = "measure") %>%
    filter(measure %in% c("cluster.number", "within.cluster.ss", "avg.silwidth", "ch", "dunn2", "sindex")) %>%
    pivot_longer(cols = starts_with("Test"), names_prefix = "Test.", names_to = "cutheight") %>%
    pivot_wider(names_from = "measure") %>%
    arrange(-ch)
  
  #- choose one height above the optimal-> gives best output
  newcutheight <- stats_LC[2,] %>% pull(cutheight)
  
  #- cut dendrogram at new height
  newLinkCommsAt(LC, newcutheight)
  
}


#- plots by # of clusters
elbow_plot <- function(statsdf, index = within.cluster.ss){
  
  ggplot(data = statsdf, ##data.frame(t(statsdf)), 
         aes(x=cluster.number, y= {{index}})) + 
    geom_jitter()+
    geom_line()+
    #labs(y = paste({{index}})) +
    theme(plot.title = element_text(hjust = 0.5))
  
}

#- plots by slice height
elbow_plot_2 <- function(statsdf, index = within.cluster.ss){
  
  ggplot(data = statsdf, 
         aes(x = {{index}}, y = cut.height)) + 
    geom_point()+
    #labs(y = paste({{index}})) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(angle = 90)) +
    scale_y_continuous(n.breaks = 20)
}


#-- custom functions for clsuter stability testing

#- functions to create summary stats tables for cluster stability testing, adapted from FPC package functions
cstats.table <- function(dist, tree, k) {
  clust.assess <- c("cluster.number","n","within.cluster.ss","average.within","average.between",
                    "wb.ratio", "ch", "dunn2","avg.silwidth", "sindex")
  clust.size <- c("cluster.size")
  stats.names <- c()
  row.clust <- c()
  output.stats <- matrix(ncol = k, nrow = length(clust.assess))
  cluster.sizes <- matrix(ncol = k, nrow = k)
  for(i in c(1:k)){
    row.clust[i] <- paste("Cluster-", i, " size")
  }
  for(i in c(2:k)){
    stats.names[i] <- paste("Test", i-1)
    
    for(j in seq_along(clust.assess)){
      output.stats[j, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.assess])[j]
      
    }
    
    for(d in 1:k) {
      cluster.sizes[d, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.size])[d]
      dim(cluster.sizes[d, i]) <- c(length(cluster.sizes[i]), 1)
      cluster.sizes[d, i]
      
    }
  }
  output.stats.df <- data.frame(output.stats)
  cluster.sizes <- data.frame(cluster.sizes)
  cluster.sizes[is.na(cluster.sizes)] <- 0
  rows.all <- c(clust.assess, row.clust)
  # rownames(output.stats.df) <- clust.assess
  output <- rbind(output.stats.df, cluster.sizes)[ ,-1]
  colnames(output) <- stats.names[2:k]
  rownames(output) <- rows.all
  is.num <- sapply(output, is.numeric)
  output[is.num] <- lapply(output[is.num], round, 2)
  output
}

cstats.table.2 <- function(dist, tree, k, height) {
  clust.assess <- c("cluster.number","n","within.cluster.ss","average.within","average.between",
                    "wb.ratio", "ch", "dunn2","avg.silwidth", "sindex")
  clust.size <- c("cluster.size")
  stats.names <- c()
  row.clust <- c()
  output.stats <- matrix(ncol = k, nrow = length(clust.assess))
  cluster.sizes <- matrix(ncol = k, nrow = k)
  for(i in c(1:k)){
    row.clust[i] <- paste("Cluster-", i, " size")
  }
  for(i in c(2:k)){
    stats.names[i] <- paste("Test", height[i])
    
    for(j in seq_along(clust.assess)){
      output.stats[j, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, h = height[i]))[clust.assess])[j]
      
    }
    
    for(d in 1:k) {
      cluster.sizes[d, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, h = height[i]))[clust.size])[d]
      dim(cluster.sizes[d, i]) <- c(length(cluster.sizes[i]), 1)
      cluster.sizes[d, i]
      
    }
  }
  output.stats.df <- data.frame(output.stats)
  cluster.sizes <- data.frame(cluster.sizes)
  cluster.sizes[is.na(cluster.sizes)] <- 0
  rows.all <- c(clust.assess, row.clust)
  # rownames(output.stats.df) <- clust.assess
  output <- rbind(output.stats.df, cluster.sizes)[ ,-1]
  colnames(output) <- stats.names[2:k]
  rownames(output) <- rows.all
  is.num <- sapply(output, is.numeric)
  output[is.num] <- lapply(output[is.num], round, 2)
  output
}

#- bootstrap stability testing function, adapted from fpc package
clusterboot_linkcomm <- function (data, B = 100, distances = (inherits(data, "dist")), 
                                  bootmethod = "boot", bscompare = TRUE, multipleboot = FALSE, 
                                  jittertuning = 0.05, noisetuning = c(0.05, 4), subtuning = floor(nrow(data)/2), 
                                  clustermethod, 
                                  noisemethod = FALSE, count = TRUE, showplots = FALSE, 
                                  dissolution = 0.5, recover = 0.75, seed = NULL, datatomatrix = TRUE, 
                                  ...) 
  
  #-- Note: in this edited function,  the jaccard testing does not take a subsampled portion of the clustered output, 
  #-- because for our use case with linkcomm testing, the original bootstrapped dataframe is aggregated into network links
  #-- so the output is always the same length for jaccard testing. 
  
{
  sumlogic <- function(x, y, relation = "eq") switch(relation, 
                                                     eq = sum(x == y, na.rm = TRUE), s = sum(x < y, na.rm = TRUE), 
                                                     l = sum(x > y, na.rm = TRUE), se = sum(x <= y, na.rm = TRUE), 
                                                     le = sum(x >= y, na.rm = TRUE))
  if (!is.null(seed)) 
    set.seed(seed)
  invisible(distances)
  if (datatomatrix) 
    data <- as.matrix(data)
  if (distances & showplots & datatomatrix) 
    dpoints <- cmdscale(data)
  n <- nrow(data)
  p <- ncol(data)
  if (datatomatrix & !distances) {
    cod <- cov(data)
    md <- colMeans(data)
  }
  lb <- length(bootmethod)
  if (distances) 
    c1 <- clustermethod(as.dist(data), input_first_net = T, ...)
  else c1 <- clustermethod(data, input_first_net = T, ...)
  if (noisemethod) {
    if (c1$nccl == 0) 
      stop("No clusters, only noise estimated!")
  }
  else c1$nccl <- c1$nc
  bootresult <- jitterresult <- noiseresult <- bojitresult <- subsetresult <- matrix(0, 
                                                                                     nrow = c1$nc, ncol = B)
  if (("jitter" %in% bootmethod) | ("bojit" %in% bootmethod)) {
    if (!datatomatrix | distances) 
      stop("datatomatrix=FALSE and distances require boot or subset as bootmethod.")
    jsd <- numeric(0)
    ecd <- eigen(cod, symmetric = TRUE)
    ecd$values[ecd$values < 0] <- 0
    ecd$values[is.na(ecd$values)] <- 0
    rotdata <- data %*% solve(t(ecd$vectors))
    for (i in 1:p) {
      sx <- sort(rotdata[, i])
      dx <- sx[2:n] - sx[1:(n - 1)]
      dx <- dx[dx > 0]
      jsd[i] <- quantile(dx, jittertuning)
    }
  }
  if ("noise" %in% bootmethod) {
    if (!datatomatrix | distances) 
      stop("datatomatrix=FALSE and distances require boot or subset as bootmethod.")
    ecd <- eigen(cod, symmetric = TRUE)
    ecd$values[ecd$values < 0] <- 0
  }
  if (showplots & datatomatrix) {
    if (distances) 
      plot(dpoints, pch = sapply(c1$partition, toString), 
           col = c1$partition)
    else plot(data, pch = sapply(c1$partition, toString), 
              col = c1$partition)
  }
  for (l in 1:lb) {
    for (i in 1:B) {
      if (count) 
        cat(bootmethod[l], i, "\n")
      if (bootmethod[l] == "boot") {
        bsamp <- sample(n, n, replace = TRUE)
        if (!multipleboot) 
          bsamp <- unique(bsamp)
        if (distances) 
          mdata <- data[bsamp, bsamp]
        else mdata <- data[bsamp, ]
      }
      if (bootmethod[l] == "subset") {
        bsamp <- sample(n, subtuning, replace = FALSE)
        if (distances) 
          mdata <- data[bsamp, bsamp]
        else mdata <- data[bsamp, ]
      }
      if (bootmethod[l] == "jitter") {
        jnoise <- matrix(0, ncol = p, nrow = n)
        for (j in 1:p) jnoise[, j] <- rnorm(n, sd = jsd[j])
        jnoise <- jnoise %*% t(ecd$vectors)
        mdata <- data + jnoise
        bsamp <- 1:n
      }
      if (bootmethod[l] == "bojit") {
        bsamp <- sample(n, n, replace = TRUE)
        jnoise <- matrix(0, ncol = p, nrow = n)
        for (j in 1:p) jnoise[, j] <- rnorm(n, sd = jsd[j])
        jnoise <- jnoise %*% t(ecd$vectors)
        mdata <- data[bsamp, ] + jnoise
      }
      if (bootmethod[l] == "noise") {
        noiseind <- as.logical(rbinom(n, 1, noisetuning[1]))
        nn <- sum(noiseind)
        jnoise <- matrix(0, ncol = p, nrow = nn)
        for (j in 1:p) jnoise[, j] <- runif(nn, min = -noisetuning[2] * 
                                              sqrt(ecd$values[j]), max = noisetuning[2] * 
                                              sqrt(ecd$values[j]))
        jnoise <- t(t(jnoise %*% t(ecd$vectors)) + md)
        mdata <- data
        mdata[noiseind, ] <- jnoise
        bsamp <- (1:n)[!noiseind]
      }
      if ("diss" %in% names(formals(clustermethod)) & distances) 
        bc1 <- clustermethod(mdata, diss = TRUE, ...)
      else bc1 <- clustermethod(mdata, ...)
      if (showplots & datatomatrix) {
        if (distances) 
          plot(dpoints[bsamp, ], pch = sapply(bc1$partition, 
                                              toString), col = bc1$partition)
        else plot(mdata, pch = sapply(bc1$partition, 
                                      toString), col = bc1$partition)
      }
      if (noisemethod) {
        effnc1 <- c1$nccl
        effnb1 <- bc1$nccl
      }
      else {
        effnc1 <- c1$nc
        effnb1 <- bc1$nc
      }
      for (j in 1:effnc1) {
        maxgamma <- 0
        if (effnb1 > 0) {
          for (k in 1:effnb1) {
            if (multipleboot) {
              if (bscompare) 
                ncases <- 1:n
              else {
                ncases <- 1
                m <- 2
                if (m <= n) {
                  if (!(bsamp[m] %in% bsamp[1:(m - 1)])) 
                    ncases <- c(ncases, m)
                  m <- m + 1
                }
              }
            }
            else ncases <- 1:length(bsamp)
            
            #--- modifications here: c1$clusterlist[[j]][bsamp][ncases] = c1$clusterlist[[j]]
            #----                    bc1$clusterlist[[k]][ncases] = bc1$clusterlist[[k]]
            cg <- switch(bootmethod[l], boot = clujaccard(c1$clusterlist[[j]], 
                                                          bc1$clusterlist[[k]], zerobyzero = 0), 
                         bojit = clujaccard(c1$clusterlist[[j]], 
                                            bc1$clusterlist[[k]], zerobyzero = 0), 
                         subset = clujaccard(c1$clusterlist[[j]], 
                                             bc1$clusterlist[[k]], zerobyzero = 0), 
                         jitter = clujaccard(c1$clusterlist[[j]], 
                                             bc1$clusterlist[[k]], zerobyzero = 0), 
                         noise = clujaccard(c1$clusterlist[[j]][!noiseind], 
                                            bc1$clusterlist[[k]][!noiseind], zerobyzero = 0))
            if (cg > maxgamma) 
              maxgamma <- cg
          }
        }
        if (bootmethod[l] == "boot") 
          bootresult[j, i] <- maxgamma
        if (bootmethod[l] == "subset") 
          subsetresult[j, i] <- maxgamma
        if (bootmethod[l] == "bojit") 
          bojitresult[j, i] <- maxgamma
        if (bootmethod[l] == "jitter") 
          jitterresult[j, i] <- maxgamma
        if (bootmethod[l] == "noise") 
          noiseresult[j, i] <- maxgamma
      }
      if (noisemethod) {
        if (c1$nc > c1$nccl) {
          j <- c1$nc
          if (bc1$nc > bc1$nccl) 
            
            #--- modifications here: c1$clusterlist[[j]][bsamp][ncases] = c1$clusterlist[[j]]
            #----                    bc1$clusterlist[[k]][ncases] = bc1$clusterlist[[k]]
            
            maxgamma <- switch(bootmethod[l], boot = clujaccard(c1$clusterlist[[c1$nc]], 
                                                                bc1$clusterlist[[bc1$nc]], zerobyzero = 0), 
                               bojit = clujaccard(c1$clusterlist[[c1$nc]][bsamp], 
                                                  bc1$clusterlist[[bc1$nc]], zerobyzero = 0), 
                               subset = clujaccard(c1$clusterlist[[c1$nc]], 
                                                   bc1$clusterlist[[bc1$nc]], zerobyzero = 0), 
                               jitter = clujaccard(c1$clusterlist[[c1$nc]], 
                                                   bc1$clusterlist[[bc1$nc]], zerobyzero = 0), 
                               noise = clujaccard(c1$clusterlist[[c1$nc]][!noiseind], 
                                                  bc1$clusterlist[[bc1$nc]][!noiseind], 
                                                  zerobyzero = 0))
          else maxgamma <- 0
          if (bootmethod[l] == "boot") 
            bootresult[j, i] <- maxgamma
          if (bootmethod[l] == "subset") 
            subsetresult[j, i] <- maxgamma
          if (bootmethod[l] == "bojit") 
            bojitresult[j, i] <- maxgamma
          if (bootmethod[l] == "jitter") 
            jitterresult[j, i] <- maxgamma
          if (bootmethod[l] == "noise") 
            noiseresult[j, i] <- maxgamma
        }
      }
    }
  }
  if (!("boot" %in% bootmethod)) 
    bootresult <- bootmean <- bootbrd <- bootrecover <- NULL
  else {
    bootmean = apply(bootresult, 1, mean, na.rm = TRUE)
    bootbrd = apply(bootresult, 1, sumlogic, y = dissolution, 
                    relation = "se")
    bootrecover = apply(bootresult, 1, sumlogic, y = recover, 
                        relation = "l")
  }
  if (!("jitter" %in% bootmethod)) 
    jitterresult <- jittermean <- jitterbrd <- jitterrecover <- NULL
  else {
    jittermean = apply(jitterresult, 1, mean, na.rm = TRUE)
    jitterbrd = apply(jitterresult, 1, sumlogic, y = dissolution, 
                      relation = "se")
    jitterrecover = apply(jitterresult, 1, sumlogic, y = recover, 
                          relation = "l")
  }
  if (!("subset" %in% bootmethod)) 
    subsetresult <- subsetmean <- subsetbrd <- subsetrecover <- NULL
  else {
    subsetmean = apply(subsetresult, 1, mean, na.rm = TRUE)
    subsetbrd = apply(subsetresult, 1, sumlogic, y = dissolution, 
                      relation = "se")
    subsetrecover = apply(subsetresult, 1, sumlogic, y = recover, 
                          relation = "l")
  }
  if (!("noise" %in% bootmethod)) 
    noiseresult <- noisemean <- noisebrd <- noiserecover <- NULL
  else {
    noisemean = apply(noiseresult, 1, mean, na.rm = TRUE)
    noisebrd = apply(noiseresult, 1, sumlogic, y = dissolution, 
                     relation = "se")
    noiserecover = apply(noiseresult, 1, sumlogic, y = recover, 
                         relation = "l")
  }
  if (!("bojit" %in% bootmethod)) 
    bojitresult <- bojitmean <- bojitbrd <- bojitrecover <- NULL
  else {
    bojitmean = apply(bojitresult, 1, mean, na.rm = TRUE)
    bojitbrd = apply(bojitresult, 1, sumlogic, y = dissolution, 
                     relation = "se")
    bojitrecover = apply(bojitresult, 1, sumlogic, y = recover, 
                         relation = "l")
  }
  if (showplots & datatomatrix) {
    if (distances) 
      plot(dpoints, pch = sapply(c1$partition, toString), 
           col = c1$partition)
    else plot(data, pch = sapply(c1$partition, toString), 
              col = c1$partition)
  }
  out <- list(result = c1, partition = c1$partition, nc = c1$nc, 
              nccl = c1$nccl, clustermethod = c1$clustermethod, B = B, 
              noisemethod = noisemethod, bootmethod = bootmethod, multipleboot = multipleboot, 
              dissolution = dissolution, recover = recover, bootresult = bootresult, 
              bootmean = bootmean, bootbrd = bootbrd, bootrecover = bootrecover, 
              jitterresult = jitterresult, jittermean = jittermean, 
              jitterbrd = jitterbrd, jitterrecover = jitterrecover, 
              subsetresult = subsetresult, subsetmean = subsetmean, 
              subsetbrd = subsetbrd, subsetrecover = subsetrecover, 
              bojitresult = bojitresult, bojitmean = bojitmean, bojitbrd = bojitbrd, 
              bojitrecover = bojitrecover, noiseresult = noiseresult, 
              noisemean = noisemean, noisebrd = noisebrd, noiserecover = noiserecover)
  class(out) <- "clboot"
  out
}

#- clustering function for bootstrapped network
cluster_network <- function(mydf, input_first_net = F, net_input_name = LC_net){
  #- input_first_net: manually input the first network (if dendrogram cut-height is already determined, and therefore cut-height testing doesn't need to be performed)
  #-- because the cutheight for this was chosen manually based on review of multiple fit statistics
  #-- but for each comparative network, the cutheight that maximized CH is used instead
  
  LC_net = {{net_input_name}}
  
  if(input_first_net == F){
    links_net <- calc_network_metrics(mydf) %>%
      filter(weight >= 1.0) %>%
      select(dis1, dis2, weight) %>%
      rename(from = dis1, to = dis2)
    
    LC_net <- linkcomm_optimal_cutheight(links_net, ncuts = 60)
  }
  
  c1 = LC_net
  nc = length(LC_net$clusters)
  
  #- compare links in df with links_list, and assign 0 cluster to any missing links
  links_list_1 = arrange(c1$edges, node1, node2, .locale = "en")
  
  links_list = links_key %>% left_join(links_list_1) %>%
    mutate(cluster = replace_na(cluster, "0"))
  
  partition = links_list$cluster
  
  cl = list()
  for (i in 1:nc) cl[[i]] <- partition == i
  
  out <- list(result = c1, nc = nc, clusterlist = cl, partition = partition,
              clustermethod = "linkcomm")
  
  return(out)
}


#-- find clusters with anchor disease PLUS at least 2 other diseases OR 2/3 diseases for non-anchor cluster
find_cluster_match_jacc <- function(cluster_name){
  #-- for each cluster in solution...
  #- extract cluster diseases
  clust_diseases <- clusters_stable_wide %>% filter(stable_cluster_id == {{cluster_name}}) %>% select(-stable_cluster_id) %>% select_if(~.x == 1) %>% names
  
  anchor_disease = case_when(cluster_name == 2 ~ "IHD",
                             cluster_name == 3 ~ "SCZ",
                             cluster_name == 4 ~ "MAD",
                             cluster_name == 5 ~ "HF",
                             cluster_name == 6 ~ "dementia"
  )
  
  if (!is.na(anchor_disease)){
    
    #--- find candidate clusters to test, 
    meta_dis_candidate_clusters <- meta_dis %>% 
      filter(if_all(clust_diseases, ~!is.na(.x))) %>% #- for studies that actually included all cluster diseases
      filter(if_all(anchor_disease, ~.x == 1)) %>% #-- keep only clusters that contained anchor disease
      arrange(clustid)
    
    #- remove anchor disease from other cluster diseases to determine test diseases
    clust_diseases_match <- setdiff(clust_diseases, anchor_disease)
    
    #-- for each cluster in candidate clusters, calculate Jaccard index against SR cluster 
    cluster_results_jaccard_clusterwise <- data.frame(stable_cluster_id = numeric(0), clustid = character(0), jacc = numeric(0), cluster_match = character(0))
    
    jaccard_candidate_clusters <- function(test_cluster_num){
      #- select corresponding row in candidate disease dataframe, and pull clusterid
      test_dis <- meta_dis_candidate_clusters[test_cluster_num,] 
      test_clustid = pull(test_dis, clustid)
      
      #- extract cluster diseases
      test_dis <- test_dis %>%
        #select(all_of(clust_diseases)) %>% #- only diseases  in the target cluster
        select(all_of(diseases)) %>% #- all diseases
        select_if(~.x == 1) %>% names
      
      #-- calculate jaccard index for each cluster
      add_row(cluster_results_jaccard_clusterwise, 
              stable_cluster_id = cluster_name, 
              clustid = test_clustid, 
              jacc = jaccard(clust_diseases, test_dis),
              cluster_match = ifelse(length(intersect(clust_diseases_match, test_dis)) >= 2, "y", "n"))
    }
    
    cluster_results_jaccard_clusterwise <- map(1:nrow(meta_dis_candidate_clusters), jaccard_candidate_clusters) %>% list_rbind() 
    
    #-- add back to SR cluster dataframe
    left_join(meta_dis_candidate_clusters, cluster_results_jaccard_clusterwise)
    
  } else {
    
    #--- find candidate clusters to test, 
    meta_dis_candidate_clusters <- meta_dis %>% 
      filter(if_all(clust_diseases, ~!is.na(.x))) %>% #-- not missing any cluster diseases
      filter(if_any(clust_diseases, ~.x == 1)) #-- need to have at least one disease from cluster (b/c there is no anchor)
    
    #-- for each cluster in candidate clusters, calculate Jaccard index against SR cluster 
    cluster_results_jaccard_clusterwise <- data.frame(stable_cluster_id = numeric(0), clustid = character(0), jacc = numeric(0), cluster_match = character(0))
    
    jaccard_candidate_clusters <- function(test_cluster_num){
      #- select corresponding row in candidate disease dataframe, and pull clusterid
      test_dis <- meta_dis_candidate_clusters[test_cluster_num,] 
      test_clustid = pull(test_dis, clustid)
      
      #- extract cluster diseases
      test_dis <- test_dis %>%
        #select(all_of(clust_diseases)) %>% #- only diseases  in the target cluster
        select(all_of(diseases)) %>% #- all diseases
        select_if(~.x == 1) %>% names
      
      #-- calculate jaccard index for each cluster
      add_row(cluster_results_jaccard_clusterwise, 
              stable_cluster_id = cluster_name, 
              clustid = test_clustid, 
              jacc = jaccard(clust_diseases, test_dis),
              cluster_match = ifelse(length(intersect(clust_diseases, test_dis)) >= 2, "y", "n"))
    }
    
    cluster_results_jaccard_clusterwise <- map(1:nrow(meta_dis_candidate_clusters), jaccard_candidate_clusters) %>% list_rbind() 
    
    #-- add back to SR cluster dataframe
    left_join(meta_dis_candidate_clusters, cluster_results_jaccard_clusterwise)
  }
  
}


#- Use jaccard index to find the most similar clusters in sensitivity analysis. 
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
