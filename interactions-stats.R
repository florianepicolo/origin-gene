#!/usr/bin/Rscript --vanilla


#### LIBRARY ####
packages_to_install <- c("XML", "plyr", "igraph", "dplyr", "tidyr", "magrittr", "ggplot2", "visNetwork", 
                         "paletteer", "scales", "ggforce", "ggthemes", "cluster", "pheatmap", "ggraph")
for (package in packages_to_install) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}



#### FUNCTIONS ####
recup_networks <- function(fxml){ 
  # initialisation
  ids <- factor(xpathSApply(fxml, "//entry", xmlGetAttr, "id")) # on récupère les attribues "id" des balises "entry"
  pathtitle <- strsplit(as.character(xpathApply(fxml, "//pathway", xmlGetAttr, "title")), " signaling")[[1]][1]
  names <- xpathSApply(fxml, "//entry", xmlGetAttr, "name")
  types <- factor(xpathSApply(fxml, "//entry", xmlGetAttr, "type"))
  items <- data.frame(ids, names, types) #on créé un dataframe avec ces infos récupérer
  
  # identifications des groupes
  group_ids <- ids[types=="group"] #récupère les ids des "groupes"
  i = 1 ; groups = list(NA,NA)
  for(idgroup in group_ids) {
    nodes <- getNodeSet(fxml, paste0("//entry[@id='",idgroup,"']//component"))
    for(indice in 1:length(nodes)){
      idgene <- xmlGetAttr(nodes[[indice]], "id")
      groups[[1]][[i]] = idgroup; groups[[2]][[i]] = idgene
      i = i + 1
    }
  }
  groups <- data.frame(idgroup = groups[[1]], idgene = groups[[2]])
  
  # items id and names (hsa) #récupérer tous les noms de chaque gène
  items_names <- list() 
  for(i in 1:length(names)){items_names[[i]] <- unlist(strsplit(names[i], "\\s+"))[1]} # un seul crochet si on récupère tout # récupérer les différents nom des gènes
  
  # les ids voies et leur noms (gestions des étiquettes avec des noms multiples)
  names(items_names) <- ids  # mettre l'id en face du nom
  allnames <- unlist(items_names) #on récupère tous les noms
  vectnames <- NA #on déclare vectnames
  for(i in 1:length(items_names)){vectnames <- c(vectnames, names(items_names)[i])} # on concatène à chaque fois pr chaque id, la répétition du nombre de nom qu'à chaque id ce nombre de fois
  vectnames <- vectnames[-1] # on retire le premier NA
  namesfromnetwork <- data.frame(vectnames, allnames)  # tous les noms récupérer dans la voie, on en fait un df avec les noms et les id correspondant
  
  # relations constitution
  entry1 <- factor(xpathSApply(fxml, "//relation", xmlGetAttr, "entry1")) 
  entry2 <- factor(xpathSApply(fxml, "//relation", xmlGetAttr, "entry2"))
  relations <- data.frame(entry1, entry2) %>% semi_join(., items, by=c("entry1"="ids"))
  relations <- semi_join(relations, items, by=c("entry2"="ids"))
  # on ajoute les relations des groups
  if(nrow(groups) > 1 & !is.na(groups[1,1])){
    interactions <- split(groups$idgene, groups$idgroup) %>% lapply(function(x) t(combn(x, 2, simplify = FALSE))) %>% unlist(recursive = FALSE)
    interactionstable <- do.call(rbind, interactions)
    colnames(interactionstable) <- c("entry1", "entry2")
    relations <- rbind(relations, interactionstable)
  }
  res <- list(items=items, relations=relations, namesfromnetwork=namesfromnetwork, groups=groups) #on fait une liste avec item, nom, relation  ### 
  return(list(res, pathtitle))
}
set_ages <- function(GRN){
  GRN$items$times <- NA
  number_of_match <- id_match <- list()
  
  for(j in 1:nrow(GRN$items)) {
    # on s'occupe des genes
    if(GRN$items$types[j]=="gene"){ #on a les dates que pr les gènes donc on ne prend qu'eux
      id_match[[j]] <- which(data_KEGG$kegg_id %in% unlist(strsplit(GRN$items$names[j], " ")[[1]])) #quand on a qu'un vectnames
      number_of_match[[j]] <- length(id_match[[j]][!is.na(id_match[[j]])])
      if(number_of_match[[j]] == 1) {GRN$items$times[j] = data_KEGG$num_clade[id_match[[j]]]}
      if(number_of_match[[j]] > 1) {GRN$items$times[j] = min(as.numeric(names(table(data_KEGG$num_clade[id_match[[j]][!is.na(id_match[[j]])]]))))}
    }
    # et maintenant des groups ! 
    if(GRN$items$types[j]=="group"){ # groupe formé de gène-gène !
      result <- GRN$groups %>% filter(idgroup == GRN$items$ids[j]) %>% left_join(GRN$items, by = c("idgene" = "ids")) %>%
        left_join(GRN$namesfromnetwork, by = c("idgene" = "vectnames")) %>% left_join(data_KEGG, by = c("allnames" = "kegg_id"))
      GRN$items$times[j] = max(result$num_clade)
    }
  }
  
  GRN$items <- GRN$items %>% dplyr::mutate(times = case_when(
    times %in% c(2, 3) ~ 2,
    times %in% c(4, 5, 6) ~ 3,
    times %in% c(7, 8) ~ 4,
    times %in% c(9) ~ 5,
    times %in% c(10) ~ 6,
    times %in% c(11) ~ 7,
    times %in% c(12, 13) ~ 8,
    times %in% c(14) ~ 9,
    times %in% c(15) ~ 10,
    times %in% c(16) ~ 11,
    times %in% c(17) ~ 12,
    times %in% c(18, 19, 20, 21, 22) ~ 13,
    times %in% c(23) ~ 14,
    times %in% c(24, 25) ~ 15,
    TRUE ~ times))
  
  return(GRN)
}
extract_first_element <- function(x) {unlist(strsplit(x, " "))[1]}
set_names <- function(GRN){
  GRN$items$first_element <- sapply(GRN$items$names, extract_first_element)
  corresponding_genes <- match(GRN$items$first_element, data_KEGG$kegg_id)
  GRN$items$gene_name <- ifelse(is.na(corresponding_genes), GRN$items$first_element, gsub(" ", "", data_KEGG$gene_name[corresponding_genes]))
  GRN$items <- subset(GRN$items, select = -c(first_element))
  # head(GRN$items)
  grouped_gene_names <- GRN$groups %>%
    left_join(GRN$items, by = c("idgene" = "ids")) %>%
    group_by(idgroup) %>%
    dplyr::summarise(grouped_gene_name = paste(unique(gene_name), collapse = ";"))
  GRN$items$gene_name <- ifelse(GRN$items$types == "group", grouped_gene_names$grouped_gene_name[match(GRN$items$ids, grouped_gene_names$idgroup)], GRN$items$gene_name)
  return(GRN)
}
visual_graph <- function(GRN){
  graph <- graph_from_data_frame(GRN[[2]], directed=TRUE, vertices=GRN[[1]])
  times_palette = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 25)
  times_colors = times_palette[GRN$items$times]
  
  plot(graph,
       layout = layout_with_fr(graph),
       vertex.label.cex = 0.8,  # Taille du texte
       vertex.size = 10,        # Taille des nœuds
       vertex.label.font = 2,   # Épaisseur de la police des étiquettes
       edge.arrow.size = 0.5,    # Taille des flèches des arêtes pour les graphes dirigés
       vertex.color = times_colors   # Utiliser les couleurs définies pour les nœuds
  )
}
is_included <- function(path, other_paths) {
  for (other_path in other_paths) {
    if (length(path) < length(other_path) && all(path == other_path[1:length(path)])) {
      return(TRUE)
    }
  }
  return(FALSE)
}
get_subgraph <- function(GRN){
  graph <- graph_from_data_frame(GRN[[2]], directed=TRUE, vertices=GRN[[1]])
  subgraphs <- Filter(function(subgraph) { vcount(subgraph) > 1 }, decompose(graph)) # on sort l'ensemble des sous-graph de plus d'un noeud du graph
  
  interactions <- get.edgelist(graph)
  df_interactions <- data.frame(from = interactions[,1],                 # premier élément
                                to = interactions[,2] ) %>%              # deuxième élement
    mutate(positions = NA, delta_age = NA, direction = NA) 
  
  for (subgraph in subgraphs){
    up_node <- which(degree(subgraph, mode = "in") == 0 & degree(subgraph, mode = "out") > 0) # on récupère l'ensemble des noeuds de début de voie 
    for (node in up_node){
      # node = 19
      all_paths <- all_simple_paths(subgraph, from = node, to = V(subgraph), mode = "out") # on récupère toutes les sous-voies commençant par un noeud spécifique (début de voie)
      
      filtered_subpaths <- all_paths[!sapply(all_paths, is_included, other_paths = all_paths)] # on retire les voies qui sont déjà dans d'autres sous-voies (donc on prend que les plus grandes)
      paths_with_ids <- lapply(filtered_subpaths, function(path) as_ids(path))
      
      # on va donner une position à chacun des éléments de chaque sous-graph
      for (path in paths_with_ids){
        position_in_path = 1
        for(i in 1:(length(path) - 1)){   # i = les différentes positions
          from = path[i]
          to = path[i+1]
          relation = which(df_interactions$from == from & df_interactions$to == to)
          
          if (!is.na(df_interactions$positions[relation]) && df_interactions$positions[relation]!=position_in_path){
            new_row <- data.frame(
              from = from,
              to = to,
              positions = position_in_path,  # nouvelle position
              delta_age = GRN$items$times[GRN$items$ids == to] - GRN$items$times[GRN$items$ids == from],
              direction = ifelse((GRN$items$times[GRN$items$ids == to] - GRN$items$times[GRN$items$ids == from]) > 0, "forward", ifelse((GRN$items$times[GRN$items$ids == to] - GRN$items$times[GRN$items$ids == from]) < 0, "backward", "simultaneous"))
            )
            df_interactions <- rbind(df_interactions, new_row)

          }else{
            # on fait delta de l'age de B - A
            df_interactions$delta_age[relation] = GRN$items$times[GRN$items$ids==to] - GRN$items$times[GRN$items$ids==from]
            df_interactions$direction[relation] = ifelse(df_interactions$delta_age[relation]>0, "forward", ifelse(df_interactions$delta_age[relation]<0, "backward", "simultaneous"))
            
            df_interactions$positions[relation] = position_in_path
          }
          position_in_path = position_in_path + 1
        }
      }
    }
  }
  df_interactions <- distinct(df_interactions)
  df_interactions <- df_interactions %>% group_by(from) %>% dplyr::mutate(weights = n_distinct(positions)) 
  
  return(df_interactions)
}
wtd.cor <-function(x,y,weights){
  wtd.mean <- function (x, weights = NULL, normwt = "ignored", na.rm = TRUE) {
    if (!length(weights)) 
      return(mean(x, na.rm = na.rm))
    if (na.rm) {
      s <- !is.na(x + weights)
      x <- x[s]
      weights <- weights[s]
    }
    sum(weights * x)/sum(weights)
  }
  wtd.var <- function (x, weights = NULL, normwt = FALSE, na.rm = TRUE, method = c("unbiased","ML")) {
    method <- match.arg(method)
    if (!length(weights)) {
      if (na.rm) 
        x <- x[!is.na(x)]
      return(var(x))
    }
    if (na.rm) {
      s <- !is.na(x + weights)
      x <- x[s]
      weights <- weights[s]
    }
    if (normwt) 
      weights <- weights * length(x)/sum(weights)
    if (normwt || method == "ML") 
      return(as.numeric(stats::cov.wt(cbind(x), weights, method = method)$cov))
    sw <- sum(weights)
    if (sw <= 1) 
      warning("only one effective observation; variance estimate undefined")
    xbar <- sum(weights * x)/sw
    sum(weights * ((x - xbar)^2))/(sw - 1)
  }
  stdz <- function (x, weight = NULL) {
    if (is.null(weight)) {
      weight <- rep(1, length(x))
    }
    x <- x - wtd.mean(x, weight, na.rm = TRUE)
    x <- x/sqrt(wtd.var(x, weight, na.rm = TRUE))
    x
  }
  
  coef(summary(lm(stdz(y, weight = weights) ~ stdz(x, weight = weights), weights = weights)))[2, ] 
}
permutation <- function(GRN, var = "times"){
  tempGRN <- GRN
  if (var == "relations"){
      # dans les relations je dois garder les bons âges
      tempGRN$relations$entry2 <- sample(tempGRN$relations$entry2)
      tempGRN$relations$entry1 <- sample(tempGRN$relations$entry1)
    }else{
      # dans les âges, je dois permuter uniquement sur les genes et group pour retirer path/cpd
      typevalide <- which(tempGRN$items$types=="gene" | tempGRN$items$types=="group")
      tempGRN$items$times[typevalide] <- sample(tempGRN$items$times[typevalide])
  }
  return(tempGRN)
}
get_degrees <- function(GRN){
  graph <- graph_from_data_frame(GRN[[2]], directed=TRUE, vertices=GRN[[1]])
  degrees <- data.frame(
    ids = V(graph)$name,                     # id des relations
    degree_in = degree(graph, mode = "in"),  # nombre de relation où l'id est activé
    degree_out = degree(graph, mode = "out") # nombre de relation où l'id active -> nécessaire pour le "weight" des correlations
  )
  return(degrees)
}
pvalue_comp <- function(x, distrib, tail=c("one", "two")){ 
  # x = la moyenne de nos tailles de cluster de données brutes
  pvalue <- NA
  med <- median(distrib, na.rm = TRUE) # médiane de la distribution des 1000 permutations
  if(med>x) { # si med > x, le tail one va être différents ! 
    if(tail=="one"){
      pvalue <- (length(which(distrib<x))+1)/length(distrib) # pvalue = nombre de fois ou la distribution est inférieur à x et on divise ça par la longuuer des simulations qu'on a fait
      results <- "-" # résultats == - car la variable brut est plus petite que la médiane
    } else {
      pvalue <- (length(which(distrib<x))+length(which(distrib>(x+2*abs(x-med))))+1)/length(distrib) # x + 2 fois la valeur de la médiane !  (valeur plus extreme que la valeur distri)
      results <- "-"
    }
  }
  if(med<=x) { # on peut changer en retirant le = ! 
    if(tail=="one"){
      pvalue <- (length(which(distrib>x))+1)/length(distrib)
      results <- "+" # la valeur brut est plus grande que la médiane
    } else {
      pvalue <- (length(which(distrib>x))+length(which(distrib<(x-2*abs(x-med))))+1)/length(distrib)
      results <- "+"
    }
  }
  diff <- x-mean(distrib)
  return(list(results, pvalue, diff)) #on retourne une liste avec les résutlats +/- et les pvalue
}
subgraphs_given_times_condition <- function (GRN, condition) {
  graph <- graph_from_data_frame(GRN[[2]], directed=TRUE, vertices=GRN[[1]]) # je transforme les relations en graph, il est dirigé, et il a besoin des items
  V(graph)$label <- GRN$items$gene_name
  subgraphs <- Filter(function(subgraph) { vcount(subgraph) > 1 }, decompose(graph))
  
  all_components <- NA # on va garder NA == TRUE 
  # NA = TRUE va augmenter surement nos tailles de cluster, mais comme on fait la même chose pour les permutées, it's ok
  
  for(subgraph in subgraphs){
    tempg <- subgraph
    v_edges <- get.edges(tempg, E(tempg)) #on récup-re toutes les branches du graph qu'on traite pr les avoir sous forme d'objet dataframe
    times_edges <- mapvalues(v_edges, V(tempg), vertex_attr(tempg, "times")) #on remplace dans un vecteur, on remplace les attribut par les temps
    edges_cond <- apply(times_edges, 1, condition) #on applique la fonction sur times_edges # ça dit si la relation entre les 2 noeuds à chauqe ligne respecte la condition donné
    edges_cond[is.na(edges_cond)] <- FALSE # ça dit ce qu'on fait des NA ! 
    # les entrées du vecteurs qui sont NA, on les met TRUE donc on continue la relation

    if (sum(edges_cond)>0){
      suba <- subgraph.edges(tempg , E(tempg)[edges_cond]) # on extrait les sous graph qu'on traite dont les relations en filtrant les relations de manière à ce qu'elle valide la condition. 
      all_components <- c(all_components, components(suba)$csize) # et la je concatène pr chaque graph avec la taille de la composante de subnaf, on a une distri des longuuers de subgraph
    }else{all_components <- c(all_components, 0)}
  }
  all_components <- all_components[-1] # on vire les NA au fur et à mesure comme on les a ajouté au début
  mean_distri = mean(all_components)
  return(mean_distri) # on sort une lsite de taille de subgrapph qui respecte la condition donné
}
calc_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
} # calcul de la valeur la plus haute parmis les valeurs les plus fréquentes

#### PLOTS #### 
plot_posids <- function(df_genes, name_path){
  p <- df_genes %>% drop_na() %>% filter(path == name_path) %>% group_by(ids, name, pos, date) %>% dplyr::summarise(n=n(), min_rank=min(pos), max_rank=max(pos)) %>% 
    arrange(min_rank, max_rank, -n, ids) %>%  mutate(ids=factor(name, levels=name %>% unique())) %>% 
    ggplot(aes(x= reorder(ids, max_rank), y=pos, size=n, col=factor(date, levels=1:15))) + geom_point() +
    theme(axis.text.x = element_text(angle=90)) + ggtitle(name_path) + labs(x = "gene name", y = "pathway posititon", col= "branch of origin", size = "n") +
    scale_color_manual(values = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 15), drop = FALSE)
  # print(p)
  plot(p, main = paste("figure", name_path))
} # figure ggplot-generankbirthpath
plot_relations <- function(df_relations){
  df_relations$dir <- factor(df_relations$direction, levels = c("forward", "simultaneous", "backward"))
  p_direction <- df_relations %>% ungroup %>% drop_na() %>% 
    select(dir, path) %>% group_by(path, dir) %>% dplyr::summarise(n = dplyr::n())  %>% group_by(path) %>%
    dplyr::mutate(percentage = n / sum(n) * 100) %>% filter(dir == "backward") %>% arrange(-percentage)  %>% pull(path)
  
  p <- df_relations %>%
    drop_na() %>%
    group_by(path, dir) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    group_by(path) %>%
    dplyr::mutate(percentage = count / sum(count) * 100) %>%
    ggplot(aes(x = percentage, y = factor(path, levels = p_direction), fill = dir)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = paste0(round(percentage), "%")), position = position_stack(vjust = 0.5), 
              size = 3, color = "black") +
    scale_fill_manual(values = c("forward" = "#5cb032", "simultaneous" = "#f8ac32", "backward" = "#e73631"),
                      breaks = c("forward", "simultaneous", "backward")) +
    theme_minimal() +
    ggtitle("Distribution of appearance order of interaction by pathway (%)") +
    labs(x = "percentage of appearance order of interaction", y = "pathway", fill="appearance order")

  print(p)
} # figure 4
plot_cumuldistribirth <- function(df_genes){
  color_palette <- paletteer_c("ggthemes::Sunset-Sunrise Diverging", 15)[1:13]
  p <- df_genes %>% ungroup() %>% 
    filter(!grepl(";", name), !grepl("cpd", name)) %>%
    drop_na() %>% 
    select(name, date) %>% 
    unique() %>% 
    group_by(date) %>%
    dplyr::summarise(count_gene = dplyr::n()) %>%
    arrange(date) %>%
    dplyr::mutate(
      cumulative_count = cumsum(count_gene),
      cumulative_percentage = cumsum(count_gene) / sum(count_gene) * 100
    ) %>%
    ggplot(aes(x = factor(date))) +
    geom_bar(aes(y = cumulative_count, fill = date), stat = "identity") +
    geom_text(
      aes(y = cumulative_count, label = paste0(round(cumulative_percentage, 1), "%")),
      vjust = -0.5, size = 3, color = "black"
    ) +
    geom_col(aes(y = count_gene), fill = "lightgray", alpha = 0.6, width = 0.5) +
    geom_text(
      aes(y = min(count_gene), label = count_gene),
      vjust = -0.5, size = 3, color = "black"
    ) +
    labs(
      x = "age rank of appearance", 
      y = "number of gene / cumulative percentage of gene", 
      title = "Distribution of gene (cumulative percentage and raw data)"
    ) +
    scale_fill_gradientn(colours = color_palette, name = "Ajout", guide = "none") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
    scale_x_discrete(labels = date_clade$label, drop = TRUE)
  print(p)  
} # figure 3
plot_delta <- function(df_relations){
  p_delta <- df_relations %>% drop_na() %>% select(delta_age, pos_from, pos_to, path) %>% 
    group_by(path, pos_from, pos_to) %>% unique() %>% dplyr::count(delta_age)
  p_median_delta <- p_delta %>% group_by(path) %>% dplyr::summarise(median = median(delta_age)) %>% arrange(median) %>% pull(path)
  p_delta <- p_delta %>% dplyr::mutate(path = factor(path, levels=p_median_delta))
  p <- p_delta %>% ggplot(.,  aes(x = delta_age, y = path)) +
    geom_boxplot() + labs(x = "age rank difference (delta)", y = "pathway") +
    ggtitle("Distribution of age rank difference by pathway") +
    scale_x_continuous(breaks = seq(-25, 25, by = 5)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red")
  plot(p)
  tab_plot <- layer_data(p) %>% select(xmin, xlower, xmiddle, xupper, xmax)
  tab_plot$path <- rep(p_median_delta, length.out = nrow(tab_plot))
  return(tab_plot)
} # figure suppl1
plot_birth <- function(df_genes){
  p_birth <- df_genes %>% ungroup() %>% drop_na() %>% select(name, date, path) %>% 
    group_by(path) %>% unique() %>% dplyr::count(date)
  p_median_birth <- p_birth %>% group_by(path) %>% dplyr::summarise(median = median(date)) %>% arrange(median) %>% pull(path)
  p_birth <- p_birth %>% dplyr::mutate(path = factor(path, levels=p_median_birth))
  p <- p_birth %>% ggplot(.,  aes(x = date, y = path)) +
    geom_boxplot() + labs(x = "age rank of appearance", y = "pathway") +
    ggtitle("Distribution of age rank of appearance by pathway") +
    scale_x_continuous(breaks = seq(-15, 15, by = 2.5))
  plot(p)
  tab_plot <- layer_data(p) %>% select(xmin, xlower, xmiddle, xupper, xmax, xmin_final, xmax_final)
  tab_plot$path <- rep(p_median_birth, length.out = nrow(tab_plot))
  return(tab_plot)
}
plot_distrirankbirth <- function(df_genes){
  order_correlation <- df_correlation %>% filter(pvalue<0.05) %>% arrange(r_correlation, path) %>% rbind(., df_correlation %>% filter(pvalue>=0.05) %>% arrange(r_correlation, path))
  p <- df_genes %>% drop_na() %>% mutate(path = factor(path, levels = order_correlation$path)) %>% ggplot(., aes(y = "",  fill = factor(date, levels=1:15))) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 15), drop = FALSE) +
    coord_cartesian() +
    theme_void() +
    facet_grid(path ~ pos, switch = "both") +
    theme(strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0)) +
    labs(x = "path", y = "pathway position", fill = "age rank of appearance") +
    ggtitle("Distribution of age rank of appearance by pathway and pathway position")
  print(p)
} #figure ggplot-distributionbirthrank
plot_distribirthdelta <- function(df_relations, name_path){
  p <- df_relations %>% drop_na() %>% dplyr::group_by(path) %>% unique() %>% dplyr::mutate(label = paste(date_from,":", date_to, sep = "")) %>% 
    dplyr::group_by(delta_age, label, date_from, date_to) %>% dplyr::filter(path==name_path) %>% dplyr::summarise(freq=dplyr::n()) %>% 
    ggplot(., aes(delta_age, freq, fill = factor(date_from, levels=1:15))) +
    geom_bar(position = "stack", stat = "identity") + theme_minimal() +
    scale_fill_manual(values = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 15), drop = FALSE) + 
    labs(x = "delta", y = "distribution of age rank appearance", fill = "age rank of appearance") +
    ggtitle(paste("Distribution of age rank of appearance by delta for", name_path))
  plot(p, main = paste("figure", name_path))
  # print(p)
} #figure suppl2


#### FILES ####
file_KEGG <- "/home/fpicolo/Desktop/Pathways/birth-animals/nouveau-run/p-allinfos-KEGG.csv"
doss_paths <- "/home/fpicolo/Desktop/Pathways/birth-animals/paths/"

data_KEGG <- read.csv(file_KEGG, sep=";")
files_paths <- dir(doss_paths) # récupère tous les fichiers d'un dossier



#### INITIALISATION ####
load("/home/fpicolo/Desktop/Pathways/birth-animals/2024/session_16012025.RData")

npermut = 1000
# les conditions : backward ou simult, backward(A plus jeune que B), simult, forward (A plus vieux que B), forward ou simult
conditions <- c(function(x){x[1]>=x[2]}, function(x){x[1]>x[2]}, function(x){x[1]==x[2]}, function(x){x[1]<x[2]}, function(x){x[1]<=x[2]})
df_pvalue = data.frame(matrix(NA, nrow = length(files_paths), ncol = 1+length(conditions)*3))
df_fbs = data.frame(matrix(NA, nrow = length(files_paths), ncol = 1+length(conditions[2:4])*3))
df_relations_allpaths <- data.frame(matrix(NA, nrow = 0, ncol = 12))
df_relations_allpaths_permutees <- data.frame(matrix(NA, nrow = 0, ncol = 12))
df_genes_allpaths <- data.frame(matrix(NA, nrow = 0, ncol = 7))
df_correlation <- data.frame(matrix(NA, nrow = length(files_paths), ncol = 3))
df_distri_fbs <- data.frame(matrix(NA, nrow = 0, ncol = 3)) %>% setNames(c("backward", "simultaneous", "forward"))
palette <- tibble(color = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 15), clade=1:15)

date_clade <- data.frame(
  num = 1:15,
  age = c(1300, 765, 743, 708, 635, 558, 570, 563, 462, 429, 415, 350, 330, 166, 98),
  name = c("Opisthokonta", "Metazoa", "Eumetazoa", "Bilateria", "Deuterostomia", "Chordata",
           "Olfactora", "Vetebrata", "Gnathostomata", "Osteichthyes", "Sarcopterygii", "Tetrapoda",
           "Amniota", "Mammalia", "Theria"),
  label = paste(date_clade$name, " (", date_clade$age, "My)", sep = "")
)


#### PROGRAMM ####
for(i in 1:length(files_paths)){
  class_path <- xmlParse(paste(doss_paths, files_paths[[i]], sep=""))
  xml_path <- xmlRoot(class_path) # le fichier xml interprété par R
  list_result <- recup_networks(fxml = xml_path) # on sépare la voie en plusieurs informations dans des tableaux
  GRN <- list_result[[1]]; name_path <- list_result[[2]]
  distri_fbs <- data.frame(matrix(NA, nrow = 0, ncol = 3)) %>% setNames(c("forward", "simultaneous", "backward"))
  
  print(paste("###", i, name_path, "###", sep=" "))
  
  GRN <- set_ages(GRN = GRN) # on associe les âges aux gènes (si paralogues : le plus ancien, si groupe : le plus récent)
  GRN <- set_names(GRN = GRN) # on ajoute les noms des gènes pour plus tard la méga voie
  # visual_graph(GRN = GRN) ## plot du graph avec les noms

  df_relations <- get_subgraph(GRN) %>% mutate(path = name_path) # donne l'ensemble des relations d'un graph, leur ordre et leur direction
  df_relations_allinfos <- df_relations %>% dplyr::rename(pos_from = positions) %>% mutate(pos_to = pos_from +1) %>% 
    left_join(GRN$items, by = c("from" = "ids")) %>% dplyr::rename(date_from = times, name_from = gene_name) %>% 
    left_join(GRN$items, by = c("to" = "ids")) %>% dplyr::rename(date_to = times, name_to = gene_name) %>% 
    select(from, to, name_from, name_to, pos_from, pos_to, date_from, date_to, delta_age, direction, weights, path)
  df_relations_allpaths <- rbind(df_relations_allpaths, df_relations_allinfos)
  
  df_degrees <- get_degrees(GRN = GRN)
  
  ##### CORRELATION POSITION VS AGE #####
  # calcul de la correlation entre les delta age et les positions avec l'ajout du poids)
  correlation_path <- wtd.cor(x=df_relations$delta_age, y=df_relations$positions, weights=(1/df_relations$weights))
  df_correlation[i,] = c(name_path, correlation_path[1], correlation_path[4])
  
  ##### PERMUTATION POUR LES TAILLES DE CLUSTERS #####
  # moyenne des tailles des clusters validant condition pour données brutes
  meancluster = c()
  for (indx_cond in 1:length(conditions)){
    cond = conditions[[indx_cond]]
    meancluster <- c(meancluster, subgraphs_given_times_condition(GRN = GRN, condition = cond))
  }
  
  df_meancluster = data.frame(matrix(NA, nrow = npermut, ncol = length(conditions)))
  progressbar <- txtProgressBar(min = 0, max = npermut, style = 3)
  for (n in 1:npermut){
    tempGRN <- permutation(GRN = GRN, var = "times") # "relation"/"times" possible aussi
    # visual_graph(tempGRN) 
    for (indx_cond in 1:length(conditions)){
      cond = conditions[[indx_cond]]
      df_meancluster[n, indx_cond] = subgraphs_given_times_condition(GRN = tempGRN, condition = cond)
    }

    # on se prépare un petit jeu de data permutée
    df_relations_allinfos_permutees <- df_relations_allinfos %>% select(-date_from, -date_to, -delta_age, -direction) %>% 
      left_join(tempGRN$items %>% select(ids, times), by = c("from" = "ids")) %>% dplyr::rename(date_from = times) %>% 
      left_join(tempGRN$items %>% select(ids, times), by = c("to" = "ids")) %>% dplyr::rename(date_to = times) %>% 
      mutate(delta_age = date_to - date_from, direction = ifelse(delta_age>0, "forward", ifelse(delta_age<0, "backward", "simultaneous"))) %>%
      select(from, to, name_from, name_to, pos_from, pos_to, date_from, date_to, delta_age, direction, weights, path)
    fbs <- df_relations_allinfos_permutees %>% ungroup() %>% drop_na() %>% select(direction) %>% group_by(direction) %>% count() %>% 
      full_join(data.frame(direction = c("backward", "simultaneous", "forward"), join_by = "direction")) %>%
      mutate(n = ifelse(is.na(n), 0, n)) %>% 
      pivot_wider(., names_from = direction, values_from = n) %>% select(backward, simultaneous, forward)
    distri_fbs <- rbind(distri_fbs, fbs)
    
    setTxtProgressBar(progressbar, n)
  }
  close(progressbar)
  
  ## et là on va comparer fbs aux permutés
  df_fbs[i, 1] = name_path
  for (indx_cond in 1:length(conditions[2:4])){ # pour ne prendre que les conditions strictes
    df_brute_fbs <- df_relations_allinfos %>% ungroup() %>% drop_na() %>% # mettre dans le bon sens
      dplyr::select(direction, path) %>% group_by(path, direction) %>% dplyr::summarise(n = dplyr::n())  %>% group_by(path) %>% 
      pivot_wider(names_from = direction, values_from = n) %>% select(path, backward, simultaneous, forward) %>% filter(path == name_path) 
    result_fbs <- pvalue_comp(x = df_brute_fbs[indx_cond+1][[1]], distrib = as.list(distri_fbs[,indx_cond])[[1]], tail = "two")
    df_fbs[i, (1+indx_cond*3)-2] = result_fbs[[1]]; df_fbs[i, (1+indx_cond*3)-1] = result_fbs[[3]]; df_fbs[i, 1+indx_cond*3] = result_fbs[[2]]
  }
  df_distri_fbs <- rbind(df_distri_fbs, t(colMeans(distri_fbs, na.rm = TRUE)))
    
  # # maintenant on va comparer les données brutes aux données permutées
  df_pvalue[i, 1] = name_path
  for (indx_cond in 1:length(conditions)){
    result_comp <- pvalue_comp(x=meancluster[indx_cond], distri=df_meancluster[,indx_cond], tail = "two")
    df_pvalue[i, (1+indx_cond*3)-2] = result_comp[[1]]; df_pvalue[i, (1+indx_cond*3)-1] = result_comp[[3]]; df_pvalue[i, 1+indx_cond*3] = result_comp[[2]]
  }
  # tableau total concernant les informations sur les gènes
  df_from <- df_relations_allinfos %>%
    select(from, name_from, pos_from, date_from) %>%
    dplyr::rename(ids = from, name = name_from, pos = pos_from, date = date_from)
  df_to <- df_relations_allinfos %>% ungroup() %>%
    select(to, name_to, pos_to, date_to) %>%
    dplyr::rename(ids = to, name = name_to, pos = pos_to, date = date_to)
  df_genes <- bind_rows(df_from, df_to) %>% distinct() %>%
    select(-contains("delta"), -contains("direction"), -contains("weights")) %>%
    left_join(df_degrees, by=c("ids"="ids")) %>% mutate(path = name_path)
  df_genes_allpaths <- rbind(df_genes_allpaths, df_genes)
}

# on remet bien les noms des colonnes parce qu'il peut y avoir eu des pertes avec rbind()
names(df_pvalue) <- c("path", "backsimult", "backsimult_sizeeffect", "p_backsimult", "backstrict", 
                      "back_sizeeffect", "p_backstrict", "simult", "simult_sizeeffect", "p_simult", "forstrict", "for_sizeeffect",
                      "p_forstrict", "forsimult", "forsimult_sizeeffect","p_forsimult")
names(df_correlation) <- c("path", "r_correlation", "pvalue")
names(df_fbs) <- c("path", "backward_signe", "backward_sizeeffect", "p_backward", "simult_signe", 
                   "simult_sizeeffect", "p_simult", "forward_signe", "forward_sizeeffect", "p_forward")
df_correlation$r_correlation = as.numeric(df_correlation$r_correlation)
df_correlation$pvalue = as.numeric(df_correlation$pvalue)


#### SAUVEGARDE ####
save.image("/home/fpicolo/Desktop/Pathways/birth-animals/2024/session_16012025.RData")


#### INFOS VOIES ####
## reparer les problèmes de mismatch
df_diff_genes <- df_genes_allpaths %>% drop_na() %>% group_by(name) %>% dplyr::summarise(n_date = n_distinct(date)) %>% 
  filter(n_date > 1) %>% ungroup()
df_genes_diff <- df_genes_allpaths %>% semi_join(df_diff_genes, by = "name") %>% 
  ungroup() %>% select(name, date, path) %>% arrange(name)
df_dates_mode <- df_genes_allpaths %>% group_by(name) %>%
  dplyr::summarise(date_mode = min(date), .groups = "drop")
df_dates_with_modes <- df_genes_allpaths %>% left_join(df_dates_mode, by = "name") %>% 
  select(-date) %>% dplyr::rename(date = date_mode)
df_genes_allpaths <- df_dates_with_modes
# permet d'identifier les problèmes de direction : 
df_diff_directions <- df_relations_allpaths %>% drop_na() %>%
  group_by(name_from, name_to) %>% dplyr::summarise(n_directions = n_distinct(direction)) %>%
  filter(n_directions > 1) %>% ungroup()
df_relations_diff <- df_relations_allpaths %>% semi_join(df_diff_directions, by = c("name_from", "name_to")) %>% 
  ungroup() %>% select(name_from, name_to, date_from, date_to, delta_age, direction, path) %>% arrange(name_from)

## réparer les bêtises pour les mauvaises directions le temps de voir ce qui ne fonctionne pas avant
df_mode_dates <- df_relations_allpaths %>% group_by(name_from, name_to) %>%
  dplyr::summarise(date_from_mode = min(date_from), date_to_mode = calc_mode(date_to), .groups = 'drop')
df_relations_with_modes <- df_relations_allpaths %>% left_join(df_mode_dates, by = c("name_from", "name_to")) %>% 
  select(-date_from, -date_to, -direction, -delta_age) %>% rename(date_from = date_from_mode, date_to = date_to_mode) %>% 
  mutate(delta_age = date_to-date_from, 
         direction = ifelse(delta_age<0, "backward", ifelse(delta_age>0, "forward", "simultaneous")))

## puis remplacer les p_sum_dir et t_sum_dir. 
df_relations_allpaths <- df_relations_with_modes


# compter le nombre d'interaction par voies
df_nbinteraction <- df_relations_allpaths %>% ungroup() %>% drop_na() %>% group_by(path) %>% count() %>% rename(nb_interactions = n)
# compter le nombre de gènes par voies
df_nbgenes <- df_genes_allpaths %>% ungroup() %>% filter(!grepl(";", name), !grepl("cpd", name)) %>% drop_na() %>% group_by(path) %>% count() %>% rename(nb_genes = n)
# compter le nombre de gènes + groupes par voies
df_genes_allpaths %>% ungroup() %>% drop_na() %>% group_by(path) %>% count()

#### INFOS GENES ####
# somme des genes par clades pour l'ensemble des voies # infos avec figure 3
# t_gene_by_date <- df_genes_allpaths %>% ungroup %>% drop_na() %>% select(name, date) %>% unique() %>% group_by(date) %>% dplyr::summarise(n = dplyr::n()) %>% arrange(-n)
t_gene_by_date <- df_genes_allpaths %>% ungroup() %>% filter(!grepl(";", name), !grepl("cpd", name)) %>% drop_na() %>% ungroup() %>% select(name, date, -ids) %>% unique() %>% group_by(date) %>% dplyr::summarise(n = dplyr::n()) %>% arrange(date)
# somme des gènes par clades par voie
p_gene_by_date <- df_genes_allpaths %>% drop_na() %>% select(name, date, path) %>% unique() %>% group_by(path, date) %>% dplyr::summarise(n = dplyr::n()) %>% arrange(-n)

## Genes de l'immunité
# test du chi2 pour montrer que les gènes impliqués dans l'immunité sont arrivés + après le clade des vertébrés
list_immun_path <- c("B cell receptor", "C type.lectin receptor", "Chemokine", "FC epsilon RI", "IL 17", 
                "NOD like receptor", "RIG I like receptor", "T cell receptor", "Toll like receptor")
chi2_immune <- df_genes_allpaths %>% ungroup() %>% filter(!grepl(";", name), !grepl("cpd", name)) %>% drop_na() %>% select(name, date, path) %>% 
  dplyr::mutate(immune = ifelse(path %in% list_immun_path, "yes", "no"), vertebrate = ifelse(date >= 12, "yes", "no")) %>% 
  select(name, immune, vertebrate) %>% unique() %>% group_by(immune, vertebrate) %>% dplyr::summarise(n = dplyr::n()) %>% 
  xtabs(n ~ immune + vertebrate, .) %>% chisq.test(.)
## peut encore être améliorer ?? là plusieurs voies qui sont immune et non sont comptées 2 fois.


## Couleurs pour KEGG
df_recup <- data.frame(matrix(NA, nrow = 0, ncol = 7))
for(i in 1:length(files_paths)){
  class_path <- xmlParse(paste(doss_paths, files_paths[[i]], sep=""))
  xml_path <- xmlRoot(class_path) # le fichier xml interprété par R
  list_result <- recup_networks(fxml = xml_path) # on sépare la voie en plusieurs informations dans des tableaux
  GRN <- list_result[[1]]; name_path <- list_result[[2]]
  
  print(paste("###", i, name_path, "###", sep=" "))
  
  GRN <- set_ages(GRN = GRN) # on associe les âges aux gènes (si paralogues : le plus ancien, si groupe : le plus récent)
  GRN <- set_names(GRN = GRN) # on ajoute les noms des gènes pour plus tard la méga voie
  subdata <- df_genes_allpaths %>% filter(path == name_path) %>% 
    mutate(color = palette$color[match(date,palette$clade)],
           kegg_id = ifelse(!grepl(";", name), sapply(GRN$items$names[match(ids, GRN$items$ids)], extract_first_element), name)) %>%
    drop_na()
  df_recup <- rbind(df_recup, subdata)
}
write.table(df_recup, "df_color.csv", sep = ";", row.names = FALSE)

#### INFOS RELATIONS ####
## Direction (forward, backward, simultaneous)
# somme des directions par voie
p_sum_direction <- df_relations_allpaths %>% ungroup %>% drop_na() %>% group_by(path) %>% dplyr::count(direction) %>% 
  pivot_wider(names_from = direction, values_from = n)
# somme des directions toute voie confondue
t_sum_direction <- df_relations_allpaths %>% ungroup() %>% drop_na() %>% select(name_from, name_to, direction) %>% 
  unique() %>% dplyr::count(direction) %>% pivot_wider(names_from = direction, values_from = n)
# compter les % de for/back/simul par voies
df_perc_fbs <- df_relations_allpaths %>% ungroup() %>% drop_na() %>% 
  dplyr::select(direction, path) %>% group_by(path, direction) %>% dplyr::summarise(n = dplyr::n())  %>% group_by(path) %>%
  dplyr::mutate(percentage = n / sum(n) * 100) %>% select(-n) %>% 
  pivot_wider(names_from = direction, values_from = percentage)


## Delta
# # petites recherches pour les deltas
df_relations_allpaths %>% drop_na() %>% group_by(path) %>% unique() %>% dplyr::mutate(label = paste(date_from,":", date_to, sep = "")) %>%
   group_by(delta_age, label, date_from, date_to) %>% filter(path=="Hippo") %>% dplyr::summarise(freq=dplyr::n()) %>% filter(delta_age==0)
# # somme de chaque delta
df_relations_allpaths %>% drop_na() %>% group_by(path) %>% unique() %>% dplyr::mutate(label = paste(date_from,":", date_to, sep = "")) %>%
   group_by(delta_age, label, date_from, date_to) %>% filter(path=="Hippo") %>% dplyr::count(delta_age)
# # compter le nombre d'interactions / type de delta
df_relations_allpaths %>% filter(path=="Hippo") %>% dplyr::mutate(label = paste(date_from, ":", date_to, sep="")) %>% group_by(delta_age) %>%
   dplyr::count(delta_age) %>% drop_na()


#### PLOTS ####
pdf("figure3-16012025.pdf", width = 10, height = 7)
plot_cumuldistribirth(df_genes = df_genes_allpaths)       # figure 3 
dev.off()
pdf("figure4-16012025.pdf", width = 10, height = 7)
plot_relations(df_relations = df_relations_allpaths)      # figure 4
dev.off()
pdf("figure5-16012025.pdf", width = 10, height = 7.8)
plot_distrirankbirth(df_genes = df_genes_allpaths)        # figure5
dev.off()
pdf("suppldata1-16012025.pdf", width = 10, height = 7)
p <- plot_delta(df_relations = df_relations_allpaths)     # suppl 1
dev.off()
list_path = unique(df_relations_allpaths$path)
pdf("suppldata2-16012025.pdf")                            # suppl 2
for(voie in list_path){
  plot_distribirthdelta(df_relations = df_relations_allpaths, name_path = voie)
}
dev.off()

names(p)[1:5] <- c("delta_min", "delta_Q1", "delta_med", "delta_Q3", "delta_max")
p
q <- plot_birth(df_genes = df_genes_allpaths)             # infos suppl 4
q$range <- paste("[", q$xmin_final, ";", q$xmax_final, "]", sep = "")
q <- q[, !(names(q) %in% c("xmin_final", "xmax_final"))]
names(q)[1:7] <- c("birth_min", "birth_Q1", "birth_med", "birth_Q3", "birth_max", "path", "birth_range")
q


### tableau suppl 4 
merged_table <- merge(df_nbinteraction, df_nbgenes, by = "path")
merged_table <- merge(merged_table, p_sum_direction, by = "path")
merged_table <- left_join(merged_table, q, by = "path")
merged_table <- merge(merged_table, df_correlation, by = "path")
merged_table <- merge(merged_table, p, by = "path")
merged_table <- merge(merged_table, df_pvalue, by = "path")
factorsvoies <-read.csv(file.choose(),sep=";") ## le fichier est dans 2024, factor
merged_table <- merge(merged_table, factorsvoies, by = "path") 
merged_table <- merge(merged_table, df_fbs, by = "path")
write.table(merged_table, "suppldata4-16012025.csv", row.names = FALSE, sep = ";")

###### KO ######
### Pour aller regarder les KO 
list_hasard_path <- merged_table %>% filter(chi2_pv > 0.05 & pvalue > 0.05 & p_forstrict > 0.05 & p_forsimult > 0.05 & p_simult > 0.05 & p_backstrict > 0.05 & p_backsimult > 0.05) %>% .$path
tab_hasard <- df_genes_allpaths %>% ungroup() %>% filter(!grepl(";", name), !grepl("cpd", name)) %>% drop_na() %>% select(name, date, path) %>% 
  dplyr::mutate(hasard = ifelse(path %in% list_hasard_path, "yes", "no"))
# liste des gènes qui sont dans des voies hasards et dans des voies pas hasard ! 
tab_mixte <- tab_hasard %>% group_by(name) %>% summarise(hasard_yes = any(hasard == "yes"), hasard_no = any(hasard == "no")) %>% ungroup() %>% 
  filter(hasard_yes & hasard_no) # %>% select(name)
# liste des gènes qui sont uniquement dans des voies hasards
tab_fullH <- tab_hasard %>% group_by(name) %>% summarise(hasard_yes = any(hasard == "yes"), hasard_no = any(hasard == "no")) %>% ungroup() %>% 
  filter(hasard_yes == TRUE & hasard_no == FALSE)
# liste des gènes qui sont uniquement dans des voies nécessité
tab_fullN <- tab_hasard %>% group_by(name) %>% summarise(hasard_yes = any(hasard == "yes"), hasard_no = any(hasard == "no")) %>% ungroup() %>% 
  filter(hasard_yes == FALSE & hasard_no == TRUE)
tab_hasard %>% filter(name == "WAS")




