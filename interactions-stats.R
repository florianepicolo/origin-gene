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
  for(i in 1:length(names)){items_names[i] <- strsplit(names[i], " ")[[1]]} # un seul crochet si on récupère tout # récupérer les différents nom des gènes
  
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
  # V(graph)$label <- GRN$items$gene_name
  times_palette = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 25)
  times_colors = times_palette[GRN$items$times]
  # plot(graph,
  #      layout = layout_as_tree(graph, circular = TRUE),
  #      vertex.label.cex = 0.8,  # Taille du texte
  #      vertex.size = 10,        # Taille des nœuds
  #      vertex.label.font = 2,   # Épaisseur de la police des étiquettes
  #      edge.arrow.size = 0.5,    # Taille des flèches des arêtes pour les graphes dirigés
  #      vertex.color = times_colors   # Utiliser les couleurs définies pour les nœuds
  # )
  plot(graph,
       layout = layout_with_fr(graph),
       vertex.label.cex = 0.8,  # Taille du texte
       vertex.size = 10,        # Taille des nœuds
       vertex.label.font = 2,   # Épaisseur de la police des étiquettes
       edge.arrow.size = 0.5,    # Taille des flèches des arêtes pour les graphes dirigés
       vertex.color = times_colors   # Utiliser les couleurs définies pour les nœuds
  )
  # plot(graph,
  #      layout = layout_as_tree(graph),
  #      vertex.label.cex = 0.8,  # Taille du texte
  #      vertex.size = 10,        # Taille des nœuds
  #      vertex.label.font = 1,   # Épaisseur de la police des étiquettes
  #      edge.arrow.size = 0.5,    # Taille des flèches des arêtes pour les graphes dirigés
  #      vertex.color = times_colors
  # )
  # plot(graph,
  #      layout = layout_with_kk(graph),
  #      vertex.label.cex = 0.8,  # Taille du texte
  #      vertex.size = 10,        # Taille des nœuds
  #      vertex.label.font = 2,   # Épaisseur de la police des étiquettes
  #      edge.arrow.size = 0.5,    # Taille des flèches des arêtes pour les graphes dirigés
  #      vertex.color = times_colors   # Utiliser les couleurs définies pour les nœuds
  # )

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
  df_interactions <- data.frame(
    from = interactions[,1],                 # premier élément
    to = interactions[,2]                    # deuxième élement
  ) %>% # merge(., degrees[, c("ids", "degree_out")], by.x = "from", by.y = "ids", all.x = TRUE)  %>%
    mutate(positions = NA, delta_age = NA, direction = NA) 
  
  for (subgraph in subgraphs){
    up_node <- which(degree(subgraph, mode = "in") == 0 & degree(subgraph, mode = "out") > 0) # on récupère l'ensemble des noeuds de début de voie 
    # down_node <- which(degree(subgraph, mode = "in") > 0 & degree(subgraph, mode = "out") == 0)
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
              # degree_out = degrees$degree_out[degrees$ids==from],
              positions = position_in_path,  # nouvelle position
              delta_age = GRN$items$times[GRN$items$ids == to] - GRN$items$times[GRN$items$ids == from],
              direction = ifelse((GRN$items$times[GRN$items$ids == to] - GRN$items$times[GRN$items$ids == from]) > 0, "backward", ifelse((GRN$items$times[GRN$items$ids == to] - GRN$items$times[GRN$items$ids == from]) < 0, "forward", "simultaneous"))
            )
            df_interactions <- rbind(df_interactions, new_row)

          }else{
            # on fait delta de l'age de B - A
            df_interactions$delta_age[relation] = GRN$items$times[GRN$items$ids==to] - GRN$items$times[GRN$items$ids==from]
            df_interactions$direction[relation] = ifelse(df_interactions$delta_age[relation]>0, "backward", ifelse(df_interactions$delta_age[relation]<0, "forward", "simultaneous"))
            
            df_interactions$positions[relation] = position_in_path
          }
          position_in_path = position_in_path + 1
        }
      }
    }
  }
  df_interactions <- distinct(df_interactions)
  df_interactions <- df_interactions %>% group_by(from) %>% dplyr::mutate(weights = dplyr::n()) 
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
  med <- median(distrib, na.rm=T) # médiane de la distribution des 1000 permutations
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
      # plot(suba,
      #      layout = layout_with_fr(suba),
      #      vertex.label.cex = 0.8,  # Taille du texte
      #      vertex.size = 10,        # Taille des nœuds
      #      vertex.label.font = 2,   # Épaisseur de la police des étiquettes
      #      edge.arrow.size = 0.5,    # Taille des flèches des arêtes pour les graphes dirigés
      #      vertex.color = components(suba)$membership   # Utiliser les couleurs définies pour les nœuds
      # )
    }else{all_components <- c(all_components, 0)}
  }
  all_components <- all_components[-1] # on vire les NA au fur et à mesure comme on les a ajouté au début
  mean_distri = mean(all_components)
  return(mean_distri) # on sort une lsite de taille de subgrapph qui respecte la condition donné
}

#### PLOTS ####
plot_posids <- function(df_genes, name_path){
  p <- df_genes %>% drop_na() %>% filter(path == name_path) %>% group_by(ids, name, pos, date) %>% dplyr::summarise(n=n(), min_rank=min(pos), max_rank=max(pos)) %>% 
    arrange(min_rank, max_rank, -n, ids) %>%  mutate(ids=factor(name, levels=name %>% unique())) %>% 
    ggplot(aes(x= reorder(ids, max_rank), y=pos, size=n, col=factor(date, levels=1:25))) + geom_point() +
    theme(axis.text.x = element_text(angle=90)) + ggtitle(name_path) + labs(x = "gene name", y = "pathway position", col= "branch of origin", size = "n") +
    scale_color_manual(values = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 25), drop = FALSE)
  print(p)
} # figure ggplot-generankbirthpath
plot_relations <- function(df_relations){
  # p <- df_relations %>% drop_na() %>% group_by(path) %>% count(direction) %>% ungroup() %>% ggplot(., aes(x=n, y=path, fill=direction)) + geom_bar(stat = "identity") +
  #   scale_fill_manual(values = c("#66C2A5","#E6F598", "#FDAE61")) + theme_minimal() +
  #   ggtitle("Distribution of direction of interaction by pathway") +
  #   labs(x = "number of interactions", y = "pathway")
  # 
  
  df_relations$dir <- factor(df_relations$direction, levels = c("backward", "simultaneous", "forward"))
  p_direction <- df_relations %>% drop_na() %>% 
    select(dir, path) %>% group_by(path, dir) %>% dplyr::summarise(n = dplyr::n())  %>% group_by(path) %>%
    dplyr::mutate(percentage = n / sum(n) * 100) %>% filter(dir == "forward") %>% arrange(-percentage)  %>% pull(path)
  
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
    scale_fill_manual(values = c("backward" = "#66C2A5", "simultaneous" = "#FDAE61", "forward" = "#E6F598"), 
                      breaks = c("backward", "simultaneous", "forward")) +
    theme_minimal() +
    ggtitle("Distribution of direction of interaction by pathway (%)") +
    labs(x = "percentage of direction of interaction", y = "path", fill="direction")

  print(p)
} # figure ggplot-distributiondirection
plot_distribirth <- function(df_genes){
  # p <- df_genes %>% drop_na() %>% select(name, date) %>% unique() %>%
  #   ggplot(., aes(x=factor(date, levels=1:25))) + geom_bar(position="stack") +
  #   labs(x = "birth clade", y = "number of gene") + ggtitle("Distribution of birth gene") +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
  #   scale_x_discrete(labels = date_clade$name, drop = FALSE)
  
  d <- df_genes %>% drop_na() %>% mutate(nwclade = case_when(
    date %in% c(2, 3) ~ 2,
    date %in% c(4, 5, 6) ~ 4,
    date %in% c(7, 8) ~ 7,
    date %in% c(12, 13) ~ 12,
    date %in% c(18, 19, 20, 21, 22) ~ 18,
    TRUE ~ date))
  d$nwclade <- as.integer(d$nwclade)

  new_date_clade <- data.frame(
    num = c(1, 2, 4, 7, 9, 10, 11, 12, 14, 15, 16, 17, 18, 23, 24),
    age = c(1300, 765, 743, 708, 635, 635, 570, 563, 462, 429, 415, 350, 330, 166, 98),
    name = c("Opisthokonta", "Metazoa", "Eumetazoa", "Bilateria", "Deuterostomia", "Chordata",
             "Olfactora", "Vetebrata", "Gnathostomata", "Osteichthyes", "Sarcopterygii", "Tetrapoda",
             "Amniota", "Mammalia", "Theria")
  )
  new_date_clade$label <- paste(new_date_clade$name, " (", new_date_clade$age, "My)", sep = "")
  p <- d %>% ungroup() %>% filter(!grepl(";", name), !grepl("cpd", name)) %>%
    drop_na() %>% select(name, nwclade) %>% unique() %>%
    ggplot(., aes(x=factor(nwclade))) + geom_bar(position="stack") +
    labs(x = "branch of origin", y = "number of gene") + ggtitle("Distribution of gene") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
    scale_x_discrete(labels = new_date_clade$label, drop = FALSE)
  print(p)  
} #figure ggplot-distributionbirth #figure 3
plot_cumuldistribirth <- function(df_genes){
  d <- df_genes_allpaths %>% drop_na() %>% dplyr::mutate(nwclade = case_when(
    date %in% c(2, 3) ~ 2,
    date %in% c(4, 5, 6) ~ 4,
    date %in% c(7, 8) ~ 7,
    date %in% c(12, 13) ~ 12,
    date %in% c(18, 19, 20, 21, 22) ~ 18,
    TRUE ~ date))
  d$nwclade <- as.integer(d$nwclade)
  
  new_date_clade <- data.frame(
    num = c(1, 2, 4, 7, 9, 10, 11, 12, 14, 15, 16, 17, 18, 23, 24),
    age = c(1300, 765, 743, 708, 635, 635, 570, 563, 462, 429, 415, 350, 330, 166, 98),
    name = c("Opisthokonta", "Metazoa", "Eumetazoa", "Bilateria", "Deuterostomia", "Chordata",
             "Olfactora", "Vetebrata", "Gnathostomata", "Osteichthyes", "Sarcopterygii", "Tetrapoda",
             "Amniota", "Mammalia", "Theria")
  )
  new_date_clade$label <- paste(new_date_clade$name, " (", new_date_clade$age, "My)", sep = "")
  
  color_palette <- paletteer_c("ggthemes::Sunset-Sunrise Diverging", 25)
  p <- d %>% ungroup() %>% filter(!grepl(";", name), !grepl("cpd", name)) %>%
    drop_na() %>% select(name, nwclade) %>% unique() %>% 
    group_by(nwclade) %>%
    dplyr::summarise(count_gene = dplyr::n()) %>%
    arrange(nwclade) %>%
    dplyr::mutate(
      cumulative_count = cumsum(count_gene),
      cumulative_percentage = cumsum(count_gene) / sum(count_gene) * 100
    ) %>%
    ggplot(aes(x = factor(nwclade))) +
    geom_bar(aes(y = cumulative_count, fill = nwclade), stat = "identity") +
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
      x = "branch of origin", 
      y = "number of gene / cumulative percentage of gene", 
      title = "Distribution of birth gene (cumulative percentage and raw data)"
    ) +
    scale_fill_gradientn(colours = color_palette, name = "Ajout", guide = "none") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
    scale_x_discrete(labels = new_date_clade$label, drop = FALSE)
  
  # p <- d %>% drop_na() %>% select(name, nwclade) %>% unique() %>%
  #   group_by(nwclade) %>% summarise(count = n()) %>% arrange(nwclade) %>%
  #   mutate(cumulative_count = cumsum(count), 
  #          cumult_perc = cumsum(count)/sum(count)*100,
  #          diff_count = cumulative_count - lag(cumulative_count, default = 0)) %>% 
  #   ggplot(., aes(x = factor(nwclade), y = cumulative_count, fill = diff_count)) +
  #   geom_bar(stat = "identity") +
  #   geom_text(aes(label = ifelse(diff_count > 0, paste("+", round(diff_count,2), "%"), "")),
  #             vjust = -0.5, size = 3, color = "black") +
  #   scale_fill_gradientn(colours = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 25), name = "Ajout") +
  #   labs(x = "birth clade", y = "number of gene", 
  #        title = "Cumulative distribution of birth gene") +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))+
  #   scale_x_discrete(labels = new_date_clade$label, drop = FALSE)
  print(p)  
}
plot_delta <- function(df_relations){
  p_delta <- df_relations %>% drop_na() %>% select(delta_age, pos_from, pos_to, path) %>% 
    group_by(path, pos_from, pos_to) %>% unique() %>% dplyr::count(delta_age)
  p_median_delta <- p_delta %>% group_by(path) %>% dplyr::summarise(median = median(delta_age)) %>% arrange(median) %>% pull(path)
  p_delta <- p_delta %>% dplyr::mutate(path = factor(path, levels=p_median_delta))
  p <- p_delta %>% ggplot(.,  aes(x = delta_age, y = path)) +
    geom_boxplot() + labs(x = "delta", y = "path") +
    ggtitle("Distribution des deltas par path") +
    scale_x_continuous(breaks = seq(-25, 25, by = 5)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red")
  plot(p)
  tab_plot <- layer_data(p) %>% select(xmin, xlower, xmiddle, xupper, xmax) 
  tab_plot$path <- rep(p_median_delta, length.out = nrow(tab_plot))
  return(tab_plot)
} #figure ggplot-distributiondelta
plot_distrirankbirth <- function(df_genes){
  order_correlation <- df_correlation %>% filter(pvalue<0.05) %>% arrange(r_correlation, path) %>% rbind(., df_correlation %>% filter(pvalue>=0.05) %>% arrange(r_correlation, path))
  p <- df_genes %>% drop_na() %>% mutate(path = factor(path, levels = order_correlation$path)) %>% ggplot(., aes(y = "",  fill = factor(date, levels=1:25))) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 25), drop = FALSE) +
    coord_cartesian() +
    theme_void() +
    facet_grid(path ~ pos, switch = "both") +
    theme(strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0)) +
    labs(x = "path", y = "pathway position", fill = "branch of origin") +
    ggtitle("Distribution of branch of origin by path and pathway position")
  print(p)
} #figure ggplot-distributionbirthrank
plot_distribirthdelta <- function(df_relations, name_path){
  p <- df_relations %>% drop_na() %>% group_by(path) %>% unique() %>% mutate(label = paste(date_from,":", date_to, sep = "")) %>% 
    group_by(delta_age, label, date_from, date_to) %>% filter(path==name_path) %>% summarise(freq=n()) %>% 
    ggplot(., aes(delta_age, freq, fill = factor(date_from, levels=1:25))) +
    geom_bar(position = "stack", stat = "identity") + theme_minimal() +
    scale_fill_manual(values = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 25), drop = FALSE) + 
    labs(x = "delta", y = "distribution of branch of origin", fill = "branch of origin") +
    ggtitle(paste("Distribution of branch of origin by delta for", name_path))
  print(p)
} #figure ggplot-distridelta


#### FILES ####
file_KEGG <- "/home/fpicolo/Desktop/Pathways/birth-animals/nouveau-run/p-allinfos-KEGG.csv"
doss_paths <- "/home/fpicolo/Desktop/Pathways/birth-animals/paths/"

data_KEGG <- read.csv(file_KEGG, sep=";")
files_paths <- dir(doss_paths) # récupère tous les fichiers d'un dossier



#### INITIALISATION ####
npermut = 1000
conditions <- c(function(x){x[1]>=x[2]}, function(x){x[1]>x[2]}, function(x){x[1]==x[2]}, function(x){x[1]<x[2]}, function(x){x[1]<=x[2]})
df_pvalue = data.frame(matrix(NA, nrow = length(files_paths), ncol = 1+length(conditions)*3))
df_relations_allpaths <- data.frame(matrix(NA, nrow = 0, ncol = 12))
df_relations_allpaths_permutees <- data.frame(matrix(NA, nrow = 0, ncol = 12))
df_genes_allpaths <- data.frame(matrix(NA, nrow = 0, ncol = 7))
df_correlation <- data.frame(matrix(NA, nrow = length(files_paths), ncol = 3))
date_clade <- data.frame(num = 1:25,
                         name = c("Yeast", "Porifera", "Placozoa", "Ctenophora", "Cnidaria", "Xenacoelomorpha", "Spiralia",
                                  "Ecdysozoa", "Echinodermata", "Cephalochordata", "Tunicata", "Petromyzontidae", "Myxini",
                                  "Chondrichthyes", "Actinopterygii", "Coelacanthidae", "Amphibia", "Testudines", "Aves",
                                  "Crocodylia", "Squamata", "Sphenodontia", "Prototheria", "Eutheria", "Metatheria"),
                         age = c(1110, 725, 550, 109, 588, 550, 610, 446, 596, 520, 446, 416, 360, 413, 396, 350, 319, 194,
                                 109, 87, 189, 200, 166, 98, 66))
palette <- tibble(color = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 25), clade=1:25)

#### PROGRAMM ####
for(i in 1:length(files_paths)){
  class_path <- xmlParse(paste(doss_paths, files_paths[[i]], sep=""))
  xml_path <- xmlRoot(class_path) # le fichier xml interprété par R
  list_result <- recup_networks(fxml = xml_path) # on sépare la voie en plusieurs informations dans des tableaux
  GRN <- list_result[[1]]; name_path <- list_result[[2]]
  
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
    # tempGRN <- permutation(GRN = GRN, var = "relations") # "relation"/"times" possible aussi
    tempGRN <- permutation(GRN = GRN, var = "times")
    # visual_graph(tempGRN) 
    for (indx_cond in 1:length(conditions)){
      cond = conditions[[indx_cond]]
      df_meancluster[n, indx_cond] = subgraphs_given_times_condition(GRN = tempGRN, condition = cond)
    }
    if(n %% 100 == 0){setTxtProgressBar(progressbar, n)}
  }
  close(progressbar)
  
  # on se prépare un petit jeu de data permutées
  df_relations_permutees <- get_subgraph(tempGRN) %>% mutate(path = name_path)
  df_relations_allinfos_permutees <- df_relations_permutees %>% dplyr::rename(pos_from = positions) %>% mutate(pos_to = pos_from +1) %>% 
    left_join(GRN$items, by = c("from" = "ids")) %>% dplyr::rename(date_from = times, name_from = gene_name) %>% 
    left_join(GRN$items, by = c("to" = "ids")) %>% dplyr::rename(date_to = times, name_to = gene_name) %>% 
    select(from, to, name_from, name_to, pos_from, pos_to, date_from, date_to, delta_age, direction, weights, path)
  df_relations_allpaths_permutees <- rbind(df_relations_allpaths_permutees, df_relations_allinfos_permutees)
  
  # maintenant on va comparer les données brutes aux données permutées
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
names(df_pvalue) <- c("path", "forsimult", "forsimult_sizeeffect", "p_forsimult", "forstrict", 
                      "for_sizeeffect", "p_forstrict", "simult", "simult_sizeeffect", "p_simult", "backstrict", "back_sizeeffect",
                      "p_backstrict", "backsimult", "backsimult_sizeeffect","p_backsimult")
names(df_correlation) <- c("path", "r_correlation", "pvalue")
df_correlation$r_correlation = as.numeric(df_correlation$r_correlation)
df_correlation$pvalue = as.numeric(df_correlation$pvalue)


#### SAUVEGARDE ####
# save.image("session.RData")

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

# Couleurs pour KEGG
df_color <- df_genes_allpaths %>% mutate(color = palette$color[match(date, palette$clade)], 
                                         kegg_id = sapply(GRN$items$names[match(ids, ids)], extract_first_element))


#### INFOS RELATIONS ####
## Direction (forward, backward, simultaneous)
# somme des directions par voie
p_sum_direction <- df_relations_allpaths %>% ungroup %>% drop_na() %>% group_by(path) %>% dplyr::count(direction) %>% 
  pivot_wider(names_from = direction, values_from = n)
# somme des directions toute voie confondue
t_sum_direction <- df_relations_allpaths %>% ungroup() %>% drop_na() %>% select(name_from, name_to, direction) %>% 
  unique() %>% dplyr::count(direction) %>% pivot_wider(names_from = direction, values_from = n)

# permet d'identifier les problèmes de direction : 
df_diff_directions <- df_relations_allpaths %>% drop_na() %>%
  group_by(name_from, name_to) %>% dplyr::summarise(n_directions = n_distinct(direction)) %>%
  filter(n_directions > 1) %>% ungroup()
df_relations_diff <- df_relations_allpaths %>% semi_join(df_diff_directions, by = c("name_from", "name_to")) %>% 
  ungroup() %>% select(name_from, name_to, date_from, date_to, delta_age, direction, path) %>% arrange(name_from)

## réparer les bêtises pour les mauvaises directions le temps de voir ce qui ne fonctionne pas avant
calc_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
df_mode_dates <- df_relations_allpaths %>% group_by(name_from, name_to) %>%
  summarise(date_from_mode = calc_mode(date_from), date_to_mode = calc_mode(date_to), .groups = 'drop')
df_relations_with_modes <- df_relations_allpaths %>% left_join(df_mode_dates, by = c("name_from", "name_to")) %>% 
  select(-date_from, -date_to, -direction, -delta_age) %>% rename(date_from = date_from_mode, date_to = date_to_mode) %>% 
  mutate(delta_age = date_from-date_to, 
         direction = ifelse(delta_age>0, "forward", ifelse(delta_age<0, "backward", "simultaneous")))

## puis remplacer les p_sum_dir et t_sum_dir. 
df_relations_allpaths <- df_relations_with_modes
# ici il faudra faire un chisq.test() entre les % des directions brute et celles des permutées. 


## Delta
# petites recherches pour les deltas
df_relations_allpaths %>% drop_na() %>% group_by(path) %>% unique() %>% dplyr::mutate(label = paste(date_from,":", date_to, sep = "")) %>%
  group_by(delta_age, label, date_from, date_to) %>% filter(path=="MAPK") %>% dplyr::summarise(freq=dplyr::n()) %>% filter(delta_age==0)
# somme de chaque delta
df_relations_allpaths %>% drop_na() %>% group_by(path) %>% unique() %>% dplyr::mutate(label = paste(date_from,":", date_to, sep = "")) %>% 
  group_by(delta_age, label, date_from, date_to) %>% filter(path=="MAPK") %>% dplyr::count(delta_age)

## Permutations
# on va regarder si les tailles des clusters respectant une condition sont du à l'aléatoire ou non : 
df_pvalue %>% filter(p_forsimult <0.05) %>% select(path, p_forsimult, forsimult) %>% arrange(forsimult)
df_pvalue %>% filter(p_forstrict <0.05) %>% select(path, p_forstrict, forstrict) %>% arrange(forstrict)
df_pvalue %>% filter(p_simult <0.05) %>% select(path, p_simult, simult) %>% arrange(simult)
df_pvalue %>% filter(p_backstrict <0.05) %>% select(path, p_backstrict, backstrict) %>% arrange(backstrict)
df_pvalue %>% filter(p_backsimult <0.05) %>% select(path, p_backsimult, backsimult) %>% arrange(backsimult)




#### PLOTS ####
plot_relations(df_relations = df_relations_allpaths) #figure4
plot_distribirth(df_genes = df_genes_allpaths) #figure3
plot_cumuldistribirth(df_genes = df_genes_allpaths) # new figure 3 ??? 
plot_delta(df_relations = df_relations_allpaths) #figure5
plot_distrirankbirth(df_genes = df_genes_allpaths) #figure7

list_path = unique(df_relations_allpaths$path)
for(voie in list_path){
  # il faut les sauvegarder
  plot_posids(df_genes = df_genes_allpaths, name_path = voie) 
  plot_distribirthdelta(df_relations = df_relations_allpaths, name_path = voie)
}





