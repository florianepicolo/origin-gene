#!/usr/bin/Rscript --vanilla

#### LIBRARY ####
packages_to_install <- c("XML", "plyr", "igraph", "dplyr", "tidyr", "magrittr", "ggplot2", "paletteer", "scales", "ggforce", "ggthemes", "cluster", "pheatmap")
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
  # nnames <- sapply(items_names, length) ## si on a tout récupérer 
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
check_ages <- function(GRN){
  GRN$items$times <- NA
  number_of_match <- id_match <- list()
  
  for(j in 1:nrow(GRN$items)) {
    ## on s'occupe des genes
    
    if(GRN$items$types[j]=="gene"){ #on a les dates que pr les gènes donc on ne prend qu'eux
      id_match[[j]] <- which(dat$kegg_id %in% unlist(strsplit(GRN$items$names[j], " ")[[1]])) #quand on a qu'un vectnames
      # id_match[[j]] <- which(dat$kegg_id %in% GRN$namesfromnetwork$allnames[GRN$namesfromnetwork$vectnames==GRN$items$ids[j]]) # on récupère les numéros de ligne dans dat qui correspondent à ces ids
      number_of_match[[j]] <- length(id_match[[j]][!is.na(id_match[[j]])])
      if(number_of_match[[j]] == 1) {GRN$items$times[j] = dat$num_clade[id_match[[j]]]}
      if(number_of_match[[j]] > 1) {GRN$items$times[j] = min(as.numeric(names(table(dat$num_clade[id_match[[j]][!is.na(id_match[[j]])]]))))}
    }
    ## et maintenant des groups ! 
    if(GRN$items$types[j]=="group"){ # groupe formé de gène-gène !
      result <- GRN$groups %>% filter(idgroup == GRN$items$ids[j]) %>% left_join(GRN$items, by = c("idgene" = "ids")) %>%
        left_join(GRN$namesfromnetwork, by = c("idgene" = "vectnames")) %>% left_join(dat, by = c("allnames" = "kegg_id"))
      GRN$items$times[j] = max(result$num_clade)
    }
  }
  return(GRN)
}
subgraph_stats <- function(GRN, pathname){
  graph <- graph_from_data_frame(GRN[[2]], directed=TRUE, vertices=GRN[[1]])
  filtered_names <- NULL
  for (node_id in V(graph)){
    for (n in node_id){
      incoming_nodes <- neighbors(graph, n, mode = "in")
      names <- V(graph)$name[incoming_nodes]
      if (length(names)>0){
        filtered_names <- union(names, filtered_names)
      }else{break}
    }
  }
  ## ici on a tous les gènes qui sont "in" d'une interaction mais pas uniquement les initiaux des voies
  filteredbis_names <- NULL 
  for (n in filtered_names){
    incoming_nodes <- neighbors(graph, n, mode = "in")
    names <- V(graph)$name[incoming_nodes]
    if (length(names)<=0){
      filteredbis_names <- union(n, filteredbis_names)
    }
  }
  
  ### on commence la boucle des stats
  corrpath <- statspath <- dates <- positions <- name <- genename <- date_stats <- NULL
  noms_colonnes <- c("subpath", "nbinteract", "mediane", "moyenne", "ecarttype", "date-mini", "date-max", "correlation", "corpvalue")
  date_stats <- data.frame(matrix(ncol = length(noms_colonnes)))
  inc_element = 1 # compter les élements
  
  for(innode in filteredbis_names){
    all_paths <- all_simple_paths(graph, from = innode, to = V(graph), mode = "out") 
    
    # Fonction pour vérifier si un chemin est inclus dans un autre chemin
    is_included <- function(path, other_paths) {
      for (other_path in other_paths) {
        if (length(path) < length(other_path) && all(path == other_path[1:length(path)])) {
          return(TRUE)
        }
      }
      return(FALSE)
    }
    
    # Filtrer les voies qui ne sont pas incluses dans d'autres voies
    filtered_paths <- all_paths[!sapply(all_paths, is_included, other_paths = all_paths)]
    paths_with_ids <- lapply(filtered_paths, function(path) as_ids(path))
    
    for(path in paths_with_ids){
      position_voie = 1
      for(i in 1:length(path)){ # i = les différentes positions
        id = path[i]
        if(GRN$items$types[GRN$items$ids==id] %in% c("gene", "group") & !is.na(GRN$items$times[GRN$items$ids==id])){
          n = GRN$items$names[GRN$items$ids==id]
          name[inc_element] = ifelse(GRN$items$types[GRN$items$ids==id] != "group", sapply(strsplit(as.character(n), " "), "[[", 1), paste0("group_", id))
          genename[inc_element] = ifelse(grepl("group", name[inc_element]), name[inc_element], dat$gene_name[dat$kegg_id==name[inc_element]])
          dates[inc_element] = GRN$items$times[GRN$items$ids==id]
          positions[inc_element] = position_voie
          inc_element = inc_element + 1
          position_voie = position_voie + 1
        }
      }
      
      # Calculez les statistiques des dates (médiane, moyenne, écart-type, minimum, maximum)
      d = dates[(inc_element-position_voie+1):(inc_element-1)] 
      pos = positions[(inc_element-position_voie+1):(inc_element-1)]
      
      ## les stats pr la sous-voie
      if(length(d)>2){
        correlation <- cor.test(d, pos, method = "pearson")
        newligne <- c(paste(path, collapse = " "), length(path), median(d), mean(d), sd(d), min(d), max(d), correlation$estimate, correlation$p.value)
      }else{
        newligne <- c(paste(path, collapse = " "), length(path), median(d), mean(d), sd(d), min(d), max(d))
      }
      date_stats <- rbind(date_stats, newligne)
    }
  }
  #### faire les stats pr la voie totale !
  dates_na <- which(is.na(dates))
  if(length(dates_na)>0){
    positions <- positions[-dates_na]
    dates <- dates[-dates_na]
    genename <- genename[-dates_na]
  }
  
  statspath <- list(
    median = median(dates),
    mean = mean(dates),
    sd = sd(dates),
    min = min(dates),
    max = max(dates))
  
  # # Calculez la corrélation entre les dates et les positions
  corrpath <- cor.test(dates, positions, method="pearson")
  
  colnames(date_stats) <- noms_colonnes
  # write.csv(date_stats, file = paste(pathname, "-correlation.csv", sep=""), row.names = FALSE)
  df <- data.frame(id=genename, rank=positions, birth=dates)
  return(list(statspath, corrpath, df))
}
compare <- function(x){ # x = relation time
  rescomp<-rep(NA, nrow(x)) # resultats sera sous forme de vecteur de taille nombre de relation
  rescomp[levels(x[,1]) [x[,1]] > levels(x[,2]) [x[,2]]] <- "foreward" # on regarde si A plus grand que B 
  rescomp[levels(x[,1]) [x[,1]] == levels(x[,2]) [x[,2]]] <- "simult" # si A == B 
  rescomp[levels(x[,1]) [x[,1]] < levels(x[,2]) [x[,2]]] <- "backward" # si A < B
  return(rescomp) 
}
pvalue_comp <- function(x, distrib, tail=c("one", "two")){ # x = nombre dans les données brutes, distrib = la distri de la simulation, tail pr pouvoir changer de test
  pvalue <- NA
  med <- median(distrib, na.rm=T) # médiane de la distribution
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
  return(list(results, pvalue)) #on retourne une liste avec les résutlats +/- et les pvalue
}
permutation <- function(GRN, namepath, nrep=1000){
  resmeantable <- matrix(NA, 5)
  condition <- function(x){x[1]>=x[2]}
  components_sizes_fore <- subgraphs_given_times_condition(GRN, condition)
  condition <- function(x){x[1]>x[2]}
  components_sizes_forestrict <- subgraphs_given_times_condition(GRN, condition)
  condition <- function(x){x[1]<=x[2]}
  components_sizes_back <- subgraphs_given_times_condition(GRN,condition)
  condition <- function(x){x[1]<x[2]}
  components_sizes_backstrict <- subgraphs_given_times_condition(GRN, condition)
  condition <- function(x){x[1]==x[2]}
  components_sizes_simu <- subgraphs_given_times_condition(GRN, condition)
  
  runs_foreward <- runs_forestrict <- runs_backward <- runs_backstrict <- runs_simult <- NA
  runs_foreward_p <- runs_forestrict_p <- runs_backward_p <- runs_backstrict_p <- runs_simult_p <- NA    
  
  mean_foreward <- mean_fstrict <- mean_backward <- mean_bstrict <- mean_simult <- NA
  for(j in 1:nrep) {
    tempGNR <- GRN # pour copier la structure notamment
    tempGNR$relations$entry2 <- tempGNR$relations$entry2[sample(1:length(tempGNR$relations$entry2))]
    # get components size under condition (add simultaneous apparitions)
    condition <- function(x){x[1]>=x[2]}
    temp_components_sizes_fore <- subgraphs_given_times_condition(tempGNR, condition)
    condition <- function(x){x[1]>x[2]}
    temp_components_sizes_forestrict <- subgraphs_given_times_condition(tempGNR, condition)
    condition <- function(x){x[1]<=x[2]}
    temp_components_sizes_back <- subgraphs_given_times_condition(tempGNR, condition)
    condition <- function(x){x[1]<x[2]}
    temp_components_sizes_backstrict <- subgraphs_given_times_condition(tempGNR, condition)
    condition <- function(x){x[1]==x[2]}
    temp_components_sizes_simu <- subgraphs_given_times_condition(tempGNR, condition)
    
    ## là on fait la moyenne, mais on peut changer avec median
    res <- permutation_test_2groups(components_sizes_fore, temp_components_sizes_fore, mean, nrep)
    runs_foreward[j] <- res[[1]][[1]] ; runs_foreward_p[j] <- res[[1]][[2]] ; mean_foreward[j] <- mean(res[[2]])
    res <- permutation_test_2groups(components_sizes_forestrict, temp_components_sizes_forestrict, mean, nrep)
    runs_forestrict[j] <- res[[1]][[1]] ; runs_forestrict_p[j] <- res[[1]][[2]] ; mean_fstrict[j] <- mean(res[[2]])
    res <- permutation_test_2groups(components_sizes_back, temp_components_sizes_back, mean, nrep)
    runs_backward[j] <- res[[1]][[1]] ; runs_backward_p[j] <- res[[1]][[2]] ; mean_backward[j] <- mean(res[[2]])
    res <- permutation_test_2groups(components_sizes_backstrict, temp_components_sizes_backstrict, mean, nrep)
    runs_backstrict[j] <- res[[1]][[1]] ; runs_backstrict_p[j] <- res[[1]][[2]] ; mean_bstrict[j] <- mean(res[[2]])
    res <- permutation_test_2groups(components_sizes_simu, temp_components_sizes_simu, mean, nrep)
    runs_simult[j] <- res[[1]][[1]] ; runs_simult_p[j] <- res[[1]][[2]] ; mean_simult[j] <- mean(res[[2]])
  }
  resmeantable[1,1] = mean(mean_foreward)
  resmeantable[2,1] = mean(mean_fstrict)
  resmeantable[3,1] = mean(mean_backward)
  resmeantable[4,1] = mean(mean_bstrict)
  resmeantable[5,1] = mean(mean_simult)
  # print(resmeantable)
  # plot
  # pdf(paste(pathname, ".pdf", sep=""), width=12, height=6)
  par(mfrow=c(1,3)) # on sort 3 éléments
  seuil = .05
  labels <- c("minus", "minus pval<.05", "plus", "plus pval<.05")
  # foreward
  counts <- c(length(which(runs_foreward=="-" & runs_foreward_p>seuil)),length(which(runs_foreward=="-" & runs_foreward_p<=seuil)),length(which(runs_foreward=="+" & runs_foreward_p>seuil)),length(which(runs_foreward=="+" & runs_foreward_p<=seuil)))
  print(paste("foreward :", counts))
  labels_to_show <- labels[counts>0]
  pie(counts, labels_to_show, main="runs_foreward", col=c("#CCFFFF","#3399FF","#FFCC33","#FF6600")) 
  legend("bottom", legend = labels, horiz = TRUE, cex = 0.8, inset = -0.2)
  title(main = pathname, line = -1.5)
  # backward
  counts <- c(length(which(runs_backward=="-" & runs_backward_p>seuil)),length(which(runs_backward=="-" & runs_backward_p<=seuil)),length(which(runs_backward=="+" & runs_backward_p>seuil)),length(which(runs_backward=="+" & runs_backward_p<=seuil)))
  print(paste("backward :", counts))
  labels_to_show <- labels[counts>0]
  pie(counts, labels_to_show, main="runs_backward", col=c("#CCFFFF","#3399FF","#FFCC33","#FF6600"))
  #simultaneous
  counts <- c(length(which(runs_simult=="-" & runs_simult_p>seuil)),length(which(runs_simult=="-" & runs_simult_p<=seuil)),length(which(runs_simult=="+" & runs_simult_p>seuil)),length(which(runs_simult=="+" & runs_simult_p<=seuil)))
  print(paste("simultaneous :", counts))
  labels_to_show <- labels[counts>0]
  pie(counts, labels_to_show, main="runs_simult", col=c("#CCFFFF","#3399FF","#FFCC33","#FF6600"))
  # dev.off()
}
permutation_test_2groups <- function(x, y, func, nrep=1000){ # x et y mes deux distri, 
  rawdiff <- func(x) - func(y) # la moyenne de x - la moyenne d ey
  nx <- length(x)
  ny <- length(y)
  diffb <- NA # différence bootstrap
  for(i in 1:nrep){ # on le fait 1000 fois
    temp <- c(x,y)[sample(1:(nx+ny))] # on mélange avec un tirage aléatoire avec le max la taille initaile des 2 groupes
    xb <- temp[1:nx] # on prends les valeurs de 1 à nx
    yb <- temp[(nx+1):(nx+ny)] # on prend le reste ! 
    diffb[i] <- func(xb) - func(yb) # on fait la différence entre les nouveaux xb et yb 
  }
  return(list(pvalue_comp(rawdiff, diffb, tail="two"), diffb)) # on compare la diff originale avec la diff des groupes mélangés ! 
}
subgraphs_given_times_condition <- function (GRN, condition) {
  g <- graph_from_data_frame(GRN[[2]], directed=TRUE, vertices=GRN[[1]]) # je transforme les relations en graph, il est dirigé, et il a besoin des items
  glist <- decompose(g) # on décompe le graph en liste de R et il fait des sous-graph 
  glist <- glist[sapply(glist, gorder)>1] # on prend le sous-graph le plus grand et ainsi de suite ! et on vire les elements seuls
  all_components <- NA # on va garder NA == FALSE 
  
  for(i in 1:length(glist)){
    tempg <- glist[[i]] #graph temporaire
    v_edges <- get.edges(tempg, E(tempg)) #on récup-re toutes les branches du graph qu'on traite pr les avoir sous forme d'objet dataframe
    times_edges <- mapvalues(v_edges, V(tempg), vertex_attr(tempg, "times")) #on remplace dans un vecteur, on remplace les attribut par les temps
    edges_cond <- apply(times_edges, 1, condition) #on applique la fonction sur times_edges # ça dit si la relation entre les 2 noeuds à chauqe ligne respecte la condition donné
    edges_cond[is.na(edges_cond)] <- FALSE # ça dit ce qu'on fait des NA ! 
    # quel sont les entrées du vecteurs qui sont NA, on les met FALSE donc on arrête la relation

      if (sum(edges_cond)>0){
      sub <- subgraph.edges(tempg , E(tempg)[edges_cond]) # on extrait les sous graph qu'on traite dont les relations en filtrant les relations de manière à ce qu'elle valide la condition. 
      all_components <- c(all_components, components(sub)$csize) # et la je concatène pr chaque graph avec la taille de la composante de subnaf, on a une distri des longuuers de subgraph
    }else{all_components <- c(0, 0)}
  }
  all_components <- all_components[-1] # on vire les NA au fur et à mesure comme on les a ajouté au début
  return(all_components) # on sort une lsite de taille de subgrapph qui respecte la condition donné
}

#### ARGUMENTS ####
args <- commandArgs(trailingOnly = TRUE) # recupère les arguments
allinfosKEGG <- args[1] ## "allinfos-KEGG.csv"
paths <- args[2]        ## "paths/"

dat <- read.csv(allinfosKEGG, sep=";")
files_name <- dir(where) # récupère tous les fichiers d'un dossier

#### INITIALISATION ####
noms_colonnes <- c("path", "nbgene", "median", "mean", "sd", "min", "max", "correlation", "corr pvalue")
date_stats_paths <- as.data.frame(matrix(nrow = 0, ncol = length(noms_colonnes)))
nbgenes <- tabtotal <- NA
df_generankbirthpath <- as.data.frame(matrix(nrow = 0, ncol = 4))
tab_interact <- as.data.frame(matrix(nrow = 0, ncol = 5))
restable <- matrix(NA, length(files_name), 6)
resmeantable <- matrix(NA, length(files_name), 5)  
allpath <- c()


##### CALL FUNCTION #####
pdf("ggplot-generankbirthpath.pdf")
for(i in 1:length(files_name)){
  listxml <- xmlParse(paste(paths, files_name[[i]], sep=""))
  p <- xmlRoot(listxml) # le fichier xml interprété par R
  r <- recup_networks(fxml = p)
  GRN <- r[[1]]; pathname <- r[[2]]
  print(paste("###", i, pathname, "###", sep=" ")) 
  GRN <- check_ages(GRN)
  s <- subgraph_stats(GRN, pathname)
  statspath <- s[[1]]; corrpath <- s[[2]]; tabrank <- s[[3]]
  allpath <- append(allpath, pathname)
  
  nbgenes[i] <- nrow(GRN$items %>% filter(types=="gene"))
  newligne <- c(pathname, nbgenes[i], statspath$median, statspath$mean, statspath$sd, statspath$min, statspath$max, corrpath$estimate, corrpath$p.value)
  date_stats_paths <- rbind(date_stats_paths, newligne)

  df_generankbirthpath <- rbind(df_generankbirthpath, tabrank %>% mutate(path=pathname) %>% as.data.frame(.))

  ##### FOREWARD/BACKWARD et SIMULT ! #####
  relations_times <- GRN$relations
  relations_times[,1] <- mapvalues(GRN$relations[,1], GRN$items$ids, GRN$items$times, warn_missing=FALSE) # remplace GNRrelation (première colonne) et la correspondance temps
  relations_times[,2] <- mapvalues(GRN$relations[,2], GRN$items$ids, GRN$items$times, warn_missing=FALSE)
  rescomp <- compare(relations_times) #on compare les relations
  t <- data.frame("foreward" = table(rescomp)["foreward"],
                  "simult" = table(rescomp)["simult"],
                  "backward" = table(rescomp)["backward"])
  t <- replace(t, is.na(t), 0)

  runs_foreward <- runs_backward <- runs_simult <- NA
  nrep = 100
  for(j in 1:nrep) { #on peut augmenter le nombre de run, pr solidifier le test
    relations_times_boot <- relations_times # on copie
    relations_times_boot[,2] <- relations_times_boot[sample(1:nrow(relations_times_boot)),2] # on mélange que la deuxième colonne, tirage sans remise !

    res_comp_boot <- compare(relations_times_boot) # on compare les relations

    d <- data.frame("foreward"=table(res_comp_boot)["foreward"],
                    "simult"=table(res_comp_boot)["simult"],
                    "backward"=table(res_comp_boot)["backward"])
    d <- replace(d, is.na(d), 0)
    # on compte les foreward du test (on a 10 valeurs pr chacun parce que 10 fois le test en l'état)
    runs_foreward[j] <- d$foreward ; runs_simult[j] <- d$simult ; runs_backward[j] <- d$backward
  }
  pvalues_foreward <- pvalue_comp(t$foreward, runs_foreward, tail="two");  restable[i,c(1,2)] <- unlist(pvalues_foreward)
  pvalues_simult <- pvalue_comp(t$simult, runs_simult, tail="two");  restable[i,c(3,4)] <- unlist(pvalues_simult)
  pvalues_backward <- pvalue_comp(t$backward, runs_backward, tail="two");  restable[i,c(5,6)] <- unlist(pvalues_backward)
  
  
  ##### PERMUTATION #####
  permutation(GRN, pathname, nrep=1000)

  ##### PLOT PAR VOIE ######
  ## figure pour afficher chaque gène en fonction de leur rank et de leur quantité à ce rang
  tab_generank <- tabrank %>% group_by(id, rank, birth) %>% dplyr::summarise(n=n(), min_rank=min(rank), max_rank=max(rank)) %>% arrange(min_rank, max_rank, -n, id) ### e trié comme on le veut
  plot_generank <- tab_generank %>% mutate(id=factor(id, levels=id %>% unique())) %>% ggplot(aes(x= reorder(id, max_rank), y=rank, size=n, col=factor(birth, levels=1:25))) + geom_point() +
    theme(axis.text.x = element_text(angle=90)) + ggtitle(pathname) + labs(x = "id", col= "birth", size = "n") +
    scale_color_manual(values = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 25), drop = FALSE)
  print(plot_generank) ## le plot du pdf "ggplot-generankbirthpath"


  t <- GRN$relations %>% left_join(GRN$namesfromnetwork, by = c("entry1" = "vectnames")) %>%
    left_join(GRN$namesfromnetwork, by = c("entry2" = "vectnames")) %>%
    left_join(dat %>% select(kegg_id, gene_name, num_clade), by = c("allnames.x" = "kegg_id")) %>%
    left_join(dat %>% select(num_clade, kegg_id, gene_name), by = c("allnames.y" = "kegg_id")) %>%
    select(-entry1, -entry2) %>% mutate("path" = pathname)
  tab_interact <- rbind(tab_interact, t) %>% as.data.frame(.)

  ###### TABLEAU RECAP #######
  tabrecap <- GRN$relations
  for(i in 1:nrow(tabrecap)){
    id1 = as.character(tabrecap$entry1[i])
    id2 = as.character(tabrecap$entry2[i])
    tabrecap$hsa1[i] <- GRN$items$names[GRN$items$ids==id1]
    tabrecap$hsa2[i] <- GRN$items$names[GRN$items$ids==id2]
    tabrecap$name1[i] <- ifelse(length(dat$gene_name[dat$kegg_id==strsplit(GRN$items$names[GRN$items$ids==id1], " ")[[1]][1]])>0,
                                dat$gene_name[dat$kegg_id==strsplit(GRN$items$names[GRN$items$ids==id1], " ")[[1]][1]], "")
    tabrecap$name2[i] <- ifelse(length(dat$gene_name[dat$kegg_id==strsplit(GRN$items$names[GRN$items$ids==id2], " ")[[1]][1]])>0,
                                dat$gene_name[dat$kegg_id==strsplit(GRN$items$names[GRN$items$ids==id2], " ")[[1]][1]], "")
    tabrecap$birth1[i] <- GRN$items$times[GRN$items$ids==id1]
    tabrecap$birth2[i] <- GRN$items$times[GRN$items$ids==id2]
    tabrecap$direction[i] <- ifelse(tabrecap$birth1[i]>tabrecap$birth2[i], "Foreward",
                              ifelse(tabrecap$birth1[i]<tabrecap$birth2[i], "Backward",
                                     ifelse(tabrecap$birth1[i]==tabrecap$birth2[i], "Simult", "???")))
    tabrecap$path[i] <- pathname
  }
  tabtotal <- rbind(tabtotal, tabrecap)
}
dev.off()

## les tableaux à sortir
names(date_stats_paths) <- c("path", "nbgene", "median", "mean", "sd", "min", "max", "correlation", "pvalue")
# write.csv(date_stats_paths, file = "correlation-total-new.csv", row.names = FALSE)
# write.csv(df_generankbirthpath, file = "gene-birth-rank-new.csv", row.names = FALSE)
names(tab_interact) <- c("kegg1", "kegg2", "namegene1", "birth1", "birth2", "namegene2", "path")
# write.csv(tab_interact, file = "interaction2by2.csv", row.names = FALSE)
# write.csv(tabtotal, "tableau-interaction-total.csv", row.names = FALSE)

### histogramme pour les comptes for/back/sim
# forhist <- apply(restable[,c(1,2)], 1, function(x) {as.numeric(paste(x[1], x[2], sep=""))}) # permet de récupérer les foreward et on rajoute + ou - à la pvalue
# simhist <- apply(restable[,c(3,4)], 1, function(x) {as.numeric(paste(x[1], x[2], sep=""))})
# bachist <- apply(restable[,c(5,6)], 1, function(x) {as.numeric(paste(x[1], x[2], sep=""))})
# pdf("histograms_double.pdf", width=10, height=7)
# par(mfrow=c(2,3))
# hist(forhist[forhist>=0], col="tomato3"); hist(simhist[simhist>=0], col="tomato3"); hist(bachist[bachist>=0], col="tomato3")
# hist(forhist[forhist<=0], col="steelblue1"); hist(simhist[simhist<=0], col="steelblue1"); hist(bachist[bachist<=0], col="steelblue1")
# dev.off()


##### DISTRIBUTION #####
#### distribution foreward/backward/simult
abc <- na.omit(tabtotal)
# abc <- abc %>% filter(is.finite(birth1) & is.finite(birth2)) %>% select(-c(entry1, entry2)) %>% unique() %>% pivot_wider(., names_from = path, values_from = direction)

tabuniqinteractpath <- abc %>% mutate(hsa1 = sapply(strsplit(as.character(hsa1), " "), "[[", 1), 
                                      hsa2 = sapply(strsplit(as.character(hsa2), " "), "[[", 1)) %>% 
  group_by(path) %>% select(hsa1, hsa2, birth1, birth2, path, direction) %>% unique()

t <- tabuniqinteractpath %>% count(direction) %>% ungroup()
pdf("ggplot-distributiondirection.pdf")
plot_distrifbs <- ggplot(t, aes(x=n, y=path, fill=direction)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("#66C2A5","#E6F598", "#FDAE61")) + theme_minimal() + #c("#66C2A5","#E6F598", "#FDAE61", "#D53E4F","#CCCCCC", "#F78E87","#0489B1", "#BDBDBD")
  ggtitle("Distribution of direction of interaction by pathway") + 
  labs(x = "Number of interactions", y = "Pathway") 
print(plot_distrifbs)
dev.off()
# t <- tabuniqinteractpath %>% count(direction) %>% pivot_wider(names_from = direction, values_from = n) # nombre d'interaction B/F/S par pathway
# t <- t %>% ungroup() %>% summarise(B=sum(Backward), F=sum(Foreward), S=sum(Simult)) ## sum total du nombre d'interaction B/F/S (toute pathway confondue)

#### distribution par voie / naissance
data_long <- dat %>% pivot_longer(cols = 10:last_col(), names_to = "path", values_to = "value") %>% filter(value == 1) %>%
  select(kegg_id, gene_name, birth_clade, num_clade, path)
# clade_mapping <- data.frame(birth_clade = unique(data_long$birth_clade),num_clade = unique(data_long$num_clade)) %>% drop_na()
# data_long %>% group_by(path, num_clade) %>% ggplot(., aes(x=num_clade, fill=path)) + geom_bar(position="fill") +
#   labs(x = "birth clade", y = "percentage of gene") + ggtitle("Distribution of birth gene")


## ici on compte plusieurs fois les mêmes gènes
data_long %>% group_by(path, num_clade) %>% drop_na() %>% ggplot(., aes(x=factor(num_clade, levels=1:25), fill=path)) + geom_bar(position="stack") +
  labs(x = "birth clade", y = "number of gene") + ggtitle("Distribution of birth gene") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) + 
  scale_x_discrete(labels = c(
    "1" = "Yeast", "2" = "Metazoa", "3" = "Eumetazoa", 
    "4" = "Bilateria", "5" = "Cnidaria", "6" = "Protostomia",
    "7" = "Lophotrochozoa", "8" = "Ecdysozoa", "10" = "Ciona", 
    "11" = "Cyclostomata", "12" = "Gnathostomata", "13" = "Actinopterygii", "14" = "Neopterygii", 
    "15" = "Osteoglossiformes", "16" = "Clupeocephala", 
    "17" = "Sarcopterygii", "18" = "Anura", "19" = "Lepidosauria", 
    "20" = "Archosauria", "21" = "Cryptodira", "22" = "Aves", 
    "23" = "Mammalia", "24" = "Marsupialia", "25" = "Eutheria"
  ))

## ici on compte qu'une seule fois chaque gène
pdf("ggplot-distributionbirth.pdf")
plot_distripathbirth <- dat %>% filter(num_clade!="") %>% ggplot(., aes(x=factor(num_clade, levels=1:25))) + geom_bar(position="stack") +
  labs(x = "birth clade", y = "number of gene") + ggtitle("Distribution of birth gene") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) + 
  scale_x_discrete(labels = c(
    "1" = "Yeast", "2" = "Metazoa", "3" = "Eumetazoa", 
    "4" = "Bilateria", "5" = "Cnidaria", "6" = "Protostomia",
    "7" = "Lophotrochozoa", "8" = "Ecdysozoa", "10" = "Ciona", 
    "11" = "Cyclostomata", "12" = "Gnathostomata", "13" = "Actinopterygii", "14" = "Neopterygii", 
    "15" = "Osteoglossiformes", "16" = "Clupeocephala", 
    "17" = "Sarcopterygii", "18" = "Anura", "19" = "Lepidosauria", 
    "20" = "Archosauria", "21" = "Cryptodira", "22" = "Aves", 
    "23" = "Mammalia", "24" = "Marsupialia", "25" = "Eutheria"
  ))
print(plot_distripathbirth)
dev.off()


###### DELTA #####
## toutes les interactions et leurs voies et le delta
tab_interactpath <- tab_interact %>% select(-kegg1, -kegg2) %>% group_by(path) %>% drop_na() %>% unique() %>% 
  mutate(delta = birth1-birth2, here=1) %>% pivot_wider(., names_from = path, values_from = here)
#write.csv(tab_interactpath, file = "interactionpropre.csv", row.names=F)

tab_nbdeltapath <- tab_interact %>% select(-kegg1, -kegg2) %>% group_by(path) %>% drop_na() %>% unique() %>% mutate(delta = birth1-birth2) %>% 
  group_by(birth1, birth2, path) %>% count(delta) %>% pivot_wider(., names_from = path, values_from = n)
#write.csv(tab_nbdeltapath, file = "nbdelta.csv", row.names=F)

tab_delta <- tab_interact %>% select(-kegg1, -kegg2) %>% group_by(path) %>% drop_na() %>% unique() %>% 
  mutate(delta = birth1-birth2) %>% group_by(birth1, birth2, path) %>% count(delta) 
med <- tab_delta %>% group_by(path) %>% summarise(median = median(delta))
med <- med %>% arrange(median) %>% pull(path)
tab_delta <- tab_delta %>% mutate(path = factor(path, levels = med))
pdf("ggplot-distributiondelta.pdf")
plot_delta <- ggplot(tab_delta, aes(x = delta, y = path)) +
  geom_boxplot() + labs(x = "Delta", y = "Path") +
  ggtitle("Distribution des deltas par path")
print(plot_delta) ## figure 6 !
dev.off()
# permet de récupérer le tableau de valeur
#ld <- layer_data(plot_delta)
#ld$path <- as.data.frame(date_stats_paths) %>% select(path) %>% arrange(path)
#ld <- ld %>% select(path, xmin, xlower, xmiddle, xupper, xmax) %>% arrange(desc(path))


#### FIGURES DIVERSES ####
data <- df_generankbirthpath

### avec petite distribution ## meilleure distri 
date_stats_paths$correlation <- as.numeric(date_stats_paths$correlation)
correl <- date_stats_paths %>% select(path, correlation) %>% arrange(correlation, path) 
data$path <- factor(data$path, levels=correl$path)
pdf("ggplot-distributionbirthrank.pdf")
plot_distrirankbirth <- ggplot(data, aes(y = "",  fill = factor(birth))) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 25), drop = FALSE) +
  coord_cartesian() +
  theme_void() +
  facet_grid(path ~ rank, switch = "both") +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0)) +
  labs(x = "Path", y = "Rank", fill = "Birth") +
  ggtitle("Distribution of birth moments by path and rank")
print(plot_distrirankbirth)
dev.off()

## fonction pour inverser les rank pour pouvoir comparer la fin de la voie ensemble
invert_ranks <- function(data) {
  max_rank <- max(data$rank)
  inverted_ranks <- max_rank - data$rank + 1
  data$rank <- inverted_ranks
  return(data)
}
inverted_df <- data %>% group_by(path) %>% do(invert_ranks(.)) %>% ungroup()
## à mettre dans plot_distrirankbirth après 


### distribution des rank*birth
plot_intersectrankbirth <- data %>% group_by(rank, birth) %>% unique() %>% summarise(n = n()) %>% 
  ggplot(., aes(x = factor(rank), y = factor(birth))) +
  geom_count(aes(size = n), color = "darkred") +
  scale_size_continuous(range = c(3, 10)) +  # Ajuster la taille des points
  xlab("Rank") +
  ylab("Birth") +
  ggtitle("Number of genes by intersection Birth-Rank")
# print(plot_intersectrankbirth)

### distribution des delta par voie 
pdf("ggplot-distridelta.pdf")
for(p in allpath){
  plot_distridelta <- tab_interact %>% select(-kegg1, -kegg2) %>% group_by(path) %>% drop_na() %>% unique() %>%
    mutate(delta = birth1 - birth2, birth_label = ifelse(birth1 < birth2, paste(birth1, ":", birth2), 
                                                             paste(birth2, ":", birth1))) %>% 
    group_by(delta, birth_label, birth1, birth2) %>%
    filter(path==p) %>%
    summarize(freq = n()) %>%
    ggplot(., aes(delta, freq, fill = factor(birth1, levels=1:25))) +
    geom_bar(position = "stack", stat = "identity") +
    theme_minimal() +
    scale_fill_manual(values = paletteer_c("ggthemes::Sunset-Sunrise Diverging", 25), drop = FALSE) + 
    labs(x = "Delta", y = "Proportion of birth", fill = "Birth") +
    ggtitle(paste("Distribution of birth moments by rank for", p))
  print(plot_distridelta)
}
dev.off()
  


### % sur les back/for/simult
# tab_interact %>% select(-kegg1, -kegg2) %>% group_by(path) %>% drop_na() %>% unique() %>% ungroup() %>% select(-path) %>% 
#   mutate(delta = birth1 - birth2, type = ifelse(delta<0, "backward", ifelse(delta==0, "simultaneous", "forward"))) %>%
#   unique() %>% group_by(type) # %>% summarise(n = n())


## Sélection aléatoire de lignes
# nombre_de_lignes <- 1000  # Spécifiez ici le nombre de lignes souhaité dans le tableau résultant
# tableau_resultat <- data[sample(nrow(data), nombre_de_lignes), ]

