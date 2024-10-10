#! /usr/bin/env python3
# author : Picolo Floriane

from pickle import TRUE
from Bio import Phylo
import csv
from collections import defaultdict
import re       # pour utiliser le multiple splicing
import argparse

def parser():
    parser = argparse.ArgumentParser(description = "BIRTH MOMENT OF HUMAN GENE")
    parser.add_argument("-t", "--tree", metavar = "PATH", type = str, required = True, help = "vertebrates tree (nhx file)")
    parser.add_argument("-m", "--mtree", metavar = "PATH", type = str, required = True, help = "metazoa tree (nhx file)")
    parser.add_argument("-l", "--list", metavar = "PATH", type = str, required = True, help = "interest gene list (txt file)")
    parser.add_argument("-c", "--clade", metavar = "PATH", type = str, required = True, help = "clades of species (csv file)")
    parser.add_argument("-p", "--paralogue", metavar = "PATH", type = str, required = False, help = "paralogue of interest gene list (txt file)")
    return parser

def recuperation_listgene(infile):
    with open(infile) as in_file:
        liste = []
        for row in in_file:
            liste.append(row.strip()) # add line without "\n"
    liste = list(set(liste)) # delete double variable 
    return liste

def recuperation_clade(infile):
    ## le fichier est construit de telle façon : espèces, clade, numclade
    with open(infile) as in_file:
        dico = defaultdict(list)
        for row in in_file:
            r = re.split(r',|\t|;', row.strip())
            dico[r[2]].append(r[0])
    return dico

def recuperation_paralogues(infile):
    with open(infile) as in_file:
        dico = defaultdict(list)
        for row in in_file.readlines()[1:]:
            r = row.strip().split(',')
            if bool(r[1]): dico[r[0]].append(r[1]) # vérifie que la liste n'est pas vide
    return dico

def explore_ensembl_tree(ensembltreefile, interestlist=None):
    dict_Ensembl = dict()
    dict_pivot = dict()
    dict_pivot_inv = defaultdict(list)

    trees = Phylo.parse(ensembltreefile, "newick") # récupération de tous les arbres phylo 
    for tree in trees: 
        pivot = str()
        list_sp = list()

        term_names = [term.name for term in tree.get_terminals()]

        if interestlist is None:
            interestlist = [element for element in term_names if element.startswith("ENSG")]

        if any(element in interestlist for element in term_names):  # vérification qu'un gène d'intérêt est dans les arbres.
            human_interest = set(interestlist) & set(term_names)

            term_comment = [term.comment for term in tree.get_terminals()]
            for num, comment in enumerate(term_comment):
                sp = comment[comment.index(":S=")+3:comment.index(":F=")]

                if sp == "Drosophila.melanogaster" or sp == "Caenorhabditis.elegans":
                    if not(pivot): pivot = term_names[num]
                else:
                    list_sp.append(sp) 

            for human in human_interest:
                dict_Ensembl[human] = list_sp
                dict_pivot[human] = pivot
    
    # permet d'avoir qu'un identifiant metazoa pour plusieurs gènes humain
    {dict_pivot_inv[v].append(k) for k, v in dict_pivot.items()}
    dict_pivot = dict(dict_pivot_inv)

    return dict_Ensembl, dict_pivot

def explore_metazoa_tree(metazoatreefile, dict_Ensembl, interestdict):
    trees = Phylo.parse(metazoatreefile, "newick")
    
    for tree in trees:
        list_sp = list()

        term_names = [term.name for term in tree.get_terminals()]
        if any(element in interestdict.keys() for element in term_names):
            metazoa_interest = set(interestdict.keys()) & set(term_names)
            term_comment = [term.comment for term in tree.get_terminals()]
            for comment in term_comment:
                sp = comment[comment.index(":S=")+3:comment.index(":D=")]
                list_sp.append(sp)
            for metazoa in metazoa_interest:
                for human in interestdict[metazoa]:
                    dict_Ensembl[human] += list_sp

    return dict_Ensembl

def check_clade(clades, list_sp):
    for nameclade, spinclade in clades.items():
        for sp in list_sp:
            if sp in spinclade:
                return nameclade, list(clades.keys()).index(nameclade)
    return "NA","NA"

def determination_clade_with_para(dict_Ensembl, interestlist, paralogues, clades):
    dict_clade = dict()
    for humangene in interestlist: # les gènes de nos listes d'interet
        if humangene in dict_Ensembl.keys():
            clade, numclade = check_clade(clades, dict_Ensembl[humangene])
        else: 
            clade = ""
            numclade = 9999
        if humangene in paralogues.keys():
            for paragene in paralogues:
                if paragene in dict_Ensembl.keys(): 
                    c, nc = check_clade(clades, dict_Ensembl[paragene])
                    if nc < numclade: 
                        numclade = nc
                        clade = c
        dict_clade[humangene] = [clade, numclade+1]
    return dict_clade

def determination_clade(dict_Ensembl, interestlist, clades):
    dict_clade = dict()
    for humangene in interestlist: # les gènes de nos listes d'interet
        if humangene in dict_Ensembl.keys():
            clade, numclade = check_clade(clades, dict_Ensembl[humangene])
        else: 
            clade = ""
            numclade = 9999
        dict_clade[humangene] = [clade, numclade+1]
    return dict_clade

def writecsv_nb_ortho(dict_Ensembl, outfile, clades):
    spamwriter = csv.writer(open(outfile, "w"), delimiter=';', quoting=csv.QUOTE_NONE, quotechar='"', escapechar='')
    list_species = list()
    for values in clades.values(): 
        for specie in values : 
            list_species.append(specie)
    
    for k, v in dict_Ensembl.items():
        spamwriter.writerow([k, v]) 

    spamwriter.writerow(["ensembl_id"] + list_species)
        
    for k, v in dict_Ensembl.items():
        ligne = [k]
        for specie in list_species:
            ligne.append(v.count(specie))

        spamwriter.writerow(ligne)

def writecsv_birth(dict_clades, outfile):
    spamwriter = csv.writer(open(outfile, "w"), delimiter=';', quoting=csv.QUOTE_NONE, quotechar='"', escapechar='')
    spamwriter.writerow(["ensembl_id", "birthclade"])
    for k, v in dict_clades.items(): 
        spamwriter.writerow([k, v])


if __name__ == "__main__":
    
    parser = parser() 
    args = parser.parse_args()

    # args
    f_nhxVertebrates = args.tree        # f_nhxEnsembl = "trees/protein_tree_v109.nhx"
    f_nhxMetazoa = args.mtree           # f_nhxMetazoa = "trees/protein_tree_v51.nhx"
    f_interest = args.list              # f_interestlist = "pathwaylist.txt"
    f_clade_sp = args.clade             # f_clade_sp = "clades-species.csv"
    f_paralogue = args.paralogues       # f_paralogue = "paralogues-v109.txt"

    # of_nbortho = "out_nbortho.csv"
    of_birth = "out_originmoment.csv"

    ## on récupère les différentes listes
    l_interest = recuperation_listgene(infile=f_interest)
    d_clade = recuperation_clade(infile=f_clade_sp)
    d_paralogues = recuperation_paralogues(infile=f_paralogue)    
    
    ## on recherche une fois pour tout le génome humain
    d_Ensembl, d_pivot = explore_ensembl_tree(ensembltreefile=f_nhxVertebrates, interestlist=l_interest)
    d_Ensembl = explore_metazoa_tree(metazoatreefile=f_nhxMetazoa, dict_Ensembl=d_Ensembl, interestdict=d_pivot)

    ## ici on récupère les informations pour nos listes d'intérêts 
    d_clade_pathway = determination_clade_with_para(dict_Ensembl=d_Ensembl, interestlist=l_interest, paralogues=d_paralogues, clades=d_clade)
    d_clade_pathway = determination_clade(dict_Ensembl=d_Ensembl, interestlist=l_interest, clades=d_clade)

    ## ici on écrit les résultats dans des fichiers csv
    # writecsv_nb_ortho(dict_Ensembl=d_Ensembl, outfile=of_nbortho, clades=d_clade)
    writecsv_birth(dict_clades=d_clade_pathway, outfile=of_birth)


## pour avoir les couleurs sur le graph KEGG
# https://www.kegg.jp/kegg-bin/show_pathway?map=hsa00510&multi_query=79868+pink%0d%0a56052+yellow,blue
