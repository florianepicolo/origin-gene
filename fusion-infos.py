#! /usr/bin/env python3
# author : Picolo Floriane

import csv
import argparse

def parser():
    parser = argparse.ArgumentParser(description = "FUSION OF INFORMATIONS")
    parser.add_argument("-i", "--infos", metavar = "PATH", type = str, required = True, help = "KEGG informations (csv file)")
    parser.add_argument("-b", "--birth", metavar = "PATH", type = str, required = True, help = "birth moment (csv file)")
    return parser

def get_info(filename):
    with open(filename, newline='') as in_file:

        lines_wout_first_line = in_file.readlines()[1:]

        dico = dict()
        list_pathway = list()

        for row in lines_wout_first_line:
            lign = row.strip().split(";")
            ## pathway name
            if "signaling" in lign[0]:
                pathname = " ".join(lign[0].split()).split(" signaling")[0]
            else: pathname = " ".join(lign[0].split()).split(" -")[0]

            if pathname not in list_pathway:
                list_pathway.append(pathname)

            ## id 
            hsa = lign[1]
            entrezgene = hsa[4:]
            uniprot = lign[2]
            ensembl = lign[3]

            ## infos 
            genenames = lign[4]
            genename = genenames.split(",")[0]
            genedesc = lign[5]
            
            if ensembl not in dico.keys(): 
                dico[ensembl] = [list(),list()]
                
                dico[ensembl][0] = [hsa, entrezgene, uniprot, ensembl, genename, genenames, genedesc]
                dico[ensembl][1] = [pathname] # on ajoute les noms des pathways dans lesquels le gène est présent
            else: 
                dico[ensembl][1] += [pathname]

    return dico, list_pathway

def get_clade(filename, dico):
    with open(filename) as in_file: 

        lines = in_file.readlines()[1:]

        for row in lines:
            lign = row.strip().split(",")
            ensembl = lign[0]
            birth = lign[1]
            clade = lign[2]

            dico[ensembl][0] += [birth, clade]

    return dico

def write_csv(out_file, dico, listpath):
    spamwriter = csv.writer(open(out_file, "w"), delimiter=';', quoting=csv.QUOTE_NONE, quotechar='"', escapechar='')
    fields = ["kegg_id", "entrezgene_id", "uniprot_id", "ensembl_id", "gene_name", "gene_other_name", "gene_description", "birth_clade", "num_clade"]
    [fields.append(path) for path in listpath]   

    spamwriter.writerow(fields)

    for genename, infos_pathways in dico.items():
        lign = list()
        infos = infos_pathways[0] #  [hsa, entrezgene, uniprot, ensembl, genename, genedesc, birth, clade]
        pathways = infos_pathways[1] # liste des pathways
        
        [lign.append(info) for info in infos]

        for pathway in listpath: 
            if pathway in pathways:
                lign.append(1)
            else: 
                lign.append("")
        spamwriter.writerow(lign)


if __name__ == "__main__":   

    parser = parser() 
    args = parser.parse_args()

    # args
    f_infos = args.infos        # f_infos = "infos-KEGG.csv"
    f_birth = args.birth        # f_birth = "birth-moment.csv"    
    
    of_allinfos = "allinfos-KEGG.csv"

    dict_info, path_list = get_info(filename=f_infos)
    dict_info = get_clade(filename=f_birth, dico=dict_info)
    write_csv(out_file=of_allinfos, dico=dict_info, listpath=path_list)

