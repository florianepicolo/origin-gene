#! /usr/bin/env python3
# author : Picolo Floriane

from bs4 import BeautifulSoup
import requests
import shutil
import csv

## ce programme permet de récupérer pour des mots clés donnés les fichiers html des voies
## ainsi qu'un tableau répertoriant les informations disponibles sur KEGG ! 
## les mots clés sont : signaling pathway et Ovarian en excluant "multiple species"

def extract_id_pathway():
    """
        "http://rest.kegg.jp/list/pathway/hsa"
        path:hsa00010	Glycolysis / Gluconeogenesis - Homo sapiens (human)
        path:hsa00030	Pentose phosphate pathway - Homo sapiens (human)
        On récupère la liste de tous les id des pathways comprenant notre keyword
    """
    pathway_list = []
    response = requests.get("http://rest.kegg.jp/list/pathway/hsa")
    soup = BeautifulSoup(response.content, 'html.parser')
    page = soup.get_text().strip().split("\n")
    for element in page:
        if ("signaling pathway" in element.split("\t")[1] or "Ovarian" in element.split("\t")[1]) and "multiple species" not in element.split("\t")[1]:
            pathway_list.append(element.split("\t")[0])
    return pathway_list

def extract_name_pathway(pathwayid):
    """
        "http://rest.kegg.jp/list/hsa04668" 
        path:hsa04668	TNF signaling pathway - Homo sapiens (human)
        On récupère la deuxième colonne de la ligne (nom du pathway)
    """
    response = requests.get("http://rest.kegg.jp/list/" + pathwayid)
    soup = BeautifulSoup(response.content, 'html.parser')
    return soup.get_text().split("\t")[1].strip()

def extract_idgene(pathwayid):
    """
        "http://rest.kegg.jp/link/hsa/hsa04668"
        path:hsa04668	hsa:10000
        path:hsa04668	hsa:10059
        On récupère toutes les deuxièmes colonnes des lignes (id du gène)
    """
    response = requests.get("http://rest.kegg.jp/link/hsa/" + pathwayid)
    soup = BeautifulSoup(response.content, 'html.parser')
    page = soup.get_text().strip().split("\n")
    listidgene = []
    for element in page: 
        listidgene.append(element.split("\t")[1])
    return listidgene

def extract_namegene(geneid):
    """
        "http://rest.kegg.jp/find/genes/hsa:10000"
        hsa:100008587	RNA5-8SN5, RN5-8S1, RNA5-8S5; RNA, 5.8S ribosomal N5
        hsa:10000	AKT3, MPPH, MPPH2, PKB-GAMMA, PKBG, PRKBG, RAC-PK-gamma, RAC-gamma, STK-2; AKT serine/threonine kinase 3
        On récupère la deuxième colonne pour la ligne correspondante à id donné en paramètre
    """
    response = requests.get("http://rest.kegg.jp/find/genes/" + geneid)
    soup = BeautifulSoup(response.content, 'html.parser')
    page = soup.get_text().strip().split("\n")
    for element in page:
        if element.split("\t")[0] == geneid:
            return element.split("\t")[1]

def extract_uniprotID(geneid):
    """
        "http://rest.kegg.jp/conv/uniprot/hsa:10458"
        hsa:10458	up:Q9UQB8
        On récupère la deuxième colonne
    """
    response = requests.get("http://rest.kegg.jp/conv/uniprot/" + geneid)
    soup = BeautifulSoup(response.content, 'html.parser')
    page = soup.get_text().strip().split("\n")
    for element in page:
        if element.split("\t")[0] == geneid:
            return element.split("\t")[1].split(":")[1]
    return ""

def download_pathway(pathwayid):
    """
        "http://rest.kegg.jp/get/hsa03320/kgml"
        On récupère l'image html de la pathway
    """
    path_img = requests.get('http://rest.kegg.jp/get/' + pathwayid + '/kgml', stream = True)
    path_name = idpathway.replace(':', '_')
    with open("paths/" + path_name + '.xml', 'wb') as f:
        shutil.copyfileobj(path_img.raw, f)

def extract_ensembl_id(gene_id):
    """
        "https://www.kegg.jp/entry/hsa03320"
        On récupère le numéro Ensembl s'il y en a un
    """
    response = requests.get("https://www.kegg.jp/entry/"+gene_id)
    soup = BeautifulSoup(response.content, 'html.parser')
    for link in soup.find_all('a'):
        text = link.get_text('href')
        if "ENS" in text:
            return text
    return ""

if __name__ == "__main__":   
    dico_info_gene = dict()
    list_ensemblid = list()
    with open("infos-KEGG.csv", "w", newline='') as csvinfosfile:
        spamwriter = csv.writer(csvinfosfile, delimiter=';', quoting=csv.QUOTE_NONE, quotechar='', escapechar=' ')
        spamwriter.writerow(["name_pathway", "id_gene", "uniprot_gene", "ensembl_id", "name_gene", "descr_gene"])
        all_large_pathway = extract_id_pathway()
        for idpathway in all_large_pathway:
            namepathway = extract_name_pathway(pathwayid=idpathway[5:])
            listgene = extract_idgene(pathwayid=idpathway[5:])
            download_pathway(idpathway)
            for idgene in listgene:
                if idgene in dico_info_gene.keys():
                    spamwriter.writerow([namepathway, idgene, dico_info_gene[idgene]])
                else:
                    uniprotID = extract_uniprotID(geneid=idgene)
                    ensemblid = extract_ensembl_id(gene_id=idgene)
                    namegene = extract_namegene(geneid=idgene)
                    dico_info_gene[idgene] = [uniprotID, ensemblid, namegene]
                    spamwriter.writerow([namepathway, idgene, uniprotID, ensemblid, namegene])
                list_ensemblid.append(ensemblid)
    with open("pathwaylist.txt", "w") as txtpathfile:
        for ensemblid in list_ensemblid:
            txtpathfile.write(ensemblid + "\n")

