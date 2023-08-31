# Determining the time of birth of human genes

# STEP 1
## KEGG Pathway Information Retrieval

This script, `get-info-KEGG.py`, allows you to retrieve information about signaling pathways from the KEGG database based on specified keywords. It fetches data such as pathway names, gene IDs, UniProt IDs, Ensembl IDs, and gene names. The fetched data is stored in a CSV file named `infos-KEGG.csv`. Additionally, the script generates a text file named `pathwaylist.txt`, containing the Ensembl IDs of genes associated with the selected pathway.

### Prerequisites
- Python 3.x
- Packages: requests, beautifulsoup4

You can install the required packages using the following command: `pip install requests beautifulsoup4`

### Usage
1. (optional) Open `get-info-KEGG.py` in a text editor. Modify the `keywords` variable on line 25 to specify your desired keywords.
*In our case, we use keywords like `"signaling pathway"` or `"ovarian"`. And exclude `"multiples species"`.*
3. Run the script: `python get-info-KEGG.py`
4. After the script finishes running, you will find two output files in the same directory:
- `infos-KEGG.csv`: Contains pathway information (Pathway Name, Gene ID, UniProt ID, Ensembl ID, Gene Name).
- `pathwaylist.txt`: Contains Ensembl IDs of genes associated with the selected pathway.
- `\path\xxxx.xml`: Contains all pathway in xml format.

# STEP 2
## Gene Birth Moment Determination

The script `get-bith.py` helps determine the appearance timing of genes using Genomicus trees and orthology relationships. This second program requires several mandatory input files as arguments:

- `-t`: Vertabrates tree (NHX file) -> Download from: [vertebrates tree NHX](https://ftp.bio.ens.psl.eu/pub/dyogen/genomicus/109.01/protein_tree.nhx.bz2) (requires decompression)
- `-m`: Metazoa tree (NHX file) -> Download from: [metazoa tree NHX](https://ftp.bio.ens.psl.eu/pub/dyogen/genomicus-metazoa/51.01/protein_tree.nhx.bz2) (requires decompression)
- `-l`: Interest gene list (TXT file) or the "pathwaylist.txt" obtained from the first program.
- `-c`: Clades of interest species (CSV file) with [species, clade, clade number] information.
- Optional: `-p`: Paralogue of interest gene list (TXT file)

The second program generates an output file named `birth-moment.csv` with columns [Ensembl ID, Clade, Clade Number].

### Prerequisites
- Python 3.x
- Packages: Biopython

You can install the required package using the following command: `pip install biopython`

### Usage
1. Run the script with required and optional arguments: `python get-bith.py -t vertebrates_tree.nhx -m metazoa_tree.nhx -l pathwaylist.txt -c clades.csv`
Replace the argument values with the actual paths to your input files.
2. After the script finishes running, you will find the output file `birth-moment.csv` containing gene birth moment information.

# STEP 3
## Gene Information Fusion

The script `fusion-infos.py` combines information from the two previous programs and generates an output file named `allinfos-KEGG.csv` with columns:

- ["kegg id", "entrezgene id", "uniprot id", "ensembl id", "gene name", "gene other name", "gene description", "birth clade", "num clade", "list pathway"]

The `list pathway` column represents as many columns as the pathways obtained from the keywords in program 1. If the gene is present in a pathway, it will be marked with a "1" in the corresponding cell.

### Prerequisites
- Python 3.x

### Usage
1. Run the script with required arguments: `python fusion-infos.py -i infos-KEGG.csv -b birth-moment.csv`
Replace the argument values with the actual paths to your input files.
2. After the script finishes running, you will find the output file `allinfos-KEGG.csv` containing merged gene information.
