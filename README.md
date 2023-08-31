# Stript to determine the birth moment of human genes

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

### Example
For instance, if you set the `keywords` to `"signaling pathway"`, the script will retrieve information related to signaling pathways. The resulting CSV file will contain columns like `Pathway Name`, `Gene ID`, `UniProt ID`, `Ensembl ID`, and `Gene Name`.





---
Feel free to customize this README further to fit your project structure and any additional information you'd like to provide. Make sure to include any relevant licensing information, contributors, or contact details as needed.
