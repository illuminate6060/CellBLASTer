# CellBLASTer
A universal plant scRNA-seq annotation tool inspired by cellular BLAST strategies.
CellBlaster is a cross-species cell type identification and annotation tool designed specifically for plant single-cell transcriptome (scRNA-seq). Through cross-species orthogroup (OG) mapping, symbolic percentage encoding, and multi-round correction algorithms, it accurately maps the query dataset to the reference database, achieving high-confidence automatic cell type annotation.

# Installation
Before running CellBlaster, ensure you have Python 3.8+ installed. You can install all required dependencies using pip:
```
pip install requests pandas numpy seaborn matplotlib tqdm scanpy scikit-learn
```
Open your terminal and clone the CellBlaster Repository.
The total installation time is around 1-2 mintunes. If error occuors, please upgrade pip and try again.
```
git clone https://github.com/illuminate6060/CellBLASTer.git
cd CellBlaster
pip install .
```
# 



# Usage1: Annotation by embedded dataset
## Database
**(1) CellBlaster features two built-in databases: Dicot and Monocot.**

    Dicot Database: For broad-leaf plants.
    
    Monocot Database: For grasses/grains.
    
**(2) Please select the one that matches your sc/snRNA-seq data.**

**(3) Detailed dataset specifications are available in the Other Information section below.**


## Python Example
Comment out your **h5ad** file using the CellBlaster software in your  Python program, as shown below. 
Sample code is in the "**tests**" directory.
```
import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

#Import the CellBlaster class from your installed package
from CellBlaster.CellBlaster import CellBlaster

# Change the current working directory to the test folder 
# Please change to real path!!!
os.chdir("../CellBlaster-main/tests")

# Define analysis parameters
dabase_type = "Dicot"
symbols = ["CRA008947", "CRA007122"]
output_path = './Output'
query_path = "./Demo_Data_SRP285040.h5ad"
query_symbol = "SRP285040"
filter_keywords = ["AthLNC", "Mt-", "cp"]

# Initialize the CellBlaster instance
# This sets up the configuration and maps the provided parameters to the object.
cellblaster = CellBlaster(
    output_path, 
    symbols, 
    dabase_type, 
    query_path, 
    query_symbol, 
    filter_keywords
)

# Execute the annotation pipeline
# This method handles database downloading, sequence generation, and cell type mapping.
result = cellblaster.Annotation(
    output_path, 
    symbols, 
    dabase_type, 
    query_path, 
    query_symbol, 
    filter_keywords
)
```
## Linux Example
```
python ../CellBlaster-main/CellBlaster-main/CellBlaster/CellBlaster.py \
    -t Dicot \
    -s CRA008947 CRA007122 \
    -o ../CellBlaster-main/tests/Output \
    -q ../CellBlaster-main/CellBlaster-main/tests/Demo_Data_SRP285040.h5ad \
    -qs SRP285040 \
    -f AthLNC Mt- cp &
```

## Argument Reference
| Argument | Shortcut | Description |
| :--- | :---: | :--- |
| `--dabase_type` | `-t` | **Required**. Database type: `Dicot` or `Monocot`. <br> Determines the Orthogroups and background datasets used. |
| `--symbols` | `-s` | **Required**. List of reference IDs (e.g., `-s CRA008947`). <br> Automatically downloads expression matrices, cell metadata, and DEGs to `01.DataBase`. |
| `--query_path` | `-q` | **Required**. Absolute path to the input `.h5ad` file containing the single-cell transcriptomic data for annotation. |
| `--query_symbol` | `-qs` | **Required**. Unique identifier for your query. <br> Used to name generated expression matrices and results in `02.QueryData`. |
| `--filter_keywords` | `-f` | **Optional**. List of keywords to filter out genes (case-insensitive). <br> Defaults to `LNC`. Can include `mt` or `cp` to remove organelle genes. |
| `--output_path` | `-o` | **Optional**. Root output directory (defaults to `./`). <br> Results are saved in the `03.Blast_Result` directory. |

# Usage2: Annotation by new defined dataset
If the built-in database does not meet your requirements, please refer to the examples in the Use_OwnData directory. The steps are as follows:

**(1) First, download the required data for OrthoFinder.**
Once the download is complete, add your specific isoform files to the download directory. The .fa files for the 11 species we provide can be found in the "Other Information" section at the end.
```
#--Step1: Download .fa file for orthofinder.
cd  path/CellBlaster/Use_OwnData
nohup python 1.Download_isoform.py -s T.aestivum  G.max L.japonicus M.truncatula -o ./Download &
```
**-s:** The prefix of the .fa files to be downloaded from the database. Multiple values can be declared, separated by spaces (e.g., A.thaliana T.aestivum G.max).

**-o:** The directory path to save the downloaded files; relative paths are supported.

**(2) Run OrthoFinder to reconstruct the Orthogroups.txt file.**
Ensure that the OrthoFinder path is correctly configured in your system.
```
#Step2：Orthofinder to genarate "Orthogroups.txt", check your orthofinder path.
nohup /path/you_software/orthofinder -f ./Download -t 80 -og -n result &
```
**-f:** The directory containing all prepared .fa files, including those downloaded in Step 1 and your own new species files.

**-t:** Number of threads for OrthoFinder to use.

**-n:** The directory where results will be stored. The Orthogroups.txt file will be located at: ./Download/OrthoFinder/Results_result/Orthogroups/Orthogroups.txt.

**(3) Perform annotation using the newly constructed dataset.**
Run the CellBlaster annotation pipeline based on your customized Orthogroups and datasets.
```
#Step3: CellBlaster for celltype annotation.
nohup python 2.New_CellBlaster.py  -O ./Download/OrthoFinder/Results_result/Orthogroups/Orthogroups.txt -s CRA008947   CRA007122 -o ./Output  -q   path/CellBlaster-main/tests/Demo_Data_SRP285040.h5ad   -qs SRP285040  -f AthLNC Mt- cp &
```

# Other Information
## **Isoform Datasets Information**
| Species | Abbreviation | Classification | Root | Leaf | Flower |
|---------|--------------|----------------|------|------|--------|
| *Arabidopsis thaliana* | *A. thaliana* | Dicotyledon | √ | √ | √ |
| *Glycine max* | *G. max* | Dicotyledon | √ | √ | |
| *Manihot esculenta* | *M. esculenta* | Dicotyledon | √ | √ | |
| *Medicago truncatula* | *M. truncatula* | Dicotyledon | √ | | |
| *Lotus japonicus* | *L. japonicus* | Dicotyledon | √ | | |
| *Fragaria vesca* | *F. vesca* | Dicotyledon | | √ | |
| *Catharanthus roseus* | *C. roseus* | Dicotyledon | | √ | |
| *Gossypium bickii* | *G. bickii* | Dicotyledon | | √ | |
| *Brassica rapa* | *B. rapa* | Dicotyledon | | √ | |
| *Gossypium hirsutum* | *G. hirsutum* | Dicotyledon | | | √ |
| *Bombax ceiba* | *B. ceiba* | Dicotyledon | | | √ |
| *Nepeta tenuifolia* | *N. tenuifolia* | Dicotyledon | | √ | |
| *Oryza sativa* | *O. sativa* | Monocotyledon | √ | √ | √ |
| *Triticum aestivum* | *T. aestivum* | Monocotyledon | √ | | |
| *Zea mays* | *Z. mays* | Monocotyledon | √ | √ | √ |
| *Sorghum bicolor* | *S. bicolor* | Monocotyledon | √ | | |
| *Setaria viridis* | *S. viridis* | Monocotyledon | √ | | |
| *Phyllostachys edulis* | *P. edulis* | Monocotyledon | √ | | |

## **Root Datasets Information**
| Species | Classification | Accession |
|---------|----------------|-----------|
| *Arabidopsis thaliana* | Dicotyledon | SRP267870 |
| *Arabidopsis thaliana* | Dicotyledon | SRP235541 |
| *Arabidopsis thaliana* | Dicotyledon | SRP171040 |
| *Arabidopsis thaliana* | Dicotyledon | SRP182008 |
| *Arabidopsis thaliana* | Dicotyledon | SRP166333 |
| *Arabidopsis thaliana* | Dicotyledon | SRP285817 |
| *Arabidopsis thaliana* | Dicotyledon | SRP273996 |
| *Arabidopsis thaliana* | Dicotyledon | SRP330542 |
| *Arabidopsis thaliana* | Dicotyledon | SRP173393 |
| *Arabidopsis thaliana* | Dicotyledon | SRP169576 |
| *Arabidopsis thaliana* | Dicotyledon | SRP148288 |
| *Arabidopsis thaliana* | Dicotyledon | SRP332285 |
| *Arabidopsis thaliana* | Dicotyledon | SRP285040 |
| *Arabidopsis thaliana* | Dicotyledon | SRP394711 |
| *Arabidopsis thaliana* | Dicotyledon | SRP422815 |
| *Arabidopsis thaliana* | Dicotyledon | SRP327656 |
| *Arabidopsis thaliana* | Dicotyledon | SRP363581 |
| *Arabidopsis thaliana* | Dicotyledon | SRP279055 |
| *Glycine max* | Dicotyledon | CRA007122 |
| *Glycine max* | Dicotyledon | CRA008947 |
| *Manihot esculenta* | Dicotyledon | SRP406470 |
| *Medicago truncatula* | Dicotyledon | SRP390780 |
| *Lotus japonicus* | Dicotyledon | SRP376527 |
| |  |  |
| *Oryza sativa* | Monocotyledon | SRP309176 |
| *Oryza sativa* | Monocotyledon | SRP250946 |
| *Oryza sativa* | Monocotyledon | CRA004082 |
| *Triticum aestivum* | Monocotyledon | CRA008788 |
| *Triticum aestivum* | Monocotyledon | GSE270342 |
| *Triticum aestivum* | Monocotyledon | SRP543892 |
| *Zea mays* | Monocotyledon | SRP145013 |
| *Zea mays* | Monocotyledon | GSE225118 |
| *Zea mays* | Monocotyledon | SRP335180 |
| *Sorghum bicolor* | Monocotyledon | SRP422815_1 |
| *Setaria viridis* | Monocotyledon | SRP422815_2 |
| *Phyllostachys edulis* | Monocotyledon | GSE229126 |

## **Leaf Datasets Information**
| Species | Classification | Accession |
|---------|----------------|-----------|
| *Arabidopsis thaliana* | Dicotyledon | SRP292306 |
| *Arabidopsis thaliana* | Dicotyledon | ERP132245 |
| *Arabidopsis thaliana* | Dicotyledon | SRP247828_1 |
| *Arabidopsis thaliana* | Dicotyledon | SRP247828_2 |
| *Arabidopsis thaliana* | Dicotyledon | SRP247828_3 |
| *Arabidopsis thaliana* | Dicotyledon | CRA002977_1 |
| *Arabidopsis thaliana* | Dicotyledon | SRP307169 |
| *Arabidopsis thaliana* | Dicotyledon | SRP280069 |
| *Arabidopsis thaliana* | Dicotyledon | SRP398011 |
| *Arabidopsis thaliana* | Dicotyledon | SRP338044 |
| *Arabidopsis thaliana* | Dicotyledon | E-MTAB-11006 |
| *Fragaria vesca* | Dicotyledon | CRA004848 |
| *Catharanthus roseus* | Dicotyledon | SRP335448 |
| *Gossypium bickii* | Dicotyledon | SRP424189 |
| *Glycine max* | Dicotyledon | CRA008947 |
| *Brassica rapa* | Dicotyledon | CRA006988 |
| *Nepeta tenuifolia* | Dicotyledon | SRP326816 |
| *Manihot esculenta* | Dicotyledon | CRA012723 |
| |  |  |
| *Oryza sativa* | Monocotyledon | SRP286275 |
| *Oryza sativa* | Monocotyledon | CRA004082 |
| *Zea mays* | Monocotyledon | SRP281914 |
| *Zea mays* | Monocotyledon | SRP325657 |
| *Zea mays* | Monocotyledon | SRP224648 |
| *Zea mays* | Monocotyledon | CRR923261 |
| *Zea mays* | Monocotyledon | SRP417893 |

## **Flower Datasets Information**
| Species | Classification | Accession |
|---------|----------------|-----------|
| *Arabidopsis thaliana* | Dicotyledon | SRP320285 |
| *Arabidopsis thaliana* | | EMTAB9174 |
| *Arabidopsis thaliana* | | SRP374045 |
| *Gossypium hirsutum* | Dicotyledon | SRP241596 |
| *Gossypium hirsutum* | | SRP379192 |
| *Bombax ceiba* | Dicotyledon | CRA009614 |
| |  |  |
| *Oryza sativa* | Monocotyledon | SRP386976 |
| *Zea mays* | Monocotyledon | SRP272727_23_26 |

### Contact Us
If you have any suggestions/ideas for CellBLASTer or are having issues trying to use it, please don't hesitate to reach out to us.
Lin Du,  3051095449@qq.com
