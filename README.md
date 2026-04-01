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
git clone https://github.com/illuminate6060/CellBlaster.git
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
query = "./Demo_Data_SRP285040.h5ad"
query_symbol = "SRP285040"
filter_keywords = ["AthLNC", "Mt-", "cp"]

# Initialize the CellBlaster instance
# This sets up the configuration and maps the provided parameters to the object.
cellblaster = CellBlaster(
    output_path, 
    symbols, 
    dabase_type, 
    query, 
    query_symbol, 
    filter_keywords
)

# Execute the annotation pipeline
# This method handles database downloading, sequence generation, and cell type mapping.
result = cellblaster.Annotation(
    output_path, 
    symbols, 
    dabase_type, 
    query, 
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
| `--query` | `-q` | **Required**. Absolute path to the input `.h5ad` file containing the single-cell transcriptomic data for annotation. |
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
nohup python 1.Download_isoform.py -s T.aestivum  G.max L.japonicus M.truncatula -o **./Download** &
```
**-s:** The prefix of the .fa files to be downloaded from the database. Multiple values can be declared, separated by spaces (e.g., A.thaliana T.aestivum G.max).

**-o:** The directory path to save the downloaded files; relative paths are supported.

**(2) Run OrthoFinder to reconstruct the Orthogroups.txt file.**
Ensure that the OrthoFinder path is correctly configured in your system.
```
#Step2：Orthofinder to genarate "Orthogroups.txt", check your orthofinder path.
nohup orthofinder -f ./**Download** -t 80 -og -n result &
```
**-f:** The directory containing all prepared .fa files, including those downloaded in Step 1 and your own new species files.

**-t: **Number of threads for OrthoFinder to use.

**-n: **The directory where results will be stored. The Orthogroups.txt file will be located at: ./Download/OrthoFinder/Results_result/Orthogroups/Orthogroups.txt.

**(3) Perform annotation using the newly constructed dataset.**
Run the CellBlaster annotation pipeline based on your customized Orthogroups and datasets.
```
#Step3: CellBlaster for celltype annotation.
nohup python 2.New_CellBlaster.py  -O ./Download/OrthoFinder/Results_result/Orthogroups/Orthogroups.txt -s CRA008947   CRA007122 -o ./Output  -q   path/CellBlaster-main/tests/Demo_Data_SRP285040.h5ad   -qs SRP285040  -f AthLNC Mt- cp &
```

# Other Information
## **Dicot Datasets Information**
| Species | Abbreviation | Count | SRA | Cell number |
| :--- | :--- | :---: | :--- | :---: |
| *Arabidopsis thaliana* | *A. thaliana* | 14 | SRP267870 | 188776 |
| | | | SRP235541 | 23519 |
| | | | SRP171040 | 33386 |
| | | | SRP182008 | 10353 |
| | | | SRP166333 | 16055 |
| | | | SRP285817 | 16572 |
| | | | SRP273996 | 9072 |
| | | | SRP330542 | 20193 |
| | | | SRP173393 | 10201 |
| | | | SRP169576 | 25338 |
| | | | SRP148288 | 14141 |
| | | | SRP332285 | 9837 |
| | | | SRP285040 | 973 |
| | | | SRP394711 | 310839 |
| *Glycine max* | *G. max* | 2 | CRA007122 | 14808 |
| | | | CRA008947 | 8021 |
| *Manihot esculenta* | *M. esculenta* | 1 | SRP406470 | 23765 |
| *Medicago truncatula* | *M. truncatula* | 1 | SRP390780 | 24978 |
| *Lotus japonicus* | *L. japonicus* | 1 | SRP376527 | 19921 |

## **Monocot Datasets Information**
| Species | Abbreviation | Count | SRA | Cell number |
| :--- | :--- | :---: | :--- | :---: |
| *Oryza sativa* | *O. sativa* | 2 | SRP309176 | 33720 |
| | | | SRP250946 | 28918 |
| *Triticum aestivum* | *T. aestivum* | 3 | CRA008788 | 3756 |
| | | | GSE270342 | 3846 |
| | | | SRP543892 | 35907 |
| *Zea mays* | *Z. mays* | 3 | SRP145013 | 19669 |
| | | | GSE225118 | 7303 |
| | | | SRP335180 | 4251 |
| *Sorghum bicolor* | *S. bicolor* | 1 | SRP422815_1 | 12839 |
| *Setaria viridis* | *S. viridis* | 1 | SRP422815_2 | 25609 |
| *Phyllostachys edulis* | *P. edulis* | 1 | GSE229126 | 3573 |

## **Isoform Datasets Information**
| Type | Parameter name | File name |
| :--- | :--- | :--- |
| Dicot | *A. thaliana* | *A. thaliana_isoform.fa* |
| Dicot | *G. max* | *G. max_isoform.fa* |
| Dicot | *L. japonicus* | *L. japonicus_isoform.fa* |
| Dicot | *M. esculenta* | *M. esculenta_isoform.fa* |
| Dicot | *M. truncatula* | *M. truncatula_isoform.fa* |
| Monocot | *O. sativa* | *O. sativa_isoform.fa* |
| Monocot | *P. edulis* | *P. edulis_isoform.fa* |
| Monocot | *S. bicolor* | *S. bicolor_isoform.fa* |
| Monocot | *S. viridis* | *S. viridis_isoform.fa* |
| Monocot | *T. aestivum* | *T. aestivum_isoform.A.fa* |
| Monocot | | *T. aestivum_isoform.B.fa* |
| Monocot | | *T. aestivum_isoform.D.fa* |
| Monocot | *Z. mays* | *Z. mays_isoform.fa* |

### Contact Us
If you have any suggestions/ideas for CellBlaster or are having issues trying to use it, please don't hesitate to reach out to us.
Lin Du,  3051095449@qq.com
