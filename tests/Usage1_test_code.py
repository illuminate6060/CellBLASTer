import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

#Import the CellBlaster class from your installed package
from CellBLASTer.CellBlaster import CellBlaster

# Change the current working directory to the test folder 
# Please change to real path!!!
os.chdir("../CellBLASTer-main/tests")

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
adata = cellblaster.Annotation(
    output_path, 
    symbols, 
    dabase_type, 
    query_path, 
    query_symbol, 
    filter_keywords
)






"""
other usage in Linux
python ../CellBlaster-main/CellBlaster-main/CellBlaster/CellBlaster.py \
    -t Dicot \
    -s CRA008947 CRA007122 \
    -o ../CellBlaster-main/tests/Output \
    -q ../CellBlaster-main/CellBlaster-main/tests/Demo_Data_SRP285040.h5ad \
    -qs SRP285040 \
    -f AthLNC Mt- cp &
"""