import CellBLASTer
import os

os.chdir("/tiandata2/dulin/Github/CellBLASTer/CellBLASTer_package/demo")

cellblaster = CellBLASTer.CellBlaster(
    database_type="Dicot",
    organ="Root",
    symbols=["SRP169576"],
    output_path="./Output",
    query="./Query_SRP285040.h5ad",
    query_symbol="SRP285040",
    reference_adata="./Reference_SRP285040.h5ad",
    reference_symbol="SRP285040",
    n_jobs=30,
)

result = cellblaster.Annotation()