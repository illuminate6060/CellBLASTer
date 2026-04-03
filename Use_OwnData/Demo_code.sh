#--Step1: Download data for orthofinder. After download, add your isoform to the Download directory
nohup python 1.Download_isoform.py -s T.aestivum  G.max L.japonicus M.truncatula -o ./Download &

#Step2：Orthofinder to genarate "Orthogroups.txt", check your orthofinder path.
nohup orthofinder -f ./Download -t 80 -og -n result &

#Step3: CellBlaster for celltype annotation.
nohup python 2.New_CellBlaster.py  -O ./Download/OrthoFinder/Results_result/Orthogroups/Orthogroups.txt -s CRA008947   CRA007122 -o ./Output  -q   ../tests/Demo_Data_SRP285040.h5ad   -qs SRP285040  -f AthLNC Mt- cp &