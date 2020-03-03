**Igblast command for simulated dataset**

`./bin/igblastn -query ../simulated_changed_header_R1_assemble-pass.fasta -germline_db_V ./bin/database/V.fasta 
-germline_db_D ./bin/database/D.fasta -germline_db_J ./bin/database/J.fasta -outfmt 19 -auxiliary_data 
./optional_file/btau_gl.aux -out ../simulated_igblast.tsv`

**IgSimulator command**

` ./ig_simulator.py --chain-type HC --num-bases 100000 --num-mutated 200000 --repertoire-size 500000 --vgenes ./data/bovine_ig_germline_genes/B_V.fasta --dgenes ./data/bovine_ig_germline_genes/IGHD.fasta --jgenes ./data/bovine_ig_germline_genes/IGHJ.fasta -o Bovine_IgSimulation_1 `

parameters
* For diverse repertoire (naive Bcell) - what I used
  num_base = 100000 ; num_mutation = 200000
* For polarised repertoire (classed switched)
  num_base = 20000 ; num_mutation = 100000

**Assembly before runing IgBlast or IMGT**

mofidy the header of the simulated reads using awk 

`awk '{print (NR%4 == 1) ? "@antibody" ++i : $0}' simulated_reads2.fq > simulated_changed_header_R2.fq`

**run asembly.py**

` AssemblePairs.py align -2 ./Data/simulated_changed_header_R2.fq -1 ./Data/simulated_changed_header_R1.fq --fasta --nproc 4 --coord illumina`

validate this command


