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

**NB you can assemble the reads independently using pear**

`pear -f simulated_reads1.fq -r simulated_reads2.fq -o simulated_pear`

then convert the fastq assembled file to fasta for annotation 

`seqtk seq -A simulated_pear.assembled.fastq > simulated_pear.assembled.fasta`

**Human Trial Simulation**

1. Simulate human read 1 and read 2

`./ig_simulator.py --chain-type HC --num-bases 100000 --num-mutated 200000 --repertoire-size 500000 --vgenes ./data/human_ig_germline_genes_upadated/human_IGHV.fa --dgenes ./data/human_ig_germline_genes_upadated/human_IGHD.fa --jgenes ./data/human_ig_germline_genes_upadated/human_IGHJ.fa -o Human_IgSimulation_Trial`

2. Assemble read 1 and 2 for annotation by IgBlast and IMGT

use awk to modify header ......(validate this step)

`pear -f human_simulated_changed_header_R1.fq -r human_simulated_changed_header_R2.fq -o human_simulated_pear_modified`

convert fastq to fasta 


3. Run Igblast 

` ./bin/igblastn -query ../human_simulated_pear_modified.assembled.fasta -germline_db_V ./bin/human_database/IGHV.fa -germline_db_D ./bin/human_database/IGHD.fa -germline_db_J ./bin/human_database/IGHJ.fa -outfmt 19 -auxiliary_data ./optional_file/human_gl.aux -out ../human_simulated_igblast_assembly_pear_modified.tsv`










