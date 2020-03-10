**Igblast command for simulated dataset**

Diverse repertoire 

`./bin/igblastn -query ../merged_reads.fasta -germline_db_V ./bin/database/V.fasta 
-germline_db_D ./bin/database/D.fasta -germline_db_J ./bin/database/J.fasta -outfmt 19 -auxiliary_data 
./optional_file/btau_gl.aux -out ../simulated_igblast_merged.tsv`

Polarized repertoire

` ./bin/igblastn -query ../merged_reads_polarised.fasta -germline_db_V ./bin/database/V.fasta -germline_db_D ./bin/database/D.fasta -germline_db_J ./bin/database/J.fasta -outfmt 19 -auxiliary_data ./optional_file/btau_gl.aux  -out ../simulated_igblast_polarized.tsv `

**IgSimulator command**

Run for diverse antibody repertoire

` ./ig_simulator.py --chain-type HC --num-bases 100000 --num-mutated 200000 --repertoire-size 500000 --vgenes ./data/bovine_ig_germline_genes/B_V.fasta --dgenes ./data/bovine_ig_germline_genes/IGHD.fasta --jgenes ./data/bovine_ig_germline_genes/IGHJ.fasta -o Bovine_IgSimulation_1 `

Run for polarized antibody repertoire

`./ig_simulator.py --chain-type HC --num-bases 20000 --num-mutated 100000 --repertoire-size 500000 --vgenes ./data/bovine_ig_germline_genes/B_V.fasta --dgenes ./data/bovine_ig_germline_genes/IGHD.fasta --jgenes ./data/bovine_ig_germline_genes/IGHJ.fasta -o Bovine_IgSimulation_Polarised`


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

**Assembly step is no longer important for this analysis, use merged.fastq file for annotation**

**To have a similar mixcr output as igblast & imgt**

Realign your reads, but this time use `-OsaveOriginalReads=true`  this parameter 

`mixcr align -s bovine --library imgt.201948-5.sv3.0.12.json.gz -p rna-seq -OsaveOriginalReads=true --report trial_simulated.txt ../mixcr_simulated_data/Annotation_500,000_rep/simulated_reads1.fq ../mixcr_simulated_data/Annotation_500,000_rep/simulated_reads2.fq trial_simulated.vdjca`

next export your alignment, for this case you are intrested with sequence ID, Vcall, Dcall and Jcall. In the export command `-descrsR1` flag is used. 

` mixcr exportAlignments -descrsR1 -vHit -dHit -jHit simulated_bovine.vdjca  simulated_alignments.txt `

You good to go....:smile:!!














