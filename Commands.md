**Igblast command for simulated dataset**

Diverse repertoire 

`./bin/igblastn -query ../merged_reads.fasta -germline_db_V ./bin/database/V.fasta 
-germline_db_D ./bin/database/D.fasta -germline_db_J ./bin/database/J.fasta -outfmt 19 -auxiliary_data 
./optional_file/btau_gl.aux -out ../simulated_igblast_merged.tsv`

Polarized repertoire

` ./bin/igblastn -query ../merged_reads_polarised.fasta -germline_db_V ./bin/database/V.fasta -germline_db_D ./bin/database/D.fasta -germline_db_J ./bin/database/J.fasta -outfmt 19 -auxiliary_data ./optional_file/btau_gl.aux  -out ../simulated_igblast_polarized.tsv `

IgBlast Version : 1.15.0

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
  
IgSimulator Version : 2.0
 ART_Illumina Version : 2.1.8
 
IgSimulator command for modified germline J gene

` ./ig_simulator.py --chain-type HC --num-bases 100000 --num-mutated 200000 --repertoire-size 500000 --vgenes ./data/bovine_ig_germline_genes/B_V.fasta --dgenes ./data/bovine_ig_germline_genes/IGHD.fasta --jgenes ./data/bovine_ig_germline_genes/IGHJ_Modified.fasta -o Bovine_IgSimulation_Modified_Jgene`


**NB you can assemble the reads independently using pear**

`pear -f simulated_reads1.fq -r simulated_reads2.fq -o simulated_pear`

then convert the fastq assembled file to fasta for annotation 

`seqtk seq -A simulated_pear.assembled.fastq > simulated_pear.assembled.fasta`

**Assembly step is no longer important for this analysis, use merged.fastq file for annotation**

**To have a similar mixcr output as igblast & imgt**

Realign your reads, but this time use `-OsaveOriginalReads=true`  this parameter 

`mixcr align -s bovine --library imgt.201948-5.sv3.0.12.json.gz -p rna-seq -OsaveOriginalReads=true --report trial_simulated.txt ../mixcr_simulated_data/Annotation_500,000_rep/simulated_reads1.fq ../mixcr_simulated_data/Annotation_500,000_rep/simulated_reads2.fq simulated_bovine.vdjca `

MiXCR Version : 3.0.10

mixcr run for polarized repertoire

` mixcr align -s bovine --library imgt.201948-5.sv3.0.12.json.gz -p rna-seq -OsaveOriginalReads=true --report Data/paired_polarized_report.txt Data/paired_polarized_reads1.fq Data/paired_polarized_reads2.fq Data/simulated_polarized_bovine.vdjca`


next export your alignment, for this case you are intrested with sequence ID, Vcall, Dcall and Jcall. In the export command `-descrsR1` flag is used. 

` mixcr exportAlignments -descrsR1 -vHit -dHit -jHit simulated_bovine.vdjca  simulated_alignments.txt `

**To get a full alignment for MiXCR annotation**

Add `-OallowPartialAlignments=true  -OallowNoCDR3PartAlignments=true` to the align command...this gives you a 98% alignment. 


You good to go....:smile:!!














