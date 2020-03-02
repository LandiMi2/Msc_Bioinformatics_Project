

SAMPLES, = glob_wildcards("Data/{samples}_reads1.fq")

rule all:
         input: expand("Data/{samples}_clones.txt", samples = SAMPLES )
              


rule align: 
    input: 
          "Data/{samples}_reads1.fq",
          "Data/{samples}_reads2.fq",
          "imgt.201948-5.sv3.0.12.json.gz"
    
    output: 
           "Data/{samples}_report.txt",
           "Data/{samples}_bovine.vdjca"

    shell:
           "mixcr align -s bovine --library {input[2]} -p rna-seq --report {output[0]} {input[0]} {input[1]} {output[1]}"


rule assemble: 
    input:
          "Data/{samples}_report.txt",
          "Data/{samples}_bovine.vdjca"

    output:
           "Data/{samples}_clones.clna",
           "Data/{samples}_report.txt"

    shell:
          "mixcr assemble --write-alignments --report {input[0]} {input[1]} {output[0]}" 


rule assemble_full_contig: 
    input:
          "Data/{samples}_clones.clna",
          "Data/{samples}_report.txt"

    output:
           "Data/{samples}_clones.clns"

    shell:
            "mixcr assembleContigs --report {input[1]} {input[0]} {output}"


rule export: 
    input:  
          "Data/{samples}_clones.clns"
    
    output: 
          "Data/{samples}_clones.txt"

    shell: 
           "mixcr exportClones -c IG -p fullImputed {input} {output}"

          
          




