rule bed2fasta:
    input:
        bed="resources/annotations/{reference}/{prefix}.bed",
        fasta="resources/genomes/{reference}_genome.fasta",
    output:
        fasta="
        


rule salmon_quant:
    input
