rule all:
    input:
        expand("SRR2584857_quast.{lines}", lines=range(250_000, 5_000_000, 250_000)),
        expand("SRR2584857_annot.{lines}", lines=range(250_000, 5_000_000, 250_000)),
        "n50_summary.csv",

rule subset_reads:
    input:
        "{sample}.fastq.gz",
    output:
        "{sample}.{subset,\d+}.fastq.gz"
    shell: """
        gunzip -c {input} | head -{wildcards.subset} | gzip -9c > {output} || true
    """

rule annotate:
    input:
        "SRR2584857-assembly.{subset}.fa"
    output:
        directory("SRR2584857_annot.{subset}")
    shell: """
       prokka --prefix {output} {input}                                       
    """

rule assemble:
    input:
        r1 = "SRR2584857_1.{subset}.fastq.gz",
        r2 = "SRR2584857_2.{subset}.fastq.gz"
    output:
        dir = directory("SRR2584857_assembly.{subset}"),
        assembly = "SRR2584857-assembly.{subset}.fa"
    threads: 16
    shell: """
       megahit -1 {input.r1} -2 {input.r2} -f -m 5e9 -t {threads} -o {output.dir}     
       cp {output.dir}/final.contigs.fa {output.assembly}                     
    """

rule quast:
    input:
        "SRR2584857-assembly.{subset}.fa"
    output:
        directory("SRR2584857_quast.{subset}")
    shell: """                                                                
       quast {input} -o {output}                                              
    """

rule n50_by_lines:
    input:
        "SRR2584857_quast.{subset}",
    output:
        "n50.{subset}.txt",
    shell: """
        echo {wildcards.subset},$(grep N50 {input}/report.tsv | cut -f2) > {output}
    """

rule summarize_n50:
    input:
        expand("n50.{lines}.txt", lines=range(250_000, 5_000_000, 250_000)),
    output:
        "n50_summary.csv"
    shell: """
        echo lines,n50 > {output}
        cat {input} >> {output}
    """
        
