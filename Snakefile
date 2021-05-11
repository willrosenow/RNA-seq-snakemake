configfile: "config.yaml"

SAMPLES = config['SAMPLES']
HISAT2_INDEX_PREFIX = config['HISAT2_INDEX_PREFIX']

print('Samples Included:')
print(SAMPLES)
####################################################################################################################################

#Target Rule (last output):
rule all:
     input: "differential_expression/gene_count_matrix.csv"

rule align_hisat:
    input:
        fq1= config['DATA_PATH'] + "{sample}" + config['R1_EXT'],
        fq2= config['DATA_PATH'] + "{sample}" + config['R2_EXT'],
        hisat2_index=expand(f"{HISAT2_INDEX_PREFIX}.{{ix}}.ht2", ix=range(1, 9))
    output: "align_hisat2/{sample}.bam"
    log: "align_hisat2/{sample}.alignment.summary"
    threads: config['THREADS']
    shell:
        "hisat2 2>{log} -p {threads} --dta -x {HISAT2_INDEX_PREFIX} "
        "-1 {input.fq1} -2 {input.fq2} | "
        "samtools sort -@ {threads} -o {output}"

rule stringtie_assemble:
    input:
        genome_gtf=config['GENOME_GTF'],
        bam="align_hisat2/{sample}.bam"
    output: "stringtie/assembled/{sample}.gtf"
    threads: config['THREADS']
    shell:
        "stringtie -p {threads} -G {input.genome_gtf} "
        "-o {output} -l {wildcards.sample} {input.bam}"

rule stringtie_merge_list:
    input: expand("stringtie/assembled/{sample}.gtf", sample=SAMPLES)
    output: "stringtie/merged_list.txt"
    run:
        with open(output[0], 'w') as f:
            for gtf in input:
                print(Path(gtf).resolve(), file=f)

rule stringtie_merge:
    input:
        genome_gtf=config['GENOME_GTF'],
        merged_list="stringtie/merged_list.txt",
        sample_gtfs=expand("stringtie/assembled/{sample}.gtf", sample=SAMPLES)
    output: "stringtie/merged.gtf"
    threads: config['THREADS']
    shell:
        "stringtie --merge -p {threads} -G {input.genome_gtf} "
        "-o {output} {input.merged_list}"


rule stringtie_quant:
    input:
        merged_gtf="stringtie/merged.gtf",
        sample_bam="align_hisat2/{sample}.bam"
    output:
        gtf="stringtie/quant/{sample}/{sample}.gtf",
        gene_abund_tab="stringtie/quant/{sample}/gene_abund.tab",
        ctabs=expand(
            "stringtie/quant/{{sample}}/{name}.ctab",
            name=['i2t', 'e2t', 'i_data', 'e_data', 't_data']
        )
    threads: config['THREADS']
    shell:
        "stringtie -e -B -p {threads} -G {input.merged_gtf} "
        "-A {output.gene_abund_tab} -o {output.gtf} {input.sample_bam}"

rule prep:
    input: expand("stringtie/quant/{sample}/{sample}.gtf", sample=SAMPLES)
    output: "differential_expression/prepDE_list.txt" 
    run:
        shell("ls -d $PWD/stringtie/quant/*/*.gtf > paths.txt")
        shell("ls ./stringtie/quant/ > filenames.txt")
        shell("mkdir -p differential_expression")
        shell("paste filenames.txt paths.txt > {output}")
        shell("rm filenames.txt paths.txt")

rule prepde:
    input: "differential_expression/prepDE_list.txt"
    output:
        gene="differential_expression/gene_count_matrix.csv",
        trans="differential_expression/transcript_count_matrix.csv"
    shell: 
        "python prepDE.py -i {input} "
        "-g {output.gene} -t {output.trans}"


print("Done")

