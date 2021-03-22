#import cutNrunFunctions
configfile: "ConfigParameters.json"
################### VARS

rule peak_calling_broad_EPIC:
    input:
        chipBed = "{sample}/bams/{sample}_{species}.sorted.rdup.bed"
    output:
        peaks="{sample}/epic2_peaks/{sample}_{species}_rdup_EPIC/{sample}_peaks_epic.bed",
    params:
    #    genome = getGenome
    log:
        "{sample}/logs/{sample}_{species}_rdup_broad.epic2"

    message: ""
    shell: """

    epic2 \
    --treatment {input.chipBed}  \
    --genome {wildcards.species} \
    -kd \
    -fdr 0.001 \

    -o {output.peaks} > {log}
"""

rule spikeInNormalization_wig2bed:
    input:
            norm = "{sample}/bams/{sample}_{species}.sorted.rdup.SInorm.bigwig",
    output:
            bed="{sample}/bams/{sample}_{species}.sorted.rdup.SInorm_wig2bed.bed",
            bed6="{sample}/bams/{sample}_{species}.sorted.rdup.SInorm_wig2_bed6.bed"
    log:
            "{sample}/logs/{sample}_{species}.spikeInNormalization_bedfiles.log"

    params:
        wigfile="{sample}/bams/{sample}_{species}.sorted.rdup.SInorm.wig"

    shell: """
#    {config[ToolPATH][UCSCscripts]}/bigWigToWig {input.norm} {params.wigfile}
    {config[ToolPATH][BEDOPSS_PATH]}/wig2bed --zero-indexed < {params.wigfile} > {output.bed}

    less {output.bed}| sed 's/$/\t*/' > {output.bed6}
    """


rule SInorm_peaks:
    input:
        NORMbdg="{sample}/bams/{sample}_{species}.sorted.rdup.SInorm.bdg",

    output:
        "{sample}/macs2_peaks/{sample}_{species}_rdup_SInorm/{sample}_SInorm_peaks.narrowPeak",

    params:
        genome = getGenome
    log:
        "{sample}/logs/{sample}_{species}_rdup_SInorm.macs2"

    shell: """
    {config[ToolPATH][CONDA_PATH]}/macs2 callpeak \
    -g {params.genome} \
    -f BED \
    --pvalue {config[Macs2Settings][pval]} \
    -t {input.NORMbdg} \
    -n {wildcards.sample}_SInorm \
    --keep-dup all \
    -B \
    --outdir {wildcards.sample}/macs2_peaks/{wildcards.sample}_{wildcards.species}_rdup_SInorm 2>>{log}
    """

rule TagDirectory_SInorm:
    input:
        NORMbdg="{sample}/bams/{sample}_{species}.sorted.rdup.SInorm_wig2_bed6.bed",

    output:
        "{sample}/bams/{sample}_{species}.sorted_rdup_SInorm_Tags/tagInfo.txt",

    log:
        "{sample}/logs/{sample}_{species}_rdup_SInorm_Tags"

    shell:"""
    {config[ToolPATH][HOMER_PATH]}/makeTagDirectory  {wildcards.sample}/bams/{wildcards.sample}_{wildcards.species}.sorted_rdup_SInorm_Tags -genome hg19 -format bed {input.NORMbdg}
    """
