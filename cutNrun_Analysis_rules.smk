#import cutNrunFunctions
configfile: "ConfigParameters.json"
################### VARS

##############################################################################################################################
### FUNCTIONS
def get_fastqs(wildcards):
    if(glob.glob("rawData" + "/"+ wildcards.sample + "/" + wildcards.sample + "_L00[0-9]_" + wildcards.readnum + "_001.fastq.gz")):
        return(sorted(glob.glob("rawData" + "/"+ wildcards.sample + "/" + wildcards.sample + "_L00[0-9]_" + wildcards.readnum + "_001.fastq.gz")))
    else:
        return(sorted(glob.glob("rawData" + "/"+ wildcards.sample + "/" + wildcards.sample + "_" + wildcards.readnum + "_001.fastq.gz")))

def getGenome(wildcards):
    if wildcards.species == "hg19":
        gen = "hs"
    elif wildcards.species == "mm10" or wildcards.species == "mm9":
        gen = "mm"
    else:
        print("Check input reference")
    return gen

def selectIndexFile(wildcards):
    #print(wildcards.species)
    if wildcards.gen_species == "hg19":
        IndexFile = {config["Species"]["BTindex_human"]}
    elif wildcards.gen_species == "mm10":
        IndexFile = {config["Species"]["BTindex_mouse"]}
    elif wildcards.gen_species == "mm9":
        IndexFile = {config["Species"]["BTindex_mouse_mm9"]}
    elif wildcards.gen_species == "sacCer3":
        IndexFile = {config["Species"]["BTindex_yeast"]}
    #IndexFile = 'BTindex_yeast'
    else:
        print("Check input reference")
    return IndexFile

def selectBowtieParams(wildcards):
    if wildcards.gen_species == "sacCer3":
        paramsBowtie = config["AlignmentParams"]["spikein-sacCer3"]
#        paramsBowtie = "--local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 --no-unal --no-mixed --no-discordant --phred33 --no-overlap --no-dovetail -I 10 -X 700"
    elif  wildcards.gen_species == "hg19" or  wildcards.gen_species == "mm10" or  wildcards.gen_species == "mm9":
        paramsBowtie = config["AlignmentParams"]["speciesParam"]
#        paramsBowtie = "--local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700"
    else:
        print("Check input reference")
    return paramsBowtie

def convert_integer(wildcards):
    inpFileName = wildcards.sample + "/macs2_peaks/" + wildcards.sample + "_" + wildcards.species + "_unmap/" + wildcards.sample + "_treat_pileup.bdg"
    # myfile = open("{wildcards.sample}/macs2_peaks/{wildcards.sample}_{wildcards.species}_unmap/{wildcards.sample}_treat_pileup.bdg")
    myfile = open(inpFileName,"r")
    # print("{wildcards.sample}_treat_pileup_Integer.bdg")
    outFileName = wildcards.sample + "/macs2_peaks/" + wildcards.sample + "_" + wildcards.species + "_unmap/" + wildcards.sample+"_treat_pileup_Integer.bdg"
    outFile = open(outFileName, "w+")

    for line in myfile:
        line = line.rstrip("\n")
        lineL = line.split("\t")
        ct = int(float(lineL[-1]))
        if(ct == 0): continue
        outFile.write("%s\t%s\t%s\t%d\n" % (lineL[0], lineL[1], lineL[2], ct))
    return(outFile)

def selectFragmentParam(wildcards):
    print(wildcards.fragSize)
    if wildcards.fragSize == "120":
        paramsBamSubset = "--maxFragmentLength 120"
    elif wildcards.fragSize == "150":
        paramsBamSubset = "--minFragmentLength 150"
    #    IndexFile = {config['BTindex_yeast']}
    #    IndexFile = 'BTindex_yeast'
    else:
        print("Check input reference")
    return paramsBamSubset


def returnLineCount(wildcards):
    inputName= wildcards.sample+ "/bams/" + wildcards.sample+ "_"+config["Species"]["spikeIn"]+"_rdup"+".bed"
    #f = open(inputName,"r")
    scale=10000
    #spikeCount=f.readlines()
    spikeCount=0
    with open(inputName) as f:
        spikeCount=sum(1 for line in f)

    scaleFactor=scale/spikeCount
    print(scaleFactor)
    print(spikeCount)
    return(scaleFactor)



##################################################################################################
rule trim_fastq:
    input:
        pair1 = "{sample}/fastqs/{sample}_R1.fastq.gz",
        pair2 = "{sample}/fastqs/{sample}_R2.fastq.gz"
    output:
        trimmed_pair1 = "{sample}/fastqs_trimmed/{sample}_R1_trim.fq.gz",
        trimmed_pair2 = "{sample}/fastqs_trimmed/{sample}_R2_trim.fq.gz"
    log:
        "{sample}/logs/{sample}.trim_adapters.log"

    message: "Trimming adapters using TRIMGALORE!"
    shell:
        """
        {config[ToolPATH][TRIMGALORE_PATH]}/trim_galore {input.pair1} {input.pair2} --paired -o {wildcards.sample}/fastqs_trimmed &&
        mv {wildcards.sample}/fastqs_trimmed/{wildcards.sample}_R1_val_1.fq.gz {wildcards.sample}/fastqs_trimmed/{wildcards.sample}_R1_trim.fq.gz &&
        mv {wildcards.sample}/fastqs_trimmed/{wildcards.sample}_R2_val_2.fq.gz {wildcards.sample}/fastqs_trimmed/{wildcards.sample}_R2_trim.fq.gz
        """

rule fastqc_trim:
    input:
            inputFiles = "{sample}/fastqs_trimmed/{sample}_{readnum}_trim.fq.gz"
    threads : 10
    output: "{sample}/qc/{sample}_{readnum}_trim_fastqc.html",
            "{sample}/qc/{sample}_{readnum}_trim_fastqc.zip"

    log:    "{sample}/logs/{sample}_{readnum}_trim.fq.err"
    message: "FastQC: fastqc {input}"
    shell: """
            {config[ToolPATH][FASTQC_PATH]} --outdir {wildcards.sample}/qc --thread {threads} --nogroup {input} \
            >> {wildcards.sample}/qc/fastqc_trim.log 2>>{log}
            """

rule trim_fastq_fastqc_trimmomatic:
    input:
        pair1 = "{sample}/fastqs/{sample}_R1.fastq.gz",
        pair2 = "{sample}/fastqs/{sample}_R2.fastq.gz"
    output:
       fwd_paired = "{sample}/fastqs_trimmed_trimmomatic/{sample}_trimmed_PE_R1.fq.gz",
       fwd_unpaired = "{sample}/fastqs_trimmed_trimmomatic/{sample}_trimmed_SE_R1.fq.gz",
       rev_paired = "{sample}/fastqs_trimmed_trimmomatic/{sample}_trimmed_PE_R2.fq.gz",
       rev_unpaired = "{sample}/fastqs_trimmed_trimmomatic/{sample}_trimmed_SE_R2.fq.gz"
    log:
        "{sample}/logs/{sample}.trim_adapters_trimmomatic.log"
    message: "Trimming adapters using TRIMMOMATIC"
    shell:  """
         java -jar {config[ToolPATH][TRIMMOMATIC_JAR]} \
             PE {input.pair1} {input.pair2} \
             {output.fwd_paired} {output.fwd_unpaired} \
             {output.rev_paired} {output.rev_unpaired} \
             ILLUMINACLIP:{config[ToolPATH][TRIM_ADAPTER_PATH]}/{config[AdapterTrimTrimmomatic][adapter]} LEADING:{config[AdapterTrimTrimmomatic][leading]} TRAILING:{config[AdapterTrimTrimmomatic][trailing]} SLIDINGWINDOW:{config[AdapterTrimTrimmomatic][window]} MINLEN:{config[AdapterTrimTrimmomatic][minlen]}

         {config[ToolPATH][FASTQC_PATH]} {fwd_paired} {rev_paired} -o qc
         """

rule bowtie2:
    input:
        trimmed_pair1 = "{sample}/fastqs_trimmed/{sample}_R1_trim.fq.gz",
        trimmed_pair2 = "{sample}/fastqs_trimmed/{sample}_R2_trim.fq.gz",
    output:
      # bam = "bams/{sample}_{species}.sorted.bam",
        unmap_bam = "{sample}/bams/{sample}_{gen_species}.sorted.unmap.bam",
        #bambai = "bams/{sample}_{species}.sorted.bam.bai"
        flagstatOut = "{sample}/flagstats/{sample}_{gen_species}_unmap.flagstat"
    threads:
        30
    log:
       "{sample}/logs/{sample}_{gen_species}.alignment.log"
    params:
       Index = selectIndexFile,
       bowtie2Parameters = selectBowtieParams
#       bowtie2Parameters = config

    message: "Aligning to the reference."
    shell:"""
        {config[ToolPATH][BOWTIE2_PATH]} -p {threads} -x {params.Index} {params.bowtie2Parameters} -1 {input.trimmed_pair1} -2 {input.trimmed_pair2} 2> {log} > {wildcards.sample}/bams/{wildcards.sample}_{wildcards.gen_species}.sam
        {config[ToolPATH][SAMTOOLS_PATH]} view -bS {wildcards.sample}/bams/{wildcards.sample}_{wildcards.gen_species}.sam |  {config[ToolPATH][SAMTOOLS_PATH]} sort -o - - > {wildcards.sample}/bams/{wildcards.sample}_{wildcards.gen_species}.sorted.bam
        {config[ToolPATH][SAMTOOLS_PATH]} view -h -F 4 -b {wildcards.sample}/bams/{wildcards.sample}_{wildcards.gen_species}.sorted.bam > {wildcards.sample}/bams/{wildcards.sample}_{wildcards.gen_species}.sorted.unmap.bam
        sleep 10
        {config[ToolPATH][SAMTOOLS_PATH]} index {wildcards.sample}/bams/{wildcards.sample}_{wildcards.gen_species}.sorted.unmap.bam
        {config[ToolPATH][SAMTOOLS_PATH]} flagstat -@ {threads} {output.unmap_bam} > {wildcards.sample}/flagstats/{wildcards.sample}_{wildcards.gen_species}_unmap.flagstat 1> {wildcards.sample}/logs/{wildcards.sample}_{wildcards.gen_species}_flagstat.log

        rm {wildcards.sample}/bams/{wildcards.sample}_{wildcards.gen_species}.sam
        rm {wildcards.sample}/bams/{wildcards.sample}_{wildcards.gen_species}.sorted.bam
        """


rule peak_calling_macs2:
    input:
        chip = "{sample}/bams/{sample}_{species}.sorted.unmap.bam"

    output:
        "{sample}/macs2_peaks/{sample}_{species}_unmap/{sample}_peaks.narrowPeak",
        "{sample}/macs2_peaks/{sample}_{species}_unmap/{sample}_summits.bed",
        "{sample}/macs2_peaks/{sample}_{species}_unmap/{sample}_treat_pileup.bdg"
    params:
        genome = getGenome
    log:
        "{sample}/logs/{sample}_{species}.macs2"

    message: "Peak calling using MACS2 - p-value 1e-3"

    shell:"""
   # export PATH=/tools/miniconda3/bin/:$PATH
   # export PYTHONPATH=/tools/miniconda3/lib/python3.7/site-packages/:$PYTHONPATH

    {config[ToolPATH][CONDA_PATH]}/macs2 callpeak -g {params.genome} -f BAMPE --pvalue {config[Macs2Settings][pval]} -t {input.chip} -n {wildcards.sample} --keep-dup all -B --outdir {wildcards.sample}/macs2_peaks/{wildcards.sample}_{wildcards.species}_unmap 2>>{log}
"""

rule generate_Bigwigs:
    input :
        chip = "{sample}/bams/{sample}_{species}.sorted.unmap.bam"
    output:
        bigwigOut = "{sample}/bams/{sample}_{species}.sorted.unmap.bigwig"
    log:
        "{sample}/logs/{sample}_{species}.bigwig"

    message: "Generating bigwig files"

    shell:"""
#        export PYTHONPATH='/tools/deepTools/:/usr/lib/:/usr/lib/python3.6/site-packages':$PYTHONPATH
        {config[ToolPATH][DEEPTOOLS_PATH]}/bamCoverage --bam {input.chip} -o {output} --normalizeUsing RPKM --numberOfProcessors 30 -bs 20 --smoothLength 60
    """

rule rmdup:
    input:
        "{sample}/bams/{sample}_{gen_species}.sorted.unmap.bam"
    output:
        bam = "{sample}/bams/{sample}_{gen_species}.sorted.rdup.bam"
    params:
        picardmetric = "{sample}/logs/{sample}_{gen_species}.markdups.metrics.txt"

    message: "Marking and removing duplicates"
    shell:"""
        java -jar {config[ToolPATH][PICARD_PATH]} MarkDuplicates INPUT={input} OUTPUT={output.bam} REMOVE_DUPLICATES=TRUE METRICS_FILE={params.picardmetric}
        sleep 10
        {config[ToolPATH][SAMTOOLS_PATH]} index {output.bam}

        """

rule peak_calling_macs2_rdup:
    input:
        chip = "{sample}/bams/{sample}_{species}.sorted.rdup.bam"

    output:
        "{sample}/macs2_peaks/{sample}_{species}_rdup/{sample}_peaks.narrowPeak",
        "{sample}/macs2_peaks/{sample}_{species}_rdup/{sample}_summits.bed",
        "{sample}/macs2_peaks/{sample}_{species}_rdup/{sample}_treat_pileup.bdg"
    params:
        genome = getGenome
    log:
        "{sample}/logs/{sample}_{species}_rdup.macs2"

    message: "Peak calling using MACS2 - p-value 1e-3"

    shell:"""
    #export PATH=/tools/miniconda3/bin/:$PATH
    #export PYTHONPATH=/tools/miniconda3/lib/python3.7/site-packages/:$PYTHONPATH

    {config[ToolPATH][CONDA_PATH]}/macs2 callpeak -g {params.genome} -f BAMPE --pvalue {config[Macs2Settings][pval]} -t {input.chip} -n {wildcards.sample} --keep-dup all -B --outdir {wildcards.sample}/macs2_peaks/{wildcards.sample}_{wildcards.species}_rdup 2>>{log}
"""

rule generate_Bigwigs_rdup:
    input :
        chip = "{sample}/bams/{sample}_{species}.sorted.rdup.bam"
    output:
        bigwigOut = "{sample}/bams/{sample}_{species}.sorted.rdup.bigwig"
    log:
        "{sample}/logs/{sample}_{species}_rdup.bigwig"

    message: "Generating bigwig files"

    shell:"""
  #      export PYTHONPATH='/tools/deepTools/:/usr/lib/:/usr/lib/python3.6/site-packages':$PYTHONPATH
        {config[ToolPATH][DEEPTOOLS_PATH]}/bamCoverage --bam {input.chip} -o {output} --normalizeUsing RPKM --numberOfProcessors 30 -bs 20 --smoothLength 60
"""


rule generateFragments:
        input :
            bamRdup = "{sample}/bams/{sample}_{species}.sorted.rdup.bam",
            bamUnmap = "{sample}/bams/{sample}_{species}.sorted.unmap.bam",
        output:
            bamRdup_frag="{sample}/bams/{sample}_{species}.sorted.rdup_{fragSize}.bam",
            bamUnmap_frag="{sample}/bams/{sample}_{species}.sorted.unmap_{fragSize}.bam",
            bigwigRdup_frag="{sample}/bams/{sample}_{species}.sorted.rdup_{fragSize}.bigwig",
            bigwigUnmap_frag="{sample}/bams/{sample}_{species}.sorted.unmap_{fragSize}.bigwig",
            bamRdup_frag_METRICS="{sample}/bams/{sample}_{species}.sorted.rdup_{fragSize}_metrics.txt",
            bamUnmap_frag_METRICS="{sample}/bams/{sample}_{species}.sorted.unmap_{fragSize}_metrics.txt"
        log:
            "{sample}/logs/{sample}_{species}_Fragments_{fragSize}.log"
        params:
            FilterParam = selectFragmentParam

        message: "Generating small fragment files"

        shell:"""
     #   export PYTHONPATH='/tools/deepTools/:/usr/lib/:/usr/lib/python3.6/site-packages':$PYTHONPATH

            {config[ToolPATH][DEEPTOOLS_PATH]}/alignmentSieve -b {input.bamRdup} {params.FilterParam} -o {output.bamRdup_frag} --filterMetrics {output.bamRdup_frag_METRICS}
            {config[ToolPATH][DEEPTOOLS_PATH]}/alignmentSieve -b {input.bamUnmap} {params.FilterParam} -o {output.bamUnmap_frag} --filterMetrics {output.bamUnmap_frag_METRICS}

            {config[ToolPATH][SAMTOOLS_PATH]} index {output.bamRdup_frag}
            {config[ToolPATH][SAMTOOLS_PATH]} index {output.bamUnmap_frag}

            {config[ToolPATH][DEEPTOOLS_PATH]}/bamCoverage --bam {output.bamRdup_frag} -o {output.bigwigRdup_frag} --normalizeUsing RPKM --numberOfProcessors 30 -bs 20 --smoothLength 60
            {config[ToolPATH][DEEPTOOLS_PATH]}/bamCoverage --bam {output.bamUnmap_frag} -o {output.bigwigUnmap_frag} --normalizeUsing RPKM --numberOfProcessors 30 -bs 20 --smoothLength 60

        """


rule peak_calling_smallFragments:
    input:
        bamRdup_frag="{sample}/bams/{sample}_{species}.sorted.rdup_{fragSize}.bam",
        bamUnmap_frag="{sample}/bams/{sample}_{species}.sorted.unmap_{fragSize}.bam"
    output:
        "{sample}/macs2_peaks/{sample}_{species}_unmap_{fragSize}/{sample}_peaks.narrowPeak",
        "{sample}/macs2_peaks/{sample}_{species}_unmap_{fragSize}/{sample}_summits.bed",
        "{sample}/macs2_peaks/{sample}_{species}_unmap_{fragSize}/{sample}_treat_pileup.bdg",
        "{sample}/macs2_peaks/{sample}_{species}_rdup_{fragSize}/{sample}_peaks.narrowPeak",
        "{sample}/macs2_peaks/{sample}_{species}_rdup_{fragSize}/{sample}_summits.bed",
        "{sample}/macs2_peaks/{sample}_{species}_rdup_{fragSize}/{sample}_treat_pileup.bdg"

    params:
        genome = getGenome
    log:
        "{sample}/logs/{sample}_{species}.macs2_Fragments_{fragSize}.log"

    message: "Peak calling using MACS2 - p-value 1e-3"

    shell:"""
    #export PATH=/tools/miniconda3/bin/:$PATH
    #export PYTHONPATH=/tools/miniconda3/lib/python3.7/site-packages/:$PYTHONPATH

    {config[ToolPATH][CONDA_PATH]}/macs2 callpeak -g {params.genome} -f BAMPE --pvalue {config[Macs2Settings][pval]} -t {input.bamUnmap_frag} -n {wildcards.sample} --keep-dup all -B --outdir {wildcards.sample}/macs2_peaks/{wildcards.sample}_{wildcards.species}_unmap_{wildcards.fragSize} 2>>{log}
    {config[ToolPATH][CONDA_PATH]}/macs2 callpeak -g {params.genome} -f BAMPE --pvalue {config[Macs2Settings][pval]} -t {input.bamRdup_frag} -n {wildcards.sample} --keep-dup all -B --outdir {wildcards.sample}/macs2_peaks/{wildcards.sample}_{wildcards.species}_rdup_{wildcards.fragSize} 2>>{log}
"""


rule speciesBam2Bed:
    input:
        chip = "{sample}/bams/{sample}_{gen_species}.sorted.rdup.bam"

    output:
        bed = "{sample}/bams/{sample}_{gen_species}_rdup.bed"
    log:
        "{sample}/logs/{sample}_{gen_species}_bed.log"

    shell:"""
        {config[ToolPATH][SAMTOOLS_PATH]} sort -n {input.chip}| \
        {config[ToolPATH][BEDTOOLS_PATH]}/bedtools bamtobed \
        -bedpe -i - > {wildcards.sample}/bams/{wildcards.sample}_{wildcards.gen_species}_rdup.bed

"""


rule spikeInNormalization:
    input:
            chip = "{sample}/bams/{sample}_{species}.sorted.rdup.bam",
    output:
            norm = "{sample}/bams/{sample}_{species}.sorted.rdup.SInorm.bigwig",
    log:
            "{sample}/logs/{sample}_{species}.spikeInNormalization.log"

    params:
        genome=getGenome,
        scaleFactorNum=returnLineCount

    shell: """
    {config[ToolPATH][DEEPTOOLS_PATH]}/bamCoverage -b {input.chip} --scaleFactor {params.scaleFactorNum} \
    -o {output.norm} -p 30 -of "bigwig" --binSize 20 --smoothLength 60
    """



rule spikeInNormalization_rpkm:
    input:
            chip = "{sample}/bams/{sample}_{species}.sorted.rdup.bam",
    output:
            norm = "{sample}/bams/{sample}_{species}.sorted.rdup.SInorm_RPKM.bigwig",
    log:
            "{sample}/logs/{sample}_{species}.spikeInNormalization.log"

    params:
        genome=getGenome,
        scaleFactorNum=returnLineCount

    shell: """

    {config[ToolPATH][DEEPTOOLS_PATH]}/bamCoverage \
    -b {input.chip} --scaleFactor {params.scaleFactorNum} \
    -o {output.norm} -p 30 -of "bigwig" --binSize 20 --normalizeUsing "RPKM" --smoothLength 60
    """



rule TagDirectory:
    input:
        bam="{sample}/bams/{sample}_{species}.sorted.rdup.bam",

    output:
        "{sample}/bams/{sample}_{species}.sorted_rdup_Tags/tagInfo.txt",

    log:
        "{sample}/logs/{sample}_{species}_rdup_Tags"

    shell:"""
    {config[ToolPATH][HOMER_PATH]}/makeTagDirectory {wildcards.sample}/bams/{wildcards.sample}_{wildcards.species}.sorted_rdup_Tags -genome hg19 {input.bam}
    """


rule HomerPeaks:
    input:
        tagDir="{sample}/bams/{sample}_{species}.sorted_rdup_Tags"
    #    bam="{sample}/bams/{sample}_{species}.sorted.rdup.bam",

    output:
        TXTfile="{sample}/homerPeaks/{sample}_{species}.sorted_rdup/{sample}_{species}_peaks.txt",
        Peaks="{sample}/homerPeaks/{sample}_{species}.sorted_rdup/{sample}_{species}_peaks.bed",

    log:
        "{sample}/logs/{sample}_{species}_rdup_HomerPeaks"

    shell:"""
    {config[ToolPATH][HOMER_PATH]}/findPeaks {input.tagDir} -style histone > {wildcards.sample}/homerPeaks/{wildcards.sample}_{wildcards.species}.sorted_rdup/{wildcards.sample}_{wildcards.species}_peaks.txt
    {config[ToolPATH][HOMER_PATH]}/pos2bed.pl -o {wildcards.sample}/homerPeaks/{wildcards.sample}_{wildcards.species}.sorted_rdup/{wildcards.sample}_{wildcards.species}_peaks.bed {wildcards.sample}/homerPeaks/{wildcards.sample}_{wildcards.species}.sorted_rdup/{wildcards.sample}_{wildcards.species}_peaks.txt
"""



rule peak_calling_macs2_rdup_BROAD:
    input:
        chip = "{sample}/bams/{sample}_{species}.sorted.rdup.bam"

    output:
        "{sample}/macs2_peaks/{sample}_{species}_rdup_broad/{sample}_peaks.broadPeak",
    params:
        genome = getGenome
    log:
        "{sample}/logs/{sample}_{species}_rdup_broad.macs2"

    message: "Peak calling using MACS2 - "

    shell:"""
    #export PATH=/tools/miniconda3/bin/:$PATH
    #export PYTHONPATH=/tools/miniconda3/lib/python3.7/site-packages/:$PYTHONPATH

    {config[ToolPATH][CONDA_PATH]}/macs2 callpeak -g {params.genome} --broad --broad-cutoff 0.1 -q 0.1 -f BAMPE  -t {input.chip} -n {wildcards.sample} --keep-dup all -B --outdir {wildcards.sample}/macs2_peaks/{wildcards.sample}_{wildcards.species}_rdup_broad 2>>{log}
"""
