import sys
import glob
import pandas as pd
############################
#import cutNrunFunctions

############################################################################
################### CONFIGS AND METADATA
configfile: "ConfigParameters.json"
METADATA = pd.read_csv("SampleMetadata.csv", sep="\t")
SAMPLES = []
SAMPLES = list(METADATA.SampleID_New)
READNUM = ["R1","R2"]
################### VARS
SPECIES = config["Species"]["species"]
SPIKEIN = config["Species"]["spikeIn"]

if SPIKEIN == "none":
    GEN_SPECIES = config["Species"]["species"]
else:
    GEN_SPECIES = [SPIKEIN,SPECIES]
print(GEN_SPECIES)
############################################################################
### OTHER NOTES
print("Pipeline for paired-end cut&run sequences with spike in - no normalization yet")

print("Starting analysis for ... ")
print("Species : ",SPECIES)
print("Spike in control : ",SPIKEIN)
print("SAMPLES :",SAMPLES)
print("Analysis directory : ", os.path.abspath(os.getcwd()))
############################################################################
### TRIMMING
ALL_TRIMMED_FASTQ = expand("{sample}/fastqs_trimmed/{sample}_{readnum}_trim.fq.gz", sample = SAMPLES, readnum=READNUM)
ALL_FASTQC_TRIMMED = expand("{sample}/qc/{sample}_{readnum}_trim_fastqc.html", sample=SAMPLES, readnum=READNUM)
### UNMAP BAM PEAKS BIGWIGS
ALL_BAMS = expand("{sample}/bams/{sample}_{gen_species}.sorted.unmap.bam", sample=SAMPLES, gen_species=GEN_SPECIES)
ALL_PEAKS = expand("{sample}/macs2_peaks/{sample}_{species}_unmap/{sample}_peaks.narrowPeak",sample=SAMPLES, species=SPECIES)
ALL_BIGWIGS = expand("{sample}/bams/{sample}_{gen_species}.sorted.unmap.bigwig", sample=SAMPLES, gen_species=GEN_SPECIES)
############################################################################
### RDUP BAM PEAKS BIGWIGS
BAMS_RDUP = expand("{sample}/bams/{sample}_{gen_species}.sorted.rdup.bam", sample=SAMPLES, gen_species=GEN_SPECIES)
PEAKS_RDUP = expand("{sample}/macs2_peaks/{sample}_{species}_rdup/{sample}_peaks.narrowPeak",sample=SAMPLES, species=SPECIES)
BIGWIGS_RDUP = expand("{sample}/bams/{sample}_{gen_species}.sorted.rdup.bigwig", sample=SAMPLES, gen_species=GEN_SPECIES)
############################################################################
### FRAGMENT SIZE FILTERED BAM PEAKS BIGWIGS
FragMin = config["CutNRunSpecific"]["MinFragSize"]
FragMax = config["CutNRunSpecific"]["MaxFragSize"]
#FRAGSIZE = [120,150]
FRAGSIZE = [FragMin,FragMax]
############################################################################
### FRAGMENT SIZE FILTERED PEAKS BAM BIGWIGS
ALL_BAMS_frag = expand("{sample}/bams/{sample}_{gen_species}.sorted.unmap_{fragSize}.bam", sample=SAMPLES, gen_species=SPECIES,fragSize=FRAGSIZE)
ALL_PEAKS_frag = expand("{sample}/macs2_peaks/{sample}_{species}_unmap_{fragSize}/{sample}_peaks.narrowPeak",sample=SAMPLES, species=SPECIES,fragSize=FRAGSIZE)
ALL_BIGWIGS_frag = expand("{sample}/bams/{sample}_{gen_species}.sorted.unmap_{fragSize}.bigwig", sample=SAMPLES, gen_species=SPECIES,fragSize=FRAGSIZE)
BAMS_RDUP_frag = expand("{sample}/bams/{sample}_{gen_species}.sorted.rdup_{fragSize}.bam", sample=SAMPLES, gen_species=SPECIES,fragSize=FRAGSIZE)
PEAKS_RDUP_frag = expand("{sample}/macs2_peaks/{sample}_{species}_rdup_{fragSize}/{sample}_peaks.narrowPeak",sample=SAMPLES, species=SPECIES,fragSize=FRAGSIZE)
BIGWIGS_RDUP_frag = expand("{sample}/bams/{sample}_{gen_species}.sorted.rdup_{fragSize}.bigwig", sample=SAMPLES, gen_species=SPECIES,fragSize=FRAGSIZE)
############################################################################
### SPIKE-IN NORMALIZATION
SI_beds=expand("{sample}/bams/{sample}_{spikeIn}_rdup.bed", sample=SAMPLES, spikeIn=SPIKEIN)
Species_beds=expand("{sample}/bams/{sample}_{gen_species}_rdup.bed",sample=SAMPLES, gen_species=[SPECIES])
############################################################################
###
PEAKS_NORM = expand("{sample}/macs2_peaks/{sample}_{species}_rdup_SInorm/{sample}_SInorm_peaks.narrowPeak",sample=SAMPLES, species=SPECIES)
BIGWIGnorm=expand("{sample}/bams/{sample}_{species}.sorted.rdup.SInorm.bigwig", sample=SAMPLES, species=SPECIES)
BIGWIGnorm_RPGC=expand("{sample}/bams/{sample}_{species}.sorted.rdup.SInorm_RPGC.bigwig", sample=SAMPLES, species=SPECIES)
BIGWIGnorm_RPKM=expand("{sample}/bams/{sample}_{species}.sorted.rdup.SInorm_RPKM.bigwig", sample=SAMPLES, species=SPECIES)

BEDGRAPHnorm=expand("{sample}/bams/{sample}_{species}.sorted.rdup.SInorm.bdg",sample=SAMPLES, species=SPECIES)
BEDnorm=expand("{sample}/bams/{sample}_{species}.sorted.rdup.SInorm.bed",sample=SAMPLES, species=SPECIES)
#BEDnorm_wig2bed=expand("{sample}/bams/{sample}_{species}.sorted.rdup.SInorm_wig2bed.bed",sample=SAMPLES, species=SPECIES)
BEDnorm_wig2bed_bed6=expand("{sample}/bams/{sample}_{species}.sorted.rdup.SInorm_wig2_bed6.bed",sample=SAMPLES, species=SPECIES)
############################################################################
### BROAD PEAKS
PEAKS_RDUP_BROAD = expand("{sample}/macs2_peaks/{sample}_{species}_rdup_broad/{sample}_peaks.broadPeak",sample=SAMPLES, species=SPECIES)
homerPeaks=expand("{sample}/homerPeaks/{sample}_{species}.sorted_rdup/{sample}_{species}_peaks.txt",sample=SAMPLES, species=SPECIES)
homerPeaksBED=expand("{sample}/homerPeaks/{sample}_{species}.sorted_rdup/{sample}_{species}_peaks.bed",sample=SAMPLES, species=SPECIES)
############################################################################
###
TagDIR=expand("{sample}/bams/{sample}_{species}.sorted_rdup_Tags/tagInfo.txt",sample=SAMPLES, species=SPECIES)
NORM_TagDIR=expand("{sample}/bams/{sample}_{species}.sorted_rdup_SInorm_Tags/tagInfo.txt",sample=SAMPLES, species=SPECIES)

############################################################################
###
EPIC_PEAKS=expand("{sample}/epic2_peaks/{sample}_{species}_rdup_EPIC/{sample}_peaks_epic.bed",sample=SAMPLES, species=SPECIES)
EPIC_PEAKS_SInorm=expand("{sample}/epic2_peaks/{sample}_{species}_rdup_SInorm_EPIC/{sample}_peaks_epic_SInorm.bed",sample=SAMPLES, species=SPECIES)
##############################################################################################################################
###
rule all:
    input:
        ALL_TRIMMED_FASTQ,
    #    ALL_FASTQC_TRIMMED,
    #    ALL_BAMS,
    #    BAMS_RDUP,
    #    BIGWIGS_RDUP
###############################################################################################################################

#####
include: "cutNrun_Analysis_rules.smk"
