peaks_SF=sample_peaks.bed
peaks_summits=sample_summits.bed


/tools/bedtools2/bin/slopBed -i ${peaks_summits} -g /Databases/mm10/mm10.chrom.sizes -b 500 >sample_summits_slop.bed

 create a fasta file for the motif analysis
/tools/bedtools2/bin/bedtools getfasta -fi /Databases/mm10/mm10.fa -bed sample_summits_slop.bed -fo sample_summits_slop.fa


SLOPfa=sample_summits_slop.fa


/tools/meme-5.1.0/bin/meme-chip -ccut 100 -dna -meme-mod anr -meme-minw 5 -meme-maxw 15 \
-db /Databases/motifs/JASPAR2020_CORE_vertebrates_non-redundant_pfms.meme \
-meme-nmotifs 50 -dreme-e 0.05 -oc sample_memechip_JASPAR2020vertebrates ${SLOPfa} -meme-p 50 &

