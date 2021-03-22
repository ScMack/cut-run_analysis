
     python2 /tools/DownloadedScripts/bed_to_gff_converter_modified.py ${PEAKS_rdup} ${GFFfile_rdup}
     python2 /tools/rose/ROSE_main.py -r ${BAMfile_rdup} -c ${CONTROL_BAMfile} -i ${GFFfile_rdup} -g ${genome} -o ${WORKDIR}/${sample}/${ROSEoutDIR}/${sample}_${genome}_rdup_rose -t 2500
