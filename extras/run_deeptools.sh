computeMatrix reference-point \
-S ${BIGWIGfiles} \
-R ${BEDfile} --outFileName ${FILEname}.gz \
--samplesLabel ${SAMPLEnames} \
--smartLabels --referencePoint center \
--beforeRegionStartLength 2000 --afterRegionStartLength 2000 \
--outFileNameMatrix ${FILEname}_heatMat.tab -p 30


plotHeatmap \
--matrixFile ${FILEname}.gz \
--startLabel "Start" --endLabel "End" --dpi 300 \
--outFileName ${FILEname}.pdf \
--outFileSortedRegion ${FILEname}.txt \
--colorMap ${COLORmap} \
--plotType "fill" --dpi 300 --zMax ${maxVal} --yMax ${maxVal}
