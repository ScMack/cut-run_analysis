#!/bin/bash

while getopts b:s:d: option
do
  case "${option}"
  in
  b) baseDIR=${OPTARG};;
  s) sampleSheet=${OPTARG};;
  \? ) echo "$0 -b (baseDirectory) -s (Samplesheet) -d (dryrun?)";;
  [?]) echo "I nvalid option: ${OPTARG} requires an argument";;
  esac
done



while IFS= read -r line
do
	# echo $line
	headerCheck="$(echo $line| awk '{print $1}')"

	if [ ${headerCheck} != 'SampleID' ]; then

		OldFile="$(echo $line| awk '{print $1}')"
		OldFilePATH="$(echo $line| awk '{print $3}')"
		NewFile="$(echo $line|awk '{print $2}')"
		echo $OldFile $NewFile $OldFilePATH


		Merge="$(echo $line| awk '{print $4}')"
		mkdir -p ${NewFile}
		mkdir -p ${NewFile}/fastqs

		if [ ${Merge} == FALSE ]; then
			
			ln -s ${OldFilePATH}/${OldFile}_R1*.fastq.gz ${baseDIR}/${NewFile}/fastqs/${NewFile}_R1.fastq.gz
			ln -s ${OldFilePATH}/${OldFile}_R2*.fastq.gz ${baseDIR}/${NewFile}/fastqs/${NewFile}_R2.fastq.gz

		elif [ ${Merge} == TRUE ];
		then
			LaneNUM=$(cat <(find ${OldFilePATH}/${OldFile}_*R1*fastq.gz | sort |wc -l))
#
			echo "number of lanes= "$LaneNUM
					echo 'ls ${OldFilePATH}/${OldFile}_*R1*fastq.gz | sort | tr "\n" " "| awk '{ print "cat",$0}')>${baseDIR}/${NewFile}/fastqs/${NewFile}_R1.fastq.gz'


			$(ls ${OldFilePATH}/${OldFile}_*R1*fastq.gz | sort | tr "\n" " "| awk '{ print "cat",$0}')>${baseDIR}/${NewFile}/fastqs/${NewFile}_R1.fastq.gz
			$(ls ${OldFilePATH}/${OldFile}_*R2*fastq.gz | sort | tr "\n" " "| awk '{ print "cat",$0}')>${baseDIR}/${NewFile}/fastqs/${NewFile}_R2.fastq.gz

		fi
	fi

done < "$sampleSheet"
