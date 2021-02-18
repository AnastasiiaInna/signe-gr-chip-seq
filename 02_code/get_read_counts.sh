#!/bin/bash

# Funtion to count number of reads in BAM files

# Param 1 : text file with abolute paths to BAM files
# Param 2 : output file name

if [ ${#@} -ne 2 ]; then
      printf "\nUsage: %s Not enough arguments\n" "${0}" > /dev/stderr;
      exit 1;
fi

BAM_LIST_FILENAME="${1}"
OUTPUT_FILENAME="${2}"

for FILE in $(cat "${BAM_LIST_FILENAME}") ; do
   samtools view -c ${FILE} >> $OUTPUT_FILENAME;
done
