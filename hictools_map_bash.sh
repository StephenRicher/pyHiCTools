#!/usr/bin/env bash

# NEED TO ADD BASH ARGUMENT PARSER SO I can have defaults!!
sample="${1}"
output="${2}"
bt2_index="${3}"
threads="${4}"
samtools="${5}"
bowtie2="${6}"

# Remaining two positional arguments are R1 and R2 respectively
shift 6

# Files saved to current working directory
x=1
for fastq in "${@}"; do

  if [ "${x}" -eq 1 ]; then
    flag="0x40"
    read="R1"
    ((x++))
  else
    flag="0x80"
    read="R2"
  fi

  "${bowtie2}" -x "${bt2_index}" -U "${fastq}" -p "${threads}" --very-sensitive 2> "${sample}"."${read}"_bt2_alignment_stats.txt \
    | awk -v OFS='\t' -v flag="${flag}" '!/^ *@/ && flag == "0x40" {$2 = $2+0x40} !/^ *@/ && flag == "0x80" {$2 = $2+0x80} {print}' \
    | "${samtools}" sort -n -O bam -m 1G -@ "${threads}" > "${sample}"."${read}".sorted.bam

done

# ERROR CHECK - IF BOTH FILES EXIST... THEN
# Can final -b flag to variable either h or b depending on sam_out flag from wrapper python script
"${samtools}" merge -n -@ "${threads}" - "${sample}".R1.sorted.bam "${sample}".R2.sorted.bam \
  | "${samtools}" fixmate -pr -@ ${threads} - - \
  | "${samtools}" view -u -F 8 -q 15 \
  | "${samtools}" fixmate -pm -@ ${threads} - - \
  | "${samtools}" view -b -f 1
