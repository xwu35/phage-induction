#!/usr/bin/env bash

set -euo pipefail

# set default values
bam_dir=""
extension=".bam"
output="read_counts.txt"

# usage function
usage() {
  echo "Usage: $0 -i bam_dir -o output" 
  echo "Options:"
  echo "  -h          Display this help message."
  echo "  -i          Directory containning the bam files [required]" 
  echo "  -e          extension of the bam files (default: .bam)"
  echo "  -o <file>   Specify an output (default: read_counts.txt)" 
 
  exit 1
}

# the OPTSTRING is defined as :i:e:o:h. This tells the script to expect two options, -i, -e and -o, all of them require an argument.
OPTSTRING=":i:e:o:h"

# parse the options
while getopts ${OPTSTRING} opt; do
  case ${opt} in
    i) bam_dir="${OPTARG}" ;;
    e) extension="${OPTARG}" ;;
    o) output="${OPTARG}" ;;
    h) usage ;;
    :) echo "Option -${OPTARG} requires an argument."; exit 1 ;;
    ?) echo "Invalid option: -${OPTARG}."; exit 1 ;;
  esac
done

# count primary alignments
echo "sample read_counts" > "$output"
for file in $bam_dir/*/*$extension; do
    count=$(samtools view -c -F 0x900 "$file")
    echo "$file $count" >> $output
done
