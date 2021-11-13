#!/usr/bin/env bash

deduped_bam_file="$1"
barcode_file="$2"
output_dir="$3"

[ ! -d "$output_dir" ] && mkdir -p "$output_dir"

# only get line for current bsub index
barcode="$(sed -n -e "${LSB_JOBINDEX}p" "$barcode_file" | sed 's/CB:Z://g')"

# create barcode file
echo "$barcode" > "$output_dir/cell.$barcode.txt"
barcode_file="$output_dir/cell.$barcode.txt"

subset-bam --cores 12 \
  --bam "$deduped_bam_file" \
  --cell-barcodes "$barcode_file" \
  --out-bam "$output_dir/cell_{$barcode}.bam"
