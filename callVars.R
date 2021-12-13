#!/usr/bin/env Rscript
renv::restore(prompt = F)
# renv::snapshot(prompt = F)
# renv::activate()
suppressPackageStartupMessages(library(deepSNV))
suppressPackageStartupMessages(library(stringr))

args = commandArgs(trailingOnly=TRUE)
bam = args[1]
output_directory = args[2]
# bam = "cell_AAAGAACCAATGTGGG-1.bam"
# bam = "cell_AAACCCAAGAAACCAT-1.bam"

cell_label = gsub(".*/", "", gsub("-1\\.bam$", "", bam))
# 16569 from wc -c of MT sequence in hg19 then confirmed by googling and seeing that number again
mtcalls = bam2R(bam,"MT", 1, 16569, q=24, mq=24, keepflag=0)
mtcalls = as.data.frame(mtcalls)

# sum the counts of all A,T,C,G to get the coverage for that position
mtcalls = mtcalls[,c("A","T","C","G")]
mtcalls$coverage = rowSums(mtcalls)
mtcalls$pos = as.integer(rownames(mtcalls))
mtcalls = mtcalls[mtcalls$coverage != 0,]
rownames(mtcalls) = NULL

if (nrow(mtcalls) ==0) {
  print("0 coverage in bam")
  quit(save = "no", status = 0)
}

results_df <- data.frame(Name=character(), Calls=numeric(), Coverage=logical())


for (line_nb in 1:nrow(mtcalls)) {
  line_calls = mtcalls[line_nb,c("A","T","C","G")]

  # filter out NA values and 0 values
  line_calls = line_calls[,names(line_calls)[!is.na(line_calls)], drop=F]
  nonnull_names = names(line_calls)[line_calls != 0]
  line_calls = line_calls[,nonnull_names, drop=F]

  for (colnb in 1:ncol(line_calls)) {
    name = paste0("pos",mtcalls$pos[line_nb],"_alt", colnames(line_calls)[colnb])
    calls = line_calls[,colnb]
    coverage = mtcalls$coverage[line_nb]
    line_scores = data.frame(Name = name, Calls = calls, Coverage = coverage)

    results_df = rbind(results_df, line_scores)
  }

}




coverage_df = results_df
coverage_df$Calls = NULL
rownames(coverage_df) = coverage_df$Name
coverage_df$Name = NULL
colnames(coverage_df) = cell_label

variants_df = results_df
variants_df$Coverage = NULL
rownames(variants_df) = variants_df$Name
variants_df$Name = NULL
colnames(variants_df) = cell_label


dir.create(output_directory, showWarnings = F, recursive = T)

saveRDS(object = coverage_df, file = paste0(output_directory, "/",cell_label, ".coverage.rds"))
saveRDS(object = variants_df, file = paste0(output_directory, "/",cell_label, ".calls.rds"))

