#!/usr/bin/env Rscript
# renv::restore(prompt = F)
# renv::snapshot(prompt = F)
# renv::activate()
suppressPackageStartupMessages(library(deepSNV))
suppressPackageStartupMessages(library(stringr))

args = commandArgs(trailingOnly=TRUE)
bam = args[1]
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

if (nrow(mtcalls) ==0) {
  stop("0 coverage in bam")
}

mtcalls$value_max = NA
mtcalls$value_secondmax = NA
mtcalls$secondmax = NA



for (line_nb in 1:nrow(mtcalls)) {
  line_calls = mtcalls[line_nb,c("A","T","C","G")]
  # set maximum value to NA to be removed later
  mtcalls$value_max[line_nb] = max(line_calls)
  line_calls[line_calls == max(line_calls)] = NA

  # filter out NA values and 0 values
  line_calls = line_calls[,names(line_calls)[!is.na(line_calls)], drop=F]

  nonnull_names = names(line_calls)[line_calls != 0]
  line_calls = line_calls[,nonnull_names, drop=F]

  # if any values left, take the max
  if (length(nonnull_names) > 0) {
    # get column with max value of what is left (so the second max)
    mtcalls$value_secondmax[line_nb] = max(line_calls)

    secondmax = colnames(line_calls)[which(line_calls == max(line_calls))]
    secondmax_str = paste(secondmax, collapse = "")
    mtcalls$secondmax[line_nb] = secondmax_str
  }
}

# remove rows that have no second-max
mtcalls = mtcalls[!is.na(mtcalls$secondmax),]
# remove columns we no longer need
mtcalls = mtcalls[,c("coverage","pos","value_secondmax","secondmax")]


# subset only those with more than one secondmax
multiple_secondmax = mtcalls[nchar(mtcalls$secondmax) > 1,]
# remove rows with multiple bases in second place
mtcalls = mtcalls[nchar(mtcalls$secondmax) == 1,]

# split multiple secondmax into multiple rows reflecting each base thats secondmax
if (nrow(multiple_secondmax > 0)) {
  for (row_nb in 1:nrow(multiple_secondmax)) {
    pos = multiple_secondmax$pos[row_nb]
    cov = multiple_secondmax$coverage[row_nb]
    secondmax = multiple_secondmax$value_secondmax[row_nb]

    secondmaxes = str_split(multiple_secondmax$secondmax[row_nb], pattern = "")[[1]]
    for (string in secondmaxes) {
      mtcalls = rbind(mtcalls, list(pos=pos, coverage=cov, value_secondmax = secondmax, secondmax=string))
    }
  }
}

rownames(mtcalls) = NULL
rownames(mtcalls) = paste0("pos",mtcalls$pos,"_alt",mtcalls$secondmax)
mtcalls$pos = NULL
mtcalls$secondmax = NULL

coverage_df = mtcalls
coverage_df$value_secondmax = NULL
colnames(coverage_df) = cell_label

variants_df = mtcalls
variants_df$coverage = NULL
colnames(variants_df) = cell_label

saveRDS(coverage_df, file = paste0("output_dir/rds/",cell_label, ".rds"))
saveRDS(variants_df, file = paste0("output_dir/rds/",cell_label, ".rds"))

