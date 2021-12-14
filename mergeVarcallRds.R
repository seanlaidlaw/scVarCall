#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
rds_directory = args[1]
chunk_nb = args[2]

i = 0
for (table_file in list.files(path = rds_directory, pattern = ".calls.rds$", full.names = T)) {
  if (i == 0) {
    mut_table = readRDS(file = table_file)
  }
  table = readRDS(file = table_file)
  table$rownames = rownames(table)

  mut_table$rownames = rownames(mut_table)
  mut_table = base::merge(x = mut_table, y = table, by = "rownames", all=T)
  rm(table)
  rownames(mut_table) = mut_table$rownames
  mut_table$rownames = NULL
  i = i+1
}
saveRDS(object = mut_table, file = paste0(rds_directory, "/chunk_", chunk_nb, ".calls.rds"))


i = 0
for (table_file in list.files(path = rds_directory, pattern = ".coverage.rds$", full.names = T)) {
  if (i == 0) {
    mut_table = readRDS(file = table_file)
  }
  table = readRDS(file = table_file)
  table$rownames = rownames(table)

  mut_table$rownames = rownames(mut_table)
  mut_table = base::merge(x = mut_table, y = table, by = "rownames", all=T)
  rm(table)
  rownames(mut_table) = mut_table$rownames
  mut_table$rownames = NULL
  i = i+1
}
saveRDS(object = mut_table, file = paste0(rds_directory, "/chunk_", chunk_nb, ".coverage.rds"))

