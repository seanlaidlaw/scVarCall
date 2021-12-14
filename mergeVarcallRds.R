#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
rds_directory <- args[1]
chunk_nb <- args[2]

merge_rds <- function(file_list) {
    i <- 0
    for (table_file in file_list) {
        if (i == 0) {
            mut_table <- readRDS(file = table_file)
        } else {
            table <- readRDS(file = table_file)
            table$rownames <- rownames(table)

            mut_table$rownames <- rownames(mut_table)
            mut_table <- base::merge(x = mut_table, y = table, by = "rownames", all = T)
            rm(table)
            rownames(mut_table) <- mut_table$rownames
            mut_table$rownames <- NULL
        }
        i <- i + 1
    }
    return(mut_table)
}

calls_rds_list <- list.files(path = rds_directory, pattern = ".calls.rds$", full.names = T)
merged_calls_table <- merge_rds(calls_rds_list)
saveRDS(object = merged_calls_table, file = paste0(rds_directory, "/chunk_", chunk_nb, ".calls.rds"))

coverage_rds_list <- list.files(path = rds_directory, pattern = ".coverage.rds$", full.names = T)
merged_coverage_table <- merge_rds(coverage_rds_list)
saveRDS(object = merged_coverage_table, file = paste0(rds_directory, "/chunk_", chunk_nb, ".coverage.rds"))
