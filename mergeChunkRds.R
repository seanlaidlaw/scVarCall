#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
rds_directories <- args

calls_rds <- list.files(
    path = rds_directories,
    pattern = ".calls.rds$",
    full.names = T
)

if (length(calls_rds) < 1) {
    stop("Unable to find the requested calls.rds files.
		 Verify youre providing the directories containing the calls.rds")
}

coverage_rds <- list.files(
    path = rds_directories,
    pattern = ".coverage.rds$",
    full.names = T
)

if (length(coverage_rds) < 1) {
    stop("Unable to find the requested coverage.rds files.
		 Verify youre providing the directories containing the coverage.rds")
}


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

merged_calls_table <- merge_rds(calls_rds)
saveRDS(object = merged_calls_table, file = paste0("merged_chunks_calls.rds"))

merged_coverage_table <- merge_rds(coverage_rds)
saveRDS(object = merged_coverage_table, file = paste0("merged_chunks_coverage.rds"))
