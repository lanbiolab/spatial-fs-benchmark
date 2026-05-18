#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args) >= 1) args[[1]] else "data/resources/human_tfs.tsv"
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

extract_tf_genes <- function(org_pkg, species_label) {
    suppressPackageStartupMessages(library(AnnotationDbi))
    suppressPackageStartupMessages(library(org_pkg, character.only = TRUE))
    pkg_ns <- asNamespace(org_pkg)
    dbfile_fun <- get(paste0(sub("\\.db$", "", org_pkg), "_dbfile"), envir = pkg_ns)
    annotation_db <- AnnotationDbi::loadDb(dbfile_fun())
    tf_goids <- c(
        "GO:0000981",
        "GO:0000977",
        "GO:0000978",
        "GO:0001071",
        "GO:0001077",
        "GO:0001078",
        "GO:0001216",
        "GO:0001228",
        "GO:0003700",
        "GO:0140110"
    )
    annotations <- AnnotationDbi::select(
        annotation_db,
        keys = keys(annotation_db, keytype = "ENTREZID"),
        columns = c("ENSEMBL", "SYMBOL", "GOALL", "ONTOLOGYALL"),
        keytype = "ENTREZID"
    )
    annotations <- annotations[
        !is.na(annotations$ENSEMBL) &
            !is.na(annotations$SYMBOL) &
            annotations$ONTOLOGYALL == "MF" &
            annotations$GOALL %in% tf_goids,
        c("ENSEMBL", "SYMBOL")
    ]
    annotations <- unique(annotations)
    colnames(annotations) <- c("ENSEMBL", "Gene")
    annotations$Species <- species_label
    annotations[order(annotations$Gene), , drop = FALSE]
}

human_tfs <- extract_tf_genes("org.Hs.eg.db", "Human")
mouse_tfs <- extract_tf_genes("org.Mm.eg.db", "Mouse")
tfs <- unique(rbind(human_tfs, mouse_tfs))
write.table(tfs, file = out_file, quote = FALSE, sep = "\t", row.names = FALSE)
