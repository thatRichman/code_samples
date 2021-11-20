set.seed(20201109)

bioc.pkgs <- c("DESeq2","Rsamtools","GenomicFeatures","apeglm")
pkgs <- c("magrittr","openxlsx","tidyverse")
for(pkg in bioc.pkgs){
  if(!pkg %in% rownames(installed.packages())){
    BiocManager::install(pkg, ask = FALSE)
  }
}

for(pkg in pkgs){
  if(!pkg %in% rownames(installed.packages())){
    install.packages(pkg, quiet = T)
  }
}

library(tidyverse)
library(magrittr)
library(openxlsx)
library(DESeq2)
library(Rsamtools)
library(GenomicFeatures)
library(apeglm)
select <- dplyr::select


make_counts_matrix <- function(tab.dir, design.tbl.path, strand.index){
  ff <- list.files(path = tab.dir, pattern = "*ReadsPerGene.out.tab$", full.names = TRUE)
  tab.files <- lapply(ff, read.table, skip = 4)
  counts <- as.data.frame(sapply(tab.files, function(x) x[,strand.index]))
  ff <- gsub("[.]ReadsPerGene[.]out[.]tab", "", ff)
  ff <- gsub("[.]/counts/", "", ff)
  colnames(counts) <- ff
  colnames(counts) <- str_extract(colnames(counts),"(?<=//).*(?=_qc)")
  row.names(counts) <- tab.files[[1]]$V1
  
  design.table <- read.table(file = design.tbl.path, header = T,sep=",")
  rownames(design.table) <- design.table$sample
  design.table$group <- paste(design.table$treatment, design.table$time) # treatment and time become single variable
  # Remove extraneous samples from counts table
  counts.fixed <- counts[,which(colnames(counts) %in% rownames(design.table))]
  design.table <- design.table[colnames(counts.fixed),]
  
  return(list(counts.fixed, design.table))
}


run_deseq <- function(dir, 
                      design.tbl.path, 
                      strand.index = 2, 
                      design = "~ group",
                      betaPrior = FALSE){
  design.mat <- make_counts_matrix(dir, design.tbl.path, strand.index)
  counts.matrix <- design.mat[[1]]
  design.table <- design.mat[[2]]
  
  dds <- DESeqDataSetFromMatrix(countData = counts.matrix,
                                colData = design.table,
                                design = as.formula(design))
  dds$treatment <- as.factor(dds$treatment)
  dds$group <- as.factor(dds$group)
  dds$treatment <- relevel(dds$treatment, ref = "REDACTED")
  dds$group <- relevel(dds$group, ref="REDACTED 24h 0.1")
  dds <- DESeq(dds, betaPrior = betaPrior)
  REDACTED <- unique(design.table$treatment)  # some REDACTED aren't present in every cell line, infer from RowData
  res <- list(dds, REDACTED) 
  names(res) <- c("dds","REDACTED")
  return(res)
}


# contrasts should be supplied as a list, with the key being one of:
# c("24h","48h","72h")
make_results_table <- function(dds, contrasts, goi, save.raw = FALSE, padj.thresh = 0.05){
  res.tbls <- list()
  for(tp in names(contrasts)){
    print(contrasts[[tp]])
    res <- results(dds, contrast = contrasts[[tp]])
    if(save.raw==FALSE){
      message("save.raw = FALSE, ignoring padj.thresh")
      res.df <- dplyr::select(dplyr::filter(rownames_to_column(as.data.frame(res)), 
                                            rowname %in% goi), 
                              log2FoldChange, padj) 
    } else {
      message("save.raw = TRUE, ignoring goi")
      if(is.null(padj.thresh)){
        message("NULL padj threshold, outputting all values")
        res.df <- rownames_to_column(as.data.frame(res))
      } else {
        message("constraining to padj threshold")
        res.df <- subset(rownames_to_column(as.data.frame(res)), padj < padj.thresh)
      }
    }
    colnames(res.df) <- c("gene", colnames(res.df)[-1])
    res.tbls[[tp]] <- res.df
  }
  res.out <- purrr::reduce(res.tbls, full_join, by="gene")
  return(res.out)
}

make_contrast_list <- function(dds, REDACTED, REDACTED){
  contrast.list <- list()
  for(REDACTED in REDACTED){
    contrast.list[[REDACTED]] = list()
    for(REDACTED in REDACTED){
      contrast.list[[REDACTED]][[REDACTED]] <- list("24h" = c("group",paste0(REDACTED, " 24h ", REDACTED), "REDACTED 24h 0.1"),
                                    "48h" = c("group",paste0(REDACTED, " 48h ", REDACTED), "REDACTED 48h 0.1"),
                                    "72h" = c("group", paste0(REDACTED, " 72h ", REDACTED), "REDACTED 72h 0.1") 
                                    )
    }
   
  }
  return(contrast.list)
} 


make_output_xlsx <- function(file, res.tbls, cell.line){
  wb <- createWorkbook(title = cell.line)
  for(name in names(res.tbls)){
    trunc = str_trunc(name, 31, "right")
    addWorksheet(wb, sheetName = trunc)
    writeData(wb, sheet = trunc, x = res.tbls[[name]])
  }
  saveWorkbook(wb, file=file, overwrite=FALSE)
}

# Each result_file for a cell line is an xlsx file
# Each sheet is a REDACTED, compared against the REDACTED control at appropriate timepoint
# columns are
# L2FC_24h | padj_24h | L2FC_48h | padj_48h | L2FC_72h | padj_72h

tab.root <- c("REDACTED")

cell.lines <- c("REDACTED","REDACTED", "REDACTED","REDACTED","REDACTED-18")

# Don't use REDACTED list, because not all REDACTED are in each REDACTED. 
# Instead, it will infer REDACTED automatically
# REDACTED = c('REDACTED')

genes <- c('REDACTED')


all.dds <- list()
for(line in cell.lines){
  dds_res <- run_deseq(dir = paste0(tab.root, line,"REDACTED/"), 
                       design.tbl.path = paste0(tab.root, line, "REDACTED"),
                       strand.index = 2, 
                       design = "~ REDACTED + REDACTED", # if you change the design you also have to change the make_contrast_list function!
                       betaPrior = FALSE)
  dds <- dds_res[["dds"]]
  REDACTED <- dds_res[["REDACTED"]]
  contrasts.list <- make_contrast_list(dds, REDACTED = REDACTED)
  
  results.tbls <- list()
  for(REDACTED in REDACTED){
    if(REDACTED != "REDACTED"){
      results.tbls[[REDACTED]] <- make_results_table(dds, contrasts.list[[REDACTED]], goi=genes, save.raw = TRUE, padj.thresh = NULL)
    }
  }
  
  make_output_xlsx(file = paste0(date(),"_",line,"REDACTED"), 
                   res.tbls = results.tbls, 
                   cell.line = line)
  all.dds[[line]] <- dds
}


save.image(file="REDACTED")


tab.root <- c("REDACTED")
for(line in cell.lines){
  dir = paste0(tab.root, line,"_gene_tabs/")
  design.tbl.path = paste0(tab.root, line, "REDACTED")
  strand.index = 2
  design.mat <- make_counts_matrix(dir, design.tbl.path, strand.index)
  counts.matrix <- design.mat[[1]]
  design.table <- design.mat[[2]]
  write.table(file = paste0("REDACTED",line,"REDACTED"), quote = F, sep = ",", row.names = F, x = rownames_to_column(counts.matrix))
  write.table(file = paste0("REDACTED",line,"REDACTED"), quote = F, sep = ",", row.names = F, x =design.table)
  
}
