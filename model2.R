# Model 2 Application Script
library(tidymodels)
library(tidyverse)
library(lightgbm)
library(bonsai)
library(seqinr)

# Read command line arguments
allArgs <- commandArgs(trailingOnly = FALSE)
input_genome <- allArgs[6]    # Input genome file
output_name <- allArgs[7]     # Output file name prefix
#reference_genome <- "GCF_000006945.2_ASM694v2_genomic.gbff"  # Reference genome used in training

currentPath=dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
modelPath=paste(currentPath,"database","model2.Rdata",sep="/")
diamondPath=paste(currentPath,"database","nr2.dmnd",sep="/")
reference_genome=paste(currentPath,"database","GCF_000006945.2_ASM694v2_genomic.gbff",sep="/")
# Gene prediction
prodigal_code <- sprintf("prodigal -a %s_target.pep -d %s_target.cds -g 11 -p single -i %s", 
                         output_name, output_name, input_genome)
system(prodigal_code)

# DIAMOND blastp comparison
diamond_blast <- sprintf("diamond blastp --db %s -q %s_groups.pep -o %s_groups_blast.out --outfmt 6 qseqid sseqid pident qcovhsp scovhsp --max-target-seqs 10000 -e 1e-5",
                         diamondPath, output_name, output_name)
system(diamond_blast)

# SNP detection with snippy
snippy_code <- sprintf("snippy --outdir %s_snippy --ref %s --ctgs %s --cpus 6",
                       output_name, reference_genome, input_genome)
system(snippy_code)

# Load the second model training results
load(modelPath)  # Contains workflow_best, feature_best

# Process blast results
blast_results <- read_tsv(sprintf("%s_blast.out", output_name), 
                          col_names = c("qseqid", "sseqid", "pident", "qcovhsp", "scovhsp"))

# Filter high-quality matches
high_quality_hits <- blast_results %>%
  filter(qcovhsp >= 80, scovhsp >= 80, pident >= 80) %>%
  group_by(qseqid) %>%
  filter(pident == max(pident)) %>%
  ungroup()
  
# 关键修正：直接使用blast结果中的sseqid构建特征名
detected_genes <- unique(high_quality_hits$sseqid)
detected_gene_features <- paste0("gene_", detected_genes)

# 初始化基因特征向量
gene_features <- feature_best[grepl("^gene_", feature_best)]
gene_vector <- rep(0, length(gene_features))
names(gene_vector) <- gene_features

# 标记检测到的基因 - 直接匹配
matched_genes <- intersect(names(gene_vector), detected_gene_features)
gene_vector[matched_genes] <- 1
  
# Process SNP features - adjust according to actual feature format
snp_features <- feature_best[grepl("^snp_", feature_best)]
snp_vector <- rep(0, length(snp_features))
names(snp_vector) <- snp_features

# Read snippy's snps.tab file
snp_tab_file <- sprintf("%s_snippy/snps.tab", output_name)
snp_tab <- read_tsv(snp_tab_file, col_names = TRUE, na = c("", "NA"), show_col_types = FALSE)

# Process SNPs - build according to actual feature name format
snp_positions <- snp_tab %>%
  filter(TYPE == "snp") %>%
  # Build according to actual feature name format (e.g., snp_pos_1679627)
  mutate(feature_name = paste0("snp_pos_", POS)) %>%
  pull(feature_name)

# Mark detected SNPs
detected_snps <- intersect(names(snp_vector), snp_positions)
snp_vector[detected_snps] <- 1

# Process Indel features - adjust according to actual feature format
indel_features <- feature_best[grepl("^indel_", feature_best)]
indel_vector <- rep(0, length(indel_features))
names(indel_vector) <- indel_features

indel_data <- snp_tab %>%
  filter(TYPE %in% c("del", "ins")) %>%
  mutate(feature_name = paste0("indel_", CHROM, "__", toupper(TYPE), "__", POS, "__", REF, "__", ALT)) %>%
  pull(feature_name)

detected_indels <- intersect(names(indel_vector), indel_data)
indel_vector[detected_indels] <- 1

# Combine all features
all_features <- c(gene_vector, snp_vector, indel_vector)
final_feature_matrix <- as.data.frame(t(all_features[feature_best]))
colnames(final_feature_matrix) <- feature_best
final_feature_matrix_numeric <- as.data.frame(lapply(final_feature_matrix, as.numeric))

# Model prediction
predictions <- augment(workflow_best, final_feature_matrix_numeric)

# Output prediction results
final_output <- data.frame(
  Sample = output_name,
  predicted_host = predictions$.pred_class,
  probability_chicken = ifelse(".pred_chicken" %in% colnames(predictions), 
                               predictions$.pred_chicken, NA),
  probability_pig = ifelse(".pred_pig" %in% colnames(predictions), 
                           predictions$.pred_pig, NA)
)
write_tsv(final_output, sprintf("%s_host_prediction.txt", output_name))

# Output detected feature details
feature_report <- data.frame(
  feature_type = c(
    rep("gene", sum(gene_vector == 1)),
    rep("snp", sum(snp_vector == 1)), 
    rep("indel", sum(indel_vector == 1))
  ),
  feature_name = c(
    names(gene_vector)[gene_vector == 1],
    names(snp_vector)[snp_vector == 1],
    names(indel_vector)[indel_vector == 1]
  )
)
write_tsv(feature_report, sprintf("%s_detected_features.txt", output_name))

cat("Analysis completed! Prediction results saved in", sprintf("%s_host_prediction.txt", output_name), "\n")