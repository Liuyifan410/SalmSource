# Model 1 Application Script
library(tidymodels)
library(tidyverse)
library(lightgbm)
library(bonsai)
library(seqinr)
library(ranger)

# Read command line arguments
allArgs <- commandArgs(trailingOnly = FALSE)
input_genome <- allArgs[6]    # Input genome file
output_name <- allArgs[7]     # Output file name prefix
currentPath=dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
modelPath=paste(currentPath,"database","model1.Rdata",sep="/")
diamondPath=paste(currentPath,"database","nr1.dmnd",sep="/")

# Load the first model training results
load(modelPath)  # Contains workflow_best, feature_best

# Step 1: Gene prediction
prodigal_code <- sprintf("prodigal -a %s_groups.pep -d %s_groups.cds -g 11 -p single -i %s", 
                         output_name, output_name, input_genome)
system(prodigal_code)

# Step 2: Compare with pangenome database
diamond_blast <- sprintf("diamond blastp --db %s -q %s_groups.pep -o %s_groups_blast.out --outfmt 6 qseqid sseqid pident qcovhsp scovhsp --max-target-seqs 10000 -e 1e-5",
                         diamondPath, output_name, output_name)
system(diamond_blast)

# Step 3: Process comparison results and build feature vectors
blast_results <- read_tsv(sprintf("%s_groups_blast.out", output_name), 
                          col_names = c("qseqid", "sseqid", "pident", "qcovhsp", "scovhsp"))

# Filter high-quality matches
high_quality_hits <- blast_results %>%
  filter(qcovhsp >= 80, scovhsp >= 80, pident >= 80) %>%
  distinct(sseqid)  # Keep only one representative per gene cluster

# Initialize all feature vectors - for all feature types
all_features <- feature_best
feature_vector <- rep(0, length(all_features))
names(feature_vector) <- all_features

# Mark detected genes
detected_genes <- unique(high_quality_hits$sseqid)
detected_genes_clean <- gsub("\\..*", "", detected_genes)  # Remove possible version numbers

# Match feature names - including all types of features
matched_features <- intersect(names(feature_vector), detected_genes_clean)
feature_vector[matched_features] <- 1

# Create final feature matrix
final_feature_matrix <- as.data.frame(t(feature_vector[feature_best]))
colnames(final_feature_matrix) <- feature_best

# Ensure all features are numeric
final_feature_matrix_numeric <- as.data.frame(lapply(final_feature_matrix, as.numeric))

# Model prediction
predictions <- augment(workflow_best, final_feature_matrix_numeric)

# Output prediction results - adjusted according to actual response levels
final_output <- data.frame(
  Sample = output_name,
  predicted_class = predictions$.pred_class,
  probability_Yes = ifelse(".pred_Yes" %in% colnames(predictions), 
                           predictions$.pred_Yes, NA),
  probability_No = ifelse(".pred_No" %in% colnames(predictions), 
                          predictions$.pred_No, NA)
)


write_tsv(final_output, sprintf("%s_groups_host_prediction.txt", output_name))

detected_features <- names(feature_vector)[feature_vector == 1]

detected_features_df <- data.frame(
  feature_id = detected_features,
  status = "detected"
)

write_tsv(detected_features_df, sprintf("%s_detected_features.txt", output_name))



















