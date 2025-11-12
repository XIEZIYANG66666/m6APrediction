#' DNA Sequence Encoding for 5-mer
#'
#' @param dna_strings A character vector of DNA 5-mer sequences.
#' @return A data frame with encoded nucleotide positions.
#' @examples
#' if (FALSE) {
#'   m6APrediction:::dna_encoding(c("ATGAT", "CGTAC"))
#' }
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}


#' Multiple Samples m6A Prediction
#'
#' This function predicts m6A probability and status for multiple samples using
#' a trained machine learning model. It processes DNA 5-mer features, combines
#' them with other features, and returns a data frame with predictions.
#'
#' @param ml_fit A trained machine learning model (e.g., from randomForest).
#' @param feature_df A data frame containing features for multiple samples.
#'   Required columns: gc_content, RNA_type, RNA_region, exon_length,
#'   distance_to_junction, evolutionary_conservation, DNA_5mer.
#' @param positive_threshold A numeric value (0-1, default 0.5) defining the
#'   threshold for classifying predictions as "Positive".
#' @return A data frame with original features plus two new columns:
#'   predicted_m6A_prob (predicted probability) and predicted_m6A_status
#'   ("Positive" or "Negative").
#' @examples
#'
#' library(m6APrediction)
#'
#'
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' feature_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
#'
#'
#' prediction_multiple(ml_fit, feature_df, positive_threshold = 0.5)
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  required_cols <- c("gc_content", "RNA_type", "RNA_region", "exon_length",
                     "distance_to_junction", "evolutionary_conservation", "DNA_5mer")
  stopifnot(all(required_cols %in% colnames(feature_df)))

  feature_df$RNA_type <- factor(
    feature_df$RNA_type,
    levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")
  )
  feature_df$RNA_region <- factor(
    feature_df$RNA_region,
    levels = c("CDS", "intron", "3'UTR", "5'UTR")
  )

  dna_features <- dna_encoding(feature_df$DNA_5mer)
  all_features <- cbind(feature_df[, !colnames(feature_df) %in% "DNA_5mer"], dna_features)

  feature_df$predicted_m6A_prob <- predict(ml_fit, all_features, type = "prob")[, "Positive"]
  feature_df$predicted_m6A_status <- ifelse(
    feature_df$predicted_m6A_prob >= positive_threshold,
    "Positive",
    "Negative"
  )

  return(feature_df)
}


#' Single Sample Prediction for m6A Sites
#'
#' This function takes individual feature values and a trained machine learning model to predict m6A probability and status for a single sample. It returns the predicted values as a named vector.
#'
#' @param ml_fit A trained machine learning model (e.g., randomForest model).
#' @param gc_content A numeric value for the GC content.
#' @param RNA_type A character string specifying the RNA type (e.g., "mRNA", "lincRNA").
#' @param RNA_region A character string specifying the RNA region (e.g., "CDS", "intron").
#' @param exon_length A numeric value for the exon length.
#' @param distance_to_junction A numeric value for the distance to junction.
#' @param evolutionary_conservation A numeric value for evolutionary conservation.
#' @param DNA_5mer A character string representing the DNA 5-mer sequence.
#' @param positive_threshold A numeric value between 0 and 1. If the predicted probability is greater than this threshold, the status is "Positive"; otherwise, "Negative".
#' @return A named vector with two elements: "predicted_m6A_prob" (predicted probability) and "predicted_m6A_status" (either "Positive" or "Negative").
#' @examples
#'
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#'
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA",
#'                   RNA_region = "CDS", exon_length = 12,
#'                   distance_to_junction = 50, evolutionary_conservation = 0.8,
#'                   DNA_5mer = "ATGAT", positive_threshold = 0.5)
#' @export
#' @import randomForest
#' @importFrom stats predict
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  single_feature <- data.frame(
    gc_content = as.numeric(gc_content),
    RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
    RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
    exon_length = as.numeric(exon_length),
    distance_to_junction = as.numeric(distance_to_junction),
    evolutionary_conservation = as.numeric(evolutionary_conservation),
    DNA_5mer = as.character(DNA_5mer),
    stringsAsFactors = FALSE
  )

  dna_features <- dna_encoding(single_feature$DNA_5mer)
  all_features <- cbind(single_feature[, !colnames(single_feature) %in% "DNA_5mer"], dna_features)

  predicted_prob <- predict(ml_fit, all_features, type = "prob")[, "Positive"]
  predicted_status <- ifelse(predicted_prob >= positive_threshold, "Positive", "Negative")

  returned_vector <- c(predicted_m6A_prob = predicted_prob, predicted_m6A_status = predicted_status)

  return(returned_vector) #return a named vector with values for predicted m6A prob and predicted m6A status
}

