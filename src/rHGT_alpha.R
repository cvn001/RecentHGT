#' This script is used to detect recent horizontal gene transfers.
#' @param collections_dir, should be the address of the input directory, string.
#' @param result_file, should be the address of the result file, string.
#' @param param_min, should be a value of lower limit of the distribuiton.
#' @param param_max, should be a value of upper limit of the distribuiton.
#' @export a text file with the number of rHGTs between every strain pair.

library(fitdistrplus)

# fit function use fitdist method
fit_distribution <- function(data, fit_method, fit_gof) {
  fit_wei <- fitdist(data, "weibull", method = fit_method, gof = fit_gof)
  stat <- gofstat(fit_wei)
  aic <- as.numeric(stat$aic)
  shape <- as.numeric(fit_wei$estimate['shape'])
  scale <- as.numeric(fit_wei$estimate['scale'])
  dist_pars <- c(shape, scale)
  return(list(dist_pars, stat, aic))
}
em_implementation <- function(data_file, trim_min, trim_max, last_bin) {
  # load and trim data
  my_data <- read.csv(data_file, sep = '\t', header = T)
  trim_data <- my_data$Identity[my_data$Identity >= trim_min & my_data$Identity < trim_max]
  last_bin_data <- my_data$Identity[my_data$Identity >= last_bin[1] & my_data$Identity <= last_bin[2]]
  # EM implementation
  m <- c(0)
  fit_method <- 'mge'
  all_gofs <- c('ADR', 'AD2R')
  max_steps <- 100
  estimate_data <- trim_data
  all_aics <- c()
  all_bin_nums <- c()
  for (fit_gof in all_gofs) {
    for (i in 1:max_steps) {
      estimate_result <- fit_distribution(estimate_data, fit_method, fit_gof)
      estimate_pars <- estimate_result[[1]]
      estimate_stat <- estimate_result[[2]]
      estimate_aic <- estimate_result[[3]]
      estimate_shape <- estimate_pars[1]
      estimate_scale <- estimate_pars[2]
      prob <- abs(exp(-(last_bin[1] / estimate_scale) ^ estimate_shape) - exp(-(last_bin[2] / estimate_scale) ^ estimate_shape))
      m <- c(m, prob)
      estimate_last_bin_num <- length(trim_data) * prob
      if (length(last_bin_data) < estimate_last_bin_num) {
        all_aics <- c(all_aics, estimate_aic)
        all_bin_nums <- c(all_bin_nums, 0)
        break()
      }
      estimate_last_bin_data <- last_bin_data[sample(length(last_bin_data), estimate_last_bin_num)]
      whole_data <- c(trim_data, estimate_last_bin_data)
      each_step_dif <- (abs(m[length(m)] - m[length(m) - 1])) * length(whole_data)
      if (each_step_dif < 1) {
        all_aics <- c(all_aics, estimate_aic)
        all_bin_nums <- c(all_bin_nums, ceiling(estimate_last_bin_num))
        break()
      } else {
        estimate_data <- whole_data
      }
    }
  }
  best_fit <- which.min(all_aics)
  best_fit_num <- all_bin_nums[best_fit]
  hgt_num <- length(last_bin_data) - best_fit_num
  return(list(all_gofs[best_fit], best_fit_num, hgt_num))
}
infer_HGT <- function(collections_dir, result_file, p_min = 50.0, p_max = 98.5) {
  all_files <- list.files(collections_dir)
  p_last_bin = c(p_max, 100.0)
  # run each file
  strain_pairs <- c()
  em_results <- c()
  hgt_results <- c()
  for(each_file in all_files){
    each_pair <- unlist(strsplit(each_file, split = '[.]'))[1]
    strain_pairs <- c(strain_pairs, each_pair)
    f <- file.path(collections_dir, each_file)
    em_result <- em_implementation(data_file = f, trim_min = p_min, trim_max = p_max, last_bin = p_last_bin)
    # em_results <- c(em_results, em_result[[2]])
    hgt_results <- c(hgt_results, em_result[[3]])
  }
  # all_results <- data.frame(strain_pairs, em_results, hgt_results)
  all_results <- data.frame(strain_pairs, hgt_results)
  colnames(all_results) <- c('Strain Pair', 'Recent HGT Number')
  # output
  write.table(all_results, file = result_file, row.names = F, quote = F, col.names = T, sep = "\t")
}
args <- commandArgs(trailingOnly = T)
collections_dir <- args[1]
result_file <- args[2]
param_min <- as.numeric(args[3])
param_max <- as.numeric(args[4])
infer_HGT(collections_dir, result_file, param_min, param_max)
