#' This script is used to detect recent horizontal gene transfers.
#' @param collections_dir, should be the address of the input directory, string.
#' @param result_file, should be the address of the result file, string.
#' @param param_min, should be a value of lower limit of the distribuiton.
#' @param param_max, should be a value of upper limit of the distribuiton.
#' @export result_file, a text file with the number of rHGTs between every strain pair.

library(fitdistrplus, quietly = T)

# fit function use fitdist method
fit_distribution <- function(data, dist, use_method, use_gof) {
  fit <- fitdist(data, dist, method = use_method, gof = use_gof)
  # print(fit)
  stat <- gofstat(fit)
  aic <- as.numeric(stat$aic)
  if (dist == 'lnorm'){
    a <- as.numeric(fit$estimate['meanlog'])
    b <- as.numeric(fit$estimate['sdlog'])
  } 
  else if (dist == 'weibull') {
    a <- as.numeric(fit$estimate['shape'])
    b <- as.numeric(fit$estimate['scale'])
  }
  else if (dist == 'gamma') {
    a <- as.numeric(fit$estimate['shape'])
    b <- as.numeric(fit$estimate['rate'])
  }
  dist_pars <- c(a, b)
  return(list(dist_pars, stat, aic))
}
em_implementation <- function(data_file, trim_min, trim_max, last_bin) {
  # load and trim data
  my_data <- read.csv(data_file, sep = '\t', header = T)
  trim_data <- my_data$Similarity[my_data$Similarity >= trim_min & 
                                    my_data$Similarity < trim_max]
  last_bin_data <- my_data$Similarity[my_data$Similarity >= last_bin[1] & 
                                        my_data$Similarity <= last_bin[2]]
  # EM implementation
  m <- c(0)
  l <- last_bin[1]
  u <- last_bin[2]
  fit_method <- 'mge'
  # dists <- c('weibull', 'lnorm', 'gamma')
  dists <- c('weibull')
  all_gofs <- c('AD2R', 'ADR', 'AD', 'KS', 'CvM')
  max_steps <- 100
  estimate_data <- trim_data
  all_aics <- c()
  all_bin_nums <- c()
  for (fit_dist in dists) {
    for (fit_gof in all_gofs) {
      for (i in 1:max_steps) {
        estimate_result <- fit_distribution(data = estimate_data, dist = fit_dist, 
                                            use_method = fit_method, use_gof = fit_gof)
        estimate_pars <- estimate_result[[1]]
        estimate_stat <- estimate_result[[2]]
        estimate_aic <- estimate_result[[3]]
        estimate_shape <- estimate_pars[1]
        estimate_scale <- estimate_pars[2]
        if (fit_dist == 'weibull') {
          prob <- pweibull(u, shape=estimate_shape, scale=estimate_scale) - 
            pweibull(l, shape=estimate_shape, scale=estimate_scale)
        } 
        else if (fit_dist == 'gamma') {
          prob <- pgamma(u, shape=estimate_shape, rate=estimate_scale) - 
            pgamma(l, shape=estimate_shape, rate=estimate_scale)
        }
        else if (fit_dist == 'lnorm') {
          prob <- plnorm(u, meanlog=estimate_shape, sdlog=estimate_scale) - 
            plnorm(l, meanlog=estimate_shape, sdlog=estimate_scale)
        }
        m <- c(m, prob)
        estimate_last_bin_num <- length(trim_data) * prob
        if (length(last_bin_data) < estimate_last_bin_num) {
          all_aics <- c(all_aics, estimate_aic)
          all_bin_nums <- c(all_bin_nums, 0)
          break()
        }
        estimate_last_bin_data <- last_bin_data[sample(length(last_bin_data), 
                                                estimate_last_bin_num)]
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
  }
  best_fit <- which.min(all_aics)
  best_fit_num <- all_bin_nums[best_fit]
  hgt_num <- length(last_bin_data) - best_fit_num
  return(list(dists[best_fit], all_gofs[best_fit], 
              best_fit_num, hgt_num))
}
infer_HGT <- function(collections_dir, result_file, 
                      p_min = 0.0, p_max = 98.5) {
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
    em_result <- em_implementation(data_file = f, trim_min = p_min,
                                   trim_max = p_max, last_bin = p_last_bin)
    hgt_results <- c(hgt_results, em_result[[4]])
  }
  all_results <- data.frame(strain_pairs, hgt_results)
  colnames(all_results) <- c('Strain Pair', 'Recent HGT Number')
  # output
  write.table(all_results, file = result_file, row.names = F, 
              quote = F, col.names = T, sep = "\t")
}
args <- commandArgs(trailingOnly = T)
input_dir <-  args[1]
output_file <- args[2]
# param_min <- as.numeric(args[3])
param_max <- as.numeric(args[3])
infer_HGT(collections_dir = input_dir, 
          result_file = output_file, 
          p_max = param_max)
