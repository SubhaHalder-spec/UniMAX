#' UniMAX test for testing equality means against umbrella-ordered alternatives (peak = h) in one-way ANOVA
#' @export
#' @param sample_data list
#' @param significance_level numeric
#' @param h numeric
#' @return Critical value numeric
#' @return Test statistic value numeric
#' @return Result Character
#' @details Testing of H_0:mu_1 = mu_2 = ... = mu_k vs H_1:mu_1 <=.....<= mu_(h-1)<= mu_h >= mu_(h+1)>=....>= mu_k (at least one strict inequality), where mu_i represents the population means of the i-th treatment. The input consists of three variables: sample_data, significance_level, and peak = h. The output consists of the critical value, the UniMAX test statistic value, and the result, which indicates whether to reject or not reject the null hypothesis.
#' @import stats
#' @author Subha Halder
UniMAX <- function(sample_data, significance_level, h){
  set.seed(456)
  sample_data <- lapply(sample_data, function(x) x[!is.na(x)])
  num_samples = 100000
  num_datasets <- length(sample_data)
  n <- sapply(sample_data, length)
  D_star_max <- numeric(num_samples)
  for (i in 1:num_samples) {
    bootstrap_samples <- lapply(sample_data, function(x) rnorm(n = length(x), mean = 0, sd = sqrt(var(x))))
    D_star_max[i] <- max(c(
      sapply(2:h, function(j) {
        (mean(bootstrap_samples[[j]]) - mean(bootstrap_samples[[j - 1]])) /
          sqrt(
            (var(bootstrap_samples[[j]]) / length(bootstrap_samples[[j]])) +
              (var(bootstrap_samples[[j - 1]]) / length(bootstrap_samples[[j - 1]]))
          )
      }),
      sapply((h + 1):num_datasets, function(j) {
        (-mean(bootstrap_samples[[j]]) + mean(bootstrap_samples[[j - 1]])) /
          sqrt(
            (var(bootstrap_samples[[j]]) / length(bootstrap_samples[[j]])) +
              (var(bootstrap_samples[[j - 1]]) / length(bootstrap_samples[[j - 1]]))
          )
      })
    ))

  }
  sort_D_star_max <- sort(D_star_max)
  quantile_value <- quantile(sort_D_star_max, probs = 1 - significance_level)
  UMAX <- max(c(
    sapply(2:h, function(i) {
      (mean(sample_data[[i]]) - mean(sample_data[[i - 1]])) /
        sqrt(
          (var(sample_data[[i]]) / length(sample_data[[i]])) +
            (var(sample_data[[i - 1]]) / length(sample_data[[i - 1]]))
        )
    }),
    sapply((h + 1):num_datasets, function(i) {
      (-mean(sample_data[[i]]) + mean(sample_data[[i - 1]])) /
        sqrt(
          (var(sample_data[[i]]) / length(sample_data[[i]])) +
            (var(sample_data[[i - 1]]) / length(sample_data[[i - 1]]))
        )
    })
  ))
  if (UMAX > quantile_value) {
    result <- "Reject null hypothesis"
  } else {
    result <- "Do not reject null hypothesis"
  }
  return(paste("Critical value:", quantile_value, "; UniMAX Test statistic:", UMAX, "; Result:", result))
}




