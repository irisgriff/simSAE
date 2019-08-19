sim_sae <- function(formula, data, data_pop, domains, iterations, prop = 0.05) {
  #source("EBLUP.R")
  response_lab <- all.vars(as.formula(formula))[1L]
  pred_lab <- all.vars(as.formula(formula))[-1L]
  if(mean(pred_lab %in% names(data_pop)) != 1){
    message("predictor names in sample data not congruent with population data")
  }
  #find averages for entire dataset
  dat_srs <- EBLUP(formula, data = data,
                       domains = domains, pop_data = data_pop) %>% pull(srs)
  X_bar <- pop_means(data, pred_lab, domains = domains)
  #bootstrap MSE
  dat_boot <- list()
  means_vec <- list()
  for(j in 1:7) {
    dat_boot[[j]] <- list()
    dat1 <- data.frame()
    while(nrow(dat1) == 0) { #ensure at least 1 domain is selected
      dat1 <- data %>% 
        group_by_at(domains) %>%
        mutate(n = n(), pick = ifelse(runif(1, -10, (2*j)^1.7) > 0, 1, 0)) %>%
        filter(pick != 0)
      }
    means_vec[[j]] <- dat1 %>%
      group_by_at(domains) %>%
      summarise_at(.vars = response_lab, .funs = "mean") %>%
      pull()
    X_bar1 <- pop_means(dat1, pred_lab, domains = domains)
    for(i in 1:iterations) {
      dat_samp <- dat1 %>% group_by_at(domains) %>%
        mutate(n = n(), pick = rbinom(n = n(),size = 1, prob = prop) +
                 # ensure that at least two units are sampled:
                 sample(c(rep(0, n()-2), 1,  1), replace = F)) %>%
        filter(pick != 0) %>%
        dplyr::select_at(.vars = c("n", response_lab, pred_lab, domains))
      dat_boot[[j]][[i]] <- EBLUP(formula, dat_samp,
                                      domains=domains, pop_data = X_bar1) %>%
        dplyr::select(n_i, srs, sre, greg, eblup, eblup_mse) %>%
        as.matrix()
    }
    cat(j)
  }
  sim_mse <- list()
  sim_bias <- list()
  sim_mse_est <- list()
  pct_better <- rep(NA, 7)
  pct_better_weighted <- rep(NA, 7)
  n_doms <- rep(NA, 7)
  pct_off_mse <- list()
  
  for(j in 1:7) {
    sim_bias[[j]] <- (lapply(1:(iterations),FUN = function(i)(dat_boot[[j]][[i]] -
                                                          means_vec[[j]])) %>%
                        Reduce(x = .,f = '+')/(iterations))[,-1]
    sim_mse[[j]] <- (lapply(1:(iterations), FUN = function(i){(dat_boot[[j]][[i]] -
                                                           means_vec[[j]])^2}) %>%
                       Reduce(x = .,f = '+')/(iterations))[,-1]
    avg_n_i <- lapply(1:iterations, function(i) {
      dat_boot[[j]][[i]][,1]/iterations}) %>%
      Reduce(x = ., f = '+')
    sim_mse_est <- sapply(1:(iterations),
                          FUN = function(i) dat_boot[[j]][[i]][,6]) %>%
      matrix(nrow = iterations, byrow = TRUE)
    pct_better[j] <- mean(((sim_mse[[j]][,2] - sim_mse[[j]][,4]) >0))
    pct_better_weighted[j] <- mean(((sim_mse[[j]][,2] - sim_mse[[j]][,4]) >0) *
                                     avg_n_i)/ sum(avg_n_i)
    n_doms[j] <- nrow(sim_mse[[j]])
    #create list of dfs
    pct_off_mse[[j]] <- data.frame(true_mse = sim_mse[[j]][,4],
                                   sd_mse_est = sapply(1:n_doms[j],
                                                       function(i) {
                                                         sd(sim_mse_est[,i])}),
                                   bias_mse_est = colMeans(sim_mse_est) -
                                     sim_mse[[j]][,4],
                                   pct_overest = sapply(1:n_doms[j],
                                                        function(i) {
                                                          mean(sim_mse_est[,i] >
                                                                 sim_mse[[j]][,4][i])}),
                                   avg_n_i = avg_n_i) %>%
      mutate(mse_mse = bias_mse_est^2 + sd_mse_est^2)
  }
  # collapse list of pct off mse
  mse_estimate_df <- bind_rows(pct_off_mse, .id = "j_ind") %>%
    mutate(n_doms = n_doms[as.numeric(j_ind)],
           pct_better = pct_better[as.numeric(j_ind)]) %>%
    group_by(j_ind) %>%
    mutate(avg_total_n = sum(avg_n_i)) %>%
    ungroup()
  # collapse list of bias and MSE matrices
  sim_df <- (sim_bias %>% lapply(function(x) as.data.frame(x) %>%
                                   gather(key = "estimator", value = "bias")) %>%
               bind_rows(.id = "j_ind")) %>%
    mutate(MSE = sim_mse %>% lapply(function(x) as.data.frame(x) %>%
                               gather(key = "estimator", value = "MSE")) %>%
            bind_rows() %>% pull("MSE")) %>%
    mutate(var = MSE - bias^2) %>%
    filter(estimator != "eblup_mse")
  # make a ggplot
  plot_eblup_mse <- ggplot(mse_estimate_df, aes(x = avg_total_n,
                                                col = n_doms, y = pct_better)) +
    geom_point(alpha = 0.9) +
    #geom_point(inherit.aes = F, aes(x = 0, y = 0), col = "red") +
    scale_color_viridis_c() +
    labs(x = "avg sample size", y = "percent EBLUP best", col = "domains")
  # collapse list of simulated mse
  mse_df <- bind_rows(lapply(sim_mse, as.data.frame), .id = "j_ind") %>%
    mutate(n_doms = n_doms[as.numeric(j_ind)],
           avg_total_n = mse_estimate_df$avg_total_n)
  # make ggplot
  plot_eblup_mse <- ggplot(mse_df, aes(x = eblup, y = sre, col = avg_total_n)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, lty = 3) +
    scale_color_viridis_c() +
    coord_fixed(1) +
    labs(x = "EBLUP empirical MSE", y = "SRE empirical MSE", col = "avg sample size")
  # make ggplot again: empirical bias
  plot_eblup_mse <- sim_df %>%
    ggplot(aes(y = var, x = bias, col = key)) + geom_point(alpha = 0.7) +
    labs(x = "empirical bias",y = "empirical variance") + scale_color_viridis_d()
  return(list(sim_df = sim_df,
              mse_est = mse_estimate_df,
              plot_eblup = plot_eblup_mse))
}