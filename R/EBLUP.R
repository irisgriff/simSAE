EBLUP <- function(formula, data, domains, pop_data = NULL, X_bar = NULL) {
  require(JoSAE)
  data <- data %>% arrange_at(domains)
  n_dom <- data %>% pull(domains) %>% unique() %>% length()
  # find response variable
  response_lab <- all.vars(as.formula(formula))[1L]
  pred_lab <- all.vars(as.formula(formula))[-1L]
  resp <- data %>% pull(response_lab) %>% as.vector()
  #make model matrix of sample means
  x_bar <- cbind(rep(1,n_dom), data %>%
                   dplyr::select_at(c(pred_lab, domains)) %>%
                   group_by_at(domains) %>%
                   summarize_all(funs(mean)) %>%
                   dplyr::select_at(pred_lab) %>%
                   as.matrix()) %>% unname()
  # fit OLS and LMM
  ols <- lm(formula, data)
  mod <- lme(formula, data, random = formula(paste0("~1|", domains)))
  #make matrix of pop means
  if(!is.null(pop_data)){
    X_bar <- cbind(rep(1,nrow(pop_data)),
                   pop_data %>% dplyr::select_at(pred_lab) %>%
                     as.matrix()) %>% unname()
  }
  else{
    if(!is.null(X_bar)){
      X_bar <- rep(NA, n_dom*(1 + length(pred_lab))) %>%
        matrix(nrow = n_dom)
    }
  }
  #check if population and sample domains are consistent
  if(!is.null(pop_data) && nrow(X_bar) != nrow(x_bar)){
    stop("number of sample domains not equal to number of population domains, try filtering by domains %in% sample$domains, or vice versa.")
  }
  # find synthetic estimates
  sre <- as.vector(X_bar %*% ols$coefficients) %>%
    unname()
  lme_sre <- as.vector(X_bar %*% mod$coefficients$fixed) %>%
    unname()
  # find estimates!
  #fixed = (mod$fitted %>% rowSums() %>% as.vector()),
  out <- cbind(data %>% as.data.frame(),
      data.frame(ols_resid = ols$residuals,
                 lmm_resid = (mod$residuals %>% rowSums() %>% as.vector()),
                 fixed_resid = mod$residuals[,1] %>% as.vector(),
                 response = resp)) %>%
    group_by_at(domains) %>%
    summarize(n_i = n(),
              srs = mean(response),
              srs_mse = var(response)/n_i,
              r = mean(ols_resid),
              r_random = mean(lmm_resid),
              fr = mean(fixed_resid),
              sd_nu = (mod$coefficients$random[[1]][,1]) %>% sd()) %>%
    dplyr::mutate(sre = sre,
                  lme_sre = lme_sre,
                  greg = sre + r,
                  gamma = (sd_nu^2 / (sd_nu^2 + mod$sigma^2/n_i)),
                  eblup = lme_sre + r_random*gamma) %>%
    dplyr::select(fr, domains, n_i, srs, srs_mse, sre, greg, eblup, gamma) %>%
    arrange_at(domains)
  #return(out)
  # EBLUP mse estimate things
  # make U matrix for each domain for MSE estimate
  makeU <- function(index, out){
    ni <- max(out[index,3], 0)
    g <- out[index,9]
    m1 <- diag(1, ni) %>% as.numeric()
    m2 <- matrix(rep(g/ni, ni^2), nrow = ni) %>% as.numeric()
    return(matrix((m1 - m2)/(mod$sigma^2), nrow = ni))
  }
  U_list <- lapply(1:n_dom, FUN = function(x) makeU(index = x,out = out))
  #make list of design matricies
  design_list <- out %>% pull(domains) %>%
    lapply(function(dom) model.matrix(formula,
                            data = filter_at(data, .vars = vars(domains),
                                .vars_predicate = function(x) x == dom)) %>%
             as.matrix() )
  # Find the EBLUP MSE estimate!
  annoying_thing <- lapply(1:n_dom, function(i) {
    t(design_list[[i]]) %*% U_list[[i]] %*% (design_list[[i]])
    }) %>% Reduce(f = "+", x = .) %>% solve()
  # Final
  out <- out %>%
    mutate(C1 = gamma*(mod$sigma^2/n_i),
           C2 = sapply(1:n_dom, function(i) {
             t(X_bar[i,] - gamma[i]*x_bar[i,]) %*%
               annoying_thing %*% 
               (X_bar[i,] - gamma[i]*x_bar[i,])
           }),
           C3 = eblup.mse.f.c3(mod,
                               eblup.mse.f.c3.asyvarcovarmat(mod, n_i),
                               n.i = n_i),
           C3_star = eblup.mse.f.c3.star(mod,
                                         eblup.mse.f.c3.asyvarcovarmat(mod,n_i),
                                         n.i = n_i,
                                         mean.resid.i = fr),
           eblup_mse = C1 + C2 + C3 + C3_star) %>%
    dplyr::select(domains, n_i, srs, srs_mse, sre,
                  greg, eblup, eblup_mse, gamma)
  return(out)
}

#find pop means
pop_means <- function(pop_data, predictor_labels, domains) {
  pop_data %>%
    group_by_at(.vars = domains) %>%
    arrange() %>%
    mutate(n_i = n()) %>%
    summarise_at(.vars = c("n_i", predictor_labels),
                 .funs = function(x) mean(x)) #%>% 
    #dplyr::select_at("n_i", domains, predictor_labels)
}


# find difference between two EBLUP dataframes
EBLUP_diff <- function(pop, samp){
  out <- dplyr::full_join(pop, samp, by = "DOMAIN") %>%
    within(rm("DOMAIN")) %>%
    t() %>%
    data.frame() %>%
    group_by(., id = gsub('\\..*', '', rownames(.))) %>%
    summarise_all(function(x) diff(x)) %>%
    data.frame() %>%
    column_to_rownames(var = 'id') %>%
    t() %>% data.frame() %>%
    mutate(domain = pop$DOMAIN) %>%
    select(domain, n_i, srs, srs_mse, sre, greg, eblup, eblup_mse, gamma) %>%
    arrange(domain)# %>%
    #filter(!is.na(n_i))
  return(out)
}
