# Tests for milags_functions
#
# Author: Jungwirt
###############################################################################


# tests for get_n0
context("Test the get_n0 function")


test_that("Getting initial biomass works", {

  # data
  biom <- c(27727155, 28619420, 31695970, 35766010, 40113020, 44429830, 47754910, 49705080)
  method_1 <- "first.observation"
  method_2 <- "minimal.observation"
  first_calc <- biom[1]
  minim_calc <- min(biom)

  expect_equal(get_n0(biom, method_1), first_calc)
  expect_equal(get_n0(biom, method_2), minim_calc)
})

# context("Test the calc_lag function")

# test_that("Calculating if lag of given model works", {


#  expect_equal()
#})



context("Test the get_def_pars function")

test_that("Calculating if getting default parameters works", {

  expect_equal(get_def_pars(), list(model = "logistic",
                                    n0_method = "first.observation",
                                    tangent_method = "local.regression",
                                    threshold = 10^2,
                                    curve_points = 3,
                                    init_gr_rate = NULL,
                                    init_lag = NULL,
                                    algorithm = "auto",
                                    max_iter = 100))
})


database <- data.frame(time = c(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10),
                 biomass = c(4396571.976, 3807332.496, 4165206.611, 5690282.713, 7727252.94, 19381419.82,
                             13744788.86, 18066675.15, 23651017.71, 29465323.75, 28528881.12, 29448677.51,
                             29144257.31, 32086465.47, 29732262.17, 29888494.33, 30720558.23, 31220300,
                             32074919.83))

context("Test the lag_biomass_incr function")

test_that("Calculating if fitting the lag to multiple growth curves based on the biomass increase method works", {


  test_df <- database 
  test_threshold <- 5000
  test_n0 <- 0
  data_test <- test_df %>% filter(FALSE) %>% mutate(lag = numeric(0))
  for (this_curve_id in unique(test_df$curve_id)) {
    data_this_curve <- test_df %>%
      filter(curve_id == this_curve_id) %>%
      left_join(test_n0, by = "curve_id") %>%
      mutate(incr_from_n0 = biomass - test_n0)


    #find where second derivative is maximal
    threshold_diff <- which(data_this_curve$incr_from_n0 >= test_threshold)
    first_threshold_diff <- threshold_diff[1]
    lag_this_curve <- data_this_curve$time[first_threshold_diff]
    data_this_curve <- data_this_curve %>%
      mutate(
        lag = round(lag_this_curve,1))
    data_test <- rbind(data_test, data_this_curve)
}
  expect_equal(lag_biomass_incr(test_df, test_threshold, test_n0), data_test)
})



context("Test the get_init_pars_baranyi function")

test_that("Getting initial parameters for Baranyi algorithm works", {

  # data
  test_df <- database 
  init_lag <- NULL
  init_gr_rate <- NULL
  min_b <- 0.2
  min_a <- 0.8
  this_n0 <- 0.5
  if (is.null(init_lag)) {
    init_lag <- calc_lag(test_df, method = "tangent", pars = get_def_pars()) %>% pull(lag) %>% unique() %>% as.numeric()
  }

  if (is.null(init_gr_rate)) {
    data_this_curve_exp <- test_df %>%
      mutate(
        max_biomass = max(biomass),
        min_threshold = this_n0 + min_b * (max_biomass - this_n0),
        max_threshold = this_n0 + min_a * (max_biomass - this_n0)) %>%
      # take only the points that are in between min and max
      filter(biomass <= max_threshold & biomass >= min_threshold)
    data_this_curve_exp$logdata <- log(data_this_curve_exp$biomass/this_n0)
    if (nrow(data_this_curve_exp %>% filter(!is.na(time) & !is.na(logdata))  > 0)) {
      mod <- lm(logdata ~ time, data = data_this_curve_exp)
      # this growth rate is assuming an exponential model so it will be generally underestimated
      init_gr_rate <- mod$coefficients[2] %>% unname()
      # we have real r = r(1-N/K) so let us take
      n_mid <- median(data_this_curve_exp$biomass)
      init_mumax <- init_gr_rate
    } else {
      init_mumax <- 0.1
    }
  } else {
    init_mumax <- init_gr_rate
  }
  test_init <- list(init_mumax = init_mumax, init_lag = init_lag)
  expect_equal(get_init_pars_baranyi(test_df, 0.5, NULL, NULL), test_init )
})

context("Test the cut_data function")
test_that("Cutting biomass data works", {

  # data
  test_df <- database 
  max_time <- 5
  data_short <- test_df %>% filter(time <= max_time)
  expect_equal(cut_the_data(test_df, max_time), data_short )
})

 context("Test the fit_max_infl_lag function")
test_that("Fitting maximal biomass data to lag works", {

  # data
  test_df <- database 
  if (!("curve_id" %in% names(test_df))) {
    test_df$curve_id <- "growth.curve"
  }
  data_new <- test_df %>% filter(FALSE) %>% mutate(log_biomass = numeric(0), diff = numeric(0), lag = numeric(0))
  for (this_curve_id in unique(test_df$curve_id)) {
    data_this_curve <- test_df %>%
      filter(curve_id == this_curve_id) %>%
      arrange(time) %>%
      as.data.frame() %>%
      mutate(log_biomass = log(biomass),
             time_diff = mean(c(diff(time), diff(dplyr::lag(time))), na.rm = TRUE),
             time_av = (time + dplyr::lag(time))/2) %>%
      mutate(
        #second.deriv.b = c(NA,second_deriv(time, log.biomass), NA),
        # central scheme
        second_deriv_b = (dplyr::lead(log_biomass) + dplyr::lag(log_biomass) - 2 * log_biomass)/time_diff^2,
        # we only look at second derivative if we know the first derivative is positive!
        # In empirical data we sometimes see a decrease in biomass but we don;t want to look at the second derivative there!
        biomass_incr = dplyr::lead(log_biomass)  > dplyr::lag(log_biomass)
      )

    #find where second derivative is maximal
    max_second_deriv_b <- max(data_this_curve$second_deriv_b, na.rm=TRUE)
    # take first point when max derivative if there are multiple
    ind_max_second_deriv_b <- which(data_this_curve$second_deriv_b == max_second_deriv_b)[1]
    lag_this_curve <- data_this_curve$time[ind_max_second_deriv_b]
    data_this_curve <- data_this_curve %>%
      mutate(
        lag = round(lag_this_curve, 1))
    data_new <- rbind(data_new, data_this_curve)
  }
  expect_equal(fit_max_infl_lag(test_df), data_new )
})

 context("Test the lag_biomass_incr function")
test_that("Biomass increase method works", {

  # data
  test_df <- database 
  if (!("curve_id" %in% names(test_df))) {
    test_df$curve_id <- "growth.curve"
  }
  threshold <- 29465323.75
   n0 <- test_df %>%
    group_by(curve_id) %>%
    arrange(time) %>%
    summarise(n0 = get_n0(biomass, "minimal.observation")) %>%
    ungroup() %>%
    mutate(log_n0 = log(n0))
   data_new <- test_df %>% filter(FALSE) %>% mutate(lag = numeric(0))
  for (this_curve_id in unique(test_df$curve_id)) {
    data_this_curve <- test_df %>%
      filter(curve_id == this_curve_id) %>%
      left_join(n0, by = "curve_id") %>%
      mutate(incr_from_n0 = biomass - n0)


    #find where second derivative is maximal
    threshold_diff <- which(data_this_curve$incr_from_n0 >= threshold)
    first_threshold_diff <- threshold_diff[1]
    lag_this_curve <- data_this_curve$time[first_threshold_diff]
    data_this_curve <- data_this_curve %>%
      mutate(
        lag = round(lag_this_curve,1))
    data_new <- rbind(data_new, data_this_curve)
  }
  expect_equal(lag_biomass_incr(test_df, threshold, n0 ), data_new )
})



 context("Test the fit_exp_lag function")
test_that("Fitting biomass data to lag with tangent method works", {

  # data
  test_df <- database 
  if (!("curve_id" %in% names(test_df))) {
    test_df$curve_id <- "growth.curve"
  }
    n0 <- test_df %>%
    group_by(curve_id) %>%
    arrange(time) %>%
    summarise(n0 = get_n0(biomass, "minimal.observation")) %>%
    ungroup() %>%
    mutate(log_n0 = log(n0))
  tangent_method <- "local.regression"
  curve_points <- 3
  data_new <- test_df %>% filter(FALSE) %>% mutate(lag = numeric(0), line_intercept = numeric(0), line_slope = numeric(0))
  for (this_curve_id in unique(test_df$curve_id)) {
    data_this_curve <- test_df %>% filter(curve_id == this_curve_id)
    this_n0 <- n0 %>% filter(curve_id == this_curve_id) %>% pull(n0)
    lag_obj <- fit_exp_lag_to_curve(data_this_curve, this_n0, tangent_method, curve_points)
    data_this_curve <- data_this_curve %>%
      mutate(lag = round(lag_obj$lag,1)) %>%
      mutate(line_intercept = lag_obj$line_intercept,
             line_slope = lag_obj$line_slope) %>%
      left_join(lag_obj$tangent_points)
    data_new <- rbind(data_new, data_this_curve)
  }

  data_new$predicted <- NA
  expect_equal(fit_exp_lag(test_df, tangent_method, n0, curve_points), data_new )
})

context("Test the get_init_pars_logistic function")
test_that("Getting initial parameters for logistic growth curve works", {

  # data
  test_df <- database
  init_K <- NULL
  init_lag <- NULL
  init_gr_rate <- NULL
  this_n0 <- get_n0(test_df$biomass, "minimal.observation")
  min_b <- 0.2
  min_a <- 0.8
  
if (is.null(init_K)) {
    max_this_data <- max(test_df$biomass, na.rm = TRUE)
    init_K <- max_this_data %>% as.numeric()
  }
  if (is.null(init_lag)) {
    init_lag <- calc_lag(test_df, method = "tangent", pars = get_def_pars()) %>% pull(lag) %>% unique() %>% as.numeric()
  }

  if (is.null(init_gr_rate)) {
    data_this_curve_exp <- test_df %>%
      mutate(
        max_biomass = max(test_df$biomass),
        min_threshold = this_n0 + min_b*(max_biomass - this_n0),
        max_threshold = this_n0 + min_a*(max_biomass - this_n0)) %>%
      # take only the points that are in between min and max
      filter(biomass <= max_threshold & biomass >= min_threshold)
    data_this_curve_exp$logdata = log(data_this_curve_exp$biomass/this_n0)
    if (nrow(data_this_curve_exp %>% filter(!is.na(time) & !is.na(logdata)) > 0)) {
      mod = lm(logdata ~ time, data = data_this_curve_exp)
      # this growth rate is assuming an exponential model so it will be generally underestimated
      init_gr_rate = mod$coefficients[2] %>% unname()
      # we have real r = r(1-N/K) so let us take
      n_mid = median(data_this_curve_exp$biomass)
      init_gr_rate = init_gr_rate/(1-n_mid/init_K)
    } else {
      init_gr_rate = 0.1
    }
    #data_this_curve_exp$predicted = predict(mod, data_this_curve_exp)
  }
  data_new <- list(init_K = init_K, init_lag = init_lag, init_gr_rate = init_gr_rate)
  expect_equal(get_init_pars_logistic(test_df, this_n0, init_K, init_lag, init_gr_rate, min_b, min_a), data_new )
})

 context("Test the choose_lag_fit_algorithm_baranyi function")
test_that("Choosing the best nls model to fit to algorithm baranyi works", {

  # data
  test_df <- database 
  LOG10N0 <- NULL
  init_lag <- NULL
  init_mumax <- NULL
  init_LOG10Nmax <- NULL
  max_iter <- 100
  lower_bound <- c(0,0,0, 0)
  tryCatch(
    expr =
      {nlsres_LM <- nlsLM(formula = baranyi,
                        data = test_df,
                        start = list(lag=init_lag, mumax=init_mumax, LOG10N0 = LOG10N0, LOG10Nmax = init_LOG10Nmax),
                        control = nls.control(maxiter = max_iter),
                        lower = lower_bound)
      },
    error = function(cond) {
      # this operator assigns value outside the error environment
      nlsres_LM <<- NA
    })
  tryCatch(
    expr =
      {nls_LM_no_bound <- nlsLM(formula = baranyi,
                                 data = test_df,
                                 start = list(lag=init_lag, mumax=init_mumax, LOG10N0 = LOG10N0, LOG10Nmax = init_LOG10Nmax),
                                 control = nls.control(maxiter = max_iter))
      },
    error = function(cond) {
      nls_LM_no_bound <<- NA
    })
  tryCatch(
    expr =
      {nls_PORT <- nls(baranyi,
                        data = test_df,
                        start = list(lag=init_lag, mumax=init_mumax, LOG10N0 = LOG10N0, LOG10Nmax = init_LOG10Nmax),
                        algorithm = "port",
                        control = nls.control(maxiter = max_iter),
                        lower = lower_bound)
      },
    error = function(cond) {
      nls_PORT <<- NA
    })

  nls_a <- compare_algorithms(nls_LM_no_bound, nls_PORT, nlsres_LM)
  # consider the model without lower bounds only if the resulted estimates are above 0

  expect_equal(choose_lag_fit_algorithm_baranyi(test_df, LOG10N0 = NULL, init_lag = NULL, init_mumax = NULL, init_LOG10Nmax = NULL, max_iter = 100, lower_bound = c(0,0,0, 0)), nls_a )
})

 context("Test the choose_lag_fit_algorithm_logistic function")
test_that("Choosing the best nls model to fit to logistic curve works", {

  # data
  test_df <- database 
  n0 <- get_n0(test_df$biomass, "minimal.observation")
  init_K <- NULL
  init_lag <- NULL
  init_gr_rate <- NULL
  max_iter <- 100
  lower_bound <- c(0,0,0)
 tryCatch(
    expr =
      {nlsres_LM <- nlsLM(formula = biomass ~ n0 + (time >= lag)*n0*(-1+K*exp(gr_rate*(time-lag))/(K - n0 + n0*exp(gr_rate*(time - lag)))),
                        data = gr_curve,
                        start = list(gr_rate = init_gr_rate, K=init_K, lag = init_lag),
                        control = nls.control(maxiter = max_iter),
                        lower = lower_bound)
      },
    error = function(cond) {
      nlsres_LM <<- NA
    })
  tryCatch(
    expr =
      {nls_LM_no_bound <- nlsLM(formula = biomass ~ n0 + (time >= lag)*n0*(-1+K*exp(gr_rate*(time-lag))/(K - n0 + n0*exp(gr_rate*(time - lag)))),
                                 data = gr_curve,
                                 start = list(gr_rate = init_gr_rate, K=init_K, lag = init_lag),
                                 control = nls.control(maxiter = max_iter))
      },
    error = function(cond) {
      nls_LM_no_bound <<- NA
    })
  tryCatch(
    expr =
      {nls_PORT <- nls(biomass ~ n0 + (time >= lag)*n0*(-1+K*exp(gr_rate*(time-lag))/(K - n0 + n0*exp(gr_rate*(time - lag)))), gr_curve,
                        list(gr_rate = init_gr_rate, K=init_K, lag = init_lag),
                        algorithm = "port",
                        control = nls.control(maxiter = max_iter),
                        lower = lower_bound)},
    error = function(cond) {
      nls_PORT <<- NA
    })
  nls_a <- compare_algorithms(nls_LM_no_bound, nls_PORT, nlsres_LM)
   expect_equal(choose_lag_fit_algorithm_logistic(test_df, n0, init_gr_rate, init_K, init_lag, max_iter, lower_bound), nls_a )
})


