
############################################# helper functions #########################################
#' get_n0
#'
#' Gets the initial biomass to relate to
#' @param biomass vector of biomass (chronologically ordered as in growth curve)
#' @param n0_method "first.observation" if the first point is taken as the initial biomass or
#' "minimal.observation" if the minimal biomass is taken is the initial point.
#' In "healthy" growth curves these options should be equivalent
#' but sometimes a drop in OD/biomass is observed at the beginning of a growth curve.
#' In this case it is not obvious what to assume the initial biomass is.
#' @import dplyr
#' @import ggplot2
#' @import minpack.lm
#' @import nlsMicrobio
#' @rawNamespace import(testthat, except = getNamespaceExports("testthat"))
#' @returns a value of the initial biomass (either the first observation or the minimum value depending on the parameter N0.method)
get_n0 <- function(biomass, n0_method) {
  # Get the initial biomass value
  if (n0_method == "first.observation") {
    n0_calc <- biomass[1]
  } else if (n0_method == "minimal.observation") {
    n0_calc <- min(biomass)
  }
  return(n0_calc)
}



#' fit_exp_lag_to_curve
#'
#' Fits the lag to one growth curve based on the basic tangent method
#' @param data a data frame with two required columns names: "time" and "biomass",
#' This is data from one growth curve only, one (mean) observation per time
#' @param n0 the initial biomass (a tangent line crossing N0 line will determine the lag)
#' @param tangent_method "local.regression" (if the tangent is fitted to a number of points around the maximal growth rate)
#' or "to.point" (if the tangent is fitted only to the point where the growth rate is maximal); defaults to "to.point"
#' @param curve_points if tangent_method = "local.regression" then curve_points is the number of points the line is fitted to;
#' defaults to 3 i.e. the point with the maximal uptake rate one point before and one point after
#' @returns line_slope: slope of the tangent line,
#' line_intercept: intercept of the tangent line,
#' lag: lag,
#' tangent_points: i..e a data frame of all points selected for fitting the line
fit_exp_lag_to_curve <- function(data, n0, tangent_method = "to.point", curve_points = 3) {
  data <- data %>% as.data.frame() %>%
    mutate(log_biomass = log(biomass),
           db = log_biomass - dplyr::lag(log_biomass,1),
           dt = time - dplyr::lag(time,1),
           time_diff = (time + dplyr::lag(time))/2,
           central_scheme_db = dplyr::lead(log_biomass) - dplyr::lag(log_biomass,1),
           central_scheme_dt = dplyr::lead(time) - dplyr::lag(time,1),
           db_over_dt = central_scheme_db/central_scheme_dt) #db/dt)
  log_n0 <- log(n0)
  # Get the maximum growth rate
  max_gr_rate <- max(data$db_over_dt, na.rm = TRUE)
  ind_max_diff <- which(data$db_over_dt == max_gr_rate)

  if (length(ind_max_diff) > 1) {
    warning("Multiple points with max derivative: taking the first one")
    ind_max_diff <- ind_max_diff[1]
  }
  if (tangent_method == "to.point") {
    # 1. take simply one maximal growth rate point and calculate a tangent line to it
    # by definition the slope of this line will be equal to its derivative i.e. the
    exp_gr_points <- data[ind_max_diff, ]
    time_max_diff <- exp_gr_points$time
    max_diff_log_biomass <- exp_gr_points$log_biomass
    # ax + b = y
    # a = max.growth.rate
    # we know one point (x,y) = (time.max.diff, log.biomass.at.max.diff)
    line_slope <- max_gr_rate
    line_intercept <- max_diff_log_biomass - line_slope*time_max_diff
  } else if (tangent_method == "local.regression") {
    # 2. Take N points around the maximal growth rate and fit a line,
    # linearly extrapolate
    around_points <- floor((curve_points-1)/2)
    exp_gr_points <- data[(ind_max_diff-around_points):(ind_max_diff+around_points),]
    mod <- lm(log_biomass ~ time, exp_gr_points)
    line_intercept <- unname(mod$coefficients[1])
    line_slope <- unname(mod$coefficients[2])
  }
  # at x=lag we have y = log.N0, so lag = (y-b)/a
  lag <- (log_n0 - line_intercept)/line_slope
  return(list(line_slope = line_slope, line_intercept = line_intercept, lag = lag,
              tangent_points = exp_gr_points %>% select(time, tangent_point = biomass)))
}




#' compare_algorithms
#'
#' Compares results of 3 objects obtained from running nls
#' @param nls_LM_no_bound first object resulting from running nls
#' @param nls_PORT second object resulting from running nls
#' @param nlsres_LM third object resulting from running nls
#' @returns the best fitting object (lowest Res.Sum Sq provided that all coefficients are nonnegative)
compare_algorithms <- function(nls_LM_no_bound, nls_PORT, nlsres_LM) {
  if (!all(is.na(nls_LM_no_bound))) {
    s <- summary(nls_LM_no_bound)
    if (all(as.numeric(s$coefficients[,1]) >= 0)) {
      nls_LM_no_bound <- nls_LM_no_bound
    } else {
      nls_LM_no_bound <- NA
    }
  }
  # first compare two bounded models
  if (!all(is.na(nls_PORT)) & !all(is.na(nlsres_LM))) {
    # if both bounded models are available compare them and choose the better one
    anova_res <- anova(nlsres_LM, nls_PORT)
    better_model_idx <- which(anova_res$`Res.Sum Sq` == min(anova_res$`Res.Sum Sq`))
    if (better_model_idx ==1) {nls_bounded <- nlsres_LM} else if (better_model_idx == 2) {nls_bounded <- nls_PORT}
  } else if (!all(is.na(nls_PORT))) {
    nls_bounded <- nls_PORT
  } else if (!all(is.na(nlsres_LM))) {
    nls_bounded <- nlsres_LM
  } else {
    nls_bounded <- NA
  }

  # the netter bounded model compare to unbounded
  if (!all(is.na(nls_bounded)) & !all(is.na(nls_LM_no_bound))) {
    # if both bounded models are available compare them and choose the better one
    anova_res <- anova(nls_bounded, nls_LM_no_bound)
    better_model_idx <- which(anova_res$`Res.Sum Sq` == min(anova_res$`Res.Sum Sq`))
    if (better_model_idx ==1) {nls_a <- nls_bounded} else if (better_model_idx == 2) {nls_a <- nls_LM_no_bound}
  } else if (!all(is.na(nls_bounded))) {
    nls_a <- nls_bounded
  } else if (!all(is.na(nls_LM_no_bound))) {
    nls_a <- nls_LM_no_bound
  } else {
    nls_a <- NA
  }
  return(nls_a)

}




#' choose_lag_fit_algorithm_baranyi
#'
#' Runs nlsLM/nls algorithms with three different parameter setups to fit the best Baranyi parameters to our data and chooses the best model
#' @param gr_curve data from one specific growth curve with the following columns: LOG10N, t
#' @param LOG10N0 init value for the LOG10N0 parameter
#' @param init_lag initial value for the lag
#' @param init_mumax initial value for the mumax parameter
#' @param init_LOG10Nmax initial value for the LOG10Nmax parameter
#' @param max_iter max. number of iterations
#' @param lower_bound lower bound for the bounded nls optimization;
#' @returns the best nls fitting object with parameters fitted to Baranyi model (lowest Res.Sum Sq provided that all coefficients are nonnegative)
choose_lag_fit_algorithm_baranyi <- function(gr_curve,
                                                     LOG10N0,
                                                     init_lag,
                                                     init_mumax,
                                                     init_LOG10Nmax,
                                                     max_iter,
                                                     lower_bound) {

  # choose best from LM and port
  # Sometimes the lower.bound argument makes the fit much worse.
  tryCatch(
    expr =
      {nlsres_LM <- nlsLM(formula = baranyi,
                        data = gr_curve,
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
                                 data = gr_curve,
                                 start = list(lag=init_lag, mumax=init_mumax, LOG10N0 = LOG10N0, LOG10Nmax = init_LOG10Nmax),
                                 control = nls.control(maxiter = max_iter))
      },
    error = function(cond) {
      nls_LM_no_bound <<- NA
    })
  tryCatch(
    expr =
      {nls_PORT <- nls(baranyi,
                        data = gr_curve,
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
  return(nls_a)

}



#' choose_lag_fit_algorithm_logistic
#'
#' Runs nlsLM/nls algorithms with three different parameter setups to fit the best Logistic model parameters to our data and chooses the best model
#' @param gr_curve data from one specific growth curve with the following columns: LOG10N, t
#' @param n0 the initial biomass
#' @param init_gr_rate initial value for the growth rate
#' @param init_K initial value for the saturation parameter K
#' @param init_lag initial value for the lag parameter
#' @param max_iter max. number of iterations; defaults to 100
#' @param lower_bound lower bound for the bounded nls optimization; defaults to 0
#' @returns the best nls fitting object with parameters fitted to logistic model (lowest Res.Sum Sq provided that all coefficients are nonnegative)
choose_lag_fit_algorithm_logistic <- function(gr_curve, n0,
                                                      init_gr_rate = init_gr_rate,
                                                      init_K = init_K,
                                                      init_lag = init_lag,
                                                      max_iter = 100,
                                                      lower_bound = c(0,0,0)) {
  # choose best from LM and port
  # Sometimes the lower.bound argument makes the fit much worse.
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
  # consider the model without lower bounds only if the resulted estimates are above 0
  return(nls_a)

}



#' calc_lag_fit_to_logistic_with_lag
#'
#' Runs nlsLM/nls algorithm of the user's choice to fit the  Logistic model parameters to our data
#' @param gr_curve data from one specific growth curve with these two columns: time and biomass
#' @param n0 the initial biomass
#' @param init_gr_rate initial value for the growth rate
#' @param init_K initial value for the saturation parameter K
#' @param init_lag initial value for the lag parameter
#' @param algorithm defaults to "auto" which chooses between bounded and unbounded Levenberg-Marquardt method and the bounded port method
#' @param max_iter max. number of iterations; defaults to 100
#' @param lower_bound lower bound for the bounded nls optimization; defaults to 0
#' @returns lag and the nls fitting object with parameters fitted to logistic model
calc_lag_fit_to_logistic_with_lag <- function(gr_curve, n0,
                                                      init_gr_rate = init_gr_rate,
                                                      init_K = init_K,
                                                      init_lag = init_lag,
                                                      algorithm = "auto",# Levenberg-Marquardt", # or "port" for nls
                                                      max_iter = 100,
                                                      lower_bound = c(0,0,0)) {
  tryCatch(
    expr =
      {if (algorithm == "auto") {
        nls_a <- choose_lag_fit_algorithm_logistic(gr_curve, n0,
                                                            init_gr_rate = init_gr_rate, init_K = init_K, init_lag = init_lag,
                                                            max_iter = max_iter,
                                                            lower_bound = lower_bound)

        # nlsLM( is a more robust version of nls, using  Levenberg-Marquardt algorithm
      } else if (algorithm == "Levenberg-Marquardt") {
        nls_a <- nlsLM(formula = biomass ~ n0 + (time >= lag) * n0 * (-1 + K * exp(gr_rate * (time - lag)) / (K - n0 + n0 * exp(gr_rate * (time - lag)))),
                       data = gr_curve,
                       start = list(gr_rate = init_gr_rate, K = init_K, lag = init_lag),
                       control = nls.control(maxiter = max_iter))
        #lower= lower.bound)
      } else if (algorithm == "port") {
        nls_a <- nls(biomass ~ n0 + (time >= lag) * n0 * (-1 + K * exp(gr_rate * (time - lag)) / (K - n0 + n0 * exp(gr_rate * (time - lag)))), gr_curve,
                     list(gr_rate = init_gr_rate, K = init_K, lag = init_lag),
                     algorithm = "port",
                     control = nls.control(maxiter = max_iter))
        #lower= lower.bound)
      } else {
        nls_a <- nls(biomass ~ n0 + (time >= lag) * n0 * (-1 + K * exp(gr_rate * (time - lag)) / (K - n0 + n0 * exp(gr_rate * (time - lag)))), gr_curve,
                     list(gr_rate = init_gr_rate, K = init_K, lag = init_lag),
                     algorithm = algorithm,
                     control = nls.control(maxiter = max_iter))
      }
        lag_N <- coef(nls_a)[names(coef(nls_a)) == "lag"] %>% unname()
      },
    error = function(cond) {
      lag_N <<- NA
      nls_a <<- NA
      #print(cond)
    })
  return(list(lag_N = lag_N, nls_a = nls_a))
}



#' calc_lag_fit_to_baranyi_with_lag
#'
#' Runs nlsLM/nls algorithms with three different parameter setups to fit the best Logistic model parameters to our data and chooses the best model
#' @param gr_curve data from one specific growth curve with these two columns: time and biomass
#' @param LOG10N0 the decimal logarithm of initial biomass
#' @param init_lag initial value for the lag parameter
#' @param init_mumax initial value for the mumax parameter
#' @param init_LOG10Nmax initial value for the LOG10Nmax parameter
#' @param algorithm defaults to "auto" which chooses between bounded and unbounded Levenberg-Marquardt method and the bounded port method
#' @param max_iter max. number of itertaions; defaults to 100
#' @param lower_bound lower.bound for the bounded nls optimisation; defaults to 0
#' @returns lag and the nls fitting object with parameters fitted to logistic model
calc_lag_fit_to_baranyi_with_lag <- function(gr_curve, LOG10N0 = NULL, init_lag = NULL, init_mumax = NULL, init_LOG10Nmax = NULL, algorithm = "auto", max_iter = 100, lower_bound = c(0,0,0, 0)) {
  tryCatch(
    expr =
      {if (algorithm == "auto") {
        nls_a <- choose_lag_fit_algorithm_baranyi(gr_curve,
                                                           LOG10N0 = LOG10N0,
                                                           init_lag = init_lag,
                                                           init_mumax = init_mumax,
                                                           init_LOG10Nmax = init_LOG10Nmax,
                                                           max_iter = max_iter,
                                                           lower_bound = lower_bound)

        # nlsLM( is a more robust version of nls, using  Levenberg-Marquardt algorithm
      } else if (algorithm == "Levenberg-Marquardt") {
        nls_a <- nlsLM(baranyi, gr_curve,
                       list(lag=init_lag, mumax=init_mumax, LOG10N0 = LOG10N0, LOG10Nmax = init_LOG10Nmax),
                       control = nls.control(maxiter = max_iter))
        #lower= lower.bound)
      } else if (algorithm == "port") {
        nls_a <- nls(baranyi, gr_curve,
                     list(lag=init_lag, mumax=init_mumax, LOG10N0 = LOG10N0, LOG10Nmax = init_LOG10Nmax),
                     algorithm = "port",
                     control = nls.control(maxiter = max_iter))
        #lower= lower.bound)
      } else {
        nls_a <- nls(baranyi, gr_curve,
                     list(lag=init_lag, mumax=init_mumax, LOG10N0 = LOG10N0, LOG10Nmax = init_LOG10Nmax),
                     algorithm = algorithm,
                     control = nls.control(maxiter = max_iter))
      }
        lag_N <- coef(nls_a)[names(coef(nls_a)) == "lag"] %>% unname()
      },
    error = function(cond) {
      lag_N <<- NA
      nls_a <<- NA
    })
  return(list(lag_N = lag_N, nls_a = nls_a))
}



#' get_init_pars_logistic
#'
#' Finds reasonable approximation for logistic growth curve parameters (K, lag. growth rate) based on the growth curve and some initial values
#' These approximations will be used as the initial values for the proper optimization algorithm run later.
#' @param data_this_curve data from one specific growth curve with these two columns: time and biomass
#' @param this_n0 the initial biomass
#' @param init_K initial value for the saturation parameter K
#' @param init_lag initial value for the lag parameter
#' @param init_gr_rate initial value for the growth rate
#' @param min_b defaults to 0.2; mina and minb define where to look for exponential phase: it will be where the biomass is between min + (max-min)*(min_a TO min_b)
#' @param min_a defaults to 0.8
#' @returns list of parameters: init_K, init_lag, init_gr_rate,
get_init_pars_logistic = function(data_this_curve, this_n0, init_K, init_lag, init_gr_rate, min_b = 0.2, min_a = 0.8) {
  if (is.null(init_K)) {
    max_this_data <- max(data_this_curve$biomass, na.rm = TRUE)
    init_K <- max_this_data %>% as.numeric()
  }
  if (is.null(init_lag)) {
    init_lag <- calc_lag(data_this_curve, method = "tangent",pars = get_def_pars()) %>% pull(lag) %>% unique() %>% as.numeric()
  }

  if (is.null(init_gr_rate)) {
    data_this_curve_exp <- data_this_curve %>%
      mutate(
        max_biomass = max(data_this_curve$biomass),
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
  return(list(init_K = init_K, init_lag = init_lag, init_gr_rate = init_gr_rate))
}



#' calc_lagistic_fit_lag
#'
#' Calculates lag based on fitting logistic model to data
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param n0 a data frame describing initial biomass for each of the curves, i.e. it has two obligatory columns: "curve_id", "N0"
#' @param init_gr_rate initial value for the growth rate, defaults to NULL in which case it will be approximated based on the data
#' @param init_K initial value for the saturation parameter K, defaults to NULL in which case it will be approximated based on the data
#' @param init_lag initial value for the lag parameter, defaults to NULL in which case it will be approximated based on the data
#' @param algorithm eg. "auto", "Levenberg-Marquardt", "port"
#' @param max_iter Maximum number of iterations
#' @param return_all_params defaults to FALSE, TRUE if you also want to get K and growth.rate apart from lag
#' @param min_b defaults to 0.2; mina and minb define where to look for exponential phase: it will be where the biomass is between min + (max-min)*(lower.bound.for.gr TO upper.bound.for.gr)
#' @param min_a defaults to 0.8
#' @returns growth curve data with additional columns  ('lag', and predicted biomass 'predicted'), and the fitting object if return.all.params was set to TRUE
calc_lagistic_fit_lag <- function(data, n0, init_gr_rate = NULL, init_K = NULL, init_lag = NULL, algorithm, max_iter, return_all_params = FALSE,
                                      min_b = 0.2, min_a = 0.8) {
  if (!("curve_id" %in% names(data))) {
    data$curve_id = NA
  }

  i <- 0
  fiting_list <- list()
  data_new <- data %>% filter(FALSE) %>% mutate( time = numeric(0), biomass = numeric(0),curve_id = character(0), lag = numeric(0), log_info = character(0))
  for (this_curve_id in unique(data$curve_id)) {
    i <- i + 1
    data_this_curve <- data %>% filter(curve_id == this_curve_id) %>% select(time, biomass, curve_id)
    this_n0 <- n0 %>% filter(curve_id == this_curve_id) %>% pull(n0)
    init_pars <- get_init_pars_logistic(data_this_curve, this_n0, init_K, init_lag, init_gr_rate)
    this_fit_obj <- calc_lag_fit_to_logistic_with_lag(gr_curve = data_this_curve,
                                                                     n0 = this_n0,
                                                                     init_gr_rate = init_pars$init_gr_rate,
                                                                     init_K = init_pars$init_K,
                                                                     init_lag = init_pars$init_lag,
                                                                     algorithm = algorithm,
                                                                     max_iter = max_iter)
    data_this_curve <- data_this_curve %>%
      mutate(lag = round(this_fit_obj$lag_N, 1))
    data_this_curve$predicted = if (!any(is.na(this_fit_obj$nls_a))) {predict(this_fit_obj$nls_a, data_this_curve) } else {data_this_curve$predicted = NA}
    data_new <- rbind(data_new, data_this_curve)

    fiting_list[[i]] <- this_fit_obj$nls_a
  }

  if (!return_all_params) {
    return(data_new)
  } else {
    return(list(data_new = data_new, mod_fit = fiting_list))
  }
}



#' get_init_pars_baranyi
#'
#' Finds reasonable approximation for baranyi growth curve parameters (init_mumax, lag) based on the growth curve and some initial values
#' These approximations will be used as the initial values for the proper optimization algorithm run later.
#' @param data_this_curve data from one specific growth curve with these two columns: time and biomass
#' @param this_n0 the initial biomass
#' @param init_lag initial value for the lag parameter
#' @param init_gr_rate initial value for the growth rate
#' @param min_b defaults to 0.2; mina and minb define where to look for exponential phase: it will be where the biomass is between min + (max-min)*(mina TO minb)
#' @param min_a defaults to 0.8
#' @returns list of parameters: init_mumax, init_lag
get_init_pars_baranyi <- function(data_this_curve, this_n0, init_lag, init_gr_rate, min_b = 0.2, min_a = 0.8) {
  if (is.null(init_lag)) {
    init_lag <- calc_lag(data_this_curve, method = "tangent", pars = get_def_pars()) %>% pull(lag) %>% unique() %>% as.numeric()
  }

  if (is.null(init_gr_rate)) {
    data_this_curve_exp <- data_this_curve %>%
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
  return(list(init_mumax = init_mumax, init_lag = init_lag))
}



#' calc_baranyi_fit_lag
#'
#' Calculates lag based on fitting baranyi model to data
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param n0 a data frame describing initial biomass for each of the curves, i.e. it has two obligatory columns: "curve_id", "N0"
#' @param init_lag initial value for the lag parameter, defaults to NULL in which case it will be approximated based on the data
#' @param init_gr_rate initial value for the growth rate, defaults to NULL in which case it will be approximated based on the data
#' @param algorithm eg. "auto", "Levenberg-Marquardt", "port", defaults to "auto"
#' @param max_iter Maximum number of iterations, defaults to 100
#' @returns growth curve data with additional columns  ('lag', and predicted biomass 'predicted')
calc_baranyi_fit_lag <- function(data, n0, init_lag = NULL, init_gr_rate = NULL, algorithm = "auto", max_iter = 100) {
  data_new <- data %>% filter(FALSE) %>% mutate(lag = numeric(0), predicted = numeric(0))
  for (this_curve_id in unique(data$curve_id)) {
    data_this_curve <- data %>% filter(curve_id == this_curve_id)
    this_n0 <- n0 %>% filter(curve_id == this_curve_id) %>% pull(n0)

    data_this_curve_for_model <- data_this_curve %>%
      mutate(LOG10N = log10(biomass), t = time) %>%
      select(LOG10N, t)
    init_LOG10N0 <- log10(n0 %>% filter(curve_id == this_curve_id) %>% pull(n0))
    init_LOG10Nmax <- max(data_this_curve_for_model$LOG10N)
    init_pars_baranyi <- get_init_pars_baranyi(data_this_curve, this_n0, init_lag, init_gr_rate)

    fit_obj_this_curve <- calc_lag_fit_to_baranyi_with_lag(gr_curve = data_this_curve_for_model,
                                                                          LOG10N0 = init_LOG10N0,
                                                                          init_lag = init_pars_baranyi$init_lag,
                                                                          init_mumax = init_pars_baranyi$init_mumax,
                                                                          init_LOG10Nmax = init_LOG10Nmax,
                                                                          algorithm = algorithm,
                                                                          max_iter = max_iter,
                                                                          lower_bound = c(0,0,0, 0))
    data_this_curve <- data_this_curve %>%
      mutate(lag = round(fit_obj_this_curve$lag_N,1))
    data_this_curve$predicted <- if (!any(is.na(fit_obj_this_curve$nls_a))) {10^(predict(fit_obj_this_curve$nls_a, data_this_curve)) } else {data_this_curve$predicted = NA}
    data_new = rbind(data_new, data_this_curve)
  }
  return(data_new)
}



############################################# exported functions #########################################
#' make_grwoth_curve_df
#'
#' Create a growth curve data frame that can be later passed to the lag clculation functions
#' @param time numeric vector of times when biomass was measured  (chronologically ordered as in growth curve)
#' @param biomass numeric vector of measured biomass values  (chronologically ordered as in growth curve)
#' @param curve_id character vector of growth curve identifiers (i.e. if there are multiple measurements done at the same time point, they should have different curve_id)
#' @returns a data frame representing growth curve data
#' @export
make_grwoth_curve_df <- function(time, biomass, curve_id = NULL) {

  if(!(is.numeric(time) | is.integer(time))) stop("Time must be a numeric vector")
  if(!(is.numeric(biomass) | is.integer(biomass))) stop("Biomass must be a numeric vector")
  length_time <- length(time)
  length_biomass <- length(biomass)

  if(length_time != length_biomass)
    stop( "time and biomass vectors must be of the same length")

  if(!is.null(curve_id) & length(curve_id) != length_biomass)
    stop("time and curve_id vectors must be of the same length")

  if (is.null(curve_id)) {
    df = data.frame(time = time, biomass = biomass)
  } else {
    df = data.frame(time = time, biomass = biomass, curve_id = curve_id)
  }
  return(df)
}



#' plot_data
#'
#' Plots the provided growth curve (one single growth curve) on logarithmic scale
#' @param data_new a data frame with two required columns names: "time" and "biomass"
#' @param log10_transform if to plot y axis (biomass) on log10 scale
#' @returns ggplot object with a growth curve
#' @export
plot_data <- function(data_new, log10_transform = TRUE) {
  if (log10_transform) {
  data_new <- data_new %>%
    mutate(biomass = log10(biomass))
  }
  g <- ggplot(data_new, aes(x = time, y = biomass))  +
    geom_line(col = "blue") +
    geom_point(col = "blue") +
    xlab("time") +
    theme(axis.text.y.right = element_text(colour="black"),
          axis.text.y = element_text(colour="blue"),
          axis.title.y = element_text(colour="blue"),
          axis.title.y.right = element_text(colour="black"))

  if (log10_transform) {
    g = g + ylab("Log10(biomass)")
  } else {
    g = g +ylab("Biomass")
  }
  return(g)

}



#' fit_exp_lag
#'
#' Fits the lag to multiple growth curves based on the basic tangent method
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param tangent_method "local.regression" (if the tangent is fitted to a number of points around the maximal growth rate)
#' or "to.point" (if the tangent is fitted only to the point where the growth rate is maximal); defaults to "to.point"
#' @param n0 the initial biomass (a tangent line crossing N0 line will determine the lag)
#' @param curve_points if tangent_method = "local.regression" then curve_points is the number of points the line is fitted to;
#' defaults to 3 i.e. the point with the maximal uptake rate one point before and one point aftter
#' @returns growth curve data (as input) together with additional columns: lag, line.intercept and line.slope
#' @export
fit_exp_lag <- function(data, tangent_method, n0, curve_points = 3) {
  if (!("curve_id" %in% names(data))) {
    data$curve_id = NA
  }
  data_new <- data %>% filter(FALSE) %>% mutate(lag = numeric(0),line_intercept = numeric(0), line_slope = numeric(0))
  for (this_curve_id in unique(data$curve_id)) {
    data_this_curve <- data %>% filter(curve_id == this_curve_id)
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
  # data.new has columns with lag, intercept, slope, log.N0
  return(data_new)
}



#' lag_biomass_incr
#'
#' Fits the lag to multiple growth curves based on the biomass increase method
# It Find where biomass increases by a predefined threshold from N0
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param threshold A value of the biomass increase that we can surely associate with the end of the lag phase rather than random variation during the lag
#' @param n0 the initial biomass (lag will be defined as the time point where the difference between biomass and N0 reaches a predefined threshold)
#' @returns growth curve data (as input) together with additional columns: N0, increase.from.N0, lag
#' @export
lag_biomass_incr <- function(data, threshold, n0) {
  data_new <- data %>% filter(FALSE) %>% mutate(lag = numeric(0))
  for (this_curve_id in unique(data$curve_id)) {
    data_this_curve <- data %>%
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
  return(data_new)
}



#' fit_max_infl_lag
#'
#' Fits the lag to multiple growth curves based on the max growth acceleration method
#' It finds where the second derivative is the largest
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @returns growth curve data (as input) together with additional columns: lag, log.biomass, time.diff, time.av, second.deriv.b, biomass.increase
#' @export
fit_max_infl_lag <- function(data) {
  #deriv <- function(t, x) diff(x) / diff(t)
  #middle_pts <- function(x) x[-1] - diff(x) / 2
  #second_deriv <- function(t,x) deriv(middle_pts(t), deriv(t, x))
  data_new <- data %>% filter(FALSE) %>% mutate(log_biomass = numeric(0), diff = numeric(0), lag = numeric(0))
  for (this_curve_id in unique(data$curve_id)) {
    data_this_curve <- data %>%
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
  return(data_new)
}



#' plot_lag_fit
#'
#' Plots the provided growth curve (one single growth curve) together with the calculated lag and and the rationale for lag calculation
#' @param data_new a data frame output by Calculate.Lag function: it needs to have the following columns: "time", "biomass", "tangent.point", "predicted.data", "threshold", "N0", "second.deriv.b", "line.intercept", "line.slope"
#' @param log10_transform if to plot y axis (biomass) on log10 scale
#' @param print_lag_info if set to "TRUE" prints the lag length on the graph
#' @returns ggplot object with a growth curve
#' @export
plot_lag_fit <- function(data_new, print_lag_info = TRUE, log10_transform = TRUE) {
  max_narm <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
  min_narm <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
  lag.method = unique(data_new$lag_calculation_method)

  data_new <- data_new %>%
    group_by(curve_id) %>%
    mutate(x_mid = mean(time),
           lag_info = paste0("Lag = ", round(lag, 3), "."),
           max_second_deriv_b = max_narm(second_deriv_b),
           min_second_deriv_b = min_narm(second_deriv_b)) %>%
    ungroup()

  if (log10_transform) {
    data_new <- data_new %>%
      mutate(
        tangent_point = log10(tangent_point),
        biomass = log10(biomass),
        predicted_data = log10(predicted_data),
        threshold = log10(threshold),
         n0 = log10(exp(log(n0))),
        line_intercept = line_intercept/log(10),
        line_slope = line_slope/log(10))
  }

  data_new <- data_new %>%
    mutate(
        y_max_for_curve = max(biomass),
        y_min_for_curve = min(biomass),
        text_y = 1.005*y_max_for_curve,
        second_deriv_b_scaled = (second_deriv_b - min_second_deriv_b)/(max_second_deriv_b - min_second_deriv_b)*(y_max_for_curve - y_min_for_curve) + y_min_for_curve
           #y.limit = 1.1*y.max.for.curve
    ) %>%
    ungroup() %>%
    mutate(min_N0 = min(biomass),
           max_N0 = max(biomass))


  max_time <- max(data_new$time)
  size_n0_line <- 1
  size_lag_line <- 1
  g <- ggplot(data_new, aes(x = time, y = biomass))  +
    geom_vline(aes(xintercept = lag), linewidth = size_lag_line, col = "red", linetype = "dashed") +
    geom_line(col = "blue") +
    xlab("time") +
    xlim(c(0, max_time)) +
    facet_grid(curve_id ~ lag_calculation_method, scales = "free_y") +
    theme_bw() +
    theme(axis.text.y.right = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.title.y.right = element_text(colour = "black"))

  if (print_lag_info) {
    g <- g +  geom_text(aes(x = x_mid, y = text_y, label = lag_info), size = 6, col = "red")
  }
  if (log10_transform) {
    g = g + ylab("log10(biomass)")
  } else {
    g = g + ylab("Biomass")
  }



  if (length(lag.method) == 1 & lag.method == "biomass increase") {
    g = g +
      geom_line(aes(y = threshold), col = "darkgreen", na.rm = TRUE)
  }
  if (length(lag.method) == 1 & lag.method == "tangent") {
    g = g +
      geom_hline(aes(yintercept = n0), linewidth = size_n0_line, col = "black", na.rm = TRUE) +
      geom_point(aes(y = tangent_point), col = "darkgreen", size = 2, na.rm = TRUE)
    # on non-log scale it doesn't make sense to plot the tangent line because it is tangent to the log of biomass
    if (log10_transform) {
      g = g +
        geom_abline(aes(intercept = line_intercept, slope = line_slope), col = "darkgreen", na.rm = TRUE)
    }
  }
  if (length(lag.method) == 1 & lag.method == "parameter fitting to a model" & any(!is.na(data_new$lag))) {
    g = g +
      geom_line(aes(y = predicted_data), col = "darkgreen", na.rm = TRUE)
  }
  if (length(lag.method) == 1 & lag.method == "max growth acceleration") {
    g = g +
      geom_line(aes(y = second_deriv_b_scaled), col = "darkgreen", alpha = 0.5, na.rm = TRUE)
  }

  return(g)
}



#' get_lag
#'
#' The most basic function that calculates lags based on growth curve data, selected method and parameters.
#' It uses calc_lag function and strips the results to only get lag parameter for each growth curve id.
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param method method of lag calculation, choose one of the follwoing: "exponential", "biomass increase", "max growth acceleration", "parameter fitting to a model"
#' @param pars a list of parameters. Get.default.parameters function can be used to get the default ones. Otherwise create your onwn list with the following names:
#' - model: if method = "parameter fitting to a model" , one of the following models needs to be chosen: "logistic", "baranyi"
#' - n0_method: first.observation" if the first point is taken as the initial biomass or
#' "minimal.observation" if the minimal biomass is taken is the initial point.
#' In "healthy" growth curves these options should be equivalent
#' but sometimes a drop in OD/biomass is observed at the beginning of a growth curve.
#' In this case it is not obvious what to assume the initial biomass is.
#' - tangent_method "local.regression" (if the tangent is fitted to a number of points around the maximal growth rate)
#' or "to.point" (if the tangent is fitted only to the point where the growth rate is maximal); defaults to "to.point"
#' - threshold: A value of the biomass increase that we can surely associate with the end of the lag phase rather than random variation durinh the lag. Defaults to 10^2
#' - curve_points: if tangent.method = "local.regression" then curve_points is the number of points the line is fitted to;
#' defaults to 3 i.e. the point with the maximal uptake rate one point before and one point after
#' - init_gr_rate: if logistic model is fitted. Defaults to  NULL in which case the initial value will be based on the data
#' - init_lag: if a logistic model is fitted, Defaults to NULL in which case the initial value will be based on the data
#' - algorithm: if method = "parameter fitting to a model", nls algorithm to run the model fit; defaults to "auto" which will choose the best between bounded and unbounded "Levenberg-Marquardt" and bounded "port"
#' - max_iter =  if method = "parameter fitting to a model", the maximum number of nls iterations, defaults to  100
#' @returns lag per each curve_id
#' @export
get_lag <- function(data, method, pars) {
  data_extended = calc_lag(data, method, pars)
  lags_data = data_extended %>%
    group_by(curve_id) %>%
    summarise(lag = unique(lag)) %>%
    ungroup()
  return(lags_data)
}




#' calc_lag
#'
#' The main function that calculates lags based on growth curve data, selected method and parameters and returns an extended growth rate data frame (extended by multiple columns with parameters related to lag calculation)
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param method method of lag calculation, choose one of the follwoing: "exponential", "biomass increase", "max growth acceleration", "parameter fitting to a model"
#' @param pars a list of parameters. Get.default.parameters function can be used to get the default ones. Otherwise create your onwn list with the following names:
#' - model: if method = "parameter fitting to a model" , one of the following models needs to be chosen: "logistic", "baranyi"
#' - n0_method: first.observation" if the first point is taken as the initial biomass or
#' "minimal.observation" if the minimal biomass is taken is the initial point.
#' In "healthy" growth curves these options should be equivalent
#' but sometimes a drop in OD/biomass is observed at the beginning of a growth curve.
#' In this case it is not obvious what to assume the initial biomass is.
#' - tangent_method "local.regression" (if the tangent is fitted to a number of points around the maximal growth rate)
#' or "to.point" (if the tangent is fitted only to the point where the growth rate is maximal); defaults to "to.point"
#' - threshold: A value of the biomass increase that we can surely associate with the end of the lag phase rather than random variation durinh the lag. Defaults to 10^2
#' - curve_points: if tangent.method = "local.regression" then curve_points is the number of points the line is fitted to;
#' defaults to 3 i.e. the point with the maximal uptake rate one point before and one point after
#' - init_gr_rate: if logistic model is fitted. Defaults to  NULL in which case the initial value will be based on the data
#' - init_lag: if a logistic model is fitted, Defaults to NULL in which case the initial value will be based on the data
#' - algorithm: if method = "parameter fitting to a model", nls algorithm to run the model fit; defaults to "auto" which will choose the best between bounded and unbounded "Levenberg-Marquardt" and bounded "port"
#' - max_iter =  if method = "parameter fitting to a model", the maximum number of nls iterations, defaults to  100
#' @returns growth curve data (time, biomass, curve_id) with the following additional columns:  log_biomass, lag, line_slope, line_intercept, lag_calc_method, predicted_data, diff, second_deriv_b, tangent_point, threshold
#' @export
calc_lag <- function(data, method, pars) {
  if (!("curve_id" %in% names(data))) {
    data$curve_id = "growth.curve"
  }

  n0 <- data %>%
    group_by(curve_id) %>%
    arrange(time) %>%
    summarise(n0 = get_n0(biomass, pars$n0_method)) %>%
    ungroup() %>%
    mutate(log_n0 = log(n0))


  if (method == "tangent") {
    sel_tangent_method <- pars$tangent_method
    data_new <- fit_exp_lag(data,
                            n0 = n0,
                            tangent_method = sel_tangent_method,
                            curve_points = pars$curve_points)
    data_new <- data_new %>%
      select(time, biomass, curve_id, lag, line_slope, line_intercept, tangent_point) %>%
      mutate(lag_calculation_method = "tangent",
             log_biomass = log(biomass),
             predicted_data = NA,
             diff = NA,
             second_deriv_b = NA,
             threshold = NA) %>%
      select(time, biomass, log_biomass, curve_id, lag, line_slope, line_intercept, lag_calculation_method, predicted_data, diff, second_deriv_b, tangent_point, threshold)


  } else if (method == "biomass increase") {
    sel_threshold <- pars$threshold
    data_new <- lag_biomass_incr(data,
                                 threshold = sel_threshold,
                                 n0 = n0)
    data_new <- data_new %>%
      select(time, biomass, curve_id, lag) %>%
      mutate(lag_calculation_method = "biomass increase",
             log_biomass = log(biomass),
             predicted_data = NA,
             second_deriv_b = NA,
             line_intercept = NA,
             line_slope = NA,
             tangent_point = NA,
             diff = NA) %>%
      left_join(n0, by = "curve_id") %>%
      mutate(threshold = n0 + sel_threshold) %>%
      select(time, biomass, log_biomass, curve_id, lag, line_slope, line_intercept, lag_calculation_method, predicted_data, diff, second_deriv_b, tangent_point, threshold)

  } else if (method == "max growth acceleration") {
    data_new <-  fit_max_infl_lag(data)
    data_new <- data_new %>%
      select(time, biomass, log_biomass, curve_id, lag, second_deriv_b) %>%
      mutate(lag_calculation_method = "max growth acceleration",
             line_intercept = NA,
             line_slope = NA,
             predicted_data = NA,
             diff = NA,
             tangent_point = NA,
             threshold= NA) %>%
      select(time, biomass, log_biomass, curve_id, lag, line_slope, line_intercept, lag_calculation_method, predicted_data, diff, second_deriv_b, tangent_point, threshold)

  } else if (method == "parameter fitting to a model") {
    sel_model <- pars$model
    if (sel_model == "logistic") {
      data_new <- calc_lagistic_fit_lag(data, n0,
                                        init_gr_rate = pars$init_gr_rate,
                                        init_K = pars$init_K,
                                        init_lag = pars$init_lag,
                                        algorithm = pars$algorithm,
                                        max_iter = pars$max_iter)

    } else if (sel_model == "baranyi") {
      data_new <- calc_baranyi_fit_lag(data,
                                       n0,
                                       init_lag = pars$init_lag,
                                       init_gr_rate = pars$init_gr_rate,
                                       algorithm = pars$algorithm,
                                       max_iter = pars$max_iter) %>%
        mutate(lag_calculation_method = "Fitting lagged baranyi")

    } else {
      stop("model not implemented")
    }
    data_new <- data_new %>%
      select(time, biomass, curve_id, lag, predicted_data = predicted) %>%
      mutate(lag_calculation_method = "parameter fitting to a model",
             log_biomass = log(biomass),
             diff = NA,
             second_deriv_b = NA,
             line_intercept = NA,
             line_slope = NA,
             tangent_point = NA,
             threshold = NA) %>%
      select(time, biomass, log_biomass, curve_id, lag, line_slope, line_intercept, lag_calculation_method, predicted_data, diff, second_deriv_b, tangent_point, threshold)

  } else {
    stop("method not implemented")
  }
  data_new$lag[data_new$lag < 0] = NA
  data_new <- data_new %>%
    left_join(n0)# %>%
  #data.new$lag[data.new$lag < 0] = 0
  return(data_new)
}


#' get_all_methods_lag
#'
#' Runs the main function that calculates lags based on growth curve data based on all possible methods.
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param biomass_incr_threshold A value of the biomass increase that we can surely associate with the end of the lag phase rather than random variation during the lag.
#' Needs to be set specifically to avoid unconscious use of the value set by default. If set to NULL, the value from pars will be taken
#' @param pars a list of parameters. defaults to the ones set by get_def_pars function. Otherwise create your own list with the following names:
#' - model: if method = "parameter fitting to a model" , one of the following models needs to be chosen: "logistic", "baranyi"
#' - n0_method: first.observation" if the first point is taken as the initial biomass or
#' "minimal.observation" if the minimal biomass is taken is the initial point.
#' In "healthy" growth curves these options should be equivalent
#' but sometimes a drop in OD/biomass is observed at the beginning of a growth curve.
#' In this case it is not obvious what to assume the initial biomass is.
#' - tangent.method "local.regression" (if the tangent is fitted to a number of points around the maximal growth rate)
#' or "to.point" (if the tangent is fitted only to the point where the growth rate is maximal); defaults to "to.point"
#' - threshold: A value of the biomass increase that we can surely associate with the end of the lag phase rather than random variation during the lag. Defaults to 10^2
#' - curve_points: if tangent_method = "local.regression" then curve_points is the number of points the line is fitted to;
#' defaults to 3 i.e. the point with the maximal uptake rate one point before and one point after
#' - init_growth.rate: if logistic model is fitted. Defaults to  NULL in which case the initial value will be based on the data
#' - init_lag: if a logistic model is fitted, Defaults to NULL in which case the initial value will be based on the data
#' - algorithm: if method = "parameter fitting to a model", nls algorithm to run the model fit; defaults to "auto" which will choose the best between bounded and unbounded "Levenberg-Marquardt" and bounded "port"
#' - max_iter =  if method = "parameter fitting to a model", the maximum number of nls iterations, defaults to  100
#' @returns growth curve data (time, biomass, curve_id) with the column: lag_calculation_method, and with the following additional columns:  log_biomass, lag, line_slope, line_intercept, lag_calculation_method, predicted_data, diff, second_deriv_b, tangent_point, threshold
#' Note that each growth curve will appear
#' @export
get_all_methods_lag <- function(data, biomass_incr_threshold, pars = NULL) {
  if (is.null(pars)) {
    pars <- get_def_pars()
  }
  if (is.null(biomass_incr_threshold))  {
    biomass_incr_threshold <- pars$threshold
  }
  pars_logistic <- pars
  pars_logistic$model <- "logistic"
  data_new_logistic <- calc_lag(data = data,
                                    method = "parameter fitting to a model",
                                    pars = pars_logistic) %>%
    mutate(lag_calculation_method = "par. fitting to logistic model")

  pars_baranyi <- pars
  pars_baranyi$model <- "baranyi"
  data_new_baranyi <- calc_lag(data = data,
                                   method = "parameter fitting to a model",
                                   pars = pars_baranyi) %>%
    mutate(lag_calculation_method = "par. fitting to baranyi model")

  data_new_max_infl <- calc_lag(data = data,
                                    method = "max growth acceleration",
                                    pars = pars) %>%
    mutate(lag_calculation_method = "max growth acceleration")

  pars_to_point <- pars
  pars_to_point$tangent_method <- "to.point"
  data_new_exp <- calc_lag(data = data,
                                method = "tangent",
                                pars = pars_to_point
  ) %>%
    mutate(lag_calculation_method = "tangent to max growth point")

  pars_local_regr <- pars
  pars_local_regr$tangent_method <- "local.regression"
  data_new_exp2 <- calc_lag(data = data,
                                 method = "tangent",
                                 pars = pars_local_regr)%>%
    mutate(lag_calculation_method = "tangent to max growth line")

  pars$threshold <- biomass_incr_threshold
  data_new_biominc <-  calc_lag(data = data,
                                    method = "biomass increase",
                                    pars = pars)%>%
    mutate(lag_calculation_method = "biomass increase")


  data_all_with_lag <- data_new_max_infl %>%
    rbind(data_new_exp) %>%
    rbind(data_new_exp2) %>%
    rbind(data_new_biominc) %>%
    rbind(data_new_baranyi) %>%
    rbind(data_new_logistic)

  data_all_with_lag$lag[data_all_with_lag$lag < 0] = NA
  return(data_all_with_lag)
}




#' get_def_pars
#' Set defaults parameters used by calc_lag function
#' @returns list of parameters
#' @export
get_def_pars <- function() {
  pars <- list(model = "logistic",
              n0_method = "first.observation",
              tangent_method = "local.regression",
              threshold = 10^2,
              curve_points = 3,
              init_gr_rate = NULL,#0.1,
              init_lag = NULL, #3.5,
              algorithm = "auto",#"default", "Levenberg-Marquardt",#
              max_iter = 100)
  return(pars)
}



#' smooth_data
#' Smoothens growth curves data
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param smooth_kind kind used for the smooth functions, defaults to "3RS3R"
#' @returns smoothened data
#' @export
smooth_data <- function(data, smooth_kind = "3RS3R") {
  if (!("curve_id" %in% names(data))) {
    data$curve_id <- "growth.curve"
  }
  data_smooth <- data  %>% filter(FALSE)
  for (this_curve_id in unique(data$curve_id)) {
    data_this_curve <- data %>% filter(curve_id == this_curve_id) %>%
      mutate(biomass_smooth = smooth(biomass, kind = smooth_kind)) %>%
      select(time, biomass = biomass_smooth, curve_id)
    data_smooth <- rbind(data_smooth, data_this_curve)
  }
  return(data_smooth)
}



#' cut_the_data
#' Subsets the data frame containing only the observations up to the specified maximum time
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param max_time max. time at which we want to cut the growth curve data
#' @returns cut data
#' @export
cut_the_data <- function(data, max_time) {
  data_short <- data %>% filter(time <= max_time)
  return(data_short)
}



#' get_theme
#'
#' This function sets a ggplot theme without grid.
#' The theme removes the major and minor grid lines, sets a white background with a gray border and adjusts the text size.
#' @param text_size defaults to 12
#' @returns a ggplot theme
#' @export
get_theme <- function(text_size = 12) {
  my_theme <- theme(
    panel.grid.major = element_blank(),#element_line(colour = "black", size = 0.05),
    panel.grid.minor = element_blank(),#element_line(colour = "black", size = 0.05),
    panel.background = element_rect(fill = "white",
                                    colour = "gray",
                                    size = 0.5, linetype = "solid"),
    text = element_text(size = text_size),
    panel.border = element_rect(linetype = "solid", fill = NA, colour = "black")
  )
  return(my_theme)
}

