
############################################# helper functions #########################################
#' Get.N0
#' 
#' Gets the initial biomass to relate to
#' @param biomass vector of biomass (chronologically ordered as in growth curve)
#' @param N0.method "first.observation" if the first point is taken as the initial biomass or 
#' "minimal.observation" if the minimal biomass is taken is the initial point. 
#' In "healthy" growth curves these options should be equivalent 
#' but sometimes a drop in OD/biomass is observed at the beginning of a growth curve. 
#' In this case it is not obvious what to assume the initial biomass is.
#' @returns a value of the inititial biomass (either the first observation or the minimum value depending on the parameter N0.method)
Get.N0 = function(biomass, N0.method) {
  # Get the initial biomass value
  if (N0.method == "first.observation") {
    N0.calc = biomass[1]
  } else if (N0.method == "minimal.observation") {
    N0.calc = min(biomass)
  }
  return(N0.calc)
}



#' Fit.Exponential.Lag.To.Curve
#' 
#' Fits the lag to one growth curve based on the basic tangent method
#' @param data a data frame with two required columns names: "time" and "biomass", 
#' This is data from one growth curve only, one (mean) observation per time
#' @param N0 the initial biomass (a tangent line crossing N0 line will determine the lag)
#' @param tangent.method "local.regression" (if the tangent is fitted to a number of points around the maximal growth rate) 
#' or "to.point" (if the tangent is fitted only to the point where the growth rate is maximal); defaults to "to.point"
#' @param n.points.in.curve if tangent.method = "local.regression" then n.points.in.curve is the number of points the line is fitted to; 
#' defaults to 3 i.e. the point with the maximal uptake rate one point before and one point aftter
#' @returns line.slope: slope of the tangent line, 
#' line.intercept: intercept of the tangent line, 
#' lag: lag,
#' tangent.points: i..e a data frame of all points selected for fitting the line
Fit.Exponential.Lag.To.Curve = function(data, N0, tangent.method = "to.point", n.points.in.curve = 3) {
  data = data %>% as.data.frame() %>%
    mutate(log.biomass = log(biomass), 
           db = log.biomass - dplyr::lag(log.biomass,1),
           dt = time - dplyr::lag(time,1),
           time.diff = (time + dplyr::lag(time))/2,
           central_scheme_db = dplyr::lead(log.biomass) - dplyr::lag(log.biomass,1),
           central_scheme_dt = dplyr::lead(time) - dplyr::lag(time,1),
           db_over_dt = central_scheme_db/central_scheme_dt) #db/dt)
  log.N0 = log(N0)
  # Get the maximum growth rate
  max.growth.rate = max(data$db_over_dt, na.rm = TRUE)
  ind.max.diff = which(data$db_over_dt == max.growth.rate)
  
  if (length(ind.max.diff) > 1) {
    warning("Multiple points with max derivative: taking the first one")
    ind.max.diff = ind.max.diff[1]
  }
  if (tangent.method == "to.point") {
    # 1. take simply one maximal growth rate point and calculate a tangent line to it
    # by definition the slope of this line will be equal to its derivative i.e. the 
    exponential.growth.points = data[ind.max.diff, ]
    time.max.diff = exponential.growth.points$time
    log.biomass.at.max.diff = exponential.growth.points$log.biomass
    # ax + b = y
    # a = max.growth.rate
    # we know one point (x,y) = (time.max.diff, log.biomass.at.max.diff)
    line.slope = max.growth.rate
    line.intercept = log.biomass.at.max.diff - line.slope*time.max.diff
  } else if (tangent.method == "local.regression") {
    # 2. Take N points around the maximal growth rate and fit a line,
    # linearly etrapolate
    n.points.around = floor((n.points.in.curve-1)/2)
    exponential.growth.points = data[(ind.max.diff-n.points.around):(ind.max.diff+n.points.around),]
    mod = lm(log.biomass ~ time, exponential.growth.points)
    line.intercept = unname(mod$coefficients[1])
    line.slope = unname(mod$coefficients[2])
  }
  # at x=lag we have y = log.N0, so lag = (y-b)/a
  lag = (log.N0 - line.intercept)/line.slope
  return(list(line.slope = line.slope, line.intercept = line.intercept, lag = lag, 
              tangent.points = exponential.growth.points %>% select(time, tangent.point = biomass)))
}




#' CompareAlgorithms
#' 
#' Comapres results of 3 objects obtained from running nls 
#' @param nlsresLM.no.bound first object resulting from running nls
#' @param nlsresPORT second object resulting from running nls
#' @param nlsresLM third object resulting from running nls
#' @returns the best fitting object (lowest Res.Sum Sq provided that all coefficients are nonnegative) 
CompareAlgorithms =  function(nlsresLM.no.bound, nlsresPORT, nlsresLM) {
  if (!all(is.na(nlsresLM.no.bound))) {
    s = summary(nlsresLM.no.bound)
    if (all(as.numeric(s$coefficients[,1]) >= 0)) {
      nlsresLM.no.bound = nlsresLM.no.bound
    } else {
      nlsresLM.no.bound = NA
    }
  }
  # first compare two bounded models 
  if (!all(is.na(nlsresPORT)) & !all(is.na(nlsresLM))) {
    # if both bounded models are available compare them and choose the better one
    anova.res = anova(nlsresLM, nlsresPORT)
    better.model.idx = which(anova.res$`Res.Sum Sq` == min(anova.res$`Res.Sum Sq`))
    if (better.model.idx ==1) {nlsres.bounded = nlsresLM} else if (better.model.idx == 2) {nlsres.bounded = nlsresPORT}
  } else if (!all(is.na(nlsresPORT))) {
    nlsres.bounded = nlsresPORT
  } else if (!all(is.na(nlsresLM))) {
    nlsres.bounded = nlsresLM
  } else {
    nlsres.bounded = NA
  }
  
  # the netter bounded model compare to unbounded
  if (!all(is.na(nlsres.bounded)) & !all(is.na(nlsresLM.no.bound))) {
    # if both bounded models are available compare them and choose the better one
    anova.res = anova(nlsres.bounded, nlsresLM.no.bound)
    better.model.idx = which(anova.res$`Res.Sum Sq` == min(anova.res$`Res.Sum Sq`))
    if (better.model.idx ==1) {nlsres = nlsres.bounded} else if (better.model.idx == 2) {nlsres = nlsresLM.no.bound}
  } else if (!all(is.na(nlsres.bounded))) {
    nlsres = nlsres.bounded
  } else if (!all(is.na(nlsresLM.no.bound))) {
    nlsres = nlsresLM.no.bound
  } else {
    nlsres = NA
  }
  return(nlsres)
  
}




#' Choose.Best.Lag.Fitting.Algorithm.Baranyi
#' 
#' Runs nlsLM/nls algorithms with three different parameter setups to fit the best Baranyi parameters to our data and chooses the best model
#' @param growthcurve data from one specific growth curve with the following columns: LOG10N, t
#' @param LOG10N0 init value for the LOG10N0 parameter
#' @param init.lag initial value for the lag
#' @param init.mumax initial value for the mumax parameter
#' @param init.LOG10Nmax initial value for the LOG10Nmax parameter
#' @param maxiter max. number of itertaions
#' @param lower.bound lower.bound for the bounded nls optimisation; 
#' @returns the best nls fitting object with parameters fitted to Baranyi modek (lowest Res.Sum Sq provided that all coefficients are nonnegative) 
Choose.Best.Lag.Fitting.Algorithm.Baranyi = function(growthcurve, 
                                                     LOG10N0,
                                                     init.lag, 
                                                     init.mumax,
                                                     init.LOG10Nmax, 
                                                     maxiter, 
                                                     lower.bound) {

  # choose best from LM and port 
  # Sometimes the lower.bound argument makes the fit much worse.
  tryCatch(
    expr = 
      {nlsresLM = nlsLM(formula = baranyi, 
          data = growthcurve, 
          start = list(lag=init.lag, mumax=init.mumax, LOG10N0 = LOG10N0, LOG10Nmax = init.LOG10Nmax),
          control = nls.control(maxiter = maxiter),
          lower = lower.bound)
  },
    error=function(cond) {
      # this operator assigns value outside the error environment
      nlsresLM <<- NA
    })
  tryCatch(
    expr = 
      {nlsresLM.no.bound = nlsLM(formula = baranyi,
                              data = growthcurve, 
                              start = list(lag=init.lag, mumax=init.mumax, LOG10N0 = LOG10N0, LOG10Nmax = init.LOG10Nmax), 
                              control = nls.control(maxiter = maxiter))
    },
    error=function(cond) {
      nlsresLM.no.bound <<- NA
    })
  tryCatch(
    expr = 
      {nlsresPORT = nls(baranyi, 
                     data = growthcurve, 
                     start = list(lag=init.lag, mumax=init.mumax, LOG10N0 = LOG10N0, LOG10Nmax = init.LOG10Nmax), 
                     algorithm = "port", 
                     control = nls.control(maxiter = maxiter),
                     lower= lower.bound)
    },
    error=function(cond) {
      nlsresPORT <<- NA
    })
  
  nlsres = CompareAlgorithms(nlsresLM.no.bound, nlsresPORT, nlsresLM)
    # consider the model without lower bounds only if the resulted westimates are above 0
    return(nlsres)
    
  }
  
  

#' Choose.Best.Lag.Fitting.Algorithm.Logistic
#' 
#' Runs nlsLM/nls algorithms with three different parameter setups to fit the best Logistic model parameters to our data and chooses the best model
#' @param growthcurve data from one specific growth curve with the following columns: LOG10N, t
#' @param N0 the initial biomass
#' @param init.growth.rate initial value for the growth rate
#' @param init.K initial value for the saturation parameter K
#' @param init.lag initial value for the lag parameter
#' @param maxiter max. number of itertaions; defaults to 100
#' @param lower.bound lower.bound for the bounded nls optimisation; defaults to 0
#' @returns the best nls fitting object with parameters fitted to logistic model (lowest Res.Sum Sq provided that all coefficients are nonnegative) 
Choose.Best.Lag.Fitting.Algorithm.Logistic = function(growthcurve, N0, 
                                             init.growth.rate = init.growth.rate, 
                                             init.K = init.K, 
                                             init.lag = init.lag,
                                             maxiter = 100,
                                             lower.bound = c(0,0,0)) {
  # choose best from LM and port 
  # Sometimes the lower.bound argument makes the fit much worse.
  tryCatch( 
    expr = 
      {nlsresLM = nlsLM(formula = biomass ~ N0 + (time >= lag)*N0*(-1+K*exp(growth.rate*(time-lag))/(K - N0 + N0*exp(growth.rate*(time - lag)))), 
                     data = growthcurve, 
                     start = list(growth.rate = init.growth.rate, K=init.K, lag = init.lag), 
                     control = nls.control(maxiter = maxiter),
                     lower= lower.bound)
        },
    error=function(cond) {
      nlsresLM <<- NA
    })
  tryCatch(
    expr = 
      {nlsresLM.no.bound = nlsLM(formula = biomass ~ N0 + (time >= lag)*N0*(-1+K*exp(growth.rate*(time-lag))/(K - N0 + N0*exp(growth.rate*(time - lag)))), 
                              data = growthcurve, 
                              start = list(growth.rate = init.growth.rate, K=init.K, lag = init.lag), 
                              control = nls.control(maxiter = maxiter))
      },
    error=function(cond) {
      nlsresLM.no.bound <<- NA
    })
  tryCatch(
    expr = 
      {nlsresPORT = nls(biomass ~ N0 + (time >= lag)*N0*(-1+K*exp(growth.rate*(time-lag))/(K - N0 + N0*exp(growth.rate*(time - lag)))), growthcurve, 
                     list(growth.rate = init.growth.rate, K=init.K, lag = init.lag), 
                     algorithm = "port", 
                     control = nls.control(maxiter = maxiter),
                     lower= lower.bound)},
    error=function(cond) {
      nlsresPORT <<- NA
    })
  
  nlsres = CompareAlgorithms(nlsresLM.no.bound, nlsresPORT, nlsresLM)
  # consider the model without lower bounds only if the resulted westimates are above 0
  return(nlsres)
  
}



#' Choose.Best.Lag.Fitting.Algorithm.Logistic
#' 
#' Runs nlsLM/nls algorithm of the user's choice to fit the  Logistic model parameters to our data
#' @param growthcurve data from one specific growth curve with these two columns: time and biomass
#' @param N0 the initial biomass
#' @param init.growth.rate initial value for the growth rate
#' @param init.K initial value for the saturation parameter K
#' @param init.lag initial value for the lag parameter
#' @param algorithm defaults to "auto" which chooses between bounded and unbounded Levenberg-Marquardt method and the bounded port method
#' @param maxiter max. number of itertaions; defaults to 100
#' @param lower.bound lower.bound for the bounded nls optimisation; defaults to 0
#' @returns lag and the nls fitting object with parameters fitted to logistic model
Calculate.Lag.Fitting.To.Logistic.With.Lag = function(growthcurve, N0, 
                                                      init.growth.rate = init.growth.rate, 
                                                      init.K = init.K, 
                                                      init.lag = init.lag,
                                                      algorithm = "auto",# Levenberg-Marquardt", # or "port" for nls
                                                      maxiter = 100,
                                                      lower.bound = c(0,0,0)) {
  tryCatch(
    expr = 
      {if (algorithm == "auto") {
          nlsres = Choose.Best.Lag.Fitting.Algorithm.Logistic(growthcurve, N0, 
                                                              init.growth.rate = init.growth.rate, init.K = init.K, init.lag = init.lag,
                                                              maxiter = maxiter,
                                                              lower.bound = lower.bound)
          
          # nlsLM( is a more robust version of nls, using  Levenberg-Marquardt algorithm
        } else if (algorithm == "Levenberg-Marquardt") {
          nlsres = nlsLM(formula = biomass ~ N0 + (time >= lag)*N0*(-1+K*exp(growth.rate*(time-lag))/(K - N0 + N0*exp(growth.rate*(time - lag)))), 
                         data = growthcurve, 
                         start = list(growth.rate = init.growth.rate, K=init.K, lag = init.lag), 
                         control = nls.control(maxiter = maxiter))
                         #lower= lower.bound)
        } else if (algorithm == "port") {
          nlsres = nls(biomass ~ N0 + (time >= lag)*N0*(-1+K*exp(growth.rate*(time-lag))/(K - N0 + N0*exp(growth.rate*(time - lag)))), growthcurve, 
                       list(growth.rate = init.growth.rate, K=init.K, lag = init.lag), 
                       algorithm = "port", 
                       control = nls.control(maxiter = maxiter))
                       #lower= lower.bound)
        } else {
          nlsres = nls(biomass ~ N0 + (time >= lag)*N0*(-1+K*exp(growth.rate*(time-lag))/(K - N0 + N0*exp(growth.rate*(time - lag)))), growthcurve, 
                       list(growth.rate = init.growth.rate, K=init.K, lag = init.lag), 
                       algorithm = algorithm, 
                       control = nls.control(maxiter = maxiter))
        }
          lagN = coef(nlsres)[names(coef(nlsres)) == "lag"] %>% unname()
          },
    error=function(cond) {
      lagN <<- NA
      nlsres <<- NA
      #print(cond)
    })
  return(list(lagN = lagN, nlsres = nlsres))
}



#' Calculate.Lag.Fitting.To.Baranyi.With.Lag
#' 
#' Runs nlsLM/nls algorithms with three different parameter setups to fit the best Logistic model parameters to our data and chooses the best model
#' @param growthcurve data from one specific growth curve with these two columns: time and biomass
#' @param N0 the initial biomass
#' @param init.growth.rate initial value for the growth rate
#' @param init.K initial value for the saturation parameter K
#' @param init.lag initial value for the lag parameter
#' @param algorithm defaults to "auto" which chooses between bounded and unbounded Levenberg-Marquardt method and the bounded port method
#' @param maxiter max. number of itertaions; defaults to 100
#' @param lower.bound lower.bound for the bounded nls optimisation; defaults to 0
#' @returns lag and the nls fitting object with parameters fitted to logistic model
Calculate.Lag.Fitting.To.Baranyi.With.Lag = function(growthcurve, LOG10N0 = NULL, init.lag = NULL, init.mumax = NULL, init.LOG10Nmax = NULL, algorithm = "auto", maxiter = 100, lower.bound = c(0,0,0, 0)) {
  tryCatch(
    expr = 
      {if (algorithm == "auto") {
          nlsres = Choose.Best.Lag.Fitting.Algorithm.Baranyi(growthcurve, 
                                                             LOG10N0 = LOG10N0, 
                                                             init.lag = init.lag, 
                                                             init.mumax = init.mumax, 
                                                             init.LOG10Nmax = init.LOG10Nmax, 
                                                             maxiter = maxiter, 
                                                             lower.bound = lower.bound)
          
          # nlsLM( is a more robust version of nls, using  Levenberg-Marquardt algorithm
        } else if (algorithm == "Levenberg-Marquardt") {
          nlsres = nlsLM(baranyi, growthcurve, 
                         list(lag=init.lag, mumax=init.mumax, LOG10N0 = LOG10N0, LOG10Nmax = init.LOG10Nmax),
                         control = nls.control(maxiter = maxiter))
                        #lower= lower.bound)
        } else if (algorithm == "port") {
          nlsres = nls(baranyi, growthcurve, 
                       list(lag=init.lag, mumax=init.mumax, LOG10N0 = LOG10N0, LOG10Nmax = init.LOG10Nmax),
                       algorithm = "port", 
                       control = nls.control(maxiter = maxiter))
                      #lower= lower.bound)
        } else {
          nlsres = nls(baranyi, growthcurve, 
                       list(lag=init.lag, mumax=init.mumax, LOG10N0 = LOG10N0, LOG10Nmax = init.LOG10Nmax),
                       algorithm = algorithm, 
                       control = nls.control(maxiter = maxiter))
        }
          lagN = coef(nlsres)[names(coef(nlsres)) == "lag"] %>% unname()
    },
    error=function(cond) {
      lagN <<- NA
      nlsres <<- NA
    })
  return(list(lagN = lagN, nlsres = nlsres))
}



#' Get.Initial.Parameters.Logistic
#' 
#' Finds reasonable approximation for logistic growth curve parameters (K, lag. growth rate) based on the growth curve and some initial values
#' These approximations will be used as the initial values for the proper optimisation algorithm run later.
#' @param data_this_curve data from one specific growth curve with these two columns: time and biomass
#' @param this.N0 the initial biomass
#' @param init.K initial value for the saturation parameter K
#' @param init.lag initial value for the lag parameter
#' @param init.growth.rate initial value for the growth rate
#' @param minb defaults to 0.2; mina and minb define where to look for exponential phase: it will be where the biomass is between min + (max-min)*(mina TO minb)
#' @param mina defaults to 0.8
#' @returns list of parameters: init.K, init.lag, init.growth.rate, 
Get.Initial.Parameters.Logistic = function(data_this_curve, this.N0, init.K, init.lag, init.growth.rate, minb = 0.2, maxb = 0.8) {
  if (is.null(init.K)) {
    max.this.data = max(data_this_curve$biomass, na.rm = TRUE)
    init.K = max.this.data %>% as.numeric()
  }
  if (is.null(init.lag)) {
    init.lag = Calculate.Lag(data_this_curve, method = "exponential",pars = Get.default.parameters()) %>% pull(lag) %>% unique() %>% as.numeric()
  }
  
  if (is.null(init.growth.rate)) {
    data_this_curve_exponential = data_this_curve %>%
      mutate(
        max.biomass = max(data_this_curve$biomass),
        min.threshold = this.N0 + minb*(max.biomass - this.N0),
        max.threshold = this.N0 + maxb*(max.biomass - this.N0)) %>%
      # take only the points that are in between min and max
      filter(biomass <= max.threshold & biomass >= min.threshold)
    data_this_curve_exponential$logdata = log(data_this_curve_exponential$biomass/this.N0)
    if (nrow(data_this_curve_exponential %>% filter(!is.na(time) & !is.na(logdata)) > 0)) {
      mod = lm(logdata ~ time, data = data_this_curve_exponential)
      # this growth rate is assumig an exponential model so it will be generally underestimated
      init.growth.rate = mod$coefficients[2] %>% unname()
      # we have real r = r(1-N/K) so let us take
      Nmid = median(data_this_curve_exponential$biomass)
      init.growth.rate = init.growth.rate/(1-Nmid/init.K)
    } else {
      init.growth.rate = 0.1
    }
    #data_this_curve_exponential$predicted = predict(mod, data_this_curve_exponential)
  }
  return(list(init.K = init.K, init.lag = init.lag, init.growth.rate = init.growth.rate))
}



#' Calculate.Lagistic.Fit.Lag
#' 
#' Calculates lag based on fitting logistic model to data 
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param N0 a data frame describing initial biomass for each of the curves, i.e. it has two obligatory columns: "curve_id", "N0"
#' @param init.growth.rate initial value for the growth rate, defaults to NULL in which case it will be approximated based on the data
#' @param init.K initial value for the saturation parameter K, defaults to NULL in which case it will be approximated based on the data
#' @param init.lag initial value for the lag parameter, defaults to NULL in which case it will be approximated based on the data
#' @param algorithm eg. "auto", "Levenberg-Marquardt", "port"
#' @param maxiter Maximum number of iterations
#' @param return.all.params defaults to FALSE, TRUE if you also want to get K and growth.rate apart from lag
#' @param low.bound.for.gr defaults to 0.2; mina and minb define where to look for exponential phase: it will be where the biomass is between min + (max-min)*(lower.bound.for.gr TO upper.bound.for.gr)
#' @param upper.bound.for.gr defaults to 0.8
#' @returns growth curve data with additional columns  ('lag', and predicted biomass 'predicted'), and the fitting object if return.all.params was set to TRUE
Calculate.Lagistic.Fit.Lag = function(data, N0, init.growth.rate = NULL, init.K = NULL, init.lag = NULL, algorithm, maxiter, return.all.params = FALSE,
                                      low.bound.for.gr = 0.2, upper.bound.for.gr = 0.8) {
  if (!("curve_id" %in% names(data))) {
    data$curve_id = NA
  }

  i=0
  fiting.list = list()
  data.new = data %>% filter(FALSE) %>% mutate( time = numeric(0), biomass = numeric(0),curve_id = character(0), lag = numeric(0), log.info = character(0))
  for (this_curve_id in unique(data$curve_id)) {
    i = i+1
    data_this_curve = data %>% filter(curve_id == this_curve_id) %>% select(time, biomass, curve_id)
    this.N0 = N0 %>% filter(curve_id == this_curve_id) %>% pull(N0)
    initial.parameters = Get.Initial.Parameters.Logistic(data_this_curve, this.N0, init.K, init.lag, init.growth.rate)
    this.fitting.object = Calculate.Lag.Fitting.To.Logistic.With.Lag(growthcurve = data_this_curve, 
                                                          N0=this.N0, 
                                                          init.growth.rate =initial.parameters$init.growth.rate, 
                                                          init.K = initial.parameters$init.K, 
                                                          init.lag = initial.parameters$init.lag, 
                                                          algorithm =algorithm, 
                                                          maxiter = maxiter)
    data_this_curve = data_this_curve %>% 
      mutate(lag = round(this.fitting.object$lagN,1))
    data_this_curve$predicted = if (!any(is.na(this.fitting.object$nlsres))) {predict(this.fitting.object$nlsres, data_this_curve) } else {data_this_curve$predicted = NA}
    data.new = rbind(data.new, data_this_curve)
    
    fiting.list[[i]] = this.fitting.object$nlsres
  }
 
  if (!return.all.params) {
  return(data.new)
  } else {
  return(list(data.new = data.new, modfit = fiting.list))  
  }
}
 


#' Get.Initial.Parameters.Baranyi
#' 
#' Finds reasonable approximation for baranyi growth curve parameters (init.mumax, lag) based on the growth curve and some initial values
#' These approximations will be used as the initial values for the proper optimisation algorithm run later.
#' @param data_this_curve data from one specific growth curve with these two columns: time and biomass
#' @param this.N0 the initial biomass
#' @param init.lag initial value for the lag parameter
#' @param init.growth.rate initial value for the growth rate
#' @param minb defaults to 0.2; mina and minb define where to look for exponential phase: it will be where the biomass is between min + (max-min)*(mina TO minb)
#' @param mina defaults to 0.8
#' @returns list of parameters: init.mumax, init.lag
Get.Initial.Parameters.Baranyi = function(data_this_curve, this.N0, init.lag, init.growth.rate, minb = 0.2, mina = 0.8) {
  if (is.null(init.lag)) {
    init.lag = Calculate.Lag(data_this_curve, method = "exponential",pars = Get.default.parameters()) %>% pull(lag) %>% unique() %>% as.numeric()
  }
  
  if (is.null(init.growth.rate)) {
    data_this_curve_exponential = data_this_curve %>%
      mutate(
        max.biomass = max(biomass),
        min.threshold = this.N0 + minb*(max.biomass - this.N0),
        max.threshold = this.N0 + mina*(max.biomass - this.N0)) %>%
      # take only the points that are in between min and max
      filter(biomass <= max.threshold & biomass >= min.threshold)
    data_this_curve_exponential$logdata = log(data_this_curve_exponential$biomass/this.N0)
    if (nrow(data_this_curve_exponential %>% filter(!is.na(time) & !is.na(logdata))  > 0)) {
      mod = lm(logdata ~ time, data = data_this_curve_exponential)
      # this growth rate is assumig an exponential model so it will be generally underestimated
      init.growth.rate = mod$coefficients[2] %>% unname()
      # we have real r = r(1-N/K) so let us take
      Nmid = median(data_this_curve_exponential$biomass)
      init.mumax = init.growth.rate
    } else {
      init.mumax = 0.1
    }
  } else {
    init.mumax = init.growth.rate
  }
  return(list(init.mumax = init.mumax, init.lag = init.lag))
}



#' Calculate.Baranyi.Fit.Lag
#' 
#' Calculates lag based on fitting baranyi model to data 
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param N0 a data frame describing initial biomass for each of the curves, i.e. it has two obligatory columns: "curve_id", "N0"
#' @param init.lag initial value for the lag parameter, defaults to NULL in which case it will be approximated based on the data
#' @param init.growth.rate initial value for the growth rate, defaults to NULL in which case it will be approximated based on the data
#' @param algorithm eg. "auto", "Levenberg-Marquardt", "port", defaults to "auto"
#' @param maxiter Maximum number of iterations, defaults to 100
#' @returns growth curve data with additional columns  ('lag', and predicted biomass 'predicted')
Calculate.Baranyi.Fit.Lag = function(data, N0, init.lag = NULL, init.growth.rate = NULL, algorithm = "auto", maxiter = 100) {
  data.new = data %>% filter(FALSE) %>% mutate(lag = numeric(0), predicted = numeric(0))
  for (this_curve_id in unique(data$curve_id)) {
    data_this_curve = data %>% filter(curve_id == this_curve_id)
    this.N0 = N0 %>% filter(curve_id == this_curve_id) %>% pull(N0)
    
    data_this_curve_for_model = data_this_curve %>%
      mutate(LOG10N = log10(biomass), t = time) %>%
      select(LOG10N, t)
    init.LOG10N0 = log10(N0 %>% filter(curve_id == this_curve_id) %>% pull(N0))
    init.LOG10Nmax = max(data_this_curve_for_model$LOG10N)
    initial.parameters.baranyi = Get.Initial.Parameters.Baranyi(data_this_curve, this.N0, init.lag, init.growth.rate)
    
    fitting.object.this.curve = Calculate.Lag.Fitting.To.Baranyi.With.Lag(growthcurve = data_this_curve_for_model, 
                                                                          LOG10N0 = init.LOG10N0, 
                                                                          init.lag = initial.parameters.baranyi$init.lag, 
                                                                          init.mumax = initial.parameters.baranyi$init.mumax, 
                                                                          init.LOG10Nmax = init.LOG10Nmax,
                                                                          algorithm = algorithm, 
                                                                          maxiter = maxiter, 
                                                                          lower.bound = c(0,0,0, 0))
    data_this_curve = data_this_curve %>%
      mutate(lag = round(fitting.object.this.curve$lagN,1))
    data_this_curve$predicted = if (!any(is.na(fitting.object.this.curve$nlsres))) {10^(predict(fitting.object.this.curve$nlsres, data_this_curve)) } else {data_this_curve$predicted = NA}
    data.new = rbind(data.new, data_this_curve)
  }
  return(data.new)
}



#' Plot.Data
#' 
#' Plots the provided growth curve (one single growth curve) on logaritmic scale
#' @param data.new a data frame with two required columns names: "time" and "biomass"
#' @returns ggplot object with a growth curve
Plot.Data = function(data.new) {
  data.new = data.new %>%
    mutate(log10.biomass = log10(biomass))
  g = ggplot(data.new)  + 
    geom_line(aes(x= time, y = log10.biomass), col = "blue") +
    geom_point(aes(x= time, y = log10.biomass), col = "blue") +
    xlab("time [h]") +
    ylab("Log10(biomass)") +
    theme(axis.text.y.right=element_text(colour="black"),
          axis.text.y=element_text(colour="blue"),
          axis.title.y=element_text(colour="blue"),
          axis.title.y.right=element_text(colour="black"))    
   return(g)
  
}



############################################# exported functions #########################################
#' Fit.Exponential.Lag
#' 
#' Fits the lag to multiple growth curves based on the basic tangent method 
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param tangent.method "local.regression" (if the tangent is fitted to a number of points around the maximal growth rate) 
#' or "to.point" (if the tangent is fitted only to the point where the growth rate is maximal); defaults to "to.point"
#' @param N0 the initial biomass (a tangent line crossing N0 line will determine the lag)
#' @param n.points.in.curve if tangent.method = "local.regression" then n.points.in.curve is the number of points the line is fitted to; 
#' defaults to 3 i.e. the point with the maximal uptake rate one point before and one point aftter
#' @returns growth curve data (as input) together with additional columns: lag, line.intercept and line.slope
#' @export
Fit.Exponential.Lag = function(data, tangent.method, N0, n.points.in.curve = 3) {
  if (!("curve_id" %in% names(data))) {
    data$curve_id = NA
  }
  data.new = data %>% filter(FALSE) %>% mutate(lag = numeric(0),line.intercept = numeric(0), line.slope = numeric(0))
  for (this_curve_id in unique(data$curve_id)) {
    data_this_curve = data %>% filter(curve_id == this_curve_id)
    this.N0 = N0 %>% filter(curve_id == this_curve_id) %>% pull(N0)
    lag.object = Fit.Exponential.Lag.To.Curve(data_this_curve, this.N0,tangent.method, n.points.in.curve)
    data_this_curve = data_this_curve %>% 
      mutate(lag = round(lag.object$lag,1)) %>%
      mutate(line.intercept = lag.object$line.intercept,
             line.slope = lag.object$line.slope) %>%
      left_join(lag.object$tangent.points)
    data.new = rbind(data.new, data_this_curve)
  }
  
  data.new$predicted = NA
  # data.new has columns with lag, intercept, slope, log.N0
  return(data.new)
}



#' Lag.Based.On.Biomass.Increase
#' 
#' Fits the lag to multiple growth curves based on the biomass increase method
# It Find where biomass increases by a predefined threshold from N0
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param threshold A value of the biomass increase that we can surely associate with the end of the lag phase rather than random variation durin the lag
#' @param N0 the initial biomass (lag will be defined as the time point where the difference bwteem  biomass and N0 reaches a predefined threshold)
#' @returns growth curve data (as input) together with additional columns: N0, increase.from.N0, lag
#' @export
Lag.Based.On.Biomass.Increase = function(data, threshold, N0) {
  data.new = data %>% filter(FALSE) %>% mutate(lag = numeric(0))
  for (this_curve_id in unique(data$curve_id)) {
    data_this_curve = data %>% 
      filter(curve_id == this_curve_id) %>% 
      left_join(N0, by = "curve_id") %>%
      mutate(increase.from.N0 = biomass - N0)
    
    
    #find where second derivative is maximal
    diff.above.threshold = which(data_this_curve$increase.from.N0 >= threshold)
    first.diff.above.threshold = diff.above.threshold[1]
    lag.this.curve = data_this_curve$time[first.diff.above.threshold]
    data_this_curve = data_this_curve %>%
      mutate(
        lag = round(lag.this.curve,1))
    data.new = rbind(data.new, data_this_curve)
  }
  return(data.new)
}



#' Fit.Max.Inflection.Lag
#' 
#' Fits the lag to multiple growth curves based on the max growth acceleration method
#' It finds where the second derivative is the largest
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @returns growth curve data (as input) together with additional columns: lag, log.biomass, time.diff, time.av, second.deriv.b, biomass.increase
#' @export
Fit.Max.Inflection.Lag = function(data) {
  #deriv <- function(t, x) diff(x) / diff(t)
  #middle_pts <- function(x) x[-1] - diff(x) / 2
  #second_deriv <- function(t,x) deriv(middle_pts(t), deriv(t, x))
  data.new = data %>% filter(FALSE) %>% mutate(log.biomass = numeric(0), diff = numeric(0), lag = numeric(0))
  for (this_curve_id in unique(data$curve_id)) {
    data_this_curve = data %>% 
      filter(curve_id == this_curve_id) %>%
      arrange(time) %>%
      as.data.frame() %>%
      mutate(log.biomass = log(biomass), 
             time.diff = mean(c(diff(time), diff(dplyr::lag(time))), na.rm = TRUE),
             time.av = (time + dplyr::lag(time))/2) %>%
      mutate(
        #second.deriv.b = c(NA,second_deriv(time, log.biomass), NA),
        # central scheme
        second.deriv.b = (dplyr::lead(log.biomass) + dplyr::lag(log.biomass) - 2*log.biomass)/time.diff^2,
        # we only look at second derivative if we know the first derivative is positive! 
        # In empirical data we sometimes see a decrease in biomass but we don;t want to look at the seconf derivative there!
        biomass.increase=dplyr::lead(log.biomass)  > dplyr::lag(log.biomass)
      ) 
    
    #find where second derivative is maximal
    max.second.deriv.b = max(data_this_curve$second.deriv.b, na.rm=TRUE)
    # take first point when max serivative if there are multiple
    ind.max.second.deriv.b = which(data_this_curve$second.deriv.b == max.second.deriv.b)[1]
    lag.this.curve = data_this_curve$time[ind.max.second.deriv.b]
    data_this_curve = data_this_curve %>%
      mutate(
        lag = round(lag.this.curve,1))
    data.new = rbind(data.new, data_this_curve)
  }
  return(data.new)
}



#' Plot.Lag.Fit
#' 
#' Plots the provided growth curve (one single growth curve) together with the calculated lag and and the rationale for lag calculation
#' @param data.new a data frame output by Calculate.Lag function: it needs to have the following columns: "time", "biomass", "tangent.point", "predicted.data", "threshold", "N0", "second.deriv.b", "line.intercept", "line.slope"
#' @param print.lag.info if set to "TRUE" prints the lag length on the graph
#' @returns ggplot object with a growth curve
#' @export
Plot.Lag.Fit = function(data.new, print.lag.info = TRUE) {
  data.new = data.new %>%
    group_by(curve_id) %>%
    mutate(x.mid = mean(time),
           lag.info = paste0("Lag = ", round(lag, 3), " [h]."),
           log.biomass = log(biomass),
           log.10.tangent.point = log10(tangent.point),
           log.10.biomass = log10(biomass),
           log.10.predicted = log10(predicted.data),
           log.10.threshold = log10(threshold),
           y.max.for.curve = max(log.10.biomass),
           y.min.for.curve = min(log.10.biomass),
           log10N0 = log10(exp(log(N0))),
           text.y = 1.005*y.max.for.curve,
           max.second.deriv.b = max(second.deriv.b, na.rm = TRUE),
           min.second.deriv.b = min(second.deriv.b, na.rm = TRUE),
           second.deriv.b.scaled = (second.deriv.b - min.second.deriv.b)/(max.second.deriv.b - min.second.deriv.b)*(y.max.for.curve - y.min.for.curve) + y.min.for.curve
           #y.limit = 1.1*y.max.for.curve
           ) %>%
    ungroup() %>%
    mutate(min.log10N0 = min(log.10.biomass),
           max.log10N0 = max(log.10.biomass),
           log10.intercept = line.intercept/log(10),
           log10.slope = line.slope/log(10))
  
  max.time = max(data.new$time)
  
  size.N0.line = 1
  size.lag.line = 1
  g = ggplot(data.new)  + 
    geom_vline(aes(xintercept = lag), size = size.lag.line, col = "red", linetype = "dashed") +
    geom_line(aes(x= time, y = log.10.biomass), col = "blue") +
    #geom_point(aes(x= time, y = log10.biomass), col = "blue") +
    geom_point(aes(x= time, y = log.10.tangent.point), col = "darkgreen", size = 2) +
    geom_line(aes(x= time, y = log.10.predicted), col = "darkgreen") +
    geom_line(aes(x= time, y = log.10.threshold), col = "darkgreen") +
    geom_line(aes(x=time, y = second.deriv.b.scaled), col = "darkgreen", alpha = 0.5) +
    geom_hline(aes(yintercept = log10N0), size = size.N0.line, col = "black") +
    geom_abline(aes(intercept = log10.intercept, slope = log10.slope), col = "darkgreen") + 
    xlab("time [h]") +
    xlim(c(0, max.time)) +
    ylab("Log10(biomass)") + 
    #ylim(c(min(data.new$log.10.biomass), max(data.new$log.10.biomass))) +
    facet_grid(curve_id~lag.calculation.method, scales = "free_y") +
  theme(axis.text.y.right=element_text(colour="black"),
          axis.text.y=element_text(colour="blue"),
          axis.title.y=element_text(colour="blue"),
          axis.title.y.right=element_text(colour="black"))    
  
  if (print.lag.info) {
    g = g +  geom_text(aes(x=x.mid, y = text.y, label = lag.info), size = 6, col = "red")
  }
  return(g)
}



#' Calculate.Lag
#' 
#' The main function that calculates lags based on growth curve data, selected method and parameters
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param method method of lag calculation, choose one of the follwoing: "exponential", "biomass increase", "max growth acceleration", "parameter fitting to a model" 
#' @param pars a list of parameters. Get.default.parameters function can be used to get the default ones. Otherwise create your onwn list with the following names: \n
#' - model: if method = "parameter fitting to a model" , one of the following models needs to be chosen: "logistic", "baranyi" \n
#' - N0.method: first.observation" if the first point is taken as the initial biomass or 
#' "minimal.observation" if the minimal biomass is taken is the initial point. 
#' In "healthy" growth curves these options should be equivalent 
#' but sometimes a drop in OD/biomass is observed at the beginning of a growth curve. 
#' In this case it is not obvious what to assume the initial biomass is. \n
#' - tangent.method "local.regression" (if the tangent is fitted to a number of points around the maximal growth rate) 
#' or "to.point" (if the tangent is fitted only to the point where the growth rate is maximal); defaults to "to.point"
#' - threshold: A value of the biomass increase that we can surely associate with the end of the lag phase rather than random variation durin the lag. defaults to 10^2 \n
#' - n.points.in.curve: if tangent.method = "local.regression" then n.points.in.curve is the number of points the line is fitted to; 
#' defaults to 3 i.e. the point with the maximal uptake rate one point before and one point aftter \n
#' - init.growth.rate: if logistic model is fitted. Defaults to  NULL in which case the initial value will be based on the data \n
#' - init.lag: if a logistic model is fitted, Defaults to NULL in which case the initial value will be based on the data \n
#' - algorithm: if method = "parameter fitting to a model", nls algorithm to run the model fit; defaults to "auto" which will choose the best between bounded and unbounded "Levenberg-Marquardt" and bounded "port" \n
#' - maxiter =  if method = "parameter fitting to a model", the maximum number of nls iterations, defaults to  100
#' @returns growth curve data (time, biomass, curve_id) with the following additional columns:  log.biomass, lag, line.slope, line.intercept, lag.calculation.method, predicted.data, diff, second.deriv.b, tangent.point, threshold
#' @export
Calculate.Lag = function(data, method, pars) {
  if (!("curve_id" %in% names(data))) {
    data$curve_id = "growth.curve"
  }
  
  N0 = data %>%
    group_by(curve_id) %>%
    arrange(time) %>%
    summarise(N0 = Get.N0(biomass, pars$N0.method)) %>%
    ungroup() %>%
    mutate(log.N0 = log(N0))

  
  if (method == "exponential") {
    selected.tangent.method = pars$tangent.method
    data.new = Fit.Exponential.Lag(data, 
                                   N0 = N0,
                                   tangent.method = selected.tangent.method,
                                   n.points.in.curve = pars$n.points.in.curve)
    data.new = data.new %>%
      select(time, biomass, curve_id, lag, line.slope, line.intercept, tangent.point) %>%
      mutate(lag.calculation.method = "exponential",
             log.biomass = log(biomass),
             predicted.data = NA,
             diff = NA,
             second.deriv.b = NA, 
             threshold = NA) %>%
      select(time, biomass, log.biomass, curve_id, lag, line.slope, line.intercept, lag.calculation.method,predicted.data, diff, second.deriv.b, tangent.point, threshold)
    
    
  } else if (method == "biomass increase") {
    selected.threshold = pars$threshold
    data.new = Lag.Based.On.Biomass.Increase(data,
                                             threshold = selected.threshold,
                                             N0 = N0)
    data.new = data.new %>% 
      select(time, biomass, curve_id, lag) %>%
      mutate(lag.calculation.method = "biomass increase",
             log.biomass = log(biomass),
             predicted.data = NA,
             second.deriv.b = NA,
             line.intercept = NA,
             line.slope = NA,
             tangent.point = NA,
             diff = NA) %>%
      left_join(N0, by = "curve_id") %>%
      mutate(threshold = N0 + selected.threshold) %>%
      select(time, biomass, log.biomass, curve_id, lag, line.slope, line.intercept, lag.calculation.method, predicted.data, diff, second.deriv.b, tangent.point, threshold)
    
  } else if (method == "max growth acceleration") {
    data.new =  Fit.Max.Inflection.Lag(data)
    data.new = data.new %>%
      select(time, biomass, log.biomass, curve_id, lag, second.deriv.b) %>%
      mutate(lag.calculation.method = "max growth acceleration",
             line.intercept = NA,
             line.slope = NA,
             predicted.data = NA,
             diff = NA,
             tangent.point = NA,
             threshold= NA) %>%
      select(time, biomass, log.biomass, curve_id, lag, line.slope, line.intercept, lag.calculation.method,predicted.data, diff, second.deriv.b, tangent.point, threshold)
  } else if (method == "parameter fitting to a model") {
    selected.model = pars$model
    if (selected.model == "logistic") {
      data.new = Calculate.Lagistic.Fit.Lag(data, N0, 
                                            init.growth.rate = pars$init.growth.rate, 
                                            init.K = pars$init.K, 
                                            init.lag = pars$init.lag,
                                            algorithm = pars$algorithm,
                                            maxiter = pars$maxiter) 
    } else if (selected.model == "baranyi") {
      data.new = Calculate.Baranyi.Fit.Lag(data, 
                                           N0, 
                                           init.lag = pars$init.lag,
                                           init.growth.rate = pars$init.growth.rate,
                                           algorithm = pars$algorithm,
                                           maxiter = pars$maxiter) %>%
        mutate(lag.calculation.method = "Fitting lagged baranyi")
        
    } else {
      error("model not implemented")
    }
    data.new = data.new %>% 
      select(time, biomass, curve_id, lag, predicted.data = predicted) %>%
      mutate(lag.calculation.method = "parameter fitting to a model",
             log.biomass = log(biomass),
             diff = NA,
             second.deriv.b = NA,
             line.intercept = NA,
             line.slope = NA,
             tangent.point = NA,
             threshold = NA) %>%
      select(time, biomass, log.biomass, curve_id, lag, line.slope, line.intercept, lag.calculation.method, predicted.data, diff, second.deriv.b, tangent.point, threshold)
    
  }
  data.new = data.new %>% 
    left_join(N0)# %>%
  #data.new$lag[data.new$lag < 0] = 0
  return(data.new) 
}



#' Get.Lags.Calculated.By.All.Methods
#' 
#' Runs the main function that calculates lags based on growth curve data based on all possible methods.
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param biomass.increase.threshold A value of the biomass increase that we can surely associate with the end of the lag phase rather than random variation durin the lag. 
#' Needs to be set specificly to avoid unconcious use of the value set by default. If set to NULL, the value from pars will be taken 
#' @param pars a list of parameters. defaults to the ones set by Get.default.parameters function. Otherwise create your onwn list with the following names: \n
#' - model: if method = "parameter fitting to a model" , one of the following models needs to be chosen: "logistic", "baranyi" \n
#' - N0.method: first.observation" if the first point is taken as the initial biomass or 
#' "minimal.observation" if the minimal biomass is taken is the initial point. 
#' In "healthy" growth curves these options should be equivalent 
#' but sometimes a drop in OD/biomass is observed at the beginning of a growth curve. 
#' In this case it is not obvious what to assume the initial biomass is. \n
#' - tangent.method "local.regression" (if the tangent is fitted to a number of points around the maximal growth rate) 
#' or "to.point" (if the tangent is fitted only to the point where the growth rate is maximal); defaults to "to.point"
#' - threshold: A value of the biomass increase that we can surely associate with the end of the lag phase rather than random variation durin the lag. defaults to 10^2 \n
#' - n.points.in.curve: if tangent.method = "local.regression" then n.points.in.curve is the number of points the line is fitted to; 
#' defaults to 3 i.e. the point with the maximal uptake rate one point before and one point aftter \n
#' - init.growth.rate: if logistic model is fitted. Defaults to  NULL in which case the initial value will be based on the data \n
#' - init.lag: if a logistic model is fitted, Defaults to NULL in which case the initial value will be based on the data \n
#' - algorithm: if method = "parameter fitting to a model", nls algorithm to run the model fit; defaults to "auto" which will choose the best between bounded and unbounded "Levenberg-Marquardt" and bounded "port" \n
#' - maxiter =  if method = "parameter fitting to a model", the maximum number of nls iterations, defaults to  100
#' @returns growth curve data (time, biomass, curve_id) with the column: lag.calculation.method, and with the following additional columns:  log.biomass, lag, line.slope, line.intercept, lag.calculation.method, predicted.data, diff, second.deriv.b, tangent.point, threshold \n
#' Note that each growth curve will appear 
#' @export
Get.Lags.Calculated.By.All.Methods = function(data, biomass.increase.threshold, pars = NULL) {
  if (is.null(pars)) {
    pars = Get.default.parameters()
  }
  if (is.null(biomass.increase.threshold))  {
  biomass.increase.threshold = pars$threshold
  }
  pars.logistic = pars
  pars.logistic$model = "logistic"
  data.new.logistic = Calculate.Lag(data = data, 
                                    method = "parameter fitting to a model", 
                                    pars = pars.logistic) %>%
    mutate(lag.calculation.method = "par. fitting\nto logistic model")
  
  pars.baranyi = pars
  pars.baranyi$model = "baranyi"
  data.new.baranyi = Calculate.Lag(data = data, 
                                   method = "parameter fitting to a model", 
                                   pars = pars.baranyi) %>%
    mutate(lag.calculation.method = "par. fitting\nto baranyi model")
  
  data.new.max.infl = Calculate.Lag(data = data, 
                                    method = "max growth acceleration", 
                                    pars = pars) %>%
    mutate(lag.calculation.method = "max\ngrowth acceleration")
  
  pars.to.point = pars
  pars.to.point$tangent.method = "to.point"
  data.new.expo = Calculate.Lag(data = data, 
                                method = "exponential",
                                pars= pars.to.point
                                ) %>%
    mutate(lag.calculation.method = "tangent to \nmax growth point")
  
  pars.local.regr = pars
  pars.local.regr$tangent.method = "local.regression"
  data.new.expo2 = Calculate.Lag(data = data, 
                                 method = "exponential",
                                 pars = pars.local.regr)%>%
    mutate(lag.calculation.method = "tangent to \nmax growth line")
  
  pars$threshold = biomass.increase.threshold
  data.new.biominc =  Calculate.Lag(data = data, 
                                    method = "biomass increase",
                                    pars=pars)%>%
    mutate(lag.calculation.method = "biomass \nincrease")
  
  
  data.all.with.lag = data.new.max.infl %>%
    rbind(data.new.expo) %>%
    rbind(data.new.expo2) %>%
    rbind(data.new.biominc) %>%
    rbind(data.new.baranyi) %>%
    rbind(data.new.logistic)
  
  data.all.with.lag$lag[data.all.with.lag$lag < 0] = NA
  return(data.all.with.lag)
}




#' Get.default.parameters
#' Set defaults parameters used by Calculate.Lag function
#' @returns list of parameters
#' @export
Get.default.parameters = function() {
  pars = list(model = "logistic",
              N0.method = "first.observation",
              tangent.method = "local.regression",
              threshold = 10^2,
              n.points.in.curve = 3,
              init.growth.rate = NULL,#0.1, 
              init.lag = NULL, #3.5,
              algorithm = "auto",#"default", "Levenberg-Marquardt",#
              maxiter = 100)
  return(pars)
}



#' Smooth.Data
#' Smoothens growth curves data
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param smooth.kind kind used for the smooth functions, defaulyts to "3RS3R"
#' @returns smoothened data
#' @export
Smooth.Data = function(data, smooth.kind = "3RS3R") {
  if (!("curve_id" %in% names(real.data))) {
    real.data$curve_id = "growth.curve"
  }
  data.smooth = data  %>% filter(FALSE)
  for (this_curve_id in unique(data$curve_id)) {
    data_this_curve = data %>% filter(curve_id == this_curve_id) %>% 
      mutate(biomass.smooth = smooth(biomass, kind = smooth.kind)) %>%
      select(time, biomass = biomass.smooth, curve_id)
    data.smooth = rbind(data.smooth, data_this_curve)
  }
  return(data.smooth)
}



#' Cut.The.Data
#' Smoothens growth curves data
#' @param data a data frame with two required columns names: "time" and "biomass",and one optional column: "curve_id"
#' This is data from may come from multiple growth curves
#' @param max.time max. time at which we want to cut the growth curve data
#' @returns cut data
#' @export
Cut.The.Data = function(data, max.time) {
  data.short = data %>% filter(time <= max.time)
  return(data.short)
}



#' Get.Theme (not exported)
#' 
#' This function sets a ggplot theme without grid
#' @param text.size defaults to 12
#' @returns a ggplot theme

Get.Theme = function(text.size = 12) {
  my_theme = theme(
    panel.grid.major = element_blank(),#element_line(colour = "black", size = 0.05),
    panel.grid.minor = element_blank(),#element_line(colour = "black", size = 0.05),
    panel.background = element_rect(fill = "white",
                                    colour = "gray",
                                    size = 0.5, linetype = "solid"),
    text = element_text(size=text.size),
    panel.border = element_rect(linetype = "solid", fill = NA, color = "black")
  )
  return(my_theme)
}
