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


Get.N0 = function(biomass, N0.method) {
  # Get the initial biomass value
  if (N0.method == "first.observation") {
    N0.calc = biomass[1]
  } else if (N0.method == "minimal.observation") {
    N0.calc = min(biomass)
  }
  return(N0.calc)
}

sumLeastSquaresFitLogisticToDeoptim = function(param, N0, data) {
  growth.rate = param[1]
  K = param[2]
  lag = param[3]
  simulatedN = Simulate.Logistic.With.Lag(N0, growth.rate, K, lag, data$time)
  difference = data$biomass - simulatedN
  return(sum(difference^2))
}


Fit.To.Logistic.With.Lag = function(N0, data) {
  deopticontrol = DEoptim.control(itermax = 100, reltol = 10^(-8), trace = 100)
  deoptim_out = DEoptim(fn = function(param) {sumLeastSquaresFitLogisticToDeoptim(param, N0, data)}, 
                        lower = c(0.01,1,0), 
                        upper=c(2,50,10),
                        control = deopticontrol)
  return(deoptim_out)
}


#' Fit.Exponential.Lag
#' @data required columns with names "time", "biomass", 
#' data from one growth curve only, one (mean) observation per time
#' @N0 "minimal.observation" or "first.observation"
#' @tangent.method
Fit.Exponential.Lag.To.Curve = function(data, N0, tangent.method = "to.point", n.points.in.curve = 3) {
  data = data %>%
    mutate(log.biomass = log(biomass), 
           db = log.biomass - dplyr::lag(log.biomass,1),
           dt = time - dplyr::lag(time,1),
           time.diff = (time + dplyr::lag(time))/2,
           db_over_dt = db/dt)
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
    time.max.diff = data$time[ind.max.diff]
    log.biomass.at.max.diff = data$log.biomass[ind.max.diff]
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
  return(list(line.slope = line.slope, line.intercept = line.intercept, lag = lag))
}

################# Lag fitting methods ########################


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
             line.slope = lag.object$line.slope) 
    data.new = rbind(data.new, data_this_curve)
  }
  
  data.new$predicted = NA
  data.new$lag.calculation.method = paste("Exponential ", tangent.method)
  # data.new has columns with lag, intercept, slope, log.N0
  return(data.new)
}
  


# Find where biomass increases by a predefined threshold from N0
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



# Find where biomass increases by a predefined threshold from the previous value
Lag.Based.On.Biomass.Increase.From.Previous.Observation = function(data, threshold) {
  data.new = data %>% filter(FALSE) %>% mutate(lag = numeric(0))
  for (this_curve_id in unique(data$curve_id)) {
    data_this_curve = data %>% filter(curve_id == this_curve_id)
    data_this_curve = data_this_curve %>% 
      mutate(diff = c(NA, diff(biomass))) 
    #find where second derivative is maximal
    diff.above.threshold = which(data_this_curve$diff >= threshold)
    first.diff.above.threshold = diff.above.threshold[1]
    lag.this.curve = data_this_curve$time[first.diff.above.threshold]
    data_this_curve = data_this_curve %>%
      mutate(
        lag = round(lag.this.curve,1))
    data.new = rbind(data.new, data_this_curve)
  }
  return(data.new)
}

# Find where the second derivative is the largest
Fit.Max.Inflection.Lag = function(data) {
  #deriv <- function(t, x) diff(x) / diff(t)
  #middle_pts <- function(x) x[-1] - diff(x) / 2
  #second_deriv <- function(t,x) deriv(middle_pts(t), deriv(t, x))
  
  
  data.new = data %>% filter(FALSE) %>% mutate(log.biomass = numeric(0), diff = numeric(0), lag = numeric(0))
  for (this_curve_id in unique(data$curve_id)) {
    data_this_curve = data %>% 
      filter(curve_id == this_curve_id) %>%
      arrange(time) %>%
      mutate(log.biomass = log(biomass), 
             time.diff = mean(c(diff(time), diff(lag(time))), na.rm = TRUE),
             time.av = (time + dplyr::lag(time))/2) %>%
      mutate(
             #second.deriv.b = c(NA,second_deriv(time, log.biomass), NA),
             # central scheme
             second.deriv.b = (lead(log.biomass) + lag(log.biomass) - 2*log.biomass)/time.diff^2,
             # we only look at second derivative if we know the first derivative is positive! 
             # In empirical data we sometimes see a decrease in biomass but we don;t want to look at the seconf derivative there!
             biomass.increase=lead(log.biomass)  > lag(log.biomass)
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
  data.new$lag.calculation.method = "max growth acceleration"
  return(data.new)
}









###################### Data Simulation ##############################

logistic_model = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Ndot   = r*N*(1-N/K)
    return(list(c(Ndot)))
  })
}

Simulate.Logistic.With.Lag = function(N0, growth.rate, K, lag, times) {
  #logistic model solution derived for example here
  #https://math.libretexts.org/Bookshelves/Calculus/Book%3A_Calculus_(OpenStax)/08%3A_Introduction_to_Differential_Equations/8.4%3A_The_Logistic_Equation
  N = (times >= lag)*(N0*K*exp(growth.rate*(times-lag))/(K - N0 + N0*exp(growth.rate*(times - lag)))) + 
      (times < lag)*N0
  logistic.simulated.data = data.frame(time = times, biomass = N)
  return(logistic.simulated.data)
}

Simulate.Logistic.With.Lag.Old = function(N0, growth.rate, K, lag, times) {
  pars <- c(r = growth.rate, K=K)
  inits = c(N = N0)
  lag.data = data.frame(time = times) %>% mutate(biomass = N0) %>% filter(time <= lag)
  lagged.times = times - lag
  positive.lagged.times = lagged.times[lagged.times >= 0]
  logistic.simulated.data <- as.data.frame(ode(inits, positive.lagged.times, logistic_model, pars)) %>%
    select(time = time, biomass = N) %>%
    mutate(time = time + lag) %>%
    filter(time > lag & time <= max(times))
  logistic.simulated.data = rbind(lag.data, logistic.simulated.data)
  return(logistic.simulated.data)
}


baranyi_and_roberts_model = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    f = 1-(N/K)
    alpha = Q/(1+Q)
    Qdot = v*Q
    Ndot   = mu_max*alpha*f*N
    return(list(c(Ndot, Qdot)))
  })
}


Simulate.Baranyi = function(Q0, growth.rate, K,  times) {
  v = growth.rate
  pars <- c(v=v,mu_max=growth.rate,K=K)
  inits = c(N = N0, Q = Q0)
  baranyi_and_roberts.simulated.data <- as.data.frame(ode(inits, times, baranyi_and_roberts_model, pars)) %>%
    select(time = time, biomass = N)
  # add additional lag
  #baranyi_and_roberts.simulated.data = baranyi_and_roberts.simulated.data %>%
  #  mutate(time = time + lag) %>%
  #  filter(time > lag)
  #baranyi_and_roberts.simulated.data = rbind(lag.data, baranyi_and_roberts.simulated.data)
  
  return(baranyi_and_roberts.simulated.data)
}


Monod_bacteria_growth_model = function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Jg =Vh*G/(Kh+G) 
    Gdot   = -Jg*N 
    Ndot   = a*Jg*N
    return(list(c(Gdot, Ndot)))
  })
}

Simulate.Monod.With.Lag = function(a,Vh,Kh, N0, G0, lag,  times) {
  pars <- c(a = a, Vh = Vh, Kh = Kh)
  inits = c(G=G0, N=N0)
  lag.data = data.frame(time = times) %>% mutate(biomass = N0) %>% filter(time <= lag)

  monod.simulated.data <- as.data.frame(ode(inits, time, Monod_bacteria_growth_model, pars)) %>%
    select(time = time, biomass = N) %>%
    mutate(time = time + lag) %>%
    filter(time > lag)
  monod.simulated.data = rbind(lag.data, monod.simulated.data)
  return(monod.simulated.data)
}


Simulate.Exponential.Data.With.Lag = function(N0, growth.rate, lag, times) {
  exponential.growth.simulated.data = data.frame(
    time = times) %>%
    mutate(lagged.time = if_else(time - lag < 0, 0, time - lag)) %>%
    mutate(biomass = N0*exp(growth.rate*lagged.time))
  return(exponential.growth.simulated.data)
}



# growthcurve has two columns: time and biomass
Calculate.Lag.Fitting.To.Logistic.With.Lag = function(growthcurve, N0, 
                                                      init.growth.rate =init.growth.rate, init.K = init.K, init.lag = init.lag,
                                                      algorithm = "port", max.iter = 100) {
  tryCatch(
   {

      nlsres = nls(biomass ~ N0 + (time >= lag)*N0*(-1+K*exp(growth.rate*(time-lag))/(K - N0 + N0*exp(growth.rate*(time - lag)))), growthcurve, 
                   list(growth.rate = init.growth.rate, K=init.K, lag = init.lag), 
                   algorithm = algorithm, 
                   control = nls.control(maxiter = max.iter))
      #plotfit(nls1, smooth = TRUE)
      lagN = coef(nlsres)[3] %>% unname()
      return(lagN)
    },
    error=function(cond) {
      lagN = NA
      print(cond)
      return(lagN)
    })
}


Calculate.Lag.Fitting.To.Baryani.With.Lag = function(growthcurve, LOG10N0 = NULL, init.lag = NULL, init.mumax = NULL, init.LOG10Nmax = NULL, algorithm = "port", maxiter = 100) {
  tryCatch(
    {nlsres = nls(baranyi, growthcurve, 
                   list(lag=init.lag, mumax=init.mumax, LOG10N0 = LOG10N0, LOG10Nmax = init.LOG10Nmax),
                   algorithm = algorithm, 
                   control = nls.control(maxiter = maxiter)
                  )
      
      
      #plotfit(nls1, smooth = TRUE)
      lagN = coef(nlsres)[1] %>% unname()
      return(lagN)
    },
    error=function(cond) {
      lagN = NA
      return(lagN)
    })
}


Calculate.Lagistic.Fit.Lag = function(data, N0, init.growth.rate, init.K, init.lag, algorithm, max.iter) {
  if (!("curve_id" %in% names(data))) {
    data$curve_id = NA
  }
  if (is.null(init.K)) {
  init.K = 10*mean(N0$N0)
  }
  data.new = data %>% filter(FALSE) %>% mutate( time = numeric(0), biomass = numeric(0),curve_id = character(0), lag = numeric(0), log.info = character(0))
  for (this_curve_id in unique(data$curve_id)) {
    data_this_curve = data %>% filter(curve_id == this_curve_id) %>% select(time, biomass, curve_id)
    this.N0 = N0 %>% filter(curve_id == this_curve_id) %>% pull(N0)
    this.lag = Calculate.Lag.Fitting.To.Logistic.With.Lag(data_this_curve, N0=this.N0, 
                                                          init.growth.rate =init.growth.rate, init.K = init.K, init.lag = init.lag, algorithm =algorithm, max.iter = max.iter)
    data_this_curve = data_this_curve %>% 
      mutate(lag = round(this.lag,1))
    data.new = rbind(data.new, data_this_curve)
  }
  data.new$lag.calculation.method = "Fitting lagged logistic"
  return(data.new)
}
  
Calculate.Baranyi.Fit.Lag = function(data, inits, N0, init.lag, init.growth.rate) {
  data.new = data %>% filter(FALSE) %>% mutate(lag = numeric(0))
  for (this_curve_id in unique(data$curve_id)) {
    data_this_curve = data %>% filter(curve_id == this_curve_id)
    data_this_curve_for_model = data_this_curve %>%
      mutate(LOG10N = log10(biomass), t = time) %>%
      select(LOG10N, t)
    init.LOG10N0 = log10(N0 %>% filter(curve_id == this_curve_id) %>% pull(N0))
    init.LOG10Nmax = max(data_this_curve_for_model$LOG10N)
    init.mumax = init.growth.rate
    lag.this.curve = Calculate.Lag.Fitting.To.Baryani.With.Lag(data_this_curve_for_model, init.LOG10N0, init.lag, init.mumax, init.LOG10Nmax)
    data_this_curve = data_this_curve %>%
      mutate(lag = round(lag.this.curve,1))
    data.new = rbind(data.new, data_this_curve)
  }
  return(data.new)
}


######################## Plot Fitted Lag #####################
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



Plot.Lag.Fit = function(data.new) {
  data.new = data.new %>%
    group_by(curve_id) %>%
    mutate(x.mid = mean(time),
           lag.info = paste0("Lag = ", round(lag, 3), " [h]."),
           log.biomass = log(biomass),
           log10.biomass = log10(biomass),
           y.max.for.curve = max(log10.biomass),
           log10N0 = log10(exp(log(N0))),
           text.y = 1.005*y.max.for.curve
           #y.limit = 1.1*y.max.for.curve
           ) %>%
    ungroup() %>%
    mutate(min.log10N0 = min(log10.biomass),
           max.log10N0 = max(log10.biomass),
           log10.intercept = line.intercept/log(10),
           log10.slope = line.slope/log(10))
  
  max.time = max(data.new$time)
    #mutate(curve_id = paste0(curve_id, ":\n", lag.info))
  
  #coef.diff = max(data.new$y.max, na.rm = TRUE)/max(data.new$diff, na.rm = TRUE)
  #coef.second.deriv.b = max(data.new$max.log10N0, na.rm = TRUE)/max(data.new$second.deriv.b, na.rm = TRUE)
  
  size.N0.line = 1
  size.lag.line = 1
  g = ggplot(data.new)  + 
    geom_line(aes(x= time, y = log10.biomass), col = "blue") +
    geom_text(aes(x=x.mid, y = text.y, label = lag.info), size = 6, col = "red") +
    geom_point(aes(x= time, y = log10.biomass), col = "blue") +
    geom_hline(aes(yintercept = log10N0), size = size.N0.line, col = "black") +
    geom_abline(aes(intercept = log10.intercept, slope = log10.slope), col = "darkgreen") + 
    geom_vline(aes(xintercept = lag), size = size.lag.line, col = "red", linetype = "dashed") +
    xlab("time [h]") +
    xlim(c(0, max.time)) +
    ylab("Log10(biomass)") +
    facet_grid(curve_id~lag.calculation.method, scales = "free_y") +
  theme(axis.text.y.right=element_text(colour="black"),
          axis.text.y=element_text(colour="blue"),
          axis.title.y=element_text(colour="blue"),
          axis.title.y.right=element_text(colour="black"))    
  
  #g =  g + 
  #  geom_line(aes(x = time, y = diff), alpha = 0.5, col = "black")  +
  #  geom_line(aes(x = time, y = second.deriv.b), alpha = 0.5, col = "black")
  
 # if ("second.deriv.b" %in% colnames(data.new)) {
  #  g= g +
  #  geom_line(aes(x = time, y = second.deriv.b)) +
  #  scale_y_continuous(name = "Log(biomass)",
  #                     sec.axis = sec_axis(~.,
  #                                         name = "Second derivative of log(biomass)"))
  #  
  #}
  
  #else if ("diff" %in% colnames(data.new)) {
  #  g =  g + geom_line(aes(x = time, y = diff))  +
  #    scale_y_continuous(name = "Log(biomass)",
  #                       sec.axis = sec_axis(~.,
  #                                          name = "Increase in log(biomass)"))
  # }
  
  return(g)
}




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
      select(time, biomass, curve_id, lag, line.slope, line.intercept) %>%
      mutate(lag.calculation.method = "exponential",
             log.biomass = log(biomass),
             predicted.data = NA,
             diff = NA,
             second.deriv.b = NA) %>%
      select(time, biomass, log.biomass, curve_id, lag, line.slope, line.intercept, lag.calculation.method,predicted.data, diff, second.deriv.b)
    
    
  } else if (method == "biomass increase") {
    selected.threshold = pars$threshold
    data.new = Lag.Based.On.Biomass.Increase(data,
                                             threshold = selected.threshold,
                                             N0 = N0)
    data.new = data.new %>% 
      select(time, biomass, curve_id, lag, diff = increase.from.N0) %>%
      mutate(lag.calculation.method = "biomass increase",
             log.biomass = log(biomass),
             predicted.data = NA,
             second.deriv.b = NA,
             line.intercept = NA,
             line.slope = NA) %>%
      select(time, biomass, log.biomass, curve_id, lag, line.slope, line.intercept, lag.calculation.method,predicted.data, diff, second.deriv.b)
    
  } else if (method == "max growth acceleration") {
    data.new =  Fit.Max.Inflection.Lag(data)
    data.new = data.new %>%
      select(time, biomass, log.biomass, curve_id, lag, second.deriv.b) %>%
      mutate(lag.calculation.method = "max growth acceleration",
             line.intercept = NA,
             line.slope = NA,
             predicted.data = NA,
             diff = NA) %>%
      select(time, biomass, log.biomass, curve_id, lag, line.slope, line.intercept, lag.calculation.method,predicted.data, diff, second.deriv.b)
  } else if (method == "parameter fitting to a model") {
    selected.model = pars$model
    if (selected.model == "logistic") {
      data.new = Calculate.Lagistic.Fit.Lag(data, N0, 
                                            init.growth.rate = pars$init.growth.rate, 
                                            init.K = pars$init.K, 
                                            init.lag = pars$init.lag,
                                            algorithm = pars$algorithm,
                                            max.iter = pars$max.iter) 
    } else if (selected.model == "baranyi") {
      inits = pars$inits
      data.new = Calculate.Baranyi.Fit.Lag(data, 
                                           inits, 
                                           N0, 
                                           pars$init.lag,
                                           pars$init.growth.rate)
    } else {
      error("model not implemented")
    }
    data.new = data.new %>% 
      select(time, biomass, curve_id, lag) %>%
      mutate(lag.calculation.method = "exponential",
             log.biomass = log(biomass),
             predicted.data = NA,
             diff = NA,
             second.deriv.b = NA,
             line.intercept = NA,
             line.slope = NA) %>%
      select(time, biomass, log.biomass, curve_id, lag, line.slope, line.intercept, lag.calculation.method, predicted.data, diff, second.deriv.b)
    
  }
  data.new = data.new %>% left_join(N0)
  return(data.new) 
}




Get.Lags.Calculated.By.All.Methods = function(data, biomass.increase.threshold) {
  pars = Get.default.parameters()
  
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
  
  return(data.all.with.lag)
}



Get.default.parameters = function() {
  pars = list(model = "logistic",
              N0.method = "first.observation",
              tangent.method = "local.regression",
              threshold = 10^2,
              n.points.in.curve = 3,
              init.growth.rate = 0.1, 
              init.lag = 3.5,
              algorithm = "default", 
              max.iter = 100)
  return(pars)
              
}


Smooth.Data = function(real.data, smooth.kind = "3RS3R") {
  if (!("curve_id" %in% names(real.data))) {
    real.data$curve_id = "growth.curve"
  }
  data.smooth = real.data  %>% filter(FALSE)
  for (this_curve_id in unique(real.data$curve_id)) {
    data_this_curve = real.data %>% filter(curve_id == this_curve_id) %>% 
      mutate(biomass.smooth = smooth(biomass, kind = smooth.kind)) %>%
      select(time, biomass = biomass.smooth, curve_id)
    data.smooth = rbind(data.smooth, data_this_curve)
  }
  return(data.smooth)
}

Cut.The.Data = function(real.data, max.time) {
  data.short = real.data %>% filter(time <= max.time)
  return(data.short)
}
