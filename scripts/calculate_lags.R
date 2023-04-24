# now for logistic model we use initial parameters calculates wisely:
# K = max(biomass)
# lag calculated with exponential method
# growth rate calculated from the data



library(ggplot2)
library(deSolve)
library(DEoptim)
library(nlsMicrobio)
library(minpack.lm)
library(dplyr)

GITHUB.PATH = "/Users/bognasmug/Dropbox/Projects/Quiesence/2022_Lags/GitHub/microbial_lag_calulator/"
SHINY.APP.PATH = sprintf("%sshiny_app/lag_calulator/", GITHUB.PATH)
EXPERIMENTAL.DATA.PATH = sprintf("%sscripts/data/chosen_exampless_biomass_ml.txt", GITHUB.PATH)
#EXPERIMENTAL.DATA.PATH = "/Users/bognasmug/MGG Dropbox/Bogna Smug/Projects/Quiesence/2022_Lags/data/chosen_exampless_biomass_ml.txt"

  #"/Users/bognasmug/Documents/GitHub/microbial_lag_calulator/shiny_app/lag_calulator/"
source(sprintf("%sR/lags_helper.R", SHINY.APP.PATH))
OUTPUT_FIGS_PATH = "/Users/bognasmug/MGG Dropbox/Bogna Smug/Projects/Quiesence/2022_Lags/Figures/2022_12_18/"
#"/Users/bognasmug/Documents/GitHub/microbial_lag_calulator/Figures/"
dir.create(OUTPUT_FIGS_PATH)
curve_names = paste0("curve_", 1:11)

text_size = 22
biomass.increase.threshold = 10^5
my_theme = Get.Theme(text_size)
blank.constant = 5*10^6
lag.methods = c("biomass \nincrease",           
                "max\ngrowth acceleration",
                "tangent to \nmax growth line",  
                "tangent to \nmax growth point",
                "par. fitting\nto baranyi model",
                "par. fitting\nto logistic model")
  
real.data = read.table(EXPERIMENTAL.DATA.PATH,
                  sep = "\t",
                  header = TRUE) %>%
  #mutate(biomass = OD,
  #       time = time_hours) 
  select(curve_id, time = time_hours, biomass = biomass_ml)#



real.data.with.lag = Get.Lags.Calculated.By.All.Methods(real.data, biomass.increase.threshold) 
jpeg(sprintf("%sFig2.png", OUTPUT_FIGS_PATH), width = 40, height=30*11/4, units = "cm", res = 600)
Plot.Lag.Fit(real.data.with.lag  %>%
               mutate(curve_id = factor(curve_id, levels = curve_names)) %>%
               mutate(lag.calculation.method = factor(lag.calculation.method, levels = lag.methods))) + my_theme
dev.off()


real.data.plus.constant = rbind(
  real.data %>% mutate(prop.N0.dead = 0),
  real.data %>% mutate(prop.N0.dead = 0.1),
  real.data %>% mutate(prop.N0.dead = 0.5),
  real.data %>% mutate(prop.N0.dead = 0.9)
  ) %>%
  group_by(curve_id) %>%
  arrange(time) %>%
  mutate(N0 = biomass[1]) %>%
  ungroup() %>%
  mutate(biomass = biomass - N0*prop.N0.dead) %>% 
  tidyr::unite(col = 'curve_id', curve_id, prop.N0.dead, sep = "_c_") %>%
  as.data.frame() %>%
  select(time, biomass, curve_id)

real.data.plus.constant.with.lag = Get.Lags.Calculated.By.All.Methods(real.data.plus.constant, biomass.increase.threshold)
lags.stats = real.data.plus.constant.with.lag %>%
  distinct(curve_id, lag, lag.calculation.method) %>%
  tidyr::separate(col = 'curve_id', into = c("curve_id", "prop.N0.dead"), sep = "_c_")

jpeg(sprintf("%sSupFig1_not_blank_corrected.png", OUTPUT_FIGS_PATH), width = 40, height=30, units = "cm", res = 600)
ggplot(lags.stats, aes(x=prop.N0.dead, y= lag, group = curve_id, col = curve_id)) +
  geom_point(size = 2) + geom_line() + 
  facet_grid(.~lag.calculation.method) + my_theme + theme(axis.text.x = element_text(angle = 45)) +
  xlab('Proportion of initial biomass that stays in non-replicative state')
#Plot.Lag.Fit(data.all.with.lag) + my_theme
dev.off()   


data.smooth = Smooth.Data(real.data)
real.data.smooth.with.lag = Get.Lags.Calculated.By.All.Methods(data.smooth, biomass.increase.threshold)

#data.sparse = real.data %>% filter(time %in% seq(0,24,2))
#real.data.sparse.with.lag = Get.Lags.Calculated.By.All.Methods(data.sparse, biomass.increase.threshold)

data.short = Cut.The.Data(real.data,12)
real.data.short.with.lag = Get.Lags.Calculated.By.All.Methods(data.short, biomass.increase.threshold)

data.not.blank.corrected = real.data %>% mutate(biomass = biomass + blank.constant)
real.data.not.blank.corrected.with.lag = Get.Lags.Calculated.By.All.Methods(data.not.blank.corrected, biomass.increase.threshold)
     
lag.data = 
rbind(real.data.smooth.with.lag %>%
  distinct(lag, curve_id, lag.calculation.method) %>%
  mutate(data.type = "smoothened"),
  real.data.with.lag %>%
    distinct(lag, curve_id, lag.calculation.method) %>%
    mutate(data.type = "original")) %>%
  #rbind(
  #real.data.sparse.with.lag %>%
  #distinct(lag, curve_id, lag.calculation.method) %>%
  #mutate(data.type = "sparse")) %>%
  rbind(
    real.data.short.with.lag%>%
      distinct(lag, curve_id, lag.calculation.method) %>%
      mutate(data.type = "cut")) %>%
  rbind(
    real.data.not.blank.corrected.with.lag %>%
      distinct(lag, curve_id, lag.calculation.method) %>%
      mutate(data.type = "not blank corrected")) %>%
  mutate(curve_id = factor(curve_id, levels = curve_names))

write.table(lag.data, file = paste0(OUTPUT_FIGS_PATH, "fig3data.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)

jpeg(sprintf("%sFig3new4.png", OUTPUT_FIGS_PATH), width = 40, height=15, units = "cm", res = 600)
ggplot(lag.data %>% rename(`Data type` = data.type)) +
  geom_point(aes(col = curve_id, y = lag, x = `Data type` )) + 
  my_theme +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(.~lag.calculation.method) + xlab("")
dev.off()

jpeg(sprintf("%sFig3.png", OUTPUT_FIGS_PATH), width = 40, height=15, units = "cm", res = 600)
ggplot(lag.data %>% rename(`Data type` = data.type)) +
  geom_boxplot(aes(col =`Data type`, y = lag, x = curve_id), outlier.colour = NA) + 
  geom_jitter(aes(col = `Data type`, y = lag, x = curve_id), alpha = 0.5) + 
  my_theme +
  theme(axis.text.x = element_text(angle = 90))
dev.off()
jpeg(sprintf("%sFig3new.png", OUTPUT_FIGS_PATH), width = 40, height=15, units = "cm", res = 600)
ggplot(lag.data %>% mutate(`Data type` = data.type)) +
  facet_wrap("data.type") +
  geom_point(aes(col =lag.calculation.method, y = lag, x = curve_id)) + 
  my_theme +
  theme(axis.text.x = element_text(angle = 90))+ xlab("")
dev.off()
jpeg(sprintf("%sFig3new2.png", OUTPUT_FIGS_PATH), width = 40, height=15, units = "cm", res = 600)
ggplot(lag.data %>% mutate(`Data type` = data.type),
       aes(y=lag, x = data.type)) +
  facet_grid(.~curve_id) +
  #geom_boxplot() +
  geom_point(aes(col =lag.calculation.method), size = 2) + 
  Get.Theme(12) +
  theme(axis.text.x = element_text(angle = 90)) + xlab("")
dev.off()
jpeg(sprintf("%sFig3new3.png", OUTPUT_FIGS_PATH), width = 40, height=15, units = "cm", res = 600)
ggplot(lag.data %>% mutate(`Data type` = data.type)) +
  facet_grid(lag.calculation.method~curve_id) +
  geom_col(aes(fill =data.type, y = lag, x = data.type)) + 
  Get.Theme(12) +
  theme(axis.text.x = element_text(angle = 90)) + xlab("")
dev.off()
  #geom_jitter(aes(col=data.type, y = lag, x = curve_id)) 

jpeg(sprintf("%sFig3a.png", OUTPUT_FIGS_PATH), width = 50, height=40, units = "cm", res = 600)
ggplot(lag.data) +
  geom_col(aes(fill=data.type, y = lag, x = curve_id), position = "dodge", width = 0.5) +
  facet_grid(lag.calculation.method~.) + my_theme
dev.off()




# simulate perfect data
time.interval = 0.1
times = seq(0,24,time.interval)
lag = 2.5
growth.rate = 0.1
K=5*10^6
growth.rate.logistic = 0.25
# Baranyi lag= ln(1+1/q0)/growth.rate
Q0= 1/(exp(2.5*growth.rate) - 1)
#Q0 = 3.5
#a = 0.016*10^7
#Vh = 300*10^(-7)
# those parameters give as a curve similar to the logistic one
a = 0.03*10^7 # altered from previous one so that the curve has similar K to the one by logistic growth
Vh = 150*10^(-7) # altered from previous one so that the curve has similar K to the one by logistic growth
Kh = 300
N0 = 10^6
byranayi_and_roberts.simulated.data.raw = Simulate.Baranyi(Q0, growth.rate, K, times) 
logistic.simulated.data = Simulate.Logistic.With.Lag(N0, growth.rate.logistic, K, lag, times)
monod.simulated.data = Simulate.Monod.With.Lag(a = a,Vh=Vh,Kh=Kh, N0=N0,G0=13.9, lag=lag, times=times)
exponential.growth.simulated.data = Simulate.Exponential.Data.With.Lag(N0, growth.rate, lag, times)


all.simulated.data =
  rbind(byranayi_and_roberts.simulated.data.raw %>% mutate(curve_id = "Baranyi.and.Roberts"),
        logistic.simulated.data %>% mutate(curve_id = "Logistic"),
        monod.simulated.data %>% mutate(curve_id = "Monod"),
        exponential.growth.simulated.data %>% select(time, biomass) %>% mutate(curve_id = "exponential"))


write.csv(all.simulated.data, file = sprintf("%sall.ssimulated_data.csv", OUTPUT_FIGS_PATH), row.names = FALSE)
jpeg(sprintf("%sFig1.png", OUTPUT_FIGS_PATH), width = 40, height=30, units = "cm", res = 600)
data.all.with.lag = Get.Lags.Calculated.By.All.Methods(all.simulated.data, biomass.increase.threshold)
Plot.Lag.Fit(data.all.with.lag %>%
               mutate(lag.calculation.method = factor(lag.calculation.method, levels = lag.methods))) + my_theme
dev.off()


# Simulate data rom logistic curve
set.seed(1)

# example data with noise
logistic.simulated.data.example.1 = Simulate.Logistic.With.Lag(N0, growth.rate.logistic, K, lag, times) %>%
  mutate(biomass = biomass + rnorm(n = length(times), mean = 0, sd = 0.05*N0)) %>%
  mutate(sd = 0.05, K = K, real.lag = lag)
logistic.simulated.data.example.2 = Simulate.Logistic.With.Lag(N0, growth.rate.logistic, K, lag, times) %>%
  mutate(biomass = biomass + rnorm(n = length(times), mean = 0, sd = 0.2*N0))%>%
  mutate(sd = 0.2, K = K, real.lag = lag)
logistic.simulated.data.example.3 = Simulate.Logistic.With.Lag(N0, growth.rate.logistic, K, lag, times) %>%
  mutate(biomass = biomass + rnorm(n = length(times), mean = 0, sd = 0.5*N0))%>%
  mutate(sd = 0.5, K = K, real.lag = lag)

logistic.simulated.data.example = rbind(logistic.simulated.data.example.1,
                                        rbind(logistic.simulated.data.example.2,
                                              logistic.simulated.data.example.3)) %>%
  rowwise() %>%
  mutate(sd = paste0("sd = ", sd))

jpeg(sprintf("%sFig4_example_noisy_curved.png", OUTPUT_FIGS_PATH), width = 60, height=30, units = "cm", res = 600)
ggplot(logistic.simulated.data.example, aes(x = time, y = biomass)) +
  geom_point() + geom_line() +
  facet_grid(.~sd) + my_theme +
  xlab("time [h]") +
  ylab("biomass [CFU/mL]")
dev.off()


Num.obs = 100
real.lag = 2.5
growth.rates = c(0.1, 0.25, 0.5)
#carrying.capacities = c(0.5*K, K, 10*K)
real.lags = c(0.5, 2.5, 5)
sd_range = seq(0.0, 0.5, 0.1)

################## LOGISTIC MODEL SIMULATIONS ###########################
all.lag.data = data.frame(curve_id = character(0), 
                             lag = numeric(0), 
                             lag.calculation.method = character(0), 
                             sd = numeric(0),
                          growth.rate = numeric(0),
                          carrying.capacity = numeric(0),
                          real.lag = numeric(0),
                          time.interval = numeric(0))

logistic.curves = data.frame(time = numeric(0),
                             biomass = numeric(0),
                             curve_id = character(0), 
                             sd = numeric(0),
                             growth.rate = numeric(0),
                             carrying.capacity = numeric(0),
                             real.lag = numeric(0),
                             time.interval = numeric(0))

logistic.simulated.data.basic = Simulate.Logistic.With.Lag(N0, growth.rate.logistic, K, real.lag, times)
ggplot(data = logistic.simulated.data.basic, aes(x = time, y = biomass )) + geom_point() + geom_line() +
  ylim(c(10^6, 5*10^6))




#all.lag.data.1 = all.lag.data %>% filter(FALSE)
#for (this.carrying.capacity in carrying.capacities) {
#      logistic.simulated.data.basic = Simulate.Logistic.With.Lag(N0, growth.rate.logistic, this.carrying.capacity, real.lag, times)
#      lag.df = Get.Lag.Fitting.Data.For.Noisy.Simulations(logistic.simulated.data.basic,
#                                                 sd_range = sd_range,
#                                                 biomass.increase.threshold,
#                                                 Num.obs = Num.obs) 
#      all.lag.data.1 = rbind(all.lag.data.1, lag.df %>% mutate(growth.rate = growth.rate.logistic,
#                                                             carrying.capacity = this.carrying.capacity,
#                                                             real.lag = real.lag))
#}
  


all.lag.data.2 = all.lag.data %>% filter(FALSE)
logistic.curves.2 = logistic.curves %>% filter(FALSE)
for (this.real.lag in real.lags) {
  logistic.simulated.data.basic = Simulate.Logistic.With.Lag(N0, growth.rate.logistic, K, this.real.lag, times)
  lag.df.obj = Get.Lag.Fitting.Data.For.Noisy.Simulations(logistic.simulated.data.basic,
                                                      sd_range = sd_range,
                                                      biomass.increase.threshold,
                                                      Num.obs = Num.obs) 
 
  logistic.curves.2 = rbind(logistic.curves.2, 
                          lag.df.obj$curves %>% mutate(growth.rate = growth.rate.logistic,
                                                       carrying.capacity=K,
                                    real.lag = this.real.lag,
                                    time.interval = time.interval))
  all.lag.data.2 = rbind(all.lag.data.2, 
                         lag.df.obj$lag.df %>% mutate(growth.rate = growth.rate.logistic,
                                           carrying.capacity = K,
                                           real.lag = this.real.lag,
                                           time.interval = time.interval))
}



all.lag.data.3 = all.lag.data %>% filter(FALSE)
logistic.curves.3 = logistic.curves %>% filter(FALSE)
for (this.growth.rate in growth.rates) {
  logistic.simulated.data.basic = Simulate.Logistic.With.Lag(N0, this.growth.rate, K, real.lag, times)
  lag.df.obj = Get.Lag.Fitting.Data.For.Noisy.Simulations(logistic.simulated.data.basic,
                                                      sd_range = sd_range,
                                                      biomass.increase.threshold,
                                                      Num.obs = Num.obs) 
  logistic.curves.3 = rbind(logistic.curves.3, 
                            lag.df.obj$curves %>% mutate(growth.rate = this.growth.rate,
                                                         carrying.capacity=K,
                                                         real.lag = real.lag,
                                                         time.interval = time.interval))
  all.lag.data.3 = rbind(all.lag.data.3, 
                         lag.df.obj$lag.df %>% mutate(growth.rate = this.growth.rate,
                                                       carrying.capacity = K,
                                                       real.lag = real.lag,
                                                       time.interval = time.interval))
}


all.lag.data.4 = all.lag.data %>% filter(FALSE) %>% mutate(time.interval = numeric(0))
logistic.curves.4 = logistic.curves %>% filter(FALSE)
for (this.time.interval in c(0.1, 0.5, 2)) {
  logistic.simulated.data.basic = Simulate.Logistic.With.Lag(N0, 
                                                             growth.rate.logistic,
                                                             K, 
                                                             real.lag,
                                                             seq(0,24,this.time.interval))
  lag.df.obj = Get.Lag.Fitting.Data.For.Noisy.Simulations(logistic.simulated.data.basic,
                                                      sd_range = sd_range,
                                                      biomass.increase.threshold,
                                                      Num.obs = Num.obs) 
  logistic.curves.4 = rbind(logistic.curves.4, 
                            lag.df.obj$curves %>% mutate(growth.rate = growth.rate.logistic,
                                                         carrying.capacity=K,
                                                         real.lag = real.lag,
                                                         time.interval = this.time.interval))
  all.lag.data.4 = rbind(all.lag.data.4, lag.df.obj$lag.df %>% mutate(growth.rate = growth.rate.logistic,
                                                           carrying.capacity = K,
                                                           real.lag = real.lag,
                                                           time.interval = this.time.interval))
}
    
#all.lag.data.1 = all.lag.data.1 %>% mutate(obs.minus.real.lag = lag - real.lag)
all.lag.data.2 = all.lag.data.2 %>% mutate(obs.minus.real.lag = lag - real.lag)
all.lag.data.3 = all.lag.data.3 %>% mutate(obs.minus.real.lag = lag - real.lag)
all.lag.data.4 = all.lag.data.4 %>% mutate(obs.minus.real.lag = lag - real.lag)

#saveRDS(all.lag.data.1, paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.carrying.capacities.rds"))    
saveRDS(all.lag.data.2, paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.lags.rds"))    
saveRDS(all.lag.data.3, paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.growth.rate.rds"))   
saveRDS(all.lag.data.4, paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.time.interval.rds"))   
saveRDS(logistic.curves.2, paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.lags.curves.rds"))    
saveRDS(logistic.curves.3, paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.growth.rate.curves.rds"))   
saveRDS(logistic.curves.4, paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.time.interval.curves.rds"))   



all.lag.data.2 = readRDS(paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.lags.rds"))    
all.lag.data.3 = readRDS(paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.growth.rate.rds"))   
all.lag.data.4 = readRDS(paste0(OUTPUT_FIGS_PATH,"all.lag.data.logistic.varying.time.interval.rds"))   



jpeg(sprintf("%sFig4_logistic_noise_vs_grwoth_rate.png", OUTPUT_FIGS_PATH), width = 60, height=30, units = "cm", res = 600)
ggplot(all.lag.data.3 %>% 
         mutate(growth.rate = paste0("Growth rate = ", growth.rate),
                sd = as.factor(sd)), 
       aes(x=sd, y=obs.minus.real.lag, colour=lag.calculation.method)) +
  geom_point(size = 0.5) +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_jitter(alpha = 0.5) +
  geom_boxplot() +
  facet_grid(growth.rate~lag.calculation.method) +
  ylab("Lag [h]: observed - expected") +
  my_theme + theme(legend.position = "none") + 
  ylim(c(-5,15))
dev.off()


jpeg(sprintf("%sFig4_logistic_noise_vs_reallag.png", OUTPUT_FIGS_PATH), width = 60, height=30, units = "cm", res = 600)
ggplot(all.lag.data.2 %>% 
         mutate(real.lag = paste0("Exoected lag = ", real.lag),
                sd = as.factor(sd)),  
       aes(x=sd, y=obs.minus.real.lag, colour=lag.calculation.method)) +
  geom_point(size = 0.5) +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_jitter(alpha = 0.5) +
  geom_boxplot() +
  facet_grid(real.lag~lag.calculation.method) +
  ylab("Lag [h]: observed - expected") +
  my_theme + theme(legend.position = "none") + 
  ylim(c(-5,15))
dev.off()

jpeg(sprintf("%sFig4_logistic_noise_vs_time_interval.png", OUTPUT_FIGS_PATH), width = 60, height=30, units = "cm", res = 600)
ggplot(all.lag.data.4 %>% 
         mutate(time.interval = paste0("Time interval = ", time.interval, " [h]"),
                sd = as.factor(sd)),
       aes(x=sd, y=obs.minus.real.lag, colour=lag.calculation.method, alpha = sd)) +
  geom_point(size = 0.5) +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_jitter(alpha = 0.5) +
  geom_boxplot() +
  facet_grid(time.interval~lag.calculation.method) +
  ylab("Lag [h]: observed - expected") +
  my_theme + theme(legend.position = "none") + 
  ylim(c(-5,15))
dev.off()



######################## SIMULATE MODNOS DATA ##############################
# parameters are chosen so that the curve resambles the logistic one
real.lag = 2.5
a = 0.03*10^7 #so that the curve has similar K to the one by logistic growth
Vh = 150*10^(-7) #so that the curve has similar K to the one by logistic growth
Kh = 300 
N0 = 10^6

Vhs = c(0.4*Vh, Vh, 2*Vh)
real.lags = c(0.5, 2.5, 5)
sd_range = seq(0.0, 0.5, 0.1)
time.interval = 0.1
times = seq(0,24,time.interval)


all.lag.data = data.frame(curve_id = character(0), 
                          lag = numeric(0), 
                          lag.calculation.method = character(0), 
                          sd = numeric(0),
                          Vh = numeric(0),
                          a = numeric(0),
                          real.lag = numeric(0),
                          time.interval = numeric(0))

monod.curves = data.frame(time = numeric(0),
                             biomass = numeric(0),
                             curve_id = character(0), 
                             sd = numeric(0),
                             Vh = numeric(0),
                             a = numeric(0),
                             real.lag = numeric(0),
                             time.interval = numeric(0))

monod.simulated.data.basic = Simulate.Monod.With.Lag(a = a,Vh=Vh,Kh=Kh, N0=N0,G0=13.9, lag=real.lag, times=times)
ggplot(data = monod.simulated.data.basic, aes(x = time, y = biomass )) + geom_point() + geom_line() +
ylim(c(10^6, 5*10^6))
# on data that is not noisy we get a good lag estimate by logistic method and all others.
data.monod.with.lag = Get.Lags.Calculated.By.All.Methods(monod.simulated.data.basic, biomass.increase.threshold)
Plot.Lag.Fit(data.monod.with.lag) + my_theme



all.lag.data.1.monod = all.lag.data %>% filter(FALSE)
monod.curves.1 = monod.curves %>% filter(FALSE)
for (this.Vh in Vhs) {
  monod.simulated.data.basic = Simulate.Monod.With.Lag(a = a,Vh=this.Vh,Kh=Kh, N0=N0,G0=13.9, lag=real.lag, times=times)
  
  lag.df.obj = Get.Lag.Fitting.Data.For.Noisy.Simulations(monod.simulated.data.basic,
                                                      sd_range = sd_range,
                                                      biomass.increase.threshold,
                                                      Num.obs = Num.obs) 
  all.lag.data.1.monod = rbind(all.lag.data.1.monod, 
                               lag.df.obj$lag.df %>% mutate(Vh = this.Vh,
                                                           a = a,
                                                           real.lag = real.lag,
                                                           time.interval = time.interval))
  
  monod.curves.1 = rbind(monod.curves.1,
                         lag.df.obj$curves  %>% mutate(Vh = this.Vh,
                                                        a = a,
                                                        real.lag = real.lag,
                                                        time.interval = time.interval))
}


monod.curves.2 = monod.curves %>% filter(FALSE)
all.lag.data.2.monod = monod.curves %>% filter(FALSE)
for (this.real.lag in real.lags) {
  monod.simulated.data.basic = Simulate.Monod.With.Lag(a = a,Vh=Vh,Kh=Kh, N0=N0,G0=13.9, lag=this.real.lag, times=times)
  lag.df.obj = Get.Lag.Fitting.Data.For.Noisy.Simulations(monod.simulated.data.basic,
                                                      sd_range = sd_range,
                                                      biomass.increase.threshold,
                                                      Num.obs = Num.obs) 
  all.lag.data.2.monod = rbind(all.lag.data.2.monod, lag.df.obj$lag.df %>% mutate(Vh = Vh,
                                                           a = a,
                                                           real.lag = this.real.lag,
                                                           time.interval = time.interval))
  monod.curves.2 = rbind(monod.curves.2,
                         lag.df.obj$curves  %>% mutate(Vh = Vh,
                                                       a = a,
                                                       real.lag = this.real.lag,
                                                       time.interval = time.interval))
}



all.lag.data.3.monod = all.lag.data %>% filter(FALSE)
monod.curves.3 = monod.curves %>% filter(FALSE)
for (this.time.interval in c(0.1, 0.5, 2)) {
  monod.simulated.data.basic = Simulate.Monod.With.Lag(a = a,Vh=Vh,Kh=Kh, N0=N0,G0=13.9, lag=real.lag, times=seq(0,24,this.time.interval))
  lag.df.obj = Get.Lag.Fitting.Data.For.Noisy.Simulations(monod.simulated.data.basic,
                                                      sd_range = sd_range,
                                                      biomass.increase.threshold,
                                                      Num.obs = Num.obs) 
  all.lag.data.3.monod = rbind(all.lag.data.3.monod, lag.df.obj$lag.df %>% mutate(Vh = Vh,
                                                           a = a,
                                                           real.lag = real.lag,
                                                           time.interval = this.time.interval))
  
  monod.curves.3 = rbind(monod.curves.3,
                         lag.df.obj$curves  %>% mutate(Vh = Vh,
                                                       a = a,
                                                       real.lag = real.lag,
                                                       time.interval = this.time.interval))
}




all.lag.data.1.monod = all.lag.data.1.monod %>% mutate(obs.minus.real.lag = lag - real.lag)
all.lag.data.2.monod = all.lag.data.2.monod %>% mutate(obs.minus.real.lag = lag - real.lag)
all.lag.data.3.monod = all.lag.data.3.monod %>% mutate(obs.minus.real.lag = lag - real.lag)

saveRDS(all.lag.data.1.monod, paste0(OUTPUT_FIGS_PATH,"all.lag.data.monod.varying.Vh.rds"))    
saveRDS(all.lag.data.2.monod, paste0(OUTPUT_FIGS_PATH,"all.lag.data.monod.varying.lags.rds"))    
saveRDS(all.lag.data.3.monod, paste0(OUTPUT_FIGS_PATH,"all.lag.data.monod.varying.time.interval.rds"))  
saveRDS(monod.curves.1, paste0(OUTPUT_FIGS_PATH,"all.lag.data.monod.varying.Vh.curves.rds"))    
saveRDS(monod.curves.2, paste0(OUTPUT_FIGS_PATH,"all.lag.data.monod.varying.lags.curves.rds"))    
saveRDS(monod.curves.3, paste0(OUTPUT_FIGS_PATH,"all.lag.data.monod.varying.time.interval.curves.rds"))  


jpeg(sprintf("%sFig4_Monod_noise_vs_Vh.png", OUTPUT_FIGS_PATH), width = 60, height=30, units = "cm", res = 600)
ggplot(all.lag.data.1.monod %>% 
         mutate(growth.rate = paste0("Vh = ", Vh), 
                sd = as.factor(sd)),
       aes(x=sd, y=obs.minus.real.lag, colour=lag.calculation.method, alpha = sd)) +
  geom_point(size = 0.5) +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_jitter(alpha = 0.5) +
  geom_boxplot() +
  facet_grid(growth.rate~lag.calculation.method) +
  ylab("Lag [h]: observed - expected") +
  my_theme + theme(legend.position = "none")  + 
  ylim(c(-5,15))
dev.off()


jpeg(sprintf("%sFig4_Monod_noise_vs_reallag.png", OUTPUT_FIGS_PATH), width = 60, height=30, units = "cm", res = 600)
ggplot(all.lag.data.2.monod %>% 
         mutate(real.lag = paste0("Exoected lag = ", real.lag),
                sd = as.factor(sd)),
       aes(x=sd, y=obs.minus.real.lag, colour=lag.calculation.method, alpha = sd)) +
  geom_point(size = 0.5) +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_jitter(alpha = 0.5) +
  geom_boxplot() +
  facet_grid(real.lag~lag.calculation.method) +
  ylab("Lag [h]: observed - expected") +
  my_theme + theme(legend.position = "none")  + 
  ylim(c(-5,15))
dev.off()


jpeg(sprintf("%sFig4_Monod_noise_vs_timeinterval.png", OUTPUT_FIGS_PATH), width = 60, height=30, units = "cm", res = 600)
ggplot(all.lag.data.3.monod %>% 
         mutate(time.interval = paste0("Time interval = ", time.interval, " [h]"),
                sd = as.factor(sd)),
       aes(x=sd, y=obs.minus.real.lag, colour=lag.calculation.method, alpha = sd)) +
  geom_point(size = 0.5) +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_jitter(alpha = 0.5) +
  geom_boxplot() +
  facet_grid(time.interval~lag.calculation.method) +
  ylab("Lag [h]: observed - expected") +
  my_theme + theme(legend.position = "none") + 
  ylim(c(-5,15))
dev.off()


dat = monod.simulated.data.basic %>% mutate(type = "monod") %>% rbind(logistic.simulated.data.basic %>% mutate(type = "logistic"))
ggplot(data = dat, aes(x = time, y = biomass, col = type)) + geom_point() + geom_line() +
      ylim(c(10^6, 5*10^6))
