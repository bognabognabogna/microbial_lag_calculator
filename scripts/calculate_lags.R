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
SHINY.APP.PATH = "/Users/bognasmug/Documents/GitHub/microbial_lag_calulator/shiny_app/lag_calulator/"
source(sprintf("%sR/lags_helper.R", SHINY.APP.PATH))
OUTPUT_FIGS_PATH = "~/MGG Dropbox/Bogna Smug/Projects/Quiesence/2022_Lags/Figures/new/"
curve_names = paste0("curve_", 1:11)

text_size = 22
biomass.increase.threshold = 10^5
my_theme = Get.Theme(text_size)
blank.constant = 5*10^6

real.data = read.table("/Users/bognasmug/MGG Dropbox/Bogna Smug/Projects/Quiesence/2022_Lags/data/chosen_exampless_biomass_ml.txt",
                  sep = "\t",
                  header = TRUE) %>%
  #mutate(biomass = OD,
  #       time = time_hours) 
  select(curve_id, time = time_hours, biomass = biomass_ml)#



real.data.with.lag = Get.Lags.Calculated.By.All.Methods(real.data, biomass.increase.threshold) 
jpeg(sprintf("%sFig2.png", OUTPUT_FIGS_PATH), width = 40, height=30*11/4, units = "cm", res = 600)
Plot.Lag.Fit(real.data.with.lag  %>%
               mutate(curve_id = factor(curve_id, levels = curve_names))) + my_theme
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

data.sparse = real.data %>% filter(time %in% seq(0,24,2))
real.data.sparse.with.lag = Get.Lags.Calculated.By.All.Methods(data.sparse, biomass.increase.threshold)

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
  rbind(
  real.data.sparse.with.lag %>%
  distinct(lag, curve_id, lag.calculation.method) %>%
  mutate(data.type = "sparse")) %>%
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
N0 = 10^6
lag = 2.5
growth.rate = 0.1
K=5*10^6
growth.rate.logistic = 0.25
# Baranyi lag= ln(1+1/q0)/growth.rate
Q0= 1/(exp(2.5*growth.rate) - 1)
#Q0 = 3.5
a = 0.016*10^7
Vh = 300*10^(-7)
Kh = 300
byranayi_and_roberts.simulated.data.raw = Simulate.Baranyi(Q0, growth.rate, K, times) 
logistic.simulated.data = Simulate.Logistic.With.Lag(N0, growth.rate.logistic, K, lag, times)
monod.simulated.data = Simulate.Monod.With.Lag(a = a,Vh=Vh,Kh=Kh, N0=N0,G0=13.9, lag=lag, times=times)
exponential.growth.simulated.data = Simulate.Exponential.Data.With.Lag(N0, growth.rate, lag, times)


all.simulated.data =
  rbind(byranayi_and_roberts.simulated.data.raw %>% mutate(curve_id = "Baranyi.and.Roberts"),
        logistic.simulated.data %>% mutate(curve_id = "Logistic"),
        monod.simulated.data %>% mutate(curve_id = "Monod"),
        exponential.growth.simulated.data %>% select(time, biomass) %>% mutate(curve_id = "exponential"))


jpeg(sprintf("%sFig1.png", OUTPUT_FIGS_PATH), width = 40, height=30, units = "cm", res = 600)
data.all.with.lag = Get.Lags.Calculated.By.All.Methods(all.simulated.data, biomass.increase.threshold)
Plot.Lag.Fit(data.all.with.lag) + my_theme
dev.off()


# Simulate data rom logistic curve
set.seed(1)

# example data with noise
logistic.simulated.data.example.1 = Simulate.Logistic.With.Lag(N0, growth.rate, K, lag, times) %>%
  mutate(biomass = biomass + rnorm(n = length(times), mean = 0, sd = 0.05*N0)) %>%
  mutate(sd = 0.05, K = K, real.lag = lag)
logistic.simulated.data.example.2 = Simulate.Logistic.With.Lag(N0, growth.rate, K, lag, times) %>%
  mutate(biomass = biomass + rnorm(n = length(times), mean = 0, sd = 0.2*N0))%>%
  mutate(sd = 0.2, K = K, real.lag = lag)
logistic.simulated.data.example.3 = Simulate.Logistic.With.Lag(N0, growth.rate, K, lag, times) %>%
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
growth.rate.logistic = 0.25
growth.rates = c(0.1, 0.25, 0.5)
carrying.capacities = c(0.5*K, K, 10*K)
real.lags = c(0.5, 2.5, 5)
sd_range = seq(0.0, 0.5, 0.1)

all.lag.data = data.frame(curve_id = character(0), 
                             lag = numeric(0), 
                             lag.calculation.method = character(0), 
                             sd = numeric(0),
                          growth.rate = numeric(0),
                          carrying.capacity = numeric(0),
                          real.lag = numeric(0))





all.lag.data.1 = all.lag.data %>% filter(FALSE)
for (this.carrying.capacity in carrying.capacities) {
      logistic.simulated.data.basic = Simulate.Logistic.With.Lag(N0, growth.rate.logistic, this.carrying.capacity, real.lag, times)
      lag.df = Get.Lag.Fitting.Data.For.Noisy.Simulations(logistic.simulated.data.basic,
                                                 sd_range = sd_range,
                                                 biomass.increase.threshold,
                                                 Num.obs = Num.obs) 
      all.lag.data.1 = rbind(all.lag.data.1, lag.df %>% mutate(growth.rate = growth.rate.logistic,
                                                             carrying.capacity = this.carrying.capacity,
                                                             real.lag = real.lag))
}
  


all.lag.data.2 = all.lag.data %>% filter(FALSE)
for (this.real.lag in real.lags) {
  logistic.simulated.data.basic = Simulate.Logistic.With.Lag(N0, growth.rate.logistic, K, this.real.lag, times)
  lag.df = Get.Lag.Fitting.Data.For.Noisy.Simulations(logistic.simulated.data.basic,
                                                      sd_range = sd_range,
                                                      biomass.increase.threshold,
                                                      Num.obs = Num.obs) 
  all.lag.data.2 = rbind(all.lag.data.2, lag.df %>% mutate(growth.rate = growth.rate.logistic,
                                                       carrying.capacity = K,
                                                       real.lag = this.real.lag))
}



all.lag.data.3 = all.lag.data %>% filter(FALSE)
for (this.growth.rate in growth.rates) {
  logistic.simulated.data.basic = Simulate.Logistic.With.Lag(N0, this.growth.rate, K, real.lag, times)
  lag.df = Get.Lag.Fitting.Data.For.Noisy.Simulations(logistic.simulated.data.basic,
                                                      sd_range = sd_range,
                                                      biomass.increase.threshold,
                                                      Num.obs = Num.obs) 
  all.lag.data.3 = rbind(all.lag.data.3, lag.df %>% mutate(growth.rate = this.growth.rate,
                                                       carrying.capacity = K,
                                                       real.lag = real.lag))
}


all.lag.data.4 = all.lag.data %>% filter(FALSE) %>% mutate(time.interval = numeric(0))
for (this.time.internal in c(0.1, 0.5, 2)) {
  logistic.simulated.data.basic = Simulate.Logistic.With.Lag(N0, growth.rate.logistic, K, real.lag,seq(0,24,this.time.internal))
  lag.df = Get.Lag.Fitting.Data.For.Noisy.Simulations(logistic.simulated.data.basic,
                                                      sd_range = sd_range,
                                                      biomass.increase.threshold,
                                                      Num.obs = Num.obs) 
  all.lag.data.4 = rbind(all.lag.data.4, lag.df %>% mutate(growth.rate = growth.rate.logistic,
                                                           carrying.capacity = K,
                                                           real.lag = real.lag,
                                                           time.internal = this.time.internal))
}
    
all.lag.data.1 = all.lag.data.1 %>% mutate(obs.minus.real.lag = lag - real.lag)
all.lag.data.2 = all.lag.data.2 %>% mutate(obs.minus.real.lag = lag - real.lag)
all.lag.data.3 = all.lag.data.3 %>% mutate(obs.minus.real.lag = lag - real.lag)
all.lag.data.4 = all.lag.data.4 %>% mutate(obs.minus.real.lag = lag - real.lag)

saveRDS(all.lag.data.1, paste0(OUTPUT_FIGS_PATH,"all.lag.data.varying.carrying.capacities.rds"))    
saveRDS(all.lag.data.2, paste0(OUTPUT_FIGS_PATH,"all.lag.data.varying.lags.rds"))    
saveRDS(all.lag.data.3, paste0(OUTPUT_FIGS_PATH,"all.lag.data.varying.growth.rate.rds"))   
saveRDS(all.lag.data.4, paste0(OUTPUT_FIGS_PATH,"all.lag.data.varying.time.interval.rds"))   


#ggplot(all.lag.data) + geom_boxplot(aes(col=lag.calculation.method, x=lag.calculation.method, y = obs.minus.real.lag)) + 
#  facet_grid(.~sd) +
#  Get.Theme(10) 


all.lag.data.1 = readRDS(paste0(OUTPUT_FIGS_PATH,"all.lag.data.varying.carrying.capacities.rds"))    
all.lag.data.2 = readRDS(paste0(OUTPUT_FIGS_PATH,"all.lag.data.varying.lags.rds"))    
all.lag.data.3 = readRDS(paste0(OUTPUT_FIGS_PATH,"all.lag.data.varying.growth.rate.rds"))   
all.lag.data.4 = readRDS(paste0(OUTPUT_FIGS_PATH,"all.lag.data.varying.time.interval.rds"))   


jpeg(sprintf("%sFig4_noise_vs_grwoth_rate.png", OUTPUT_FIGS_PATH), width = 60, height=30, units = "cm", res = 600)
ggplot(all.lag.data.3 %>% 
         filter(carrying.capacity == K & real.lag == 2.5) %>%
         rowwise() %>%
         mutate(growth.rate = paste0("Growth rate = ", growth.rate)), 
       aes(x=sd, y=obs.minus.real.lag, colour=lag.calculation.method)) +
       geom_point(size = 0.5) +
       geom_hline(aes(yintercept = 0), col = "black") +
       geom_quantile(quantiles = c(0.25, 0.5, 0.75), size = 2, aes(alpha = ..quantile..)) +
       facet_grid(growth.rate~lag.calculation.method) +
  ylab("Lag [h]: observed - expected") +
       #scale_alpha_manual(values = c(0.75, 0.5, 0.25), breaks = c(0.3,0.9, 0.7)) +
  my_theme + theme(legend.position = "none") #+
  #ylim(c(-10,10))
dev.off()

jpeg(sprintf("%sFig4_noise_vs_K.png", OUTPUT_FIGS_PATH), width = 60, height=30, units = "cm", res = 600)
ggplot(all.lag.data.1 %>% 
         filter(growth.rate == growth.rate.logistic & real.lag == 2.5)  %>%
         mutate(carrying.capacity = plyr::mapvalues(carrying.capacity, 
                                                    from = carrying.capacities, to = c("low", "medium", "high"))) %>%
         rowwise() %>%
         mutate(carrying.capacity = paste0("Carrying capacity = ", carrying.capacity)),  
       aes(x=sd, y=obs.minus.real.lag, colour=lag.calculation.method)) +
  geom_point(size = 0.5) +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_quantile(quantiles = c(0.25, 0.5, 0.75), size = 2, aes(alpha = ..quantile..)) +
  facet_grid(carrying.capacity~lag.calculation.method) +
  #scale_alpha_manual(values = c(0.75, 0.5, 0.25), breaks = c(0.3,0.9, 0.7)) +
  ylab("Lag [h]: observed - expected") +
#scale_alpha_manual(values = c(0.75, 0.5, 0.25), breaks = c(0.3,0.9, 0.7)) +
my_theme + theme(legend.position = "none") #+
 # ylim(c(-10,10))
dev.off()


jpeg(sprintf("%sFig4_noise_vs_reallag.png", OUTPUT_FIGS_PATH), width = 60, height=30, units = "cm", res = 600)
ggplot(all.lag.data.2 %>% 
         filter(growth.rate == growth.rate.logistic & carrying.capacity == K) %>%
         rowwise() %>%
         mutate(real.lag = paste0("Exoected lag = ", real.lag)),  
       aes(x=sd, y=obs.minus.real.lag, colour=lag.calculation.method)) +
  geom_point(size = 0.5) +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_quantile(quantiles = c(0.25, 0.5, 0.75), size = 2, aes(alpha = ..quantile..)) +
  facet_grid(real.lag~lag.calculation.method) +
  #scale_alpha_manual(values = c(0.75, 0.5, 0.25), breaks = c(0.3,0.9, 0.7)) +
  ylab("Lag [h]: observed - expected") +
  #scale_alpha_manual(values = c(0.75, 0.5, 0.25), breaks = c(0.3,0.9, 0.7)) +
  my_theme + theme(legend.position = "none") #+
  #ylim(c(-10,10))
dev.off()

jpeg(sprintf("%sFig4_noise_vs_time_interval.png", OUTPUT_FIGS_PATH), width = 60, height=30, units = "cm", res = 600)
ggplot(all.lag.data.4 %>%
         filter(growth.rate == growth.rate.logistic & carrying.capacity == K & real.lag == 2.5) %>%
         rowwise() %>%
         mutate(time.internal = paste0("Time interval = ", time.internal, " [h]")),  
       aes(x=sd, y=obs.minus.real.lag, colour=lag.calculation.method)) +
  geom_point(size = 0.5) +
  geom_hline(aes(yintercept = 0), col = "black") +
  geom_quantile(quantiles = c(0.25, 0.5, 0.75), size = 2, aes(alpha = ..quantile..)) +
  facet_grid(time.internal~lag.calculation.method) +
  #scale_alpha_manual(values = c(0.75, 0.5, 0.25), breaks = c(0.3,0.9, 0.7)) +
  ylab("Lag [h]: observed - expected") +
  #scale_alpha_manual(values = c(0.75, 0.5, 0.25), breaks = c(0.3,0.9, 0.7)) +
  my_theme + theme(legend.position = "none") #+
  #ylim(c(-10,10))
dev.off()

# data sparsity



#ggplot(all.lag.data) + 
#  geom_boxplot(aes(x=lag.calculation.method, col = lag.calculation.method, y = obs.minus.real.lag)) + 
#  my_theme +
#  facet_grid(.~sd) +        geom_hline(aes(yintercept = 0), col = "black")  + 
#  ylim(c(-10,10))

#params = Fit.To.Logistic.With.Lag(N0, logistic.simulated.data)
#mu.fitted = params$optim$bestmem[1] %>% as.numeric()
#k.fitted = params$optim$bestmem[2] %>% as.numeric()
#lag.fitted = params$optim$bestmem[3] %>% as.numeric()

#logistic model solution derived for example here
#https://math.libretexts.org/Bookshelves/Calculus/Book%3A_Calculus_(OpenStax)/08%3A_Introduction_to_Differential_Equations/8.4%3A_The_Logistic_Equation


