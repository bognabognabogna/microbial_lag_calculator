library(dplyr)
library(ggplot2)
library(deSolve)
library(DEoptim)
library(nlsMicrobio)
source("/Users/bognasmug/MGG Dropbox/Bogna Smug/Projects/Quiesence/2022_mdpi_Lags/shiny_app/lag_calulator/R/lags_helper.R")
OUTPUT_FIGS_PATH = "~/MGG Dropbox/Bogna Smug/Projects/Quiesence/2022_mdpi_Lags/Figures/"
curve_names = paste0("curve_", 1:11)

text_size = 22
biomass.increase.threshold = 10^2
my_theme = Get.Theme(text_size)


real.data = read.table("/Users/bognasmug/MGG Dropbox/Bogna Smug/Projects/Quiesence/2022_mdpi_Lags/data/chosen_exampless_biomass_ml.txt",
                  sep = "\t",
                  header = TRUE) %>%
  #mutate(biomass = OD,
  #       time = time_hours) 
  select(curve_id, time = time_hours, biomass = biomass_ml)#


jpeg(sprintf("%sFig2a.png", OUTPUT_FIGS_PATH), width = 40, height=30*11/4, units = "cm", res = 600)
real.data.with.lag = Get.Lags.Calculated.By.All.Methods(real.data, biomass.increase.threshold)
Plot.Lag.Fit(real.data.with.lag  %>%
               mutate(curve_id = factor(curve_id, levels = curve_names))) + my_theme
dev.off()


data.smooth = Smooth.Data(real.data)
jpeg(sprintf("%sFig2b_smooth.png", OUTPUT_FIGS_PATH), width = 30, height=30*11/4, units = "cm", res = 600)
real.data.smooth.with.lag = Get.Lags.Calculated.By.All.Methods(data.smooth, biomass.increase.threshold)
Plot.Lag.Fit(real.data.smooth.with.lag %>%
               mutate(curve_id = factor(curve_id, levels = curve_names))) + my_theme
dev.off()


data.sparse = real.data %>% filter(time %in% seq(0,24,2))
real.data.sparse.with.lag = Get.Lags.Calculated.By.All.Methods(data.sparse, biomass.increase.threshold)

data.short = Cut.The.Data(real.data,12)
real.data.short.with.lag = Get.Lags.Calculated.By.All.Methods(data.short, biomass.increase.threshold)

     
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
  mutate(curve_id = factor(curve_id, levels = curve_names))

jpeg(sprintf("%sFig3.png", OUTPUT_FIGS_PATH), width = 40, height=15, units = "cm", res = 600)
ggplot(lag.data %>% rename(`Data type` = data.type)) +
  geom_boxplot(aes(col =`Data type`, y = lag, x = curve_id), outlier.colour = NA) + 
  geom_jitter(aes(col = `Data type`, y = lag, x = curve_id), alpha = 0.5) + 
  my_theme +
  theme(axis.text.x = element_text(angle = 90))
dev.off()
  #geom_jitter(aes(col=data.type, y = lag, x = curve_id)) 

jpeg(sprintf("%sFig3a.png", OUTPUT_FIGS_PATH), width = 50, height=40, units = "cm", res = 600)
ggplot(lag.data) +
  geom_col(aes(fill=data.type, y = lag, x = curve_id), position = "dodge", width = 0.5) +
  facet_grid(lag.calculation.method~.) + my_theme
dev.off()
#dd = data1 %>% mutate(biomass = OD)
#lo <- loess(dd$biomass~dd$time)
#dd$smoothOD = predict(lo)
#Fit.Max.Inflection.Lag(dd %>% mutate(biomass = smoothOD))
#sm=ksmooth(data$time, data$OD)
#Fit.Max.Inflection.Lag(data.frame(biomass = sm$y,time = sm$x))

# try also rollmean as a smoother


#TODO add methods
# Byranyi and Roberts 1994



# simulate perfect data
time.interval = 0.5
N0 = 10^6
time = seq(0,24,time.interval)
lag = 2.5
growth.rate = 0.1
K=1*10^7
times = seq(0,24, 0.5)
# Baranyi lag= ln(1+1/q0)/growth.rate
#Q0= 1/(exp(2.5*growth.rate) - 1)
Q0 = 3.5
a = 0.016*10^7
Vh = 300*10^(-7)
Kh = 300
byranayi_and_roberts.simulated.data.raw = Simulate.Baranyi(Q0, growth.rate, K, times) 
logistic.simulated.data = Simulate.Logistic.With.Lag(N0, growth.rate, K, lag, times)
monod.simulated.data = Simulate.Monod.With.Lag(a = a,Vh=Vh,Kh=Kh, N0,G0=13.9, lag, times)
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




#params = Fit.To.Logistic.With.Lag(N0, logistic.simulated.data)
#mu.fitted = params$optim$bestmem[1] %>% as.numeric()
#k.fitted = params$optim$bestmem[2] %>% as.numeric()
#lag.fitted = params$optim$bestmem[3] %>% as.numeric()

#logistic model solution derived for example here
#https://math.libretexts.org/Bookshelves/Calculus/Book%3A_Calculus_(OpenStax)/08%3A_Introduction_to_Differential_Equations/8.4%3A_The_Logistic_Equation


