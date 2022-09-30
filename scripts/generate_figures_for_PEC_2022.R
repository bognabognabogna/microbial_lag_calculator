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
OUTPUT_FIGS_PATH = "/Users/bognasmug/MGG Dropbox/Bogna Smug/Presentations/2022_09_21_Torun"
curve_names = paste0("curve_", 1:11)

text_size = 10
biomass.increase.threshold = 10^5
my_theme = Get.Theme(text_size)
blank.constant = 5*10^6

real.data = read.table("/Users/bognasmug/MGG Dropbox/Bogna Smug/Projects/Quiesence/2022_Lags/data/chosen_exampless_biomass_ml.txt",
                  sep = "\t",
                  header = TRUE) %>%
  #mutate(biomass = OD,
  #       time = time_hours) 
  select(curve_id, time = time_hours, biomass = biomass_ml)



# simulate perfect data
time.interval = 0.5
times = seq(0,24,time.interval)
N0 = 10^6
lag = 5
growth.rate = 0.1
K=5*10^6
growth.rate.logistic = 0.25
# Baranyi lag= ln(1+1/q0)/growth.rate
Q0= 1/(exp(2.5*growth.rate) - 1)
#Q0 = 3.5
a = 0.016*10^7
Vh = 300*10^(-7)
Kh = 300
monod.simulated.data = Simulate.Monod.With.Lag(a = a,Vh=Vh,Kh=Kh, N0=N0,G0=13.9, lag=lag, times=times)




p=ggplot(monod.simulated.data,
  aes(x = time, y = biomass)) +
  geom_line() +
  geom_point() + 
  xlab("Time [h]") +
  ylab("Population size [CFU/mL]") +
  scale_y_log10() +
  my_theme 
ggsave(filename = "Fig1", plot = p, path =OUTPUT_FIGS_PATH, device = "png")



p=ggplot(real.data %>% filter(curve_id == "curve_11"), 
         aes(x = time, y = biomass)) +
  geom_line() +
  geom_point() + 
  xlab("Time [h]") +
  ylab("Population size [CFU/mL]") +
  scale_y_log10() +
  my_theme 
ggsave(filename = "Fig2_11", plot = p, path =OUTPUT_FIGS_PATH, device = "png")

p=ggplot(real.data %>% filter(curve_id == "curve_10"), 
         aes(x = time, y = biomass)) +
  geom_line() +
  geom_point() + 
  xlab("Time [h]") +
  ylab("Population size [CFU/mL]") +
  scale_y_log10() +
  my_theme 
ggsave(filename = "Fig2_10", plot = p, path =OUTPUT_FIGS_PATH, device = "png")


p=ggplot(real.data %>% filter(curve_id == "curve_7"), 
         aes(x = time, y = biomass)) +
  geom_line() +
  geom_point() + 
  xlab("Time [h]") +
  ylab("Population size [CFU/mL]") +
  scale_y_log10() +
  my_theme 
ggsave(filename = "Fig2_7", plot = p, path =OUTPUT_FIGS_PATH, device = "png")



curve11lag = Get.Lags.Calculated.By.All.Methods(real.data %>% filter(curve_id == "curve_11"), 
                                               biomass.increase.threshold)  

curve11lag1 = curve11lag %>%
  filter(lag.calculation.method == "tangent to \nmax growth point") 
Plot.Lag.Fit(curve11lag1, print.lag.info = FALSE) + my_theme + ylim(c(6.5,7.5))
ggsave(filename = "Fig3_tangent", plot = last_plot(), path =OUTPUT_FIGS_PATH, device = "png")

curve11lag2 = curve11lag %>%
  filter(lag.calculation.method == "max\ngrowth acceleration" ) %>%
  mutate(N0 = NA)
Plot.Lag.Fit(curve11lag2) + my_theme + ylim(c(6.5,7.5))
ggsave(filename = "Fig3_acceleration", plot = last_plot(), path =OUTPUT_FIGS_PATH, device = "png")


curve11lag3 = curve11lag %>%
  filter(lag.calculation.method == "biomass \nincrease" ) %>%
  mutate(N0 = NA)
Plot.Lag.Fit(curve11lag3) + my_theme + ylim(c(6.5,7.5))
ggsave(filename = "Fig3_biomassincrease", plot = last_plot(), path =OUTPUT_FIGS_PATH, device = "png")



monod.simulated.data.lags = Get.Lags.Calculated.By.All.Methods(monod.simulated.data,
                                                biomass.increase.threshold)  

monod.simulated.data.lag1 = monod.simulated.data.lags %>%
  filter(lag.calculation.method == "tangent to \nmax growth point")
Plot.Lag.Fit(monod.simulated.data.lag1, print.lag.info = FALSE) + my_theme + ylim(c(5.7,6.6))
ggsave(filename = "Fig3_tangent", plot = last_plot(), path =OUTPUT_FIGS_PATH, device = "png")

monod.simulated.data.lag2 = monod.simulated.data.lags %>%
  filter(lag.calculation.method == "max\ngrowth acceleration" ) %>%
  mutate(N0 = NA)
Plot.Lag.Fit(monod.simulated.data.lag2, print.lag.info = FALSE) + my_theme + ylim(c(5.7,6.6))
ggsave(filename = "Fig3_acceleration", plot = last_plot(), path =OUTPUT_FIGS_PATH, device = "png")

monod.simulated.data.lag1 = monod.simulated.data.lags %>%
  filter(lag.calculation.method == "biomass \nincrease")
Plot.Lag.Fit(monod.simulated.data.lag1, print.lag.info = FALSE) + my_theme + ylim(c(5.7,6.6))
ggsave(filename = "Fig3_biomass", plot = last_plot(), path =OUTPUT_FIGS_PATH, device = "png")

monod.simulated.data.lag1 = monod.simulated.data.lags %>%
  filter(lag.calculation.method == "par. fitting\nto logistic model")
Plot.Lag.Fit(monod.simulated.data.lag1, print.lag.info = FALSE) + my_theme + ylim(c(5.7,6.6))
ggsave(filename = "Fig3_logistic_fit", plot = last_plot(), path =OUTPUT_FIGS_PATH, device = "png")


 N0 = monod.simulated.data$biomass[1] 
 noise = rnorm(59, mean = 0, sd = 0.05*N0)
 #noise = c(rnorm(10, mean = 0, sd = 0.15*N0), rnorm(49, mean = 0, sd = 0.05*N0))
 monod.simulated.data.with.noise = monod.simulated.data %>%
  mutate(biomass = biomass + noise)
 
 Plot.Lag.Fit(monod.simulated.data.with.noise.lags %>%
                filter(lag.calculation.method == "par. fitting\nto logistic model") %>%
                mutate(N0 = NA), print.lag.info = FALSE) + my_theme + ylim(c(5.7,6.6))
 ggsave(filename = "Fig3_noise_model_fitting", plot = last_plot(), path =OUTPUT_FIGS_PATH, device = "png")
 
 
 
 
 monod.simulated.data.with.noise.lags = Get.Lags.Calculated.By.All.Methods(monod.simulated.data.with.noise,
                                                                5*biomass.increase.threshold)  

 monod.simulated.data.with.noise.lag2 = monod.simulated.data.with.noise.lags %>%
   filter(lag.calculation.method == "max\ngrowth acceleration" ) %>%
   mutate(N0 = NA)
 Plot.Lag.Fit(monod.simulated.data.with.noise.lag2) + my_theme + ylim(c(5.7,6.6))
 ggsave(filename = "Fig3_noise_acceleration", plot = last_plot(), path =OUTPUT_FIGS_PATH, device = "png")
 
 Plot.Lag.Fit(monod.simulated.data.with.noise.lags %>%
                filter(lag.calculation.method == "par. fitting\nto logistic model") %>%
                mutate(N0 = NA), print.lag.info = FALSE) + my_theme + ylim(c(5.7,6.6))
 ggsave(filename = "Fig3_noise_model_fitting", plot = last_plot(), path =OUTPUT_FIGS_PATH, device = "png")
 
 
 Plot.Lag.Fit(monod.simulated.data.with.noise.lags %>%
                filter(lag.calculation.method == "biomass \nincrease"), print.lag.info = FALSE) + my_theme + ylim(c(5.7,6.6))
 ggsave(filename = "Fig3_noise_biomass_increase", plot = last_plot(), path =OUTPUT_FIGS_PATH, device = "png",
        height = 10, width = 10, units = "cm")
 
 
 real.data.with.lag = Get.Lags.Calculated.By.All.Methods(real.data, biomass.increase.threshold) 
 
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
   mutate(curve_id = factor(curve_id, levels = curve_names)) %>%
   mutate(original.data.type = data.type == "original")
 
 jpeg(sprintf("%s/Fig_dataPpeproc.png", OUTPUT_FIGS_PATH), width = 40, height=15, units = "cm", res = 600)
 ggplot(lag.data %>% rename(`Data type` = data.type)) +
   geom_boxplot(aes(y=lag, x = `Data type`, fill = original.data.type)) +
   geom_point(aes(col = curve_id, y = lag, x = `Data type` )) + 
   Get.Theme(14) +
   theme(axis.text.x = element_text(angle = 90)) +
   facet_grid(.~lag.calculation.method) + xlab("") +
   scale_fill_manual(values = c("TRUE" = "darkblue", "FALSE" = "white")) +
   guides(fill = "none")
 dev.off()
 
 original.data = lag.data %>% filter(data.type == "original")
 ggplot(original.data, aes(x = curve_id, y = lag, fill = lag.calculation.method)) +
geom_col(position = "dodge", width =0.7) + 
   theme_bw() +theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + xlab("") + ylab("Lag [h]")
 ggsave(filename = "Fig4_methods_differ", plot = last_plot(), path =OUTPUT_FIGS_PATH, device = "png",
        heigh = 10, width = 20, units = "cm")
 