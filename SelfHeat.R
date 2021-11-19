# This is an R-script which implements the sensor self heating correction for Li-7500 open path IRGA during cold conditions. 
# The procedure was originally developed by Burba et al. (2008) and applied to arctic conditions by Oechel et al. (2014).
# 
# The script:
# 1. imports testdata
# 2. performs assessment of the optimal parameters for the Burba correction after Oechel et al. (2014). 
# 3. performs sonsor heating correction using the optimal parameters
# 4. exports the results as .png and .csv files
# 
# NOTE: It is necessary to assign the optimal weighting parameter (WEIGHT_NEW) manually. 
# 
# Felix Nieberding (f.nieberding@tu-bs.de)
# 2019-12-07


rm(list=ls())
Sys.setenv(TZ='UTC', LANG = "en_US.UTF-8")

library(tidyverse)
library(lubridate)
library(data.table)
library(ggpmisc)
library(gridExtra)

setwd(dir = "~/TransTiP/_NamCo_sync/4_Texte/SelfHeatR/")

# custom variables --------------------------------------------------------
# (If needed) Change this number to get different threshold of what we consider day and night (solar radiantion in W/m2)
RAD_THR = 100

# Change if needed
S_THR = 100

# Set this number to the temperature below which you do not expect a change in respiration anymore, e.g. -35 °C in Oechel et al. (2014)
TEMP_THR = -35

# Set this number to the tempearature where the daily Burba correction (here HC) should yield non-negative values. Positive values are unplausible because it would implicate that the sensor is colder than ambient air.
TEMP_THR2 = -10

# This is the temperature threshold under which to apply the Burba correction
TEMP_THR3 = 0

# This is the month of the year for which we assume well frozen conditions (e.g. January to April in Oechel et al. (2014))
FROZEN_START = 1 
FROZEN_END = 4

# Import data -------------------------------------------------------------
# set working directory to the folder where your EddyPro output fiels are, you need the full_output and fluxnet files
df_input <- fread("Burba_testdata.csv") %>%
  mutate(DATE = as.Date(DOY - 1, origin = "2000-01-01"),
         DATETIME = as.POSIXct(paste(as.Date(DOY - 1, origin = "2000-01-01"), structure(as.integer(TIME_DEC*60*60), class="ITime"), dec = " "))) %>%
  rename(F_CO2 = FCO2_final_op,
         RAD = Rnet,
         TEMP_C = Air_T,
         RHOCP = RhoCp,
         MEAN_U = Mean_U)


# Burba correction --------------------------------------------------------
# function: calculations for Burba correction
Burba_correction <- function(df, WEIGHT, RAD_THR, S_THR) {
  df %>%
    mutate(TEMP_K = TEMP_C + 273.16, # Air Temperature in K
           K_AIR = 0.000067 * TEMP_C + 0.024343, # (W/mK)
           TI_BOT_NEB = ifelse(RAD > RAD_THR,  # Bottom sensor housing from NEBraska (Burba 2008) temperature in K
                               0.944 * TEMP_C + 2.57, 
                               0.883 * TEMP_C + 2.17) + 273.16,
           TI_TOP_NEB = ifelse(RAD > RAD_THR,  # TOP sensor housing from NEBraska (Burba 2008) temperature in K
                               1.005 * TEMP_C + 0.24, 
                               1.008 * TEMP_C - 0.4) + 273.16,
           T_BOT_NEW = 0.01 * (WEIGHT * TI_BOT_NEB + (100 - WEIGHT) * TI_TOP_NEB), # Bottom sensor housing from NEW sensor (Burba 2008) temperature in K
           TI_SPAR = ifelse(RAD > RAD_THR,
                            1.01 * TEMP_C + 0.36,
                            1.01 * TEMP_C - 0.17) + 273.16,
           SI_TOP = K_AIR * (TI_TOP_NEB - TEMP_K) * (0.0225 + (0.0028 * sqrt(0.045 / MEAN_U) + 0.00025 / MEAN_U + 0.0045)) / (0.0225 * (0.0028 * sqrt(0.045 / MEAN_U) + 0.00025 / MEAN_U + 0.0045)),
           SI_BOT = K_AIR * (T_BOT_NEW - TEMP_K) / (0.004 * sqrt(0.065 / MEAN_U) + 0.004), # (W/m²) NEBRASKA
           SI_BOT_NEW = 0.01 * (S_THR * SI_BOT + (100 - S_THR) * SI_TOP), # (W/m²) NEW
           SIP_SPAR = K_AIR * (TI_SPAR - TEMP_K) / (0.0025 * log((0.0025 + 0.0058 * sqrt(0.005 / MEAN_U)) / 0.0025)) * 0.15, # (W/m²)
           PD = 44.6 * 28.97 * P / 101.3 * 273.16 / TEMP_K, # (g/m³)
           HC = (SI_BOT_NEW + SI_TOP + SIP_SPAR) / RHOCP * CO2 / TEMP_K * (1 + 1.6077 * H2O / PD), # (mg/m²s) heating correction
           FC = F_CO2 * 0.044, # Fc before heating correction 
           FC_HC = ifelse(TEMP_C < TEMP_THR3, F_CO2 * 0.044 + HC, F_CO2 * 0.044), # (µmol/m²s) Fc corrected, apply only if Temp_C below 0 °C, above no correction needed
           FC_HC_DIV = HC / F_CO2, # ??
           FC_CUM = (cumsum(coalesce(FC, 0)) + FC * 0) * 30 * 60 * (12 / 44) / 1000, #  (µmol/m²s) Cumulative FC before heating correction
           FC_HC_CUM = (cumsum(coalesce(FC_HC, 0)) + FC_HC * 0) * 30 * 60 * (12 / 44) / 1000) #  (µmol/m²s) Cumulative FC after heating correction
}

# Find Optimal Weight for inclined sensor -------------------------------------------
# Build empty containers of desired size
df_WEIGHT = data.table(Tbot_NE = seq(0,100,by=1)) %>%
  mutate(Ttop_NE = 100 - Tbot_NE,
         Slope = rep(0,length(Tbot_NE)),
         N_neg = rep(0,length(Tbot_NE)))

SLOPE = list()
NEG = vector(length = nrow(df_WEIGHT))

# apply drift correction for every T_bot_NE 
for (i in 1:nrow(df_WEIGHT)) {
  WEIGHT = df_WEIGHT$Tbot_NE[i]

  df_slope <- df_input %>%
    filter(TEMP_C < TEMP_THR & month(DATE) >= FROZEN_START & month(DATE) <= FROZEN_END) %>%
    Burba_correction(df = ., WEIGHT = WEIGHT, RAD_THR = RAD_THR, S_THR = S_THR) 
  
  SLOPE[[i]] <- unname(coef(lm(df_slope$FC_HC ~ df_slope$TEMP_C))[2]) * 1000
  
  df_neg <- df_input %>%
    filter(TEMP_C < TEMP_THR2 & month(DATE) >= FROZEN_START & month(DATE) <= FROZEN_END) %>%
    Burba_correction(df = ., WEIGHT = WEIGHT, RAD_THR = RAD_THR, S_THR = S_THR) %>%
    group_by(DATE) %>% 
    summarise_at(.vars = c("HC"), .funs = c("mean"), na.rm=T)
  
  NEG[[i]] <- sum(df_neg$HC < 0, na.rm = T)
} 

df_WEIGHT$Slope <- do.call("rbind", SLOPE)
df_WEIGHT$N_neg <- NEG

# display table for Adjustemt for inclined sensor 
df_WEIGHT

stop("Set optimal weight before continuing")

# find weighting where slope of Fc vs T is closest to zero in cold temperatures, yet number of uptake hours is smallest  below TEMP_THR2, 
# set this optimal value as new weight for the following Burba correction
WEIGHT_NEW <- 63

# Example from Oechel 2014 ------------------------------------------------
# Fill out this table for each weighting and find weighting where slope of Fc vs T is closest to zero in cold temperatures, yet number of uptake hours is smallest  below -10C
# Weighing the Tbot regression				
#                                   Tbot-NE 		
#                                       Ttop-NE
#                                                Slope of flux vs T at Ta < -35 (*1000)						
#                                                            Number of negative daily correction at Ta < -10 (implausible)
#                                    0	50		   0.90        161
# Adjustemt for inclined sensor			50	50			 0.080				48
#                                   55	45			-0.002				15
#                                   60	40			-0.085	 			 6
#                                   61	39			-0.102				 3
#                                   62	38			-0.119				 3
#                     Optimal --->	63	37			-0.135				 1
#                                   64	36			-0.152				 1
#                                   65	35			-0.169				 1
#                                   70	30			-0.250				 0
#                                   80	20			-0.420				 0
#                                   90	10			-0.580				 0
# Vertical sensor in Nebraska			 100	 0			-0.750				 0

# Burba correction --------------------------------------------------------
# set weight according to Optimal value, if needed adjust RAD_THR and S_THR
df_Burba <- Burba_correction(df = df_input, WEIGHT = WEIGHT_NEW, RAD_THR = RAD_THR, S_THR = S_THR)

# 24-hour binned average and standard deviation
df_Burba_daily <-
  df_Burba %>% 
  group_by(DATE) %>% 
  summarise_at(.vars = c("TEMP_C", "HC", "FC", "FC_HC"), .funs = c("mean", "sd"), na.rm=T)

# Figures for Tuning ------------------------------------------------------
# Plot cummulative corrected and uncorrected F_CO2
p1 <- df_Burba %>%
  select(DATETIME, FC_CUM, FC_HC_CUM) %>%
  gather(FC_CUM, FC_HC_CUM, key = series, value = FCO2) %>%
  ggplot(aes(DATETIME, FCO2, color = series)) +
  geom_line(na.rm=T, size = 2) +
  labs(x = "TIME", y = expression('Cummulative CO'[2]*' flux [mg C m'^-2*']')) +
  scale_color_manual(labels = c("uncorrected", "corrected"),
                     values = c("royalblue1", "orangered1")) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
  theme_light() +
  theme(axis.title.x = element_blank())

# Plot daily mean corrected and uncorrected F_CO2
p2 <- df_Burba_daily %>%
  gather(FC_mean,FC_HC_mean, key = series, value = FCO2) %>%
  ggplot(aes(DATE, FCO2, color = series)) +
  geom_line(na.rm=T, size = 1) +
  labs(x = "Time", y = expression('30-min CO'[2]*' flux ['*mu*'mol m'^-2*' s'^-1*']')) +
  scale_color_manual(labels = c("corrected", "uncorrected"),
                     values = c("orangered1", "royalblue1")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  theme_light() +
  theme(axis.title.x = element_blank())

# Mean daily heating correction during cold conditions
p3 <- df_Burba_daily %>%
  filter(TEMP_C_mean < TEMP_THR2 & month(DATE) >= FROZEN_START & month(DATE) <= FROZEN_END) %>%
  ggplot(aes(TEMP_C_mean, HC_mean)) +
  geom_abline(slope = 0, intercept = 0, color = "grey60", size =1.5) +
  geom_point(na.rm=T) +
  geom_smooth(method = "lm", se = FALSE, na.rm = T) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, na.rm = T) +
  labs(x = 'Temperature [°C]', y = expression('Daily mean heating correction [mg m'^-2*']')) +
  scale_y_continuous(limits = c(-0.03,0.05)) +
  theme_light()

# 30-min CO2 Flux against TEMP_C during cold conditions
p4 <- df_Burba %>%   
  select(DATE, TEMP_C, DOY_TIME, FC, FC_HC) %>%
  filter(TEMP_C < TEMP_THR & month(DATE) >= FROZEN_START & month(DATE) <= FROZEN_END) %>% 
  gather(FC, FC_HC, key = series, value = flux) %>%
  ggplot(aes(TEMP_C, flux, color = series)) +
  geom_point(na.rm=T) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, na.rm = T) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, na.rm = T) + 
  labs(x = 'Temperature [°C]', y = expression('30-min CO'[2]*' flux ['*mu*'mol m'^-2*' s'^-1*']')) +
  scale_color_manual(labels = c("uncorrected", "corrected"),
                     values = c("royalblue1", "orangered1")) +
  theme_light()

# daily CO2 Flux against T_SOIL during cold conditions
p5 <- df_Burba_daily %>%   
  select(DATE, TEMP_C_mean, FC_mean, FC_HC_mean) %>%
  filter(TEMP_C_mean < TEMP_THR & month(DATE) >= FROZEN_START & month(DATE) <= FROZEN_END) %>% 
  gather(FC_mean, FC_HC_mean, key = series, value = flux) %>%
  ggplot(aes(TEMP_C_mean,flux, color=series)) +
  geom_point(na.rm=T) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, na.rm = T) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, na.rm = T) +
  labs(x = 'Temperature [°C]', y = expression('Mean daily CO'[2]*' flux ['*mu*'mol m'^-2*' s'^-1*']')) +
  scale_color_manual(labels = c("corrected", "uncorrected"),
                     values = c("orangered1", "royalblue1")) +
  theme_light()

# Daily CO2 FLUX against T_SOIL during cold conditions 
p6 <- df_Burba_daily %>%   
  select(DATE, TEMP_C_mean, FC_mean, FC_HC_mean, FC_sd, FC_HC_sd) %>%
  filter(TEMP_C_mean < TEMP_THR & month(DATE) >= FROZEN_START & month(DATE) <= FROZEN_END) %>% 
  ggplot(aes(x = TEMP_C_mean, y = FC_mean, fill = 'Uncorrected')) +
  geom_errorbar( aes(ymin = FC_mean - FC_sd, ymax = FC_mean + FC_sd), width = 0.4, color = "royalblue1", alpha = 0.6, size = 1, na.rm = TRUE) +
  geom_errorbar( aes(ymin = FC_HC_mean - FC_HC_sd, ymax = FC_HC_mean + FC_HC_sd), width = 0.4, color = "orangered1", alpha = 0.6, size = 1, na.rm = TRUE) +
  geom_point(aes(x = TEMP_C_mean, y = FC_HC_mean, fill = 'corrected'), na.rm = TRUE, size = 2, shape = "triangle filled") +
  geom_point(na.rm = TRUE, size = 2, shape = "triangle down filled") +
  scale_color_manual(values = c("black", "black")) +
  scale_fill_manual(values = c("orangered1", "royalblue1")) +
  labs(x = 'Temperature [°C]', y = expression('Mean daily CO'[2]*' flux ['*mu*'mol m'^-2*' s'^-1*']')) +
  geom_abline(slope = 0, intercept = 0, color = "grey60", size =1, alpha = 0.6) +
  theme_light()

# matrix to set heights
lay <- rbind(c(1,1), c(2,2), c(3,4), c(5,6))
# arrange plots
plots <- arrangeGrob(grobs = list(p1, p2, p3, p4, p5, p6), layout_matrix = lay)

# save plots
ggsave("heatingCorrection.png", plots, width = 24, height = 12)

# export results ----------------------------------------------------------
df_export <- df_Burba %>%
  select(HC, FC, FC_HC, FC_CUM, FC_HC_CUM) 

cbind.data.frame(df, df_export) %>%
  mutate(DATETIME = as.character(DATETIME)) %>%
  fwrite(file = "heatingCorrection.csv")

