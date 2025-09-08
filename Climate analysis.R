# Climate analysis for Gentry tree project

# Set your working directory

# Load libraries
library(dplyr)
library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(MuMIn)
library(stringr)
library(ggpubr)
library(PairedData)
library(raster)
library(terra)
library(reshape2)
library(zoo)

###analyze climate data. 
#read CRU data to calculate mcwd
cru_temp <- read_excel("CEDA_climate.xlsx", 
                                sheet = "temperature")

cru_temp <- cru_temp %>%
  pivot_longer(cols = -c(1), names_to = "month", values_to = "temp") %>%
  mutate(month = str_replace(month, "X", "")) %>%
  mutate(month = as.numeric(month))

cru_tmin <- read_excel("CEDA_climate.xlsx", 
                       sheet = "tmin")

cru_tmin <- cru_tmin %>%
  pivot_longer(cols = -c(1), names_to = "month", values_to = "tmin") %>%
  mutate(month = str_replace(month, "X", "")) %>%
  mutate(month = as.numeric(month))

cru_tmax <- read_excel("CEDA_climate.xlsx", 
                       sheet = "tmax")

cru_tmax <- cru_tmax %>%
  pivot_longer(cols = -c(1), names_to = "month", values_to = "tmax") %>%
  mutate(month = str_replace(month, "X", "")) %>%
  mutate(month = as.numeric(month))

cru_cld <- read_excel("CEDA_climate.xlsx", 
                      sheet = "cloud cover")

cru_cld <- cru_cld %>%
  pivot_longer(cols = -c(1), names_to = "month", values_to = "cloud_cover") %>%
  mutate(month = str_replace(month, "X", "")) %>%
  mutate(month = as.numeric(month))

cru_vpr <- read_excel("CEDA_climate.xlsx", 
                      sheet = "vapor pressure")

cru_vpr <- cru_vpr %>%
  pivot_longer(cols = -c(1), names_to = "month", values_to = "vapor_pressure") %>%
  mutate(month = str_replace(month, "X", "")) %>%
  mutate(month = as.numeric(month))

cru_precip <- read_excel("CEDA_climate.xlsx", 
                      sheet = "precipitation")

cru_precip <- cru_precip %>%
  pivot_longer(cols = -c(1), names_to = "month", values_to = "precip") %>%
  mutate(month = str_replace(month, "X", "")) %>%
  mutate(month = as.numeric(month))

cru_ET <- read_excel("CEDA_climate.xlsx", 
                         sheet = "evapotranspiration")

cru_ET <- cru_ET %>%
  pivot_longer(cols = -c(1), names_to = "month", values_to = "ET") %>%
  mutate(month = str_replace(month, "X", "")) %>%
  mutate(month = as.numeric(month))

cru_climate <- left_join(cru_temp, cru_cld, by = c("Year", "month")) %>%
  left_join(cru_vpr, by = c("Year", "month")) %>%
  left_join(cru_precip, by = c("Year", "month")) %>%
  left_join(cru_tmin, by = c("Year", "month")) %>%
  left_join(cru_tmax, by = c("Year", "month")) %>%
  left_join(cru_ET, by = c("Year", "month")) 

cru_climate$vapor_pressure <- cru_climate$vapor_pressure / 10 #convert to kilopascals

#calculate evapotranspiration using the Penman-Monteith equation

calculate_ET <- function(temp, cloud_cover, vapor_pressure, tmax, tmin) {
  
  #constants
  lambda <- 2.45
  gamma <- (0.00163*(98.4/2.45)) #psychometric constant
  
  # Saturation vapor pressure at T_max and T_min
  e_s_max <- 0.6108 * exp((17.27 * tmax) / (tmax + 237.3))
  e_s_min <- 0.6108 * exp((17.27 * tmin) / (tmin + 237.3))
  
  # Mean saturation vapor pressure
  e_s <- (e_s_max + e_s_min) / 2
  
  # Actual vapor pressure
  e_a <- vapor_pressure
  
  # Slope of vapor pressure curve
  delta <- (4098 * e_s) / (temp + 237.3)^2
  
  # Net radiation (estimated using cloud cover percentage). 35.5 is the average extraterrestrial radiation for Madre de Dios.
  Rn <- 0.77 * (35.5 * (1 - cloud_cover/100))
  
  # Wind speed
  U <- 0.8 #average wind speed in Puerto Maldonado in m/s
  
  #Penman-Monteith equation
  ET <- ((0.408 * delta * (Rn)) + (gamma * (900 / (temp + 273)) * U * (e_s - e_a))) / 
    (delta + gamma * (1 + 0.34 * U))
  
  return(ET)
}

cru_climate$evap_day <- calculate_ET(cru_climate$temp, cru_climate$cloud_cover, cru_climate$vapor_pressure, cru_climate$tmax, cru_climate$tmin)


#calculate monthly totals
thirty <- c(4, 6, 9, 11)
thirtyone <- c(1, 3, 5, 7, 8, 10, 12)

cru_climate$evap <- ifelse(cru_climate$month %in% thirty, cru_climate$evap_day * 30, ifelse(cru_climate$month %in% thirtyone, cru_climate$evap_day * 31, cru_climate$evap_day * 28))

cru_climate$wd <- pmax(0, cru_climate$evap - cru_climate$precip)

# calculate saturation vapor pressure and relative humidity
cru_climate$sat_vp <- 0.6113 * exp((2500000 / 461)*((1/273.15) - (1/ (cru_climate$temp + 273.15))))
cru_climate$RH <- (cru_climate$vapor_pressure / cru_climate$sat_vp) * 100

#write.csv(cru_climate, "cru_monthly_climate.csv")

#summarize to calculate mcwd for each year
climate_year <- cru_climate %>%
  group_by(Year) %>%
  summarise(precip = sum(precip), 
            tmax = mean(tmax), 
            tmin = mean(tmin), 
            tmean = mean(temp),
            mcwd = max(cumsum(wd)) * -1)

summary(lm(mcwd ~ Year, data = climate_year))

#make the year column in a date format
climate_year$Year <- as.Date(paste0(climate_year$Year, "-01-01"))

ggplot(climate_year, aes(x = Year, y = mcwd)) +
  annotate("rect", xmin = as.Date("2005-01-01"), xmax = as.Date("2006-01-01"), ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#C4961A") +
  annotate("rect", xmin = as.Date("1997-01-01"), xmax = as.Date("1998-01-01"), ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#C4961A") +
  annotate("rect", xmin = as.Date("2010-01-01"), xmax = as.Date("2011-01-01"), ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#C4961A") +
  annotate("rect", xmin = as.Date("2015-01-01"), xmax = as.Date("2017-01-01"), ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#C4961A") +
  annotate("rect", xmin = as.Date("2023-01-01"), xmax = as.Date("2025-01-01"), ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#C4961A") +
  geom_vline(xintercept = as.Date(c("1984-01-01", "1985-01-01", "1986-01-01", "1987-01-01", "2003-01-01", "2008-01-01", "2011-01-01", "2017-01-01", "2023-01-01")), linetype = "dashed", color = "gray45", linewidth = 1) +
  geom_smooth(linetype = "solid", color = "black", method = "lm", se = FALSE) +
  geom_line(linewidth = 1.5, color = "#D16103") +
   labs(title = "Maximum cumulative water deficit by year",
       x = "Date",
       y = "MCWD (mm)") +
  
  theme_bw()

#input co2 data
co2  <- read.csv(file = "co2_mm_mlo.csv")
co2$date <- as.Date(paste0(co2$year, "-", co2$month, "-01"))
co2 <- co2[co2$year >= 1983,]

ggplot(co2, aes(x = date, y = average)) +
  geom_vline(xintercept = as.Date(c("1984-01-01", "1985-01-01", "1986-01-01", "1987-01-01", "2003-01-01", "2008-01-01", "2011-01-01", "2017-01-01", "2023-01-01")), linetype = "dashed", color = "gray45", linewidth = 1) +
  geom_line(linewidth = 1, color= "#4E84C4") +
  labs(title = "Atmospheric [CO2]",
       x = "Date",
       y = "CO2 mole fraction (ppm)") +
  #scale_x_continuous(breaks = seq(1980, 2020, by = 5)) +
  theme_bw()

co22 <- co2 %>%
  group_by(year) %>%
  summarize(mean = mean(average))

summary(lm(mean ~ year, data = co22))

#make a new dataframe with mcwd from climate_year and average from co2
new <- merge(climate_year, co2, by.x = "Year", by.y = "date", all = TRUE)
#now make it long keeping only average and mcwd
new <- new[,c("Year", "average", "mcwd")]
new <- melt(new, id.vars = "Year")

#interpolate missing values
new$value <- na.approx(new$value, na.rm = FALSE)

labels <- c("average" = "CO2 mole fraction (ppm)", "mcwd" = "MCWD (mm)", "mat" = "Mean annual temperature (C)")
cols <- c("#4E84C4", "#D16103", "darkred")

#make vector of year and month 
supp_table1 <- read.csv("supp_table1.csv")
dates <- as.Date(paste0(supp_table1$year, "-", supp_table1$Month, "-01"))

#make new dataframe for mean temperature from 1983 to 2023, with temp starting at 26.7 and increasing by 0.036 per year
#mat <- data.frame(Year = seq(as.Date("1983-01-01"), as.Date("2024-07-01"), by = "month"))
#mat$variable <- "mat"
#mat$value <- seq(from = 26.7, by = 0.003, length.out = nrow(mat))

#add mat rows to new as a new variable
#new <- rbind(new, mat)


#plot
ggplot(new, aes(x = Year, y = value, color = variable)) +
  annotate("rect", xmin = as.Date("2005-01-01"), xmax = as.Date("2006-01-01"), ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#C4961A") +
  annotate("rect", xmin = as.Date("1997-01-01"), xmax = as.Date("1998-01-01"), ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#C4961A") +
  annotate("rect", xmin = as.Date("2010-01-01"), xmax = as.Date("2011-01-01"), ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#C4961A") +
  annotate("rect", xmin = as.Date("2015-01-01"), xmax = as.Date("2017-01-01"), ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#C4961A") +
  annotate("rect", xmin = as.Date("2023-01-01"), xmax = as.Date("2025-01-01"), ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "#C4961A") +
  geom_vline(xintercept = dates, linetype = "dashed", color = "gray45", linewidth = 1) +
  geom_vline(xintercept = as.Date("2023-09-01"), linetype = "dashed", color = "gray45", linewidth = 1) +
  geom_smooth(data = subset(new, variable == "mcwd"), linetype = "solid", color = "black", method = "lm", se = FALSE) +
  geom_line(linewidth = 1.2) +
  #geom_smooth(linetype = "solid", color = "black", method = "lm", se = FALSE) +
  labs(title = "Changes in climate",
       x = "Date",
       y = NULL) +
  scale_color_manual(values = cols) +
  facet_grid(variable ~ ., scales = "free_y", labeller = labeller(variable = labels), switch = "y") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0, hjust = 0),
        strip.placement = "outside",
        strip.background = element_blank(),
        legend.position = "none") 

