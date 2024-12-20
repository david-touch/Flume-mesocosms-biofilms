
#Flume-mesocosms-biofilms

#Graph of flow and temperature 

library(tidyverse)
library(reshape2)
library(ggpubr)
library(tidyquant)
library(zoo)
library(readr)
library(lubridate)
library(tidyr)

data_2022 <- read.csv("data_2022.csv")
data_2022 <- data_2022 %>% 
  filter(!row_number() %in% c(1:38, 9927:9930))
dates_2022 <- read.csv("dates_2022.csv")
dates_2022$Date <- as.Date(dates_2022$Date, "%Y-%m-%d")
dates_2022$Date <- format(as.POSIXct(dates_2022$Date, format = "%Y-%m-%d"), "%Y-%m-%d %H:%M:%S")
dates_2022$Time[dates_2022$Time == "2022-06-10 09:30:00"] <- "2022-06-10 00:00:00"
timepoints_2022 <- dates_2022 %>% 
  filter(Time != "NA")
dates_2022 <- dates_2022 %>% rename(Sampling = Time, Time = Date)

data_2022_temp <- data_2022 %>% 
  select(date, tankTemp0,tankTemp1) %>% 
  dplyr::rename(Time = date, warm = tankTemp1, control = tankTemp0) %>%
  melt(variable.name = "Treatment", value.name = "Temperature") %>% 
  mutate(Date = substring(Time, 1, 10))
data_2022_temp$Date <- as.Date(data_2022_temp$Date, "%Y-%m-%d")
data_2022_temp$Date <- format(data_2022_temp$Date, "%b %d")

##Plot temperature
df.sliced_temp <- data_2022_temp %>%
  filter(Treatment == "control") %>% 
  slice(round(seq(1, n(), length.out = 6)))

data_2022_temp_avg <- data_2022_temp %>%
  group_by(Date, Treatment) %>% 
  mutate(TempAv = mean(Temperature, na.rm = TRUE)) %>%
  unite("TreatDate", c(Treatment,Date), remove = FALSE) %>% 
  distinct(TreatDate, .keep_all = TRUE) %>% 
  arrange(Time) %>% 
  left_join(dates_2022 %>% select(Time, Incubation), by = "Time") %>% 
  mutate(Incubation = ifelse(is.na(Incubation), 0, Incubation))
data_2022_temp_avg$Time[data_2022_temp_avg$Time == "2022-06-10 09:30:00"] <- "2022-06-10 00:00:00"

df.sliced_temp_avg <- data_2022_temp_avg %>%
  ungroup() %>% 
  select(-TempAv, -TreatDate) %>% 
  filter(Treatment == "control") %>% 
  filter(Date == "Jun 10" | Date == "Jun 30" | Date == "Jul 20" | Date == "Aug 09" | Date == "Aug 29" | Date == "Sep 18") %>% 
  rename(Datetoadd = Date) %>% 
  mutate(Datetoadd = paste0("(", Datetoadd, ")")) %>% 
  select(-Treatment, -Temperature)

desired_incubation <- seq(0, 100, by = 5)
filtered_data <- data_2022_temp_avg %>%
  filter(Incubation %in% desired_incubation) %>%
  select(Time, Incubation) %>% 
  left_join(df.sliced_temp_avg %>% select(Incubation, Datetoadd), by = "Incubation") %>% 
  mutate(Datetoadd = ifelse(is.na(Datetoadd), "", Datetoadd)) %>% 
  unite(Label, Incubation, Datetoadd, sep = " ", remove = FALSE) %>% 
  filter(Treatment == "control") %>% 
  mutate(Label = str_replace(Label, "(\\d+)\\s", "\\1\n"))

#Figure 1C
plot_temp_avg <- data_2022_temp_avg %>% ggplot(aes(x = Time, y = TempAv)) +
  geom_line(aes(color = Treatment, group = Treatment),linewidth = 0.75, lineend = "round") +
  scale_color_manual(values = alpha(c("royalblue", "darkorange"), 1),labels = c("Control","Warm")) +
  geom_point(data = timepoints_2022, aes(x = Date, y = 0), color = "black", fill = "black", size = 2.5, shape = 25, position = position_dodge(width=0.75)) +
  scale_x_discrete(
    breaks = filtered_data$Time,
    labels = setNames(filtered_data$Label, filtered_data$Time)) +
  ylab("Temperature (ºC)") +
  xlab("Days of growth") +
  labs(color = "Temperature\nregimes") +
  theme_classic() +
  theme(legend.title = element_text(size = 16, face="bold"),
        legend.text = element_text(size = 16), 
        axis.title.y = element_text(size=16, face="bold"), 
        axis.title.x = element_text(size=16, face="bold"),
        axis.text = element_text(size = 14),
        strip.background = element_blank(),
        axis.line=element_line(),
        legend.position = "right", 
        legend.justification = "left") +
  guides(color = guide_legend(override.aes = list(linewidth = 3)))
plot_temp_avg
ggsave("plot_temp_avg.png", plot_temp_avg)


data_2022_temp_window <- data_2022_temp %>% 
  mutate(Time2 = Time)
data_2022_temp_window$Time2 <- as.Date(data_2022_temp_window$Time2, "%Y-%m-%d") 
data_2022_temp_window <- data_2022_temp_window %>% 
  filter(Time2 > '2022-06-29') %>% 
  filter(Time2 < '2022-07-04')

#Figure 1C inset
plot_temp_window <- data_2022_temp_window %>% ggplot(aes(x = Time, y = Temperature)) +
  geom_line(aes(color = Treatment, group = Treatment),linewidth = 0.75, lineend = "round") +
  scale_color_manual(values = alpha(c("royalblue", "darkorange"), 1),labels = c("Ambient","Warm")) +
  xlab("Days of growth") +
  ylim(5,20) +
  labs(color = "Temperature\nregimes") +
  theme_classic() +
  theme(legend.title = element_text(size = 16, face="bold"),
        legend.text = element_text(size = 16), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        axis.line.x=element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = "none") +
  guides(color = guide_legend(override.aes = list(linewidth = 3)))
plot_temp_window
ggsave("plot_temp_window.png", plot_temp_window, width = 3.5, height = 3, units = "cm")


##Plot Flow rate 
data_flowrate_2022 <- read.csv("20220610_20220922_flowrate.csv")
data_flowrate_2022 <- data_flowrate_2022 %>% 
  select(-Q_to_pump) %>% 
  melt(variable.name = "Treatment", value.name = "Flow") %>%
  mutate(Date = substring(time, 1, 10)) %>% 
  dplyr::rename(Time = time)
data_flowrate_2022$Date <- as.Date(data_flowrate_2022$Date, "%d/%m/%Y")
data_flowrate_2022$Date <- format(data_flowrate_2022$Date, "%b %d")
data_flowrate_2022$Time <- gsub("/","-",data_flowrate_2022$Time)
data_flowrate_2022$Time <- format(as.POSIXct(data_flowrate_2022$Time, format = '%d-%m-%Y %H:%M'), "%Y-%m-%d %H:%M:%S")

df.sliced_flow <- data_flowrate_2022 %>%
  filter(Treatment == "natural") %>% 
  slice(round(seq(1, n(),
                  length.out = 6)))

#Figure 1B
plot_flow <- data_flowrate_2022 %>% ggplot(aes(x = Time, y = Flow)) +
  geom_line(aes(color = Treatment, group = Treatment), linewidth = 0.75, lineend = "round") +
  scale_color_manual(values = alpha(c("#85929e", "#ec7063", "#52be80", "#5dade2"), 1), labels = c("Natural","Intermittent","Stochastic","Constant")) +
  geom_point(data = timepoints_2022, aes(x = Date, y = 0.02), color = "black", fill = "black", size = 2.5, shape = 25) +
  ylab("Discharge (L/s)") +
  xlab("Days of growth") +
  labs(color = "Flow regimes") +
  scale_x_discrete(
    breaks = filtered_data$Time,
    labels = setNames(filtered_data$Label, filtered_data$Time)) +
  theme_classic() +
  theme(legend.title = element_text(size = 16, face="bold"),
        legend.text = element_text(size = 16), 
        axis.title.y = element_text(size=16, face="bold"), 
        axis.title.x = element_text(size=16, face="bold"),
        axis.text = element_text(size = 14),
        strip.background = element_blank(),
        axis.line=element_line(),
        legend.position = "right",
        legend.justification = "left") +
  guides(color = guide_legend(override.aes = list(linewidth = 3)))
plot_flow

#Figure B-C combined
plot_treatments <- ggarrange(plot_flow, plot_temp_avg, ncol = 1, nrow = 2, align = "v")
plot_treatments
ggsave("plot_treatments.png", plot_treatments,  width = 30, height = 20, units = "cm")

#Mean flow and SD
mean(data_flowrate_2022$Flow,na.rm = TRUE)
sd(data_flowrate_2022$Flow,na.rm = TRUE)


##SD Difference between warm and control
# Convert "Time" to a datetime object if not already
data_2022_temp$Time <- as.POSIXct(data_2022_temp$Time, format="%Y-%m-%d %H:%M:%S")
# Split the data into control and warm
control_df <- data_2022_temp[data_2022_temp$Treatment == "control", ]
warm_df <- data_2022_temp[data_2022_temp$Treatment == "warm", ]
# Merge the control and warm dataframes on Time
merged_df <- merge(control_df, warm_df, by="Time", suffixes = c("_control", "_warm"))
# Calculate the temperature difference
merged_df$Temp_Difference <- merged_df$Temperature_warm - merged_df$Temperature_control
# Average and SD
mean(merged_df$Temp_Difference,na.rm = TRUE)
sd(merged_df$Temp_Difference,na.rm = TRUE)


#Graph of channel depth and flow velocity

# flume characteristics
# slope of the channel
slope <- 0.03 # 3%

# channel width (meters)
width_cm <- 3        # 3 cm
width <- width_cm / 100 # Convert cm to meters

# channel height (meters)
height_cm <- 4        # 4 cm
height <- height_cm / 100 # Convert cm to meters

# manning roughness coefficient
n <- 0.035

# define functions for calculations
# function to calculate water depth (y) for a given discharge (Q)
calculate_y <- function(Q, b, S, n, h_max) {
  # equation to solve: Q = (1/n) * (R)^(2/3) * S^(1/2) * A
  # for a rectangular channel: A = b * y, P = b + 2y, R = A / P
  
  # define the equation as a function of y
  eq <- function(y) {
    A <- b * y             # cross-sectional area
    P <- b + 2 * y         # wetted perimeter
    R <- A / P             # hydraulic radius
    V <- (1 / n) * (R^(2/3)) * sqrt(S) # manning's equation for velocity
    Q_calc <- V * A        # volumetric flow rate
    return(Q_calc - Q)     # the equation to be zeroed
  }
  
  # check if the discharge is feasible for the given channel height
  Q_at_h_max <- eq(h_max)
  
  if (Q_at_h_max < 0) {
    warning("Discharge Q =", Q, "m³/s exceeds channel capacity for the given height.")
    return(NA)
  }
  
  # use uniroot to find the root (water depth) within the interval (0, h_max)
  y_solution <- tryCatch({
    uniroot(eq, lower = 0.0001, upper = h_max)$root
  }, error = function(e) {
    warning(paste("Failed to find water depth for Q =", Q, "m³/s"))
    return(NA)
  })
  
  return(y_solution)
}

# function to calculate flow velocity (V) given water depth (y)
calculate_V <- function(y, b, S, n) {
  A <- b * y             # cross-sectional area
  P <- b + 2 * y         # wetted perimeter
  R <- A / P             # hydraulic radius
  V <- (1 / n) * (R^(2/3)) * sqrt(S) # manning's equation for velocity
  return(V)
}

# load and prepare data
data <- read_csv("20220610_20220922_flowrate.csv")

# discharge measurements
discharge_types <- c("natural", "intermittent", "stochastic", "constant")

results_flow <- data %>%
  # step 1: adjust the discharge measurements by dividing by 3 (sub-channels per flume) and converting to m³/s
  mutate(
    # adjust discharge by dividing by 3 for each type
    across(all_of(discharge_types), ~ .x / 3, .names = "Q_L_s_{col}"), # adjusted discharge in L/s
    # convert from L/s to m³/s
    across(starts_with("Q_L_s_"), ~ .x / 1000, .names = "Q_m3_s_{sub('Q_L_s_', '', .col)}")
  ) %>%
  # step 2: calculate y_cm and V_m_s for each discharge type
  rowwise() %>%
  mutate(
    # for 'natural' discharge
    y_m_natural = calculate_y(Q_m3_s_natural, width, slope, n, height),
    y_cm_natural = round(y_m_natural * 100, 3),                    # rounded to 3 decimals
    V_m_s_natural = round(ifelse(!is.na(y_m_natural), calculate_V(y_m_natural, width, slope, n), NA), 3),
    
    # for 'intermittent' discharge
    y_m_intermittent = calculate_y(Q_m3_s_intermittent, width, slope, n, height),
    y_cm_intermittent = round(y_m_intermittent * 100, 3),          # rounded to 3 decimals
    V_m_s_intermittent = round(ifelse(!is.na(y_m_intermittent), calculate_V(y_m_intermittent, width, slope, n), NA), 3),
    
    # for 'stochastic' discharge
    y_m_stochastic = calculate_y(Q_m3_s_stochastic, width, slope, n, height),
    y_cm_stochastic = round(y_m_stochastic * 100, 3),              # rounded to 3 decimals
    V_m_s_stochastic = round(ifelse(!is.na(y_m_stochastic), calculate_V(y_m_stochastic, width, slope, n), NA), 3),
    
    # for 'constant' discharge
    y_m_constant = calculate_y(Q_m3_s_constant, width, slope, n, height),
    y_cm_constant = round(y_m_constant * 100, 3),                  # rounded to 3 decimals
    V_m_s_constant = round(ifelse(!is.na(y_m_constant), calculate_V(y_m_constant, width, slope, n), NA), 3)
  ) %>%
  ungroup() %>%
  # select only the useful columns
  select(
    time,
    natural, y_cm_natural, V_m_s_natural,
    intermittent, y_cm_intermittent, V_m_s_intermittent,
    stochastic, y_cm_stochastic, V_m_s_stochastic,
    constant, y_cm_constant, V_m_s_constant
  ) %>% 
  replace(is.na(.), 0)


##Plot depth
results_depth <- results_flow %>%
  rename(Time = time, Natural =y_cm_natural, Intermittent = y_cm_intermittent, Stochastic =y_cm_stochastic, Constant = y_cm_constant) %>% 
  select(Time, Natural, Intermittent, Stochastic, Constant) %>% 
  melt(variable.name = "Treatment", value.name = "Depth")

results_depth$Time <- gsub("/","-",results_depth$Time)
results_depth$Time <- format(as.POSIXct(results_depth$Time, format = '%d-%m-%Y %H:%M'), "%Y-%m-%d %H:%M:%S")

#Figure S1A
plot_depth <- results_depth %>% ggplot(aes(x = Time, y = Depth)) +
  geom_line(aes(color = Treatment, group = Treatment), linewidth = 0.75, lineend = "round") +
  scale_color_manual(values = alpha(c("#85929e", "#ec7063", "#52be80", "#5dade2"), 1), labels = c("Natural","Intermittent","Stochastic","Constant")) +
  geom_point(data = timepoints_2022, aes(x = Date, y = 0.02), color = "black", fill = "black", size = 2.5, shape = 25) +
  ylab("Channel depth (cm)") +
  xlab("Days of growth") +
  labs(color = "Flow regimes") +
  scale_x_discrete(
    breaks = filtered_data$Time,
    labels = setNames(filtered_data$Label, filtered_data$Time)) +
  theme_classic() +
  theme(legend.title = element_text(size = 16, face="bold"),
        legend.text = element_text(size = 16), 
        axis.title.y = element_text(size=16, face="bold"), 
        axis.title.x = element_text(size=16, face="bold"),
        axis.text = element_text(size = 14),
        strip.background = element_blank(),
        axis.line=element_line(),
        legend.position = "right",
        legend.justification = "left") +
  guides(color = guide_legend(override.aes = list(linewidth = 3)))
plot_depth

results_velocity <- results_flow %>%
  rename(Time = time, Natural = V_m_s_natural, Intermittent = V_m_s_intermittent, Stochastic = V_m_s_stochastic, Constant = V_m_s_constant) %>% 
  select(Time, Natural, Intermittent, Stochastic, Constant) %>% 
  melt(variable.name = "Treatment", value.name = "Velocity") %>% 
  mutate(Velocity = Velocity*100)

results_velocity$Time <- gsub("/","-",results_velocity$Time)
results_velocity$Time <- format(as.POSIXct(results_velocity$Time, format = '%d-%m-%Y %H:%M'), "%Y-%m-%d %H:%M:%S")

#Figure S1A
plot_velocity <- results_velocity %>% ggplot(aes(x = Time, y = Velocity)) +
  geom_line(aes(color = Treatment, group = Treatment), linewidth = 0.75, lineend = "round") +
  scale_color_manual(values = alpha(c("#85929e", "#ec7063", "#52be80", "#5dade2"), 1), labels = c("Natural","Intermittent","Stochastic","Constant")) +
  geom_point(data = timepoints_2022, aes(x = Date, y = 0.02), color = "black", fill = "black", size = 2.5, shape = 25) +
  ylab("Velocity (cm/s)") +
  xlab("Days of growth") +
  labs(color = "Flow regimes") +
  scale_x_discrete(
    breaks = filtered_data$Time,
    labels = setNames(filtered_data$Label, filtered_data$Time)) +
  theme_classic() +
  theme(legend.title = element_text(size = 16, face="bold"),
        legend.text = element_text(size = 16), 
        axis.title.y = element_text(size=16, face="bold"), 
        axis.title.x = element_text(size=16, face="bold"),
        axis.text = element_text(size = 14),
        strip.background = element_blank(),
        axis.line=element_line(),
        legend.position = "right",
        legend.justification = "left") +
  guides(color = guide_legend(override.aes = list(linewidth = 3)))
plot_velocity

#Figure S1
FigureS1 <- ggarrange(plot_depth, plot_velocity, labels = c("A", "B"), font.label = list(size = 18), ncol = 1, nrow = 2, align = "v")
FigureS1
ggsave("FigureS1.pdf", FigureS1,  width = 30, height = 20, units = "cm")


# Calculate average depth and velocity
mean(results_depth$Depth,na.rm = TRUE)
sd(results_depth$Depth,na.rm = TRUE)
mean(results_velocity$Velocity,na.rm = TRUE)
sd(results_velocity$Velocity,na.rm = TRUE)

