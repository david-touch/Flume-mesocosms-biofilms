
#Flume-mesocosms-biofilms

#Graph of flow and temperature 

library(tidyverse)
library(reshape2)
library(ggpubr)
library(tidyquant)
library(zoo)

data_2022 <- read.csv("data_2022.csv")
data_2022 <- data_2022 %>% 
  filter(!row_number() %in% c(1:38, 9927:9930))
dates_2022 <- read.csv("dates_2022.csv")
dates_2022$Date <- as.Date(dates_2022$Date, "%Y-%m-%d")
dates_2022$Date <- format(as.POSIXct(dates_2022$Date, format = "%Y-%m-%d"), "%Y-%m-%d %H:%M:%S")

timepoints_2022 <- dates_2022 %>% 
  filter(Time != "NA")

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
  mutate(TempAv = mean(Temperature)) %>%
  drop_na() %>% 
  unite("TreatDate", c(Treatment,Date), remove = FALSE) %>% 
  distinct(TreatDate, .keep_all = TRUE)

df.sliced_temp_avg <- data_2022_temp_avg %>%
  select(-TempAv, -TreatDate) %>% 
  filter(Treatment == "control") %>% 
  filter(Date == "Jun 10" | Date == "Jun 30" | Date == "Jul 21" | Date == "Aug 11" | Date == "Aug 31" | Date == "Sep 21")

#Figure 1C
plot_temp_avg <- data_2022_temp_avg %>% ggplot(aes(x = Time, y = TempAv)) +
  geom_line(aes(color = Treatment, group = Treatment),linewidth = 0.75, lineend = "round") +
  scale_color_manual(values = alpha(c("royalblue", "darkorange"), 1),labels = c("Control","Warm")) +
  geom_point(data = timepoints_2022, aes(x = Date, y = 0), color = "black", fill = "black", size = 2.5, shape = 25, position = position_dodge(width=0.75)) +
  ylab("Temperature (ºC)") +
  xlab("Date") +
  labs(color = "Thermal regimes") +
  scale_x_discrete(breaks = df.sliced_temp_avg$Time, labels = df.sliced_temp_avg$Date, expand = expansion(mult = c(0, 0.01))) +
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
  xlab("Date") +
  ylim(5,20) +
  labs(color = "Thermal regimes") +
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
  xlab("Date") +
  labs(color = "Flow regimes") +
  scale_x_discrete(breaks = df.sliced_flow$Time, labels = df.sliced_temp$Date, expand = expansion(mult = c(0, 0.01))) +
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

