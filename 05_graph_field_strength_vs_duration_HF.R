
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
getwd() # working directory set to folder of .R file

#rm(list=ls())  #clear global environment of stuff from previous sessions

library(pacman)
pacman::p_load(dplyr, tidyr, ggplot2, ggpubr, scales, RColorBrewer, openxlsx, sqldf)


###############################################################################
# Standardizing field strengths / power densities, EMF source, insect species
###############################################################################


df <- openxlsx::read.xlsx("tables/data_table_HFLF_1.xlsx", sheet = 1, na.strings = "NA") 

table(df$EMF_source)
# E_Field in V/m,  B_Field in uT

# conditional mutate: if _ then _ else
# to standardize all E-Field values to V/m where power densities given in mW/m2
df <- mutate(df, E_Field = if_else ( is.na(E_Field), sqrt(Power_mWm2*0.377), E_Field))

df$E_Field
df$SRF <- gsub(",", ".", df$SRF)
df$SRF <- as.numeric(df$SRF)

# to standardize all E-Field values to V/m where absorbed power given as SAR, 
# using formula from Geronikolou 2014: https://doi.org/10.1371/journal.pone.0112139
# SAR = (E_Field^2) * 1.19/1000
# SAR / (1.19/1000) = (E_Field^2) 
# (E_Field^2) = SAR * 1000/1.19
# E_Field = sqrt( SAR * 1000/1.19 )

df$SAR <- gsub(",", ".", df$SAR)
df$SAR <- as.numeric(df$SAR)
df <- mutate(df, E_Field = if_else ( is.na(E_Field), sqrt(SAR*1000/1.19), E_Field))

df$E_Field
df$B_Field
str(df)

# to standardize all E-Field values to V/m where RF intensity given as uT, using formula B = E/c
# or E = B * c
# since B is in microTesla, needs to be divided by 10^6
# for c = 3*10^8
# E = B * 300

df <- mutate(df, E_Field = if_else ( is.na(E_Field) & (EMF_type == "HF"), B_Field*300, E_Field))
df$E_Field

# get power density from E-Field if power density NA, for HF-EMF only
df <- mutate(df, Power_mWm2 = if_else ( is.na(Power_mWm2) & (EMF_type == "HF"), (E_Field^2)/0.377, Power_mWm2))
df$Power_mWm2 <- round(df$Power_mWm2, digits = 5)
df$E_Field <- round(df$E_Field, digits = 5)

# calculate equivalent magnetic field strength in uT for HF-EMF, with formula c = E*B 
#df <- mutate(df, B_Field = if_else ( is.na(B_Field), E_Field/300, B_Field))

names(df)

table(df$EMF_source)
df$EMF_source <- as.factor(df$EMF_source)

table(df$Insect_type)
df$Insect_type <- as.factor(df$Insect_type)

# renaming "key" column into "study"
names(df)[names(df) == 'key'] <- 'study'

# saving updated/completed table
openxlsx::write.xlsx(df, "tables/data_table_HFLF_2.xlsx", rowNames = F, colWidths = 15, overwrite = T, firstRow = T, firstCol = T)


# making table for SAR calculations used later
keep_col <- c("study","Title","SignalGen","Wavetype","CellsAnimal","N_exposed","N_control","SRF","SAR","E_Field",
              "CummHrs","Bioeffects", "Percent.change.vs.control","Confid","Direction_of_effect","Insect_type","EMF_type","EMF_source")
newdf <- subset(df, EMF_type =='HF')
newdf <- newdf %>% select(all_of(keep_col))
newdf <- arrange(newdf, SRF, CellsAnimal)


# create subdirectory for data tables
dir.create(file.path(getwd(), "tables"))
openxlsx::write.xlsx(newdf, "tables/SAR_frequencies.xlsx", rowNames = F, colWidths = 15, overwrite = T, firstRow = T, firstCol = T)


########################################################################
# Plots and histograms of duration and field strength / power density 
########################################################################

df <- openxlsx::read.xlsx("tables/data_table_HFLF_2.xlsx", sheet = 1, na.strings = "NA") 
df <- df[grep("HF", df$EMF_type),] # select only HF experiments


p1 <- ggplot(aes(x = E_Field), data = df) + geom_histogram(binwidth = 0.1) +
  xlab("Field strength [V / m]") + ylab("Number of experiments") + theme_classic2() +
  scale_x_log10(labels = label_number(accuracy = 0.01)) +  annotation_logticks(sides="b")
p1


p2 <- ggplot(aes(x = CummHrs), data = df) + geom_histogram(binwidth = 0.1) +
  xlab("Exposure duration [hours]") + ylab("Number of experiments") + theme_classic2() +
  scale_x_log10(labels = label_number(accuracy = 0.1, )) +  annotation_logticks(sides="b")
p2

ggarrange(p1, p2, ncol = 2, nrow = 1, align = "hv",
          hjust=0, label.x=0.3, legend="bottom", common.legend = T)

#ggsave("Strahlungsintensitaet_Expositionsdauer.jpg", width = 10, height = 4)

# cleaning up table
table(df$EMF_source)
df$EMF_source[grep("Oven", df$EMF_source)] <- NA
df$EMF_source[grep("FM Radio", df$EMF_source)] <- NA
df$EMF_source[grep("Coil system", df$EMF_source)] <- NA # removing because low frequency fields difficult to compare to HF
df$EMF_source[grep("Bluetooth", df$EMF_source)] <- NA
#df$EMF_source[grep("WiFi", df$EMF_source)] <- NA
df$EMF_source[grep("PEF", df$EMF_source)] <- NA

table(df$EMF_source)
df$EMF_source <- as.factor(df$EMF_source)

myvars <- names(df) %in% c("study", "Title", "EMF_source", "Insect_type", "E_Field", "CummHrs")
df2 <- df[myvars]
names(df2)
df2 <- df2 %>% mutate_at("Title", ~replace_na(.,""))
df2 <- mutate(df2, CummHrs = if_else ( is.na(CummHrs), 8760, CummHrs)) # set exposure duration for base station studies to 1 year (8760 hours)

#View(df2)

df2 <- na.omit(df2) # keeping only data lines with both E-field and duration values
length(df2$study) # 239 experiments
distinct(df2, study) # in 48 studies

color_table <- tibble(
  EMF_source = c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system"),
  #EMF_source = c("Basisstation","DECT","Handy","Signalgenerator","WLAN,"Spulensystem"),
  Color =       c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#999999", "#FF7F00"))

p3 <- ggplot(aes(x = CummHrs, y = E_Field, color = EMF_source), data = df2) + 
  geom_point() + theme_classic() + 
  scale_x_log10(labels = label_number(accuracy = 0.1, trim=T)) + 
  scale_y_log10(labels = label_number(accuracy = 0.01, trim=T)) +
  xlab("Exposure duration [hours]") + ylab("Field strength [V / m]") +
  annotation_logticks() + 
  scale_color_manual(values = color_table$Color, drop = FALSE) +
  theme(legend.position = "bottom") + labs(color='Field source') + geom_point(data=df2 %>% group_by(EMF_source) %>% summarise_at(vars(CummHrs,E_Field), median), size=10, shape=3, stroke=0.75) #for adding median or mean of point cloud

p3

ggsave("figures/Power_density_vs_duration_HF.jpg", width = 6, height = 6)
ggsave("figures/Power_density_vs_duration_HF.png", width = 6, height = 6)


############################################
# Pie-chart HF frequencies
############################################

df <- openxlsx::read.xlsx("tables/data_table_HFLF_2.xlsx", sheet = 1, na.strings = "NA") 
df <- df[grep("HF", df$EMF_type),] # select only HF experiments

myvars <- names(df) %in% c("study", "EMF_source", "SRF")
df <- df[myvars]
names(df)

df$SRF <- as.factor(df$SRF)
frequencies <- summary(df$SRF)
frequencies <- sort(frequencies, decreasing = T)
head(frequencies)
frequencies
length(df$SRF)
length(frequencies)

table_freq <- df %>%
  group_by(SRF) %>% 
  filter(n() > 1) %>% 
  summarise(n = n(), percentage = sum(n)/length(df))

table_freq

x <- summary(df$SRF)
labels <- names(x)
x <- as.vector(x)

jpeg("figures/piechart_frequencies[MHz].jpeg", quality = 95)
pie(x,labels)
dev.off()

################################################################
# getting max and min values of E-field and exposure duration
################################################################

df <- openxlsx::read.xlsx("tables/data_table_HFLF_2.xlsx", sheet = 1, na.strings = "NA")
df <- df[grep("HF", df$EMF_type),] # select only HF experiments

myvars <- names(df) %in% c("study", "EMF_source", "E_Field", "CummHrs")
df2 <- df[myvars]
names(df2)
df2 <- df[complete.cases(df[,c("study", "EMF_source", "E_Field", "CummHrs")]),] 


boxplot(df2$E_Field)
min(df2$E_Field)
(min(df2$E_Field)^2)/0.377
max(df2$E_Field)
(max(df2$E_Field)^2)/0.377

mean(df2$E_Field)
median(df2$E_Field)
sd(df2$E_Field)

min(df2$CummHrs)*60*60  # in Sekunden
max(df2$CummHrs)/(24*30) # in Monaten

boxplot(df2$CummHrs)
mean(df2$CummHrs)
median(df2$CummHrs)

EMF <- df[complete.cases(df[,c("study", "EMF_source", "E_Field")]),]
mobile <- sqldf("SELECT * FROM EMF WHERE EMF_source = 'Mobile phone'")
#View(df)
nrow(mobile)
mean(mobile$E_Field) # 18.66
(mean(mobile$E_Field)^2)/0.377 # 924.1
median(mobile$E_Field) # 16.19
(median(mobile$E_Field)^2)/0.377 # 695.26

tower <- sqldf("SELECT * FROM EMF WHERE EMF_source = 'Base station'")
nrow(tower)
mean(tower$E_Field) # 0.69
(mean(tower$E_Field)^2)/0.377 # 1.28
median(tower$E_Field) # 0.3475
(median(tower$E_Field)^2)/0.377 # 0.32

signalz <- sqldf("SELECT * FROM df WHERE EMF_source = 'RF signal generator'")
signalz$Wavetype

pulsed <- sqldf("SELECT * FROM signalz WHERE WaveType = 'Pulsed'")
FM <- sqldf("SELECT * FROM signalz WHERE WaveType = 'Continuous, 50-kHz FM'")
continuous <- sqldf("SELECT * FROM signalz WHERE WaveType = 'Continuous'")

nrow(pulsed)/nrow(signalz) # percentage experimental groups that used pulsed signal
nrow(FM)/nrow(signalz) # percentage experimental groups that used FM signal
1-(nrow(pulsed)/nrow(signalz))-(nrow(FM)/nrow(signalz)) # percentage experimental groups that used continuous signal

