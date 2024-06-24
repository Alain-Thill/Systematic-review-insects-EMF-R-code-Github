
# run this first!

#rm(list=ls())  #clear global environment of stuff from previous sessions

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to active folder

dir.create(file.path(getwd(), "tables")) # create subdirectory for data tables

pacman::p_load(dplyr, tidyr, ggplot2, ggpubr, scales, RColorBrewer, openxlsx)


#########################################################
# Verifying exposure duration values (optional)
#########################################################

df1 <- openxlsx::read.xlsx("data_table_HF.xlsx", sheet = 1)

df2 <- openxlsx::read.xlsx("data_table_LF.xlsx", sheet = 1) 
df2$RFLow <- as.character(df2$RFLow)
df2$RFHigh <- as.character(df2$RFHigh)
df2$SAR <- as.character(df2$SAR)
df2$Output_power <- as.character(df2$Output_power)
df <- full_join(df1,df2)

# total number of individual studies
distinct(df, key)

# standardizing duration: transforming European comma notation to US dot notation 
df$mins <- gsub(",", ".", df$mins)
df$hours <- gsub(",", ".", df$hours)
df$days <- gsub(",", ".", df$days)
df$Hrperday <- gsub(",", ".", df$Hrperday)

# filling in 0 in NA fields
df <- df %>% mutate_at("mins", ~replace_na(.,"0"))
df <- df %>% mutate_at("hours", ~replace_na(.,"0"))
df <- df %>% mutate_at("days", ~replace_na(.,"0"))
df <- df %>% mutate_at("weeks", ~replace_na(.,"0"))
df <- df %>% mutate_at("years", ~replace_na(.,"0"))
df <- df %>% mutate_at("Hrperday", ~replace_na(.,"24"))
df <- df %>% mutate_at("DaysWk", ~replace_na(.,"7"))
df$mins <- as.numeric(df$mins)
df$hours <- as.numeric(df$hours)
df$days <- as.numeric(df$days)
df$weeks <- as.numeric(df$weeks)
df$years <- as.numeric(df$years)
df$Hrperday <- as.numeric(df$Hrperday)
df$DaysWk <- as.numeric(df$DaysWk)

#control: recalculating duration from minute, hour etc. fields
attach(df)
df$CummHrs2 <- (mins/60) + hours + (days*Hrperday) + (weeks*DaysWk*Hrperday) + (years*365*24)
detach(df)

control <- df[c("key","CummHrs","CummHrs2")]
control
#View(control) # No errors found. Minute -- year fields will be removed later
# 0 duration is for inventories of insects at base stations, where duration of exposure is unknown


##############################################################
# Pie-chart and bar-chart of insect species/types
##############################################################

# putting HF and LF data tables together (adding EMF_type variable):
df1 <- openxlsx::read.xlsx("data_table_HF.xlsx", sheet = 1)
df1 <- df1 %>% arrange(key) # sort by study key, i.e. author-year
df1$EMF_type <- "HF"
df1 <- mutate(df1, EMF_type = if_else ( SignalGen == '50-Hz magnetic field', 'LF', EMF_type))

df2<-openxlsx::read.xlsx("data_table_LF.xlsx", sheet = 1)
df2 <- df2 %>% arrange(key) # sort by study key
df2$EMF_type <- "LF"
df2$RFLow <- as.character(df2$RFLow)
df2$RFHigh <- as.character(df2$RFHigh)
df2$SAR <- as.character(df2$SAR)
df2$Output_power <- as.character(df2$Output_power)

df <- full_join(df1,df2)

str(df)
myvars <- names(df) %in% c("mins", "hours", "days", "weeks", "years", "Hrperday", "DaysWk")
df <- df[!myvars] # removing columns listed in "myvars" from dataframe (!myvars)

names(df)

df$EMF_source <- NA
  
# standardizing naming of insects
df$Insect_type <- df$CellsAnimal  # new column for insect types

df$Insect_type[grep("ant", df$Insect_type)] <- "Ant"
df$Insect_type[grep("Ants", df$Insect_type)] <- "Ant"
df$Insect_type[grep("Tetramorium", df$Insect_type)] <- "Ant"
df$Insect_type[grep("beetle", df$Insect_type)] <- "Beetle"
df$Insect_type[grep("Beetles", df$Insect_type)] <- "Beetle"
df$Insect_type[grep("Orius predatory bug", df$Insect_type)] <- "Beetle"
df$Insect_type[grep("Lasioderma", df$Insect_type)] <- "Beetle"
df$Insect_type[grep("cockroach", df$Insect_type)] <- "Cockroach"
df$Insect_type[grep("Drosophila", df$Insect_type)] <- "Drosophila"
df$Insect_type[grep("honey", df$Insect_type)] <- "Honeybee"
df$Insect_type[grep("Honey", df$Insect_type)] <- "Honeybee"
df$Insect_type[grep("Wasps", df$Insect_type)] <- "Wasp"
df$Insect_type[grep("Vespula", df$Insect_type)] <- "Wasp"
df$Insect_type[grep("parasitic wasp", df$Insect_type)] <- "Wasp"
df$Insect_type[grep("locust", df$Insect_type)] <- "Locust"
df$Insect_type[grep("Locusta", df$Insect_type)] <- "Locust"
df$Insect_type[grep("butterfly", df$Insect_type)] <- "Butterfly"
df$Insect_type[grep("common cutworm", df$Insect_type)] <- "Butterfly"
df$Insect_type[grep("Spodoptera", df$Insect_type)] <- "Butterfly"
df$Insect_type[grep("Bicyclus", df$Insect_type)] <- "Butterfly"
df$Insect_type[grep("springtail", df$Insect_type)] <- "Springtail"
df$Insect_type[grep("ticks", df$Insect_type)] <- "Tick" # technically not an insect, but arachnid
df$Insect_type[grep("stick insect", df$Insect_type)] <- "Stick insect"

table(df$Insect_type)
length(df$Insect_type) # 405 experiments or experimental groups in total

# standardizing naming of EMF sources
df$EMF_source <- df$SignalGen  # new column for insect types

df$EMF_source[grep("DECT", df$EMF_source)] <- "DECT"
df$EMF_source[grep("Base", df$EMF_source)] <- "Base station"
df$EMF_source[grep("station", df$EMF_source)] <- "Base station"
df$EMF_source[grep("GSM", df$EMF_source)] <- "Mobile phone"
df$EMF_source[grep("DCS", df$EMF_source)] <- "Mobile phone"
df$EMF_source[grep("mobile", df$EMF_source)] <- "Mobile phone"
df$EMF_source[grep("Mobile", df$EMF_source)] <- "Mobile phone"
df$EMF_source[grep("Magnetic", df$EMF_source)] <- "Coil system 50/60 Hz"
df$EMF_source[grep("magnetic", df$EMF_source)] <- "Coil system 50/60 Hz"
df$EMF_source[grep("line", df$EMF_source)] <- "Power line"
df$EMF_source[grep("Line", df$EMF_source)] <- "Power line"
df$EMF_source[grep("WiFi", df$EMF_source)] <- "WiFi"
df$EMF_source[grep("RF", df$EMF_source)] <- "RF signal generator"
df$EMF_source[grep("Plate capacitor", df$EMF_source)] <- "Plate capacitor 50 Hz"

table(df$EMF_source) # to verify if all names are similar, and grouping is tidy and meaningful

df <- relocate(df, EMF_source, EMF_type, .after = "SignalGen")
df <- relocate(df, Insect_type, .after = "CellsAnimal")

names(df)

# to remove newly added studies (post publication)
df <- subset(df, key != "Vilic2024")
df <- subset(df, key != "Treder2023")
df <- subset(df, key != "Thill2018")

# saving HFLF table with standardized names
openxlsx::write.xlsx(df, "tables/data_table_HFLF_1.xlsx", rowNames = F, colWidths = 15, overwrite = T, firstRow = T, firstCol = T)


###########################################
# Pie-chart and bar-chart of insect types
###########################################

df <- openxlsx::read.xlsx("tables/data_table_HFLF_1.xlsx", sheet = 1, na.strings = "NA") 

sort(table(df$Insect_type), decreasing = F)

# selecting only insect types with more than 5 studies
table_insects <- df %>%
  group_by(Insect_type) %>% 
  filter(n() > 5) %>% 
  summarise(n = n())

table_insects

# calculate missing number, i.e. insect types with 5 or less studies
other_insects <- length(df$Insect_type) - sum(table_insects$n)
other_insects

# add "Other" to summary table
table_insects <- rbind(table_insects, list("Other", other_insects))


# to make this in German
#table_insects$Insect_type <- c("Ameise","Käfer","Schabe","Drosophila","Honigbiene","Heuschrecke","Springschwanz","Wespe","Andere")

# make simple pie-chart:
amount <- table_insects$n
labels <- table_insects$Insect_type
amount <- as.vector(amount)
pie(amount,labels)

################################################
# make nicer ggplot pie-chart and bar-chart:
################################################

count.data <- data.frame(
  Insect = labels,
  n = amount,
  prop = amount*100/sum(amount))


# Add label position
count.data <- count.data %>%
  arrange(desc(Insect)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)

count.data
maximum <- max(count.data$n)
count.data$prop <- round(count.data$prop, digits = 1)

p3 <-ggplot(count.data, aes(x = "", y = prop, fill = Insect)) +
  geom_bar(width = 1, stat = "identity", color = "white",  alpha=3/4) +
  geom_text(aes(y = lab.ypos, label = percent(prop/100)), size=3) +
  coord_polar("y", start = 0) + theme_void()
p3
#ggsave("piechart_species.pdf", width = 5, height = 4)
#ggsave("piechart_species.jpg", width = 5, height = 4)


# make bar-chart, vertical and horizontal
p3b <- ggplot(count.data, aes(x = reorder(Insect, -n), y = n, fill = Insect)) +
  geom_bar(stat = "identity",  alpha=3/4) + 
  geom_text(aes(y = n/2, label = percent(prop/100)), size=3) +
  theme(axis.text.x=element_text(angle=45, hjust=0.9), legend.position = "none") +
  labs(x = "Insect", y = "Number of experiments")
p3b


p3c <- ggplot(count.data, aes(x = n, y = reorder(Insect, +n), fill = Insect)) +
  geom_bar(stat = "identity", alpha=3/4) + 
  geom_text(aes(x = n/2, label = percent(prop/100)), size=3) +
  theme_classic(base_size = 12) + theme(legend.position = "none") +
  labs(title ="A", y = "", x = "Number of experiments") + geom_vline(xintercept = c(maximum), linetype = 3) +
  scale_x_continuous(n.breaks = 8)
p3c

#ggsave("figures/barchart_species.jpg", width = 4, height = 3)
#ggsave("figures/barchart_species.png", width = 4, height = 3)


##############################################################
# Pie-chart and bar-chart of EMF-types
##############################################################

df <- openxlsx::read.xlsx("tables/data_table_HFLF_1.xlsx", sheet = 1, na.strings = "NA") 


table_EMF <- df %>%
  group_by(EMF_source) %>% 
  filter(n() > 5) %>% 
  summarise(n = n())

table_EMF

# calculate missing number, i.e. insect types with 5 or less studies
other_EMF <- length(df$EMF_source) - sum(table_EMF$n)
other_EMF
table_EMF <- rbind(table_EMF, list("Other", other_EMF))

# to make this in German
#table_EMF$EMF_source <- c("Basisstation","DECT","Helmholtzspule","Handy","Plattenkondensator","Hochspannung","RF Signalgenerator","Andere") 

amount <- table_EMF$n
labels <- table_EMF$EMF_source
amount <- as.vector(amount)
pie(amount,labels)

count.data2 <- data.frame(
  EMF = labels,
  n = amount,
  prop = amount*100/sum(amount))

count.data2

# Add label position
count.data2 <- count.data2 %>%
  arrange(desc(EMF)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)

# get maximum count, for bar-chart line
maximum2 <- max(count.data2$n)
count.data2$prop <- round(count.data2$prop, digits = 1)


p4 <- ggplot(count.data2, aes(x = "", y = prop, fill = EMF)) +
  geom_bar(width = 1, stat = "identity", color = "white", alpha=3/4) +
  geom_text(aes(y = lab.ypos, label = percent(prop/100)), size=3) +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette="Set1") +
  theme_void()
p4


p4b <- ggplot(count.data2, aes(x = reorder(EMF, -n), y = n, fill = EMF)) +
  geom_bar(stat = "identity",  alpha=3/4) + 
  geom_text(aes(y = n/2, label = percent(prop/100)), size=3) +
  theme(axis.text.x=element_text(angle=45, hjust=0.9), legend.position = "none") +
  scale_fill_brewer(palette="Set1") +
  labs(x = "EMF", y = "Number of experiments")
p4b

count.data2$EMF
# putting in same order and color as in the following duration vs EMF field strength graphs and toxicity graphics (R codes 3,4)
count.data2$EMF <- factor(count.data2$EMF, levels=c("Base station","DECT","Mobile phone","RF signal generator","Coil system 50/60 Hz","Other","Power line","Plate capacitor 50 Hz","WiFi"), ordered=TRUE)
# for German:
#count.data2$EMF <- factor(count.data2$EMF, levels=c("Basisstation","DECT","Handy","RF Signalgenerator","Helmholtzspule","Andere","Hochspannung","Plattenkondensator"), ordered=TRUE)

p4c <- ggplot(count.data2, aes(x = n, y = reorder(EMF, +n), fill = EMF)) +
  geom_bar(stat = "identity", alpha=3/4) + geom_text(aes(x = n/2, label = percent(prop/100)), size=3) +
  theme_classic(base_size = 12) + theme(legend.position = "none") +
  labs(title ="B", y = "", x = "Number of experiments") + geom_vline(xintercept = c(maximum2), linetype = 3) +
  scale_fill_brewer(palette="Set1") +
  scale_x_continuous(n.breaks = 7)
p4c

#ggsave("figures/barchart_EMF.jpg", width = 4, height = 3)
#ggsave("figures/barchart_EMF.png", width = 4, height = 3)


# both bar-charts side-by-side
ggarrange(p3c, p4c, ncol = 2, nrow = 1, align = "hv")

ggsave("figures/barcharts_Insect_EMF.jpg", width = 10, height = 4)
ggsave("figures/barcharts_Insect_EMF.png", width = 10, height = 4)

ggarrange(p3c, p4c, ncol = 1, nrow = 2, align = "hv")

# both bar-charts arranged top-bottom
ggsave("figures/barcharts_Insect_EMF(vertical).jpg", width = 5, height = 8)
ggsave("figures/barcharts_Insect_EMF(vertical).png", width = 5, height = 8)

