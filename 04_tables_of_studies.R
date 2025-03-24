

# run this first!
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to active directory

dir.create(file.path(getwd(), "tables")) # create subfolder for saving tables

#rm(list=ls())  #clear global environment of stuff from previous sessions

if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(dplyr, tidyr, gdata, ggplot2, ggpubr, scales, RColorBrewer, sqldf, xtable, stringr)


############################################################
#  table of studies by category: HF
############################################################

df <- openxlsx::read.xlsx("tables/data_table_HFLF_1.xlsx", sheet = 1, na.strings = "NA") 

df <- df[grep("HF", df$EMF_type),]  # select only HF studies
df$EMF_type
names(df)

EMF <- df[c("key","EMF_source")]
EMF <- unique(EMF)
EMF <- aggregate(EMF_source~key, EMF, paste, collapse=", ")
EMF

insect <- df[c("key","Insect_type")]
insect <- unique(insect)
insect <- aggregate(Insect_type~key, insect, paste, collapse=", ")
insect

title <- df[c("key","Title")]
title <- na.omit(title)
head(title)

df2 <- full_join(insect, EMF)
df2 <- full_join(df2, title)
str(df2)

df <- df2
papers <- as.character(df$key)

for (i in 1:length(papers)) {
  string <- papers[i]
  string1 <- substr(string, 1, nchar(string)-4)
  df$Author[i] <- string1
  df$Year[i] <- str_extract_all(string,"\\(?[0-9,.]+\\)?")[[1]]}

df$Author <- gsub("1", "", df$Author)
df$Author <- gsub("2", "", df$Author)
df$Year <- as.numeric(df$Year)
str(df)

# reorder columns
myvars <- names(df) %in% c("Author","Year","Insect_type","EMF_source","Title")
df <- df[myvars]

df <- relocate(df, Insect_type, .after = c("Year"))
df <- relocate(df, EMF_source, .after = c("Insect_type"))
df <- relocate(df, Title, .after = c("EMF_source"))

str(df)

names(df)
colnames(df) <- c("Author","Year","Insect","Signal generator","Title")
colnames(df)

#View(df)

openxlsx::write.xlsx(df, "tables/Table_of_studies_HF.xlsx", rowNames = F, colWidths = "auto", overwrite = T, firstRow = T, firstCol = T)

# Latex print-out
df <- openxlsx::read.xlsx("tables/Table_of_studies_HF.xlsx", sheet = 1, na.strings = "NA") 

df <- arrange_at(df, vars(Year, Author), desc)


# for too long lists of insect or EMF types, change to "various"
df <- mutate(df, Insect = if_else ( nchar(Insect) >= 30, "Various", Insect))
df <- mutate(df, Signal.generator = if_else ( nchar(Signal.generator) >= 30, "Various", Signal.generator))

df$Title <- as.character(df$Title)
nchar(df$Title)

for (i in 1:length(df$Title)) {
  string<- df$Title[i]
  string1 <- substr(string, 1, 77)
  if (nchar(string1) > 76) {  df$Title[i] <- paste0(string1,"..")}
}

#View(df)

xhf <- xtable::xtable(df, caption = 'Included HF studies')
print(xhf, include.rownames = F)


############################################################
#  table of studies by category: LF
############################################################


df <- openxlsx::read.xlsx("tables/data_table_HFLF_1.xlsx", sheet = 1, na.strings = "NA") 

df <- df[grep("LF", df$EMF_type),]  # select only HF studies


EMF <- df[c("key","EMF_source")]
EMF <- unique(EMF)
EMF <- aggregate(EMF_source~key, EMF, paste, collapse=", ")
EMF

insect <- df[c("key","Insect_type")]
insect <- unique(insect)
insect <- aggregate(Insect_type~key, insect, paste, collapse=", ")
insect

title <- df[c("key","Title")]
title <- na.omit(title)
head(title)

df2 <- full_join(insect, EMF)
df2 <- full_join(df2, title)
str(df2)

df <- df2
papers <- as.character(df$key)

for (i in 1:length(papers)) {
  string <- papers[i]
  string1 <- substr(string, 1, nchar(string)-4)
  df$Author[i] <- string1
  df$Year[i] <- str_extract_all(string,"\\(?[0-9,.]+\\)?")[[1]]}

df$Author <- gsub("1", "", df$Author)
df$Author <- gsub("2", "", df$Author)
df$Year <- as.numeric(df$Year)
str(df)

# reorder columns
myvars <- names(df) %in% c("Author","Year","Insect_type","EMF_source","Title")
df <- df[myvars]

df <- relocate(df, Insect_type, .after = c("Year"))
df <- relocate(df, EMF_source, .after = c("Insect_type"))
df <- relocate(df, Title, .after = c("EMF_source"))

str(df)

colnames(df)
colnames(df) <- c("Author","Year","Insect","Signal generator","Title")

#View(df)

openxlsx::write.xlsx(df, "tables/Table_of_studies_LF.xlsx", rowNames = F, colWidths = "auto", overwrite = T, firstRow = T, firstCol = T)

# Latex print-out
df <- openxlsx::read.xlsx("tables/Table_of_studies_LF.xlsx", sheet = 1, na.strings = "NA") 

df <- arrange_at(df, vars(Year,Author), desc)

df <- mutate(df, Insect = if_else ( nchar(Insect) >= 30, "Various", Insect))
df <- mutate(df, Signal.generator = if_else ( nchar(Signal.generator) >= 30, "Various", Signal.generator))

df$Title <- as.character(df$Title)
nchar(df$Title)

for (i in 1:length(df$Title)) {
  string<- df$Title[i]
  string1 <- substr(string, 1, 77)
  if (nchar(string1) > 76) {  df$Title[i] <- paste0(string1,"..")}
}

#View(df)

xhf <- xtable::xtable(df, caption = 'Included LF studies')
# Latex code can be copy-pasted into latex editor
print(xhf, include.rownames = F)

