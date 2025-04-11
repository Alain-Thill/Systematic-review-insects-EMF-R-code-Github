
# run this first!
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(dplyr, tidyr, gdata, ggplot2, ggpubr, scales, RColorBrewer, openxlsx, sqldf, data.table)


################################################
# Create "collapsed" summary table
################################################

df <- openxlsx::read.xlsx("tables/HFLF_meta_table.xlsx", sheet = 1, na.strings = "NA") 

names(df)

HF <- df[grep("HF", df$type),]
LF <- df[grep("LF", df$type),]
names(HF)

df_start <- df[c("study","Title","Effect","EMF_type")]
df_start <- na.omit(df_start)
#View(df_start)

df_EMF <- unique(df[c("study","EMF_source")])
df_EMF <- aggregate(EMF_source ~ study, df_EMF, paste, collapse=", ")
df_EMF
df_insect <- unique(df[c("study","Insect_type")])
df_insect <- aggregate(Insect_type ~ study, df_insect, paste, collapse=", ")
df_insect
df_Efield <- unique(df[c("study","E_Field")])
df_Efield <- arrange_all(df_Efield)
df_Efield <- aggregate(E_Field ~ study, df_Efield, paste, collapse=", ")
#View(df_Efield)

min_Efield <- df_Efield %>% mutate(E_Field = strsplit(E_Field, ", ")) 

min_Efield$min_E <- lapply(min_Efield$E_Field, `[[`, 1)
min_Efield$min_E <- as.numeric(min_Efield$min_E)
min_Efield <- min_Efield[,c(1,3)]
names(min_Efield) <- c("study","min_E_Field")
str(min_Efield)
#View(min_Efield)

df_bioeffects <- unique(df[c("study","Bioeffect_cat")])
df_bioeffects <- aggregate(Bioeffect_cat ~ study, df_bioeffects, paste, collapse=", ")
#View(df_bioeffects)

df_direct_effect <- unique(df[c("study","Direction_of_effect")])
df_direct_effect <- aggregate(Direction_of_effect ~ study, df_direct_effect, paste, collapse=", ")

df_merge <- full_join(df_start, df_EMF, by = c("study"), copy=T)
df_merge <- full_join(df_merge, df_insect, by = c("study"), copy=T)
df_merge <- full_join(df_merge, min_Efield, by = c("study"), copy=T)
df_merge <- full_join(df_merge, df_bioeffects, by = c("study"), copy=T)
df_merge <- full_join(df_merge, df_direct_effect, by = c("study"), copy=T)

df_merge <- df_merge %>% mutate_at("Direction_of_effect", ~replace_na(.,"none"))
str(df_merge)
View(df_merge)

openxlsx::write.xlsx(df_merge, "tables/Bioeffects_HFLF_aggregated.xlsx", rowNames = F, colWidths = "auto", overwrite = T, firstRow = T, firstCol = T)


########################################################
# Table median of all variables grouped by EMF source
########################################################

df <- openxlsx::read.xlsx("tables/HFLF_meta_table.xlsx", sheet = 1, na.strings = "NA") 

str(df)
View(df)

df <- df[grep("HF", df$EMF_type),]

table1 <- df %>%
  group_by(EMF_source) %>% 
  summarise(
    across(where(is.numeric), ~ median(.x, na.rm = TRUE)), 
    across(where(is.factor), nlevels),
    n = n(), 
  )
table1
View(table1)


table2 <- df %>%
  group_by(study, EMF_source, Bioeffect_cat) %>% 
  summarise(
    across(where(is.numeric), ~ median(.x, na.rm = TRUE)), 
    across(where(is.factor), nlevels),
    n = n(), 
  )

View(table2)



########################################################################################################
# code for selecting only the experimental group with strongest measured effect among all exp. groups
########################################################################################################

DECT <- sqldf("SELECT * FROM df WHERE EMF_source LIKE '%DECT%'")
DECT <- as.data.table(DECT)
DECT2 <- DECT[DECT[, .I[log_ROM == max(log_ROM)], by=experiment]$V1]
View(DECT2)
length(DECT$study)
