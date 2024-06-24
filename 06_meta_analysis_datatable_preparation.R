

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
getwd() # working directory set to folder of .R file

#rm(list=ls())  # clear global environment of stuff from previous sessions


##########################################################################
# Meta-Analysis preparation: deriving standard errors, collating tables
##########################################################################

pacman::p_load(dplyr, tidyr, ggplot2, ggpubr, scales, RColorBrewer, openxlsx, sqldf, stringr, dmetar, 
               metafor, meta, bayesmeta, RoBMA, data.table)

#################################################################
# Standard error from mean and standard deviation
#################################################################

# for those studies that provided means and SD, or just as figures where data points could be manually extracted from with a graph digitizer
df <- openxlsx::read.xlsx("HF_means_sd.xlsx", sheet = 1, na.strings = "NA") 

#?metacont # deriving log ratio of means and log standard error from means + SD
meta_df <- metacont(n_e, mean_e, sdv_e, n_c, mean_c, sdv_c,
                    data = df, studlab = paste(key), sm = "ROM")

str(meta_df)

meta_df$pval # precise p value
meta_df$TE # log ROM value
exp(meta_df$TE) # Ratio of means (ROM) value

df$p <- meta_df$pval
df$log_ROM <- meta_df$TE
df$log_SE <- meta_df$seTE

# backtransform percent changes from more precise ROM values derived with "metacont" from extracted plot values
df <- mutate(df, Percent.change.vs.control = exp(log_ROM)-1)
df$Percent.change.vs.control

# Adding an "experiment_id" index
setDT(df)[, experiment_id := paste0(key, "_", rowid(key))]

#View(df)
openxlsx::write.xlsx(df, "tables/HF_se_from_means_sd.xlsx", rowNames = F, colWidths = 15, overwrite = T, firstRow = T, firstCol = T)

################################
# Collating tables
################################

df1 <- openxlsx::read.xlsx("tables/data_table_HFLF_2.xlsx", sheet = 1, na.strings = "NA") 
names(df1)

# swapping N_exposed and N_control columns to back
df1 <- relocate(df1, N_exposed, N_control, .after = last_col())

# Adding an "experiment_id" index to comfortably transfer new values to dataframes
setDT(df1)[, experiment_id := paste0(study, "_", rowid(study))]

df2 <- read.xlsx("tables/HF_se_from_means_sd.xlsx", sheet = 1, na.strings = "NA") 

#duplicated(df1$experiment_id) # to verify that index key is unique
#duplicated(df2$experiment_id) # to verify that index key is unique

df2 <- df2[c("p","n_e","n_c","mean_e","mean_c","sdv_e","sdv_c","log_ROM","log_SE","experiment_id")]

# A left_join(x, y) keeps all observations in x.
# A right_join(x, y) keeps all observations in y.
# Here, we want to join the smaller data table (with effect size and SE values derived from means and SD) to the full-length data table.
# To keep the more precise values from the smaller table, concerned columns need to have different names or else the old values are kept.

df <- left_join(df1, df2, by="experiment_id")

remove_col <- c("Output_power","Organs","Stats_method","RFLow","RFHigh")
df <- df %>% select(!all_of(remove_col))
names(df)
#View(df)

# Standardizing p and n

df <- mutate(df, Confid = if_else ( is.na(Confid), p, Confid)) # if Confid NA, input p value derived from means + SD
df$p <- df$Confid # replace confid with p
df <- df %>% select(!Confid) # get rid of confid
names(df)

# remove sample type designing words (hive / adult / larva) from sample size columns of exposed and control group
sample_size <- df$N_exposed
for (i in 1:length(sample_size)) {
  string <- sample_size[i]
  df$n_e[i] <- str_extract_all(string,"\\(?[0-9,.]+\\)?")[[1]]}

sample_size <- df$N_control
for (i in 1:length(sample_size)) {
  string <- sample_size[i]
  df$n_c[i] <- str_extract_all(string,"\\(?[0-9,.]+\\)?")[[1]]}

df$n_e <- as.numeric(df$n_e)
df$n_c <- as.numeric(df$n_c)


# Conversion of the relevant variables given in studies as a percentage deviation to the control
# into the respective multiplier that exists between irradiated and control
# e.g. increase by 100% (i.e. from 100 to 200) = doubling = factor 2
# or reduction by 50% (i.e. from 100 to 50) = halving = factor 0.5

#df$Ratio_of_means <- df$Percent.change.vs.control + 1  # this is sufficient if percent are given as an integer, e.g. +100% = 200% = 2, +0% = 100% = 1, -50% = 100-50% = 0.5

# deriving precise ratio of means values from the means + SD values
df$Ratio_of_means <- exp(df$log_ROM)

# for experiments lacking log_ROM derived from means + SD, get ratio of means from percent change
df <- mutate(df, Ratio_of_means = if_else ( is.na(Ratio_of_means), df$Percent.change.vs.control + 1, Ratio_of_means))

df <- relocate(df, Ratio_of_means, .after = "Percent.change.vs.control")

openxlsx::write.xlsx(df, "tables/data_table_HFLF_3.xlsx", rowNames = F, colWidths = 15, overwrite = T, firstRow = T, firstCol = T)


#############################################
# Determining standard error from p-values
#############################################

df1 <- openxlsx::read.xlsx("tables/data_table_HFLF_3.xlsx", sheet = 1, na.strings = "NA") 

names(df1)


# assigning p = 0.5 for experiments declaring "no statistical significance" and not giving precise p value
# this is done since if only those experiments that mention an exact p value < 0.05 are included,
# but all experiments not giving p value are excluded, this increases bias
df2 <- mutate(df1, p = if_else ( is.na(p), 0.5, p))

df2 <- df2[complete.cases(df2[,c("Ratio_of_means","p","n_e")]),] # se.from.p needs an effect size, p and n

#?se.from.p
df3 <- se.from.p(effect.size = df2$Ratio_of_means, p = df2$p, N = df2$n_e,
                 effect.size.type= "ratio", calculate.g = F)


df4 <- cbind(df2, df3)
names(df4)

df4 <- df4[c("logEffectSize","logStandardError","experiment_id","p")]  # grabbing the columns with newly computed values from se.from.p

df <- left_join(df1, df4, by="experiment_id")
#View(df)  # inspecting the results (scroll to right side) reveals the standard errors derived from p values are bigger than those derived from means + SD
names(df)

# # filling in log_ROM and log_SE derived from p, if none already present, so the most precise value is conserved
df <- mutate(df, log_ROM = if_else ( is.na(log_ROM), logEffectSize, log_ROM))
df <- mutate(df, log_SE = if_else ( is.na(log_SE), logStandardError, log_SE))
df <- mutate(df, p.x = if_else ( is.na(p.x), p.y, p.x)) # to indicate about 20 "no effect" experiments where p was assigned as 0.5

#View(df)

# manually assigning similar log_SE to Myan2015 no effect experiment since for an exactly zero effect size (or 1 as ROM), 
# standard error also is 0, leading to zero contribution in the meta-analysis, which seems unbalanced
df <- mutate(df, log_SE = if_else ( log_ROM == 0, 0.02, log_SE))

names(df)
df <- subset(df, select = -c(logEffectSize, logStandardError, p.y))
# renaming "p.x" column into "p"
names(df)[names(df) == 'p.x'] <- 'p'


# list studies with p = 0.5
pval_0.5 <- sqldf("SELECT * FROM df WHERE p == 0.5")
pval_0.5$study

# list studies with p = 0.05
pval_0.05 <- sqldf("SELECT * FROM df WHERE p == 0.05")
pval_0.05$study


#############################################
# Additional formatting
#############################################

# nicer Author Year label
for (i in 1:length(df$study)) {
  string <- df$study[i]
  df$Author[i] <- unlist(str_split(string, "\\(?[0-9]+\\)?"))[1]
  df$an_index[i] <- unlist(str_split(string, "\\(?[0-9]+\\)?"))[2]
  df$Year[i] <- str_extract_all(string,"\\(?[0-9,.]+\\)?")[[1]]
  df$Year[i] <- paste0(df$Year[i], df$an_index[i])
  df$Author_Year[i] <- paste0(df$Author[i][[1]]," ",df$Year[i])}

df <- subset(df, select = -c(Author, an_index, Year))


# making subgroups based on experimental groups sharing the same control group
df$experiment <- df$study

# making subgroup for the 2 different experiments in Panagopoulos2010a
df <- df %>% mutate(df, experiment = if_else ( SRF == "900" & study == "Panagopoulos2010a", "Panagopoulos2010a_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( SRF == "1800" & study == "Panagopoulos2010a", "Panagopoulos2010a_2", experiment))

# making subgroup for the 2 different experiments in Panagopoulos2010c
df <- df %>% mutate(df, experiment = if_else ( SRF == "900" & study == "Panagopoulos2010c", "Panagopoulos2010c_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( SRF == "1800" & study == "Panagopoulos2010c", "Panagopoulos2010c_2", experiment))

# making subgroup for the 4 different experiments in Sagioglou2014
df <- df %>% mutate(df, experiment = if_else ( SRF == "100" & study == "Sagioglou2014", "Sagioglou2014_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( SRF == "395" & study == "Sagioglou2014", "Sagioglou2014_2", experiment))
df <- df %>% mutate(df, experiment = if_else ( SRF == "682" & study == "Sagioglou2014", "Sagioglou2014_3", experiment))
df <- df %>% mutate(df, experiment = if_else ( SRF == "900" & study == "Sagioglou2014", "Sagioglou2014_4", experiment))

# making subgroup for the 3 different experiments in Manta2013
df <- df %>% mutate(df, experiment = if_else ( X == "Male Flies" & study == "Manta2013", "Manta2013_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( X == "Female Flies Bodies" & study == "Manta2013", "Manta2013_2", experiment))
df <- df %>% mutate(df, experiment = if_else ( X == "Female Flies Ovaries" & study == "Manta2013", "Manta2013_3", experiment))

# making subgroup for the 4 different experiments in Vijver2014
df <- df %>% mutate(df, experiment = if_else ( Insect_type == "Springtail" & study == "Vijver2014", "Vijver2014_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( Insect_type == "Beetle" & study == "Vijver2014", "Vijver2014_2", experiment))
df <- df %>% mutate(df, experiment = if_else ( Insect_type == "Wasp" & study == "Vijver2014", "Vijver2014_3", experiment))
df <- df %>% mutate(df, experiment = if_else ( Insect_type == "Drosophila" & study == "Vijver2014", "Vijver2014_4", experiment))

# making subgroup for the 4 different experiments in Vilic2024
df <- df %>% mutate(df, experiment = if_else ( Bioeffects == "Glutathione S-transferase (GST) activity" & study == "Vilic2024", "Vilic2024_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( Bioeffects == "Catalase (CAT) activity" & study == "Vilic2024", "Vilic2024_2", experiment))
df <- df %>% mutate(df, experiment = if_else ( Bioeffects == "Superoxide dismutase (SOD) activity" & study == "Vilic2024", "Vilic2024_3", experiment))
df <- df %>% mutate(df, experiment = if_else ( Bioeffects == "Thiobarbituric acid reactive substance (TBARS) concentration" & study == "Vilic2024", "Vilic2024_4", experiment))

# making subgroup for the 3 different experiments in Thill2018
df <- df %>% mutate(df, experiment = if_else ( Insect_type == "Honeybee" & study == "Thill2018", "Thill2018_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( EMF_source == "Base station" & study == "Thill2018", "Thill2018_2", experiment))
df <- df %>% mutate(df, experiment = if_else ( Bioeffects == "Increased mortality compared to control group" & study == "Thill2018", "Thill2018_3", experiment))

df$experiment

# save expanded table
openxlsx::write.xlsx(df, "tables/data_table_HFLF_4.xlsx", rowNames = F, colWidths = 15, overwrite = T, firstRow = T, firstCol = T)

