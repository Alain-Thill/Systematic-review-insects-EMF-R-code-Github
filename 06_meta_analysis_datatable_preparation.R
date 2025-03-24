
#run this first!

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
getwd() # working directory set to folder of .R file

#rm(list=ls())  # clear global environment of stuff from previous sessions

if (!require("pacman")) {
  install.packages("pacman")
}

if (!require("remotes")) {
  install.packages("remotes")
}
# run once to istall package "dmetar"
#remotes::install_github("MathiasHarrer/dmetar")

pacman::p_load(dplyr, tidyr, ggplot2, ggpubr, scales, RColorBrewer, openxlsx, sqldf, stringr, dmetar, 
               metafor, meta, bayesmeta, RoBMA, data.table)


##########################################################################
# Meta-Analysis preparation: deriving standard errors, collating tables
##########################################################################

#################################################################
# Standard error from mean and standard deviation
#################################################################

# for those studies that provided means and SD, or just as figures where data points could be manually extracted from with a graph digitizer
df <- openxlsx::read.xlsx("HF_means_sd.xlsx", sheet = 1, na.strings = "NA") 

#?metacont # deriving log ratio of means and log standard error from means + SD
meta_df <- metacont(n_e, mean_e, sdv_e, n_c, mean_c, sdv_c,
                    data = df, studlab = paste(key), sm = "ROM", cluster=key)

str(meta_df)

#meta_df$pval # precise p value
meta_df$TE # log ROM value
exp(meta_df$TE) # Ratio of means (ROM) value
#View(cbind(meta_df$studlab, meta_df$TE, exp(meta_df$TE), exp(meta_df$seTE)))

#df$p <- meta_df$pval
df$log_ROM <- meta_df$TE
df$log_SE <- meta_df$seTE

# Adding an "experiment_id" index
setDT(df)[, experiment_id := paste0(key, "_", rowid(key))]

#View(df)
openxlsx::write.xlsx(df, "tables/HF_se_from_means_sd.xlsx", rowNames = F, colWidths = 15, overwrite = T, firstRow = T, firstCol = T)

################################
# Collating tables
################################

df1 <- openxlsx::read.xlsx("tables/data_table_HFLF_2.xlsx", sheet = 1, na.strings = "NA") 
names(df1)
df1$Confid <- as.numeric(df1$Confid)

# swapping N_exposed and N_control columns to back
df1 <- relocate(df1, N_exposed, N_control, .after = last_col())

# Adding an "experiment_id" index to comfortably transfer new values to dataframes
setDT(df1)[, experiment_id := paste0(study, "_", rowid(study))]

df2 <- read.xlsx("tables/HF_se_from_means_sd.xlsx", sheet = 1, na.strings = "NA") 

#duplicated(df1$experiment_id) # to verify that index key is unique
#duplicated(df2$experiment_id) # to verify that index key is unique

df2 <- df2[c("n_e","n_c","mean_e","mean_c","sdv_e","sdv_c","log_ROM","log_SE","experiment_id")]

# A left_join(x, y) keeps all observations in x.
# A right_join(x, y) keeps all observations in y.
# Here, we want to join the smaller data table (with effect size and SE values derived from means and SD) to the full-length data table.
# To keep the more precise values from the smaller table, concerned columns need to have different names or else the old values are kept.

df <- left_join(df1, df2, by="experiment_id")

remove_col <- c("Output_power","Organs","Stats_method","RFLow","RFHigh")
df <- df %>% select(!all_of(remove_col))
names(df)

#View(df)
#View(cbind(df$experiment_id, df$log_ROM, df$log_SE, df$Confid))

# Standardizing p and n
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

# Conversion of the relevant variables given in studies as a percentage deviation from the control
# into the respective multiplier that exists between irradiated and control
# e.g. increase by 100% (i.e. from 100 to 200) = doubling = factor 2
# or reduction by 50% (i.e. from 100 to 50) = halving = factor 0.5

# filling in missing log_ROM values from percent changes
df <- mutate(df, log_ROM = if_else ( is.na(log_ROM), log(df$Percent.change.vs.control + 1), log_ROM))

# backtransform ROM values from log_ROM derived with "metacont" from extracted plot values
df <- mutate(df, Ratio_of_means = exp(log_ROM))
df <- relocate(df, Ratio_of_means, .after = "Percent.change.vs.control")

# backtransform percent changes from more precise log_ROM values derived with "metacont" from extracted plot values
df <- mutate(df, Percent.change.vs.control = exp(log_ROM)-1)
df$Percent.change.vs.control <- percent(df$Percent.change.vs.control)

openxlsx::write.xlsx(df, "tables/data_table_HFLF_3.xlsx", rowNames = F, colWidths = 15, overwrite = T, firstRow = T, firstCol = T)


#############################################
# Determining standard error from p-values
#############################################

df1 <- openxlsx::read.xlsx("tables/data_table_HFLF_3.xlsx", sheet = 1, na.strings = "NA") 

#View(df1)
names(df1)
str(df1$p)

# assigning p = 0.5 for experiments declaring "no statistical significance" and not giving precise p value
# this is done since if only those experiments that mention an exact p value < 0.05 are included,
# but all experiments not giving p value are excluded, this increases bias
# assigning an imprecise p value to calculate standard errors has the disadvantage that it distorts meta-analysis funnel plots
df2 <- mutate(df1, p = if_else ( is.na(p) & is.na(log_SE), 0.5, p))

df2 <- df2[complete.cases(df2[,c("Ratio_of_means","p","n_e")]),] # se.from.p needs an effect size, p and number of samples

# calculate total samples per experiment
df2 <- mutate(df2, N_total = if_else ( is.na(n_c), n_e, n_e+n_c))

#View(df2)

#?se.from.p
df3 <- se.from.p(effect.size = df2$Ratio_of_means, p = df2$p, N = df2$N_total,
                 effect.size.type= "ratio", calculate.g = F)

df4 <- cbind(df2, df3)
#View(df4)
names(df4)

df4 <- df4[c("logEffectSize","logStandardError","experiment_id","p")]  # grabbing the columns with newly computed values from se.from.p

df <- left_join(df1, df4, by="experiment_id")

names(df)

# filling in log_ROM and log_SE derived from p, if none already present, so the most precise value is conserved
df <- mutate(df, log_ROM = if_else ( is.na(log_ROM), logEffectSize, log_ROM))
df <- mutate(df, log_SE = if_else ( is.na(log_SE), logStandardError, log_SE))
df <- mutate(df, p.x = if_else ( is.na(p.x), p.y, p.x)) # to indicate about 20 "no effect" experiments where p was assigned as 0.5


View(cbind(df$experiment_id, df$log_ROM, df$log_SE, df$logStandardError, df$p.x))
# inspecting the results reveals the standard errors derived from p values are often bigger than those derived from means + SD


# manually assigning similar log_SE to Miyan2014 no effect experiment since for an exactly zero effect size (or 1 as ROM), 
# standard error also is 0, leading to zero contribution in the meta-analysis, which seems unbalanced
df <- mutate(df, log_ROM = if_else ( experiment_id == 'Miyan2014_3', 0.013, log_ROM))
df <- mutate(df, log_SE = if_else ( experiment_id == 'Miyan2014_3', 0.02, log_SE))


names(df)
df <- subset(df, select = -c(logEffectSize, logStandardError, p.y))
# renaming "p.x" column into "p"
names(df)[names(df) == 'p.x'] <- 'p'


# list studies with p = 0.5
pval_0.5 <- sqldf("SELECT * FROM df WHERE p == 0.5")
str(pval_0.5$study)

# list studies with p = 0.05
pval_0.05 <- sqldf("SELECT * FROM df WHERE p == 0.05")
str(pval_0.05$study)


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

# Base station:

# making subgroup for the 4 different experiments in Vijver2014
df <- df %>% mutate(df, experiment = if_else ( Insect_type == "Springtail" & study == "Vijver2014", "Vijver2014_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( Insect_type == "Beetle" & study == "Vijver2014", "Vijver2014_2", experiment))
df <- df %>% mutate(df, experiment = if_else ( Insect_type == "Wasp" & study == "Vijver2014", "Vijver2014_3", experiment))
df <- df %>% mutate(df, experiment = if_else ( Insect_type == "Drosophila" & study == "Vijver2014", "Vijver2014_4", experiment))

# DECT:
# making subgroups for the 3 different experiments in Manta2013
df <- df %>% mutate(df, experiment = if_else ( X == "Male Flies" & study == "Manta2013", "Manta2013_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( X == "Female Flies Bodies" & study == "Manta2013", "Manta2013_2", experiment))
df <- df %>% mutate(df, experiment = if_else ( X == "Female Flies Ovaries" & study == "Manta2013", "Manta2013_3", experiment))

# Mobile phone:
# making subgroup for the 3 different experiments in Chavdoula2010
df <- df %>% mutate(df, experiment = if_else ( Bioeffects == "DNA Fragmentation" & study == "Chavdoula2010", "Chavdoula2010_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( Bioeffects == "Developmental Effects (actin-cytoskeleton disorganization)" & study == "Chavdoula2010", "Chavdoula2010_2", experiment))
df <- df %>% mutate(df, experiment = if_else ( Bioeffects == "Reduced Reproductive Capacity" & study == "Chavdoula2010", "Chavdoula2010_3", experiment))

# making subgroups for the 2 different experiments in Panagopoulos2007a
df <- df %>% mutate(df, experiment = if_else ( SRF == "900" & study == "Panagopoulos2007a", "Panagopoulos2007a_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( SRF == "1800" & study == "Panagopoulos2007a", "Panagopoulos2007a_2", experiment))

# making subgroups for the 2 different experiments in Panagopoulos2007b
df <- df %>% mutate(df, experiment = if_else ( SRF == "900" & study == "Panagopoulos2007b", "Panagopoulos2007b_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( SRF == "1800" & study == "Panagopoulos2007b", "Panagopoulos2007b_2", experiment))

# making subgroups for the 2 different experiments in Panagopoulos2010a
df <- df %>% mutate(df, experiment = if_else ( SRF == "900" & study == "Panagopoulos2010a", "Panagopoulos2010a_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( SRF == "1800" & study == "Panagopoulos2010a", "Panagopoulos2010a_2", experiment))

# making subgroups for the 2 different experiments in Panagopoulos2010c
df <- df %>% mutate(df, experiment = if_else ( SRF == "900" & study == "Panagopoulos2010c", "Panagopoulos2010c_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( SRF == "1800" & study == "Panagopoulos2010c", "Panagopoulos2010c_2", experiment))

# making subgroups for the 2 + 3 (Signal gen) different experiments in Vilic2017
#df <- df %>% mutate(df, experiment = if_else ( study == "Vilic2017", experiment_id, experiment))


# Signal generator:
# making subgroups for the 3 different experiments in Gunes2021
df <- df %>% mutate(df, experiment = if_else ( SRF == "900" & study == "Gunes2021", "Gunes2021_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( SRF == "1800" & study == "Gunes2021", "Gunes2021_2", experiment))
df <- df %>% mutate(df, experiment = if_else ( SRF == "2100" & study == "Gunes2021", "Gunes2021_3", experiment))

# making subgroups for the 4 different experiments in Sagioglou2014
df <- df %>% mutate(df, experiment = if_else ( SRF == "100" & study == "Sagioglou2014", "Sagioglou2014_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( SRF == "395" & study == "Sagioglou2014", "Sagioglou2014_2", experiment))
df <- df %>% mutate(df, experiment = if_else ( SRF == "682" & study == "Sagioglou2014", "Sagioglou2014_3", experiment))
df <- df %>% mutate(df, experiment = if_else ( SRF == "900" & study == "Sagioglou2014", "Sagioglou2014_4", experiment))

# making subgroups for the 3 different experiments in Wang2021
df <- df %>% mutate(df, experiment = if_else ( Bioeffects == "Increased level of activity at night, Reduced total sleep duration" & study == "Wang2021", "Wang2021_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( Bioeffects == "Reduced number of movements" & study == "Wang2021", "Wang2021_2", experiment))
df <- df %>% mutate(df, experiment = if_else ( Bioeffects == "Reduced total level of activity" & study == "Wang2021", "Wang2021_3", experiment))


# Coil system: (all studies have multiple independent groups)
# misc

# making subgroup for the 4 different experiments in Vilic2024
df <- df %>% mutate(df, experiment = if_else ( Bioeffects == "Glutathione S-transferase (GST) activity" & study == "Vilic2024", "Vilic2024_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( Bioeffects == "Catalase (CAT) activity" & study == "Vilic2024", "Vilic2024_2", experiment))
df <- df %>% mutate(df, experiment = if_else ( Bioeffects == "Superoxide dismutase (SOD) activity" & study == "Vilic2024", "Vilic2024_3", experiment))
df <- df %>% mutate(df, experiment = if_else ( Bioeffects == "Thiobarbituric acid reactive substance (TBARS) concentration" & study == "Vilic2024", "Vilic2024_4", experiment))


# making subgroup for the 3 different experiments in Thill2018
df <- df %>% mutate(df, experiment = if_else ( Insect_type == "Honeybee" & study == "Thill2018", "Thill2018_1", experiment))
df <- df %>% mutate(df, experiment = if_else ( Insect_type == "Bumblebee" & study == "Thill2018", "Thill2018_2", experiment))
df <- df %>% mutate(df, experiment = if_else ( EMF_source == "Base station" & study == "Thill2018", "Thill2018_3", experiment))
df <- df %>% mutate(df, experiment = if_else ( Bioeffects == "Increased mortality compared to control group" & study == "Thill2018", "Thill2018_4", experiment))

df$experiment
#View(df)

# save expanded table
openxlsx::write.xlsx(df, "tables/data_table_HFLF_4.xlsx", rowNames = F, colWidths = 15, overwrite = T, firstRow = T, firstCol = T)

