

# run this first!
pacman::p_load(dplyr, tidyr, openxlsx, sqldf, data.table, meta, metafor, poolr, gridExtra, metap)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

#rm(list=ls())  #clear global environment of stuff from previous sessions


#########################################################################################
#  Table of probabilities for n "Effect" experiments among p total experiments: counting
#########################################################################################

df <- openxlsx::read.xlsx("tables/HFLF_meta_table.xlsx", sheet = 1, na.strings = "NA") 

df$Exp_Effect <- "" # adding "effect" variable for all experiments within studies 
df$Exp_Effect <- as.character(df$Exp_Effect)
df <- mutate(df, Exp_Effect = if_else ( Direction_of_effect== 'detrimental', "Effect", Exp_Effect)) # If harmful, classify as "Effect"
df <- mutate(df, Exp_Effect = if_else ( Direction_of_effect== 'beneficial', "No effect", Exp_Effect)) # If positive effect, classify as "No effect"
df <- mutate(df, Exp_Effect = if_else ( Direction_of_effect== 'none', "No effect", Exp_Effect)) # If no effect, classify as "No effect"
df <- mutate(df, Exp_Effect = if_else ( Direction_of_effect== 'uncertain', "Unknown", Exp_Effect)) # If uncertain effect, classify as "NA", for removal
df$Exp_Effect[grep("Unknown", df$Exp_Effect)] <- NA

df <- df[complete.cases(df[,"Exp_Effect"]),] # remove "uncertain" experiments

HF <- df[grep("HF", df$EMF_type),]
LF <- df[grep("LF", df$EMF_type),]

table(df$Direction_of_effect)
table(df$Exp_Effect)

# Building a table

`EMF type` <- c("All","LF","LF","LF","HF","HF","HF","HF","HF","HF","HF","Base station","DECT","Helmholtz coil","Mobile phone","Power line","Signal generator")
`E-field strength` <- c("All","All","All","All","All","All","All","< 41 V/m","< 6 V/m","< 6 V/m","< 6 V/m",
                        "All","All","All","All","All","All")
`Insect` <- c("All","All","Drosophila","Honeybee","All","Drosophila","Honeybee","All","All","Drosophila","Honeybee","All","All","All","All","All","All")
odds_exp <- tibble(`EMF type`,`E-field strength`,`Insect`)

odds_exp
odds_exp$`Number of experiments` <- 0
odds_exp$`Effect` <- 0
odds_exp$`No Effect` <- 0

#nrow(subset(HF, Direction_of_effect == "detrimental"))
#nrow(HF[grep("detrimental", HF$Direction_of_effect),]) # same as line above


# All
nrow(df)
nrow(sqldf("SELECT * FROM df WHERE Exp_Effect = 'Effect'"))
nrow(sqldf("SELECT * FROM df WHERE Exp_Effect = 'No effect'"))
odds_exp$`Number of experiments`[1] <- nrow(df)
odds_exp$`Effect`[1] <- nrow(sqldf("SELECT * FROM df WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[1] <- nrow(sqldf("SELECT * FROM df WHERE Exp_Effect = 'No effect'"))

# LF
nrow(sqldf("SELECT * FROM LF WHERE Exp_Effect = 'Effect'"))
nrow(sqldf("SELECT * FROM LF WHERE Exp_Effect = 'No effect'"))
odds_exp$`Number of experiments`[2] <- nrow(LF)
odds_exp$`Effect`[2] <- nrow(sqldf("SELECT * FROM LF WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[2] <- nrow(sqldf("SELECT * FROM LF WHERE Exp_Effect = 'No effect'"))


drosophila <- sqldf("SELECT * FROM LF WHERE Insect_type = 'Drosophila'")

nrow(sqldf("SELECT * FROM drosophila WHERE Exp_Effect = 'Effect'"))
nrow(sqldf("SELECT * FROM drosophila WHERE Exp_Effect = 'No effect'"))
odds_exp$`Number of experiments`[3] <- nrow(drosophila)
odds_exp$`Effect`[3] <- nrow(sqldf("SELECT * FROM drosophila WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[3] <- nrow(sqldf("SELECT * FROM drosophila WHERE Exp_Effect = 'No effect'"))


bees <- sqldf("SELECT * FROM LF WHERE Insect_type = 'Honeybee'")

nrow(sqldf("SELECT * FROM bees WHERE Exp_Effect = 'Effect'"))
nrow(sqldf("SELECT * FROM bees WHERE Exp_Effect = 'No effect'"))
odds_exp$`Number of experiments`[4] <- nrow(bees)
odds_exp$`Effect`[4] <- nrow(sqldf("SELECT * FROM bees WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[4] <- nrow(sqldf("SELECT * FROM bees WHERE Exp_Effect = 'No effect'"))

# HF
nrow(sqldf("SELECT * FROM HF WHERE Exp_Effect = 'Effect'"))
nrow(sqldf("SELECT * FROM HF WHERE Exp_Effect = 'No effect'"))
odds_exp$`Number of experiments`[5] <- nrow(HF)
odds_exp$`Effect`[5] <- nrow(sqldf("SELECT * FROM HF WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[5] <- nrow(sqldf("SELECT * FROM HF WHERE Exp_Effect = 'No effect'"))


drosophila <- sqldf("SELECT * FROM HF WHERE Insect_type = 'Drosophila'")

nrow(sqldf("SELECT * FROM drosophila WHERE Exp_Effect = 'Effect'"))
nrow(sqldf("SELECT * FROM drosophila WHERE Exp_Effect = 'No effect'"))
odds_exp$`Number of experiments`[6] <- nrow(drosophila)
odds_exp$`Effect`[6] <- nrow(sqldf("SELECT * FROM drosophila WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[6] <- nrow(sqldf("SELECT * FROM drosophila WHERE Exp_Effect = 'No effect'"))


drosophila %>% 
  group_by(EMF_source) %>% 
  filter(n() > 1) %>% 
  summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)), n = n())


bees <- sqldf("SELECT * FROM HF WHERE Insect_type = 'Honeybee'")

nrow(sqldf("SELECT * FROM bees WHERE Exp_Effect = 'Effect'"))
nrow(sqldf("SELECT * FROM bees WHERE Exp_Effect = 'No effect'"))
odds_exp$`Number of experiments`[7] <- nrow(bees)
odds_exp$`Effect`[7] <- nrow(sqldf("SELECT * FROM bees WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[7] <- nrow(sqldf("SELECT * FROM bees WHERE Exp_Effect = 'No effect'"))


icnirp <- sqldf("SELECT * FROM HF WHERE E_Field <= 41")
#View(icnirp)

nrow(sqldf("SELECT * FROM icnirp WHERE Exp_Effect = 'Effect'"))
nrow(sqldf("SELECT * FROM icnirp WHERE Exp_Effect = 'No effect'"))
odds_exp$`Number of experiments`[8] <- nrow(icnirp)
odds_exp$`Effect`[8] <- nrow(sqldf("SELECT * FROM icnirp WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[8] <- nrow(sqldf("SELECT * FROM icnirp WHERE Exp_Effect = 'No effect'"))


count <- sqldf("SELECT * FROM HF WHERE E_Field <= 6")

odds_exp$`Number of experiments`[9] <- nrow(count)
odds_exp$`Effect`[9] <- nrow(sqldf("SELECT * FROM count WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[9] <- nrow(sqldf("SELECT * FROM count WHERE Exp_Effect = 'No effect'"))

drosophila <- sqldf("SELECT * FROM count WHERE Insect_type = 'Drosophila'")
odds_exp$`Number of experiments`[10] <- nrow(drosophila)
odds_exp$`Effect`[10] <- nrow(sqldf("SELECT * FROM drosophila WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[10] <- nrow(sqldf("SELECT * FROM drosophila WHERE Exp_Effect = 'No effect'"))

bees <- sqldf("SELECT * FROM count WHERE Insect_type = 'Honeybee'")
odds_exp$`Number of experiments`[11] <- nrow(bees)
odds_exp$`Effect`[11] <- nrow(sqldf("SELECT * FROM bees WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[11] <- nrow(sqldf("SELECT * FROM bees WHERE Exp_Effect = 'No effect'"))

table(df$EMF_source)

# by EMF_source
base <- sqldf("SELECT * FROM HF WHERE EMF_source = 'Base station'")

nrow(sqldf("SELECT * FROM base WHERE Exp_Effect = 'Effect'"))
nrow(sqldf("SELECT * FROM base WHERE Exp_Effect = 'No effect'"))
odds_exp$`Number of experiments`[12] <- nrow(base)
odds_exp$`Effect`[12] <- nrow(sqldf("SELECT * FROM base WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[12] <- nrow(sqldf("SELECT * FROM base WHERE Exp_Effect = 'No effect'"))

DECT <- sqldf("SELECT * FROM df WHERE EMF_source = 'DECT'")
odds_exp$`Number of experiments`[13] <- nrow(DECT)
odds_exp$`Effect`[13] <- nrow(sqldf("SELECT * FROM DECT WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[13] <- nrow(sqldf("SELECT * FROM DECT WHERE Exp_Effect = 'No effect'"))

helmholtz <- sqldf("SELECT * FROM df WHERE EMF_source = 'Coil system 50/60 Hz'")
nrow(helmholtz) # 93
nrow(sqldf("SELECT * FROM helmholtz WHERE Exp_Effect = 'Effect'")) # 25
nrow(sqldf("SELECT * FROM helmholtz WHERE Exp_Effect = 'No effect'")) # 15
odds_exp$`Number of experiments`[14] <- nrow(helmholtz)
odds_exp$`Effect`[14] <- nrow(sqldf("SELECT * FROM helmholtz WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[14] <- nrow(sqldf("SELECT * FROM helmholtz WHERE Exp_Effect = 'No effect'"))

mobile <- sqldf("SELECT * FROM df WHERE EMF_source = 'Mobile phone'")
nrow(mobile) # 107
nrow(sqldf("SELECT * FROM mobile WHERE Exp_Effect = 'Effect'")) # 50
nrow(sqldf("SELECT * FROM mobile WHERE Exp_Effect = 'No effect'")) # 8
odds_exp$`Number of experiments`[15] <- nrow(mobile)
odds_exp$`Effect`[15] <- nrow(sqldf("SELECT * FROM mobile WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[15] <- nrow(sqldf("SELECT * FROM mobile WHERE Exp_Effect = 'No effect'"))

power <- sqldf("SELECT * FROM df WHERE EMF_source = 'Power line'")
odds_exp$`Number of experiments`[16] <- nrow(power)
odds_exp$`Effect`[16] <- nrow(sqldf("SELECT * FROM power WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[16] <- nrow(sqldf("SELECT * FROM power WHERE Exp_Effect = 'No effect'"))

siggen <- sqldf("SELECT * FROM df WHERE EMF_source = 'RF signal generator'")
nrow(siggen) # 68
nrow(sqldf("SELECT * FROM siggen WHERE Exp_Effect = 'Effect'")) # 27
nrow(sqldf("SELECT * FROM siggen WHERE Exp_Effect = 'No effect'")) # 11
odds_exp$`Number of experiments`[17] <- nrow(siggen)
odds_exp$`Effect`[17] <- nrow(sqldf("SELECT * FROM siggen WHERE Exp_Effect = 'Effect'"))
odds_exp$`No Effect`[17] <- nrow(sqldf("SELECT * FROM siggen WHERE Exp_Effect = 'No effect'"))

#View(odds_exp)
openxlsx::write.xlsx(odds_exp, "tables/coin_toss_odds_exp.xlsx", rowNames = F, overwrite = T)


################################################################
# Table of binomial probabilities for "effect" experiments
################################################################

df <- openxlsx::read.xlsx("tables/coin_toss_odds_exp.xlsx", sheet = 1, na.strings = "NA") 

df$Number.of.experiments <- df$No.Effect+df$Effect

df$`P-value` <- dbinom(df$Effect, df$Number.of.experiments, 0.5)
df$Percent_effect <- percent(round(df$Effect/df$Number.of.experiments, digits = 2))

#calculate 95% confidence interval
confid <- binconf(x=df$Effect, n=df$Number.of.experiments, alpha=.05, include.x=T, include.n=T, return.df = T)
confid$Lower <- round(confid$Lower, digits=2)
confid$Upper <- round(confid$Upper, digits=2)
confid$Lower <- 100*confid$Lower
confid$Upper <- 100*confid$Upper
df$`95% CI` <- paste0(confid$Lower,"--",confid$Upper," %")

# optional
#confid$Upper-confid$Lower = 2*Z_a * sqrt(Variance)  # Z alpha/2 is 1.96 for alpha 0.05
#sqrt(Variance) = (confid$Upper-confid$Lower)/2*1.96
#df$Variance = ((confid$Upper-confid$Lower)/3.92)^2

names(df) <- c("EMF type", "E-field strength", "Insect", "Experiments (n= )", "Effect", "No effect", "p-value", "Effect (%)", "95% CI") 
View(df)

##############################
# table of probabilities
##############################

options(digits = 4)
options(scipen = 999)
df$`p-value`
df$`p-value` <- round(df$`p-value`, digits=4)
df$`p-value` <- as.character(df$`p-value`)
df <- mutate(df, `p-value` = if_else ( `p-value` < 0.0001, "< 0.0001", `p-value`))
#View(df)

openxlsx::write.xlsx(df, "tables/Table_of_probabilities_(experiments).xlsx", rowNames = F, overwrite = T, colWiths = "auto")

pdf("tables/table_of_probabilities_exp.pdf", height=5, width=10) # save to pdf
grid.table(df, rows=NULL)
dev.off()

# latex code for table
latexcode <- xtable::xtable(df, caption = 'Pooling the data of many experiments', auto = T, digits = 4)
print(latexcode, include.rownames = F)

