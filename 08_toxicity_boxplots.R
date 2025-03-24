
# run this first

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
getwd() # working directory set to folder of .R file

#rm(list=ls())  #clear global environment of stuff from previous sessions

if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(dplyr, tidyr, ggplot2, ggpubr, scales, RColorBrewer, openxlsx, sqldf, data.table, metafor, RoBMA)



###########################################################################
# Toxicity: standardizing measured effects, calculating toxicity indicator
###########################################################################

df <- openxlsx::read.xlsx("tables/data_table_HFLF_5.xlsx", sheet = 1, na.strings = "NA") 

names(df)

# Normalization of log_ROM to compare "increased ovarian apoptosis" with decreased reproductive capacity
# in general, we want all uncertain and detrimental bioeffects to be counted as positive, and all beneficial effects as negative
# but there is one exception, base station studies where reduced or increased abundance was observed: here invert
df <- mutate(df, log_ROM = if_else ( EMF_source == "Base station" & Bioeffect_cat == "Altered behavior", -log_ROM, log_ROM))

# use "Not Or Not" logic to select all other studies not "base station + altered behavior"
df <- mutate(df, log_ROM = if_else ( EMF_source != "Base station" | Bioeffect_cat != "Altered behavior", abs(log_ROM), log_ROM))

# Assignment of negative effect size for all experiments that find beneficial effects
df <- mutate(df, log_ROM = if_else ( Direction_of_effect == "beneficial", -abs(log_ROM), log_ROM))

# Normalization of ROM quotients < 1 to their inverse value (normalized ROM)
df <- mutate(df, ROM_norm = if_else ( Ratio_of_means< 1, 1/Ratio_of_means, Ratio_of_means))
# Assignment of negative effect size for all experiments that find beneficial effects
df <- mutate(df, ROM_norm = if_else ( Direction_of_effect == "beneficial", 1/ROM_norm, ROM_norm))
df <- relocate(df, ROM_norm, .after = "Ratio_of_means")


# Calculate the variance
df <- escalc(measure = "ROM", ni=n_e, yi=log_ROM, sei=log_SE, data=df, slab=study)
df$sei <- sqrt(df$vi)

# Adding z values
df$z <- logOR2z(df$log_ROM)
df$se_z <- se_logOR2se_z(df$log_SE, df$log_ROM)

#View(df)

openxlsx::write.xlsx(df, "tables/HFLF_meta_table.xlsx", rowNames = F, colWidths = 15, overwrite = T, firstRow = T, firstCol = T)


###########################################################
# Counting percentages of present data among all studies
###########################################################

df <- openxlsx::read.xlsx("tables/data_table_HFLF_3.xlsx", sheet = 1, na.strings = "NA")

# percentage of experiments indicating effect size values
eff_size <- sqldf("SELECT * FROM df WHERE Ratio_of_means != 'NA'")

length(eff_size$study)/length(df$study)
length(unique(eff_size$study))/length(unique(df$study))
# 61.8 % of exp
# 55.3 % of studies

# percentage of experiments with p values BEFORE assigning p = 0.5 to "no effect" studies
pval <- sqldf("SELECT * FROM df WHERE p NOT NULL")

length(pval$study)/length(df$study)
length(unique(pval$study))/length(unique(df$study))
# 53.6 % of exp
# 58.7 % of studies

df <- openxlsx::read.xlsx("tables/data_table_HFLF_4.xlsx", sheet = 1, na.strings = "NA")

# percentage of experiments with p values AFTER assigning p = 0.5 to "no effect" studies
pval <- sqldf("SELECT * FROM df WHERE p NOT NULL")
length(pval$study)/length(df$study)
length(unique(pval$study))/length(unique(df$study))
# 58.1 % of exp
# 62.8 % of studies

# percentage of experiments with SE values (amenable to meta-analysis)
logSE <- sqldf("SELECT * FROM df WHERE log_SE NOT NULL")
length(logSE$study)/length(df$study)
length(unique(logSE$study))/length(unique(df$study))
# 58.9 % of exp
# 47.9 % of studies


################################################################
# Plots of toxicity
################################################################

pacman::p_load(dplyr, tidyr, ggplot2, ggpubr, scales, RColorBrewer, openxlsx, sqldf)

###########################
# formatting
###########################

df1 <- openxlsx::read.xlsx("tables/HFLF_meta_table.xlsx", sheet = 1, na.strings = "NA") 
dev.off()
#View(df1)

# for observational base station studies, set duration to one year (higher than most insects' lifetime)
df1 <- mutate(df1, CummHrs = if_else ( is.na(CummHrs), 365*24, CummHrs))

df1$EMF_source <- as.factor(df1$EMF_source)
sort(summary(df1$EMF_source), decreasing = F)


df1 <- subset(df1, EMF_source=="Base station" | EMF_source=="DECT" | EMF_source=="Mobile phone" | EMF_source=="RF signal generator" | EMF_source=="WiFi" | EMF_source=="Coil system 50/60 Hz" )

df1$EMF_source <- gsub("RF signal generator", "Signal generator", df1$EMF_source)
df1$EMF_source[grep("Coil", df1$EMF_source)] <- "Coil system"

# indicating plot order 
df1$EMF_source <- factor(df1$EMF_source, levels=c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system"), ordered=TRUE)
table(df1$EMF_source)

df2 <- df1[complete.cases(df1[,"ROM_norm"]),]  # select all rows that have effect size value

# Field parameter plots
ggplot(aes(y = log(ROM_norm), x = log(B_Field), col=EMF_source), data = df1) + 
  geom_point() + scale_color_brewer(palette="Set1") + labs(col="Field source")
#ggsave("figures/log_ROM~log_B-Field.jpg", width = 8, height = 8)

ggplot(aes(y = log(ROM_norm), x = log(E_Field), col=EMF_source), data = df1) + 
  geom_point() + scale_color_brewer(palette="Set1") + labs(col="Field source")
#ggsave("figures/log_ROM~log_E-Field.jpg", width = 8, height = 8)


####################################################
# Boxplots of effect size (toxicity)
####################################################

# custom color table
color_table <- tibble(
  EMF_source = c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system"),
  #EMF_source = c("Basisstation","DECT","Handy","Signalgenerator","WLAN,"Spulensystem"),
  Color =       c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#999999", "#FF7F00"))


r1<-ggplot(aes(x = EMF_source, y = (ROM_norm-1), col=EMF_source), data = df2) + 
  geom_boxplot() + xlab(NULL) + labs(col="Field source", title="A") +
  scale_y_continuous(limits = c(-0.75, 3), labels = percent) +  theme_classic2() +
  ylab(expression(paste("Percent change vs. control, norm."))) + #, symbol('\256')))) +
  theme(axis.title.y = element_text(size = 11), plot.margin = unit(c(0.1,0,0.1,0), 'lines')) + 
  scale_color_manual(values = color_table$Color, drop = FALSE)
  #scale_color_brewer(palette="Set1") 
r1

samp <- dplyr::count(df2, EMF_source) # sample sizes
samp
ggplot_build(r1)$data
ggstat <- ggplot_build(r1)$data # retrieving the data from the plot using ggplot_build

ggwhisk <- ggstat[[1]]$ymax
ggwhisk <- data.frame(samp, whisk = ggwhisk) # combining that with the sample size, and call that data in geom_text

r1 <- r1 + geom_text(data = ggwhisk, size = 3,
                   aes(x = EMF_source, y = whisk, label = paste0("n =", n), vjust = -0.2, hjust = -0.3)) + 
              geom_hline(yintercept = 0, linetype = 'dashed')
r1


r2 <- ggplot(aes(x = EMF_source, y = ROM_norm, col=EMF_source), data = df2) + 
  geom_boxplot() + xlab(NULL) + labs(col="Field source", title="A") +
  scale_y_continuous(limits = c(0.3, 4)) +  theme_classic2() +
  ylab(expression(paste("Effect size (ROM)")))+
  theme(axis.title.y = element_text(size = 11), plot.margin = unit(c(0,0,0,0), 'lines')) +
  scale_color_manual(values = color_table$Color, drop = FALSE)
  #scale_color_brewer(palette="Dark2") #+ geom_hline(yintercept = c(1,2))
r2

samp2 <- dplyr::count(df2, EMF_source) # sample sizes
samp2
ggstat2 <- ggplot_build(r2)$data # retrieving the data from the plot using ggplot_build
ggwhisk2 <- ggstat2[[1]]$ymax
ggwhisk2 <- data.frame(samp2, whisk = ggwhisk2) # combining that with the sample size, and call that data in geom_text

r2 <- r2 + geom_text(data = ggwhisk2, size = 3, aes(x = EMF_source, y = whisk, label = paste0("n =", n), vjust = -0.2, hjust = -0.3))+
           geom_hline(yintercept = 1, linetype = 'dashed')

r2

ggarrange(r1, r2, ncol = 1, nrow = 2, 
          labels = c("",""),  hjust=0, align="v",
          legend="none", common.legend = T)

################################################################
# Boxplots of toxicity by emitted or absorbed radiation dose
################################################################

# keep only experiments that found some toxicity
df3 <- subset(df1, Direction_of_effect == "detrimental" | Direction_of_effect == "uncertain")
df3 <- subset(df3, Ratio_of_means!=1)

#df3 <- df3 %>% filter(Confid < 0.051)  # keep only statistically significant findings
names(df3)

df3 <- df3[complete.cases(df3[,c("E_Field","CummHrs","ROM_norm")]),]  # keep only rows with E-field, duration and effect size 
df3 <- subset(df3, EMF_source!="Coil system")
#df3 <- subset(df3, EMF_source!="WiFi")

# various formulas to estimate toxicity:
#df3 <- mutate(df3, tox_norm = ROM_norm/(Power_mWm2)) # effect normalized to power density
df3 <- mutate(df3, tox_norm = ROM_norm/(Power_mWm2 * CummHrs)) # effect normalized to emitted radiation dose
# given as ROM / (mW/m^2 * h)


length(df1$study)
length(df2$study)
length(df3$study)

table(df3$EMF_source)

df3$EMF_source <- factor(df3$EMF_source, levels=c("Base station","DECT","Mobile phone","Signal generator","WiFi"), ordered=TRUE)

r3 <- ggplot(aes(x = EMF_source, y = log10(tox_norm), col=EMF_source), data = df3) + 
  geom_boxplot() + xlab(NULL)  + theme_classic2() +  scale_y_continuous() + #breaks = NULL
  labs(col="Field source", title="A") +
  ylab(expression(paste("log10 (Toxicity / Emitted power dose)"))) + 
  theme(axis.title.y = element_text(size = 11), plot.margin = unit(c(0.2,0.2,0.2,0.2), 'lines')) + 
  annotation_logticks(sides="l") + 
  scale_color_manual(values = color_table$Color, drop = FALSE)
  #scale_color_brewer(palette="Set1")

r3

samp3 <- dplyr::count(df3, EMF_source) # sample sizes
samp3
ggstat3 <- ggplot_build(r3)$data # retrieving the data from the plot using ggplot_build
ggwhisk3 <- ggstat3[[1]]$ymax
ggwhisk3 <- data.frame(samp3, whisk = ggwhisk3) # combining that with the sample size, and call that data in geom_text

r3 <- r3 + geom_text(data = ggwhisk3, size = 3,
                     aes(x = EMF_source, y = whisk, label = paste0("n =", n), vjust = -0.2, hjust = -0.3))
r3


########################################################
# SAR calculations
########################################################

pacman::p_load(dplyr, tidyr, ggplot2, ggpubr, scales, RColorBrewer, openxlsx, sqldf)

# opened SAR_frequencies.xlsx in Excel, added conversion factor based on Thielens and de Borre papers, resaved as SAR_calc

df4 <- openxlsx::read.xlsx("SAR_calc.xlsx", sheet = 1, na.strings = "NA")
#View(df4)
#df4 <- subset(df4, EMF_source!="WiFi")

df4$Insect_mass_mg <- NA
# assigning body mass per insect species or type
df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Drosophila", 1, Insect_mass_mg))
df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Ixodes ticks", 1, Insect_mass_mg))
df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Dermacentor ticks", 1, Insect_mass_mg))
df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Folsomia springtail", 1, Insect_mass_mg))
df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Asobara parasitic wasp", 1, Insect_mass_mg))

df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Musca domestica", 12, Insect_mass_mg))
df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Ant", 12, Insect_mass_mg))
df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Vespula maculata", 40, Insect_mass_mg))
df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Beetle", 100, Insect_mass_mg))

df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Honeybee", 900, Insect_mass_mg))
df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Wasps", 900, Insect_mass_mg))
df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Bee flies", 900, Insect_mass_mg))
df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Wild bees", 900, Insect_mass_mg))
df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Hoverflies", 900, Insect_mass_mg))

df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Bicyclus anynana", 2000, Insect_mass_mg))
df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Cockroach", 3000, Insect_mass_mg))
df4 <- mutate(df4, Insect_mass_mg = if_else ( Insect_type == "Locusta migratoria", 5000, Insect_mass_mg))

names(df4)

df4$Energy_uptake_nW_per_mWm2 <- df4$Energy_uptake_nW_per_Vm/0.377
df4$Power_mWm2 <- (df4$E_Field^2)/0.377
df4$SAR <- df4$Energy_uptake_nW_per_mWm2 * df4$Power_mWm2 * 10^-6 / (df4$Insect_mass_mg/1000)

# unit: mW / g  = W/kg

# control:
#df4$SAR2 <- df4$Energy_uptake_nW_per_mWm2 * df4$Power_mWm2 * 10^-6 / (df4$Insect_mass_mg/1000)
#df4 <- mutate(df4, SAR = if_else ( is.na(SAR), SAR2, SAR))
#View(df4[c("study","SAR","SAR2")])


# keep only experiments that found some toxicity or potential toxicity
df4 <- subset(df4, Direction_of_effect == "detrimental" | Direction_of_effect == "uncertain")

df4$Ratio_of_means

# Normierung der Quotienten < 1 auf ihren Umkehrwert
df4 <- mutate(df4, ROM_norm = if_else ( Ratio_of_means< 1, 1/Ratio_of_means, Ratio_of_means))
df4$ROM_norm
#df4$ROM_norm <- round(df4$ROM_norm, digits=3)
df4 <- subset(df4, Ratio_of_means!=1) # removing "no effect", since here only actually observed toxicity is to be quantified

df4 <- df4[complete.cases(df4[,c("SAR","CummHrs","ROM_norm")]),]  # keep only rows with SAR, duration and effect size 
table(df4$EMF_source)
#df4 <- subset(df4, EMF_source!="Bluetooth")
#df4 <- subset(df4, EMF_source!="Microwave Oven")

df4 <- mutate(df4, tox_norm = ROM_norm/(Power_mWm2 * CummHrs)) # effect normalized to emitted radiation dose
df4 <- mutate(df4, tox_norm2 = ROM_norm/(SAR * CummHrs)) # effect normalized to (absorbed) radiation dose
#View(df4)

length(df4$study) # 182 experiments with specification of B-field, duration and effect size
unique(df4$study) # out of 39 studies 

table(df4$EMF_source)
df4$EMF_source <- factor(df4$EMF_source, levels=c("Base station","DECT","Mobile phone","Signal generator","WiFi"), ordered=TRUE)

r4 <- ggplot(aes(x = EMF_source, y = log10(tox_norm2), col=EMF_source), data = df4) + 
  geom_boxplot() + xlab(NULL)  + theme_classic2() +  scale_y_continuous() + #breaks = NULL
  labs(col="Field source", title="B") +
  ylab(expression(paste("log10 (Toxicity / SAR dose)"))) + 
  theme(axis.title.y = element_text(size = 11), plot.margin = unit(c(0.2,0.2,0.2,0.2), 'lines')) + 
  annotation_logticks(sides="l") + 
  scale_color_manual(values = color_table$Color, drop = FALSE)
  #scale_color_brewer(palette="Set1")

r4

samp4 <- dplyr::count(df4, EMF_source) # sample sizes
samp4
ggstat4 <- ggplot_build(r4)$data # retrieving the data from the plot using ggplot_build
ggwhisk4 <- ggstat4[[1]]$ymax
ggwhisk4 <- data.frame(samp4, whisk = ggwhisk4) # combining that with the sample size, and call that data in geom_text

r4 <- r4 + geom_text(data = ggwhisk4, size = 3,
                     aes(x = EMF_source, y = whisk, label = paste0("n =", n), vjust = -0.2, hjust = -0.3))
r4

r3b <- ggplot(aes(x = EMF_source, y = log10(tox_norm), col=EMF_source), data = df4) + 
  geom_boxplot() + xlab(NULL)  + theme_classic2() +  scale_y_continuous() + #breaks = NULL
  labs(col="Field source", title="A") +
  ylab(expression(paste("log10 (Toxicity / Emitted power dose)"))) + 
  theme(axis.title.y = element_text(size = 11), plot.margin = unit(c(0.2,0.2,0.2,0.2), 'lines')) + 
  annotation_logticks(sides="l") + 
  scale_color_manual(values = color_table$Color, drop = FALSE)
  #scale_color_brewer(palette="Set1")

r3b

samp3b <- dplyr::count(df4, EMF_source) # sample sizes
samp3b
ggstat3b <- ggplot_build(r3b)$data # retrieving the data from the plot using ggplot_build
ggwhisk3b <- ggstat3b[[1]]$ymax
ggwhisk3b <- data.frame(samp3b, whisk = ggwhisk3b) # combining that with the sample size, and call that data in geom_text

r3b <- r3b + geom_text(data = ggwhisk3b, size = 3,
                     aes(x = EMF_source, y = whisk, label = paste0("n =", n), vjust = -0.2, hjust = -0.3))
r3b

ggarrange(r3, r4, ncol = 1, nrow = 2, 
          labels = c("",""),  hjust=0, align="v",
          legend="none", common.legend = T)

ggarrange(r3b, r4, ncol = 1, nrow = 2, 
          labels = c("",""),  hjust=0, align="v",
          legend="none", common.legend = T)

ggsave("figures/Toxicity_per_emitted_and_absorbed_dose.jpg", width = 10, height = 7)
ggsave("figures/Toxicity_per_emitted_and_absorbed_dose.png", width = 10, height = 7)


ggplot(aes(y = log(ROM_norm), x = log(Power_mWm2 * CummHrs), col=EMF_source), data = df4) + 
  geom_point() + scale_color_brewer(palette="Set1") + labs(col="Field source")

ggplot(aes(y = log(ROM_norm), x = log(SAR * CummHrs), col=EMF_source), data = df4) + 
  geom_point() + scale_color_brewer(palette="Set1") + labs(col="Field source")

ggplot(aes(y = log(ROM_norm), x = log(SAR), col=EMF_source), data = df4) + 
  geom_point() + scale_color_brewer(palette="Set1") + labs(col="Field source")

#ggsave("figures/log_ROM~log_SAR.jpg", width = 8, height = 8)
#ggsave("figures/log_ROM~log_SAR.png", width = 8, height = 8)

ggplot(aes(y = log(ROM_norm), x = log(Power_mWm2), col=EMF_source), data = df4) + 
  geom_point() + scale_color_brewer(palette="Set1") + labs(col="Field source")

ggplot(aes(y = log(ROM_norm), x = log(E_Field), col=EMF_source), data = df4) + 
  geom_point() + scale_color_brewer(palette="Set1") + labs(col="Field source")



