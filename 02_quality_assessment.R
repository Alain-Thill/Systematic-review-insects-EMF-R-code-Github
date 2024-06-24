
# run this first!
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#rm(list=ls())  # clear global environment of stuff from previous sessions

library(pacman)

pacman::p_load(dplyr, RefManageR, bibtex, sqldf, stringi, tidyverse, openxlsx)

############################################################
#  Getting list of studies older than 1980
############################################################

bib <- ReadBib("JabRef_5_HFLF.bib", check = "error")
df <- as.data.frame(bib, row.names = TRUE)
df$Identifier <- row.names(df)
str(df)

# grabbing year from date field
df$year <- substr(df$date, 1, 4)
df$year <- as.numeric(df$year)

# normalize ordering of subgroups that were previously manually assigned in JabRef
df$groups <- gsub("HF, HF_LF_exp", "HF_LF_exp, HF", df$groups)
df$groups <- gsub("LF, HF_LF_exp", "HF_LF_exp, LF", df$groups)
table(df$groups)

# using sql semantic to grab all HF studies
HF = sqldf("SELECT * FROM df WHERE [groups] LIKE 'HF_LF_exp, HF%' ")
#View(HF)
HF_list <- HF$Identifier
HF_list <- sort(HF_list)
HF_list <- as.data.frame(HF_list)
HF_list

LF = sqldf("SELECT * FROM df WHERE [groups] LIKE 'HF_LF_exp, LF%' ")
#View(LF)
LF_list <- LF$Identifier
LF_list <- sort(LF_list)
LF_list <- as.data.frame(LF_list)
LF_list

############################################################
#  Getting lists of accepted or rejected studies
############################################################

pacman::p_load(openxlsx, dplyr, sqldf, RefManageR, ggplot2, patchwork, ggpubr)

df <- openxlsx::read.xlsx("quality_HF.xlsx", sheet = 1) 

rejected_HF = sqldf("SELECT [Criteria] FROM df WHERE [Total.criteria.met] < 11 ")
rejected_HF
okay_quality_HF = sqldf("SELECT [Criteria] FROM df WHERE [Total.criteria.met] > 10 ")
okay_quality_HF


p1 <- ggplot(data = df, aes(Total.criteria.met)) + 
  labs(title = "Quality evaluation", y = "Number of studies", x = "Quality score (1--13)") + 
  theme_classic() + geom_histogram(binwidth = 0.5, fill="orange") +
  theme(legend.position = "bottom")
p1


df2 <- openxlsx::read.xlsx("quality_LF.xlsx", sheet = 1) 

rejected_LF = sqldf("SELECT [Criteria] FROM df2 WHERE [Total.criteria.met] < 11 ")
rejected_LF
okay_quality = sqldf("SELECT [Criteria] FROM df2 WHERE [Total.criteria.met] > 10 ")
okay_quality

p2 <- ggplot(data = df2, aes(Total.criteria.met)) + 
  labs(title = "Quality evaluation", y = "Number of studies", x = "Quality score (1--13)") + 
  theme_classic() + geom_histogram(binwidth = 0.5, fill="blue") +
  theme(legend.position = "bottom")
p2


# joing HF and LF groups
df$type <- "HF"
df2$type <- "LF"
df3 <- full_join(df,df2)

ggplot(df3, aes(Total.criteria.met, fill = type)) +
  geom_histogram(binwidth = 0.5, position="dodge") +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  labs(title = "", y = "Number of studies", x = "Quality score", fill = "Type") +
  #labs(title = "", y = "Anzahl an Studien", x = "Qualitätsbewertung", fill = "Typ") +
  theme_classic() +
  geom_vline(xintercept = 10.5, linetype="dashed", color = "black", size=1, ) +
  theme(legend.position = "bottom")


ggsave("figures/quality_scores.jpg", width = 4, height = 4)
ggsave("figures/quality_scores.png", width = 4, height = 4)


