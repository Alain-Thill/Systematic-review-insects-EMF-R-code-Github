
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
getwd() # working directory set to folder of .R file

#rm(list=ls())  #clear global environment of stuff from previous sessions


################################################################
# Plots of estimated toxicity
################################################################

pacman::p_load(dplyr, tidyr, ggplot2, ggpubr, scales, RColorBrewer, openxlsx, sqldf)

###########################
# formatting
###########################

df1 <- read.xlsx("tables/HFLF_meta_table.xlsx", sheet = 1, na.strings = "NA") 


# for observational base station studies, set exposure duration to one year
df1 <- mutate(df1, CummHrs = if_else ( is.na(CummHrs), 365*24, CummHrs))

df1 <- subset(df1, EMF_source=="Base station" | EMF_source=="DECT" | EMF_source=="Mobile phone" | EMF_source=="RF signal generator" | EMF_source=="WiFi" | EMF_source=="Coil system 50/60 Hz" )
df1$EMF_source <- gsub("RF signal generator", "Signal generator", df1$EMF_source)
df1$EMF_source <- gsub("Coil system 50/60 Hz", "Coil system", df1$EMF_source)

sort(table(df1$EMF_source), decreasing = F)
df1 <- subset(df1, EMF_source  != "WiFi")


# German translation
# df1$EMF_source <- gsub("Mobile phone", "Handy", df1$EMF_source)
# df1$EMF_source <- gsub("Base station", "Basisstation", df1$EMF_source)
# df1$EMF_source <- gsub("Signal generator", "Signalgenerator", df1$EMF_source)
# df1$EMF_source <- gsub("WiFi", "WLAN", df1$EMF_source)
# df1$EMF_source <- gsub("Coil system", "Spulensystem", df1$EMF_source)

# indicating plot order
df1$EMF_source <- factor(df1$EMF_source, levels=c("Base station","DECT","Mobile phone","Signal generator","Coil system"), ordered=TRUE)
#df1$EMF_source <- factor(df1$EMF_source, levels=c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system"), ordered=TRUE)
#df1$EMF_source <- factor(df1$EMF_source, levels=c("Basisstation","DECT","Handy","Signalgenerator","WLAN","Spulensystem"), ordered=TRUE)
table(df1$EMF_source)

# color_table <- tibble(
#   EMF_source = c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system"),
#   #EMF_source = c("Basisstation","DECT","Handy","Signalgenerator","WLAN,"Spulensystem"),
#   Color =       c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#999999", "#FF7F00"))

color_table <- tibble(
  EMF_source = c("Base station","DECT","Mobile phone","Signal generator","Coil system"),
  #EMF_source = c("Basisstation","DECT","Handy","Signalgenerator","WLAN,"Spulensystem"),
  Color =       c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))


df1 <- mutate(df1, ROM_norm = if_else ( Direction_of_effect == "beneficial", 1/ROM_norm, ROM_norm))
df2 <- df1[complete.cases(df1[,"ROM_norm"]),]  # select all rows that have effect size value


r1<-ggplot(aes(x = EMF_source, y = (ROM_norm-1), col=EMF_source), data = df2) + 
  geom_boxplot() + xlab(NULL) + #labs(col="EMF", title="A") +
  labs(col="EMF", title="A") +
  scale_y_continuous(limits = c(-0.75, 2.75), labels = percent) +  theme_classic2() +
  ylab(expression(paste("Percent change vs. control, norm."))) +
  #ylab(expression(paste("Proz. Veränderung ggü. Kontrolle"))) + 
  theme(axis.title.y = element_text(size = 11), plot.margin = unit(c(0.1,0,0.1,0), 'lines')) + 
  scale_color_manual(values = color_table$Color, drop = FALSE)

r1

samp <- dplyr::count(df2, EMF_source) # sample sizes
samp
ggplot_build(r1)$data
ggstat <- ggplot_build(r1)$data # retrieving the data from the plot using ggplot_build

ggwhisk <- ggstat[[1]]$ymax
ggwhisk <- data.frame(samp, whisk = ggwhisk) # combining that with the sample size, and call that data in geom_text

r1 <- r1 + geom_text(data = ggwhisk, size = 3,
                     aes(x = EMF_source, y = whisk, label = paste0("n =", n), vjust = -0.2, hjust = -0.3)) + 
                     geom_hline(yintercept = 0, linetype = 'dashed') +
                     geom_hline(yintercept = 0.5, linetype="dotted")
r1


r2 <- ggplot(aes(x = EMF_source, y = ROM_norm, col=EMF_source), data = df2) + 
  geom_boxplot() + xlab(NULL) + labs(col="EMF", title="A") +
  scale_y_continuous(limits = c(0.3, 4)) +  theme_classic2() +
  ylab(expression(paste("Effect size (ROM)")))+
  #ylab(expression(paste("Verhältnis der Mittelwerte")))+
  theme(axis.title.y = element_text(size = 11), plot.margin = unit(c(0,0,0,0), 'lines')) +
  scale_color_manual(values = color_table$Color, drop = FALSE)
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


########################################################
# adding effect sizes from meta-analysis (see part 6)
########################################################


es_data <- read.csv("tables/all_es_data.csv", header=T)
str(es_data)
es_data

es_data <- arrange(es_data, estimate, EMF_source)
es_data[es_data$EMF_source == 'DECT',]
#es_data[es_data$EMF_source == 'Base station',]

# basic rule:
# if CI doesn't overlap null-range of effect size, here 1, then effect significant at p=< 0.05

# assigning confidence stars according to p value using the "meta" value; change to "metafor" value by changing estimate to 4
es_data$sig <- ""
es_data <- mutate(es_data, sig = if_else (estimate == 4 & pval < 0.05, "*", sig))
es_data <- mutate(es_data, sig = if_else (estimate == 4 & pval < 0.01, "**", sig))
es_data <- mutate(es_data, sig = if_else (estimate == 4 & pval < 0.001, "***", sig))
es_data

es_data <- subset(es_data, EMF_source  != "WiFi")

## German translation
#es_data$EMF_source <- gsub("Mobile phone", "Handy", es_data$EMF_source)
#es_data$EMF_source <- gsub("Base station", "Basisstation", es_data$EMF_source)
#es_data$EMF_source <- gsub("Signal generator", "Signalgenerator", es_data$EMF_source)
#es_data$EMF_source <- gsub("Coil system", "Spulensystem", es_data$EMF_source)

#es_data$EMF_source <- factor(es_data$EMF_source, levels=c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system"), ordered=TRUE)
es_data$EMF_source <- factor(es_data$EMF_source, levels=c("Base station","DECT","Mobile phone","Signal generator","Coil system"), ordered=TRUE)
#es_data$EMF_source <- factor(es_data$EMF_source, levels=c("Basisstation","DECT","Handy","Signalgenerator","Spulensystem"), ordered=TRUE)

str(es_data)
#View(es_data)
es_data

# two_estimates_only <- es_data %>% filter(estimate < 2)
# two_estimates_only

five_estimates <- es_data %>% filter(estimate < 5)
five_estimates

mean_estimate <- five_estimates %>%
  group_by(EMF_source) %>% 
  summarise(k = mean(k), ROM = mean(ROM), LLCI = mean(LLCI), ULCI = mean(ULCI))

mean_estimate$estimate <- 5
mean_estimate$pval <- ""
mean_estimate$sig <- ""
mean_estimate

five_estimates_with_mean <- rbind(five_estimates, mean_estimate)
five_estimates_with_mean
mean_estimate <- five_estimates_with_mean[five_estimates_with_mean$estimate == 5,]
mean_estimate

color_table <- tibble(
  EMF_source = c("Base station","DECT","Mobile phone","Signal generator","Coil system"),
  #EMF_source = c("Basisstation","DECT","Handy","Signalgenerator","WLAN,"Spulensystem"),
  Color =       c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))

p <- ggplot(five_estimates, aes(x = EMF_source, y = ROM, col=EMF_source, group = estimate)) +
  geom_pointrange(aes(group = estimate, y = ROM, ymin = LLCI, ymax = ULCI), position = position_dodge(width = 0.5)) +
  xlab(NULL) + labs(col="EMF source", title="B") +
  scale_y_continuous(limits = c(0.7, 2.75)) +  theme_classic2() +
  ylab(expression(paste("Pooled effect size (ROM)")))+
  #ylab(expression(paste("Gepoolte Wirkung (Verh. der Mittelwerte)")))+
  theme(axis.title.y = element_text(size = 11)) +
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_hline(yintercept = 1.5, linetype="dotted")

p <- p + scale_color_manual(values = color_table$Color, drop = FALSE) 
p <- p + geom_text(data = es_data[es_data$estimate == 0,], size = 3, aes(x = EMF_source, y = ULCI, label = paste0("n =", k, " ", sig), vjust = -2.5, hjust = 1.5))
p <- p + geom_pointrange(data = mean_estimate, aes(group = estimate, y = ROM, ymin = LLCI, ymax = ULCI), linewidth = 1, position = position_nudge(x = 0.3))


p

ggarrange(r2, p, ncol = 1, nrow = 2, align="v", legend="none", common.legend = T)

ggsave("figures/Estimated_Toxicity_ROM_5_bars_with_mean.jpg", width = 12, height = 6)
ggsave("figures/Estimated_Toxicity_ROM_5_bars_with_mean.png", width = 12, height = 6)


es_data <- arrange(es_data, EMF_source, estimate)
es_data
split_base_station <- es_data[c(6:35),] # select only base station subgroups (altered behavior vs direct markers of toxicity)
split_base_station


p2 <- ggplot(split_base_station, aes(x = EMF_source, y = ROM-1, col=EMF_source, group = estimate)) +  
  geom_pointrange(aes(group = estimate, y = ROM-1, ymin = LLCI-1, ymax = ULCI-1), position = position_dodge(width = 0.5)) +
  xlab(NULL) + labs(col="Field source", title="C") +
  scale_y_continuous(limits = c(-0.5, 2.75), labels = percent) +  theme_classic2() +
  ylab(expression(paste("Pooled effect size (percent change)")))+
  #ylab(expression(paste("Gepoolte Wirkung (proz. Veränderung)")))+
  theme(axis.title.y = element_text(size = 11)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_hline(yintercept = 0.5, linetype="dotted")

p2
p2 <- p2 + scale_color_manual(values = color_table$Color, drop = FALSE) 
p2 <- p2 + geom_text(data = es_data[es_data$estimate == 0 & es_data$EMF_source != "Base station",], size = 3, aes(x = EMF_source, y = ULCI-1, label = paste0("n =", k, " ", sig), vjust = -1, hjust = 1.35))
#p2 <- p2 + geom_text(data = es_data[es_data$estimate == 0 & es_data$EMF_source != "Basisstation",], size = 3, aes(x = EMF_source, y = ULCI-1, label = paste0("n =", k, " ", sig), vjust = -1, hjust = 1.35))

p2 <- p2 + geom_text(data = es_data[es_data$estimate == 5,], size = 3, aes(x = EMF_source, y = ULCI-1, label = paste0("n =", k, " ", sig), vjust = -2, hjust = 2))
p2 <- p2 + geom_text(data = es_data[es_data$estimate == 10,], size = 3, aes(x = EMF_source, y = ULCI-1, label = paste0("n =", k, " ", sig), vjust = -2, hjust = 0))

p2

ggarrange(r1, p2, ncol = 1, nrow = 2, align="v",
          legend="none", common.legend = T)

ggsave("figures/Estimated_Toxicity_percent_base_split_5_bars.jpg", width = 12, height = 6)
ggsave("figures/Estimated_Toxicity_percent_base_split_5_bars.png", width = 12, height = 6)


es_data <- arrange(es_data, EMF_source, estimate)
es_data
base_station_totals_and_tox <- subset(es_data, k != "22") # select base station totals and "direct markers of toxicity" subgroup)
base_station_totals_and_tox


p3 <- ggplot(base_station_totals_and_tox, aes(x = EMF_source, y = ROM-1, col=EMF_source, group = estimate)) +  
  geom_pointrange(aes(group = estimate, y = ROM-1, ymin = LLCI-1, ymax = ULCI-1), position = position_dodge(width = 0.5)) +
  xlab(NULL) + labs(col="Field source", title="C") +
  scale_y_continuous(limits = c(-0.5, 1.75), labels = percent) +  theme_classic2() +
  ylab(expression(paste("Pooled effect size (percent change)")))+
  #ylab(expression(paste("Gepoolte Wirkung (proz. Veränderung)")))+
  theme(axis.title.y = element_text(size = 11)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_hline(yintercept = 0.5, linetype="dotted")

p3
p3 <- p3 + scale_color_manual(values = color_table$Color, drop = FALSE) 
p3 <- p3 + geom_text(data = es_data[es_data$estimate == 0,], size = 3, aes(x = EMF_source, y = ULCI-1, label = paste0("n =", k, " ", sig), vjust = -1.75, hjust = 1.75))
p3 <- p3 + geom_text(data = es_data[es_data$estimate == 5,], size = 3, aes(x = EMF_source, y = ULCI-1, label = paste0("n =", k, " ", sig), vjust = -2, hjust = 0))

p3

ggarrange(r1, p3, ncol = 1, nrow = 2, align="v",
          legend="none", common.legend = T)

ggsave("figures/Estimated_Toxicity_percent_base_with_tox_5_bars.jpg", width = 12, height = 6)
ggsave("figures/Estimated_Toxicity_percent_base_with_tox_5_bars.png", width = 12, height = 6)

