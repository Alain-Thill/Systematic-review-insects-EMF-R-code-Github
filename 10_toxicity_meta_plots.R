
# run this first!
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
getwd() # working directory set to folder of .R file

#rm(list=ls())  #clear global environment of stuff from previous sessions

if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(dplyr, tidyr, ggplot2, ggpubr, scales, RColorBrewer, openxlsx, sqldf)


################################################################
# Plots of estimated toxicity
################################################################

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
#df1 <- subset(df1, EMF_source  != "WiFi") # removing WiFi data

# indicating plot order
df1$EMF_source <- factor(df1$EMF_source, levels=c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system"), ordered=TRUE)

table(df1$EMF_source)

# define color table to link a specific color to each EMF type for all plots
color_table <- c("Base station" = "red", "Base station (tox)" = "red", "Base station (behavior)" = "red",
                 "DECT" = "#377EB8", 
                 "Mobile phone" = "#4DAF4A", 
                 "Signal generator" = "#984EA3",
                 "Coil system" = "#FF7F00", 
                 "WiFi" = "#999999")


df2 <- df1[complete.cases(df1[,"ROM_norm"]),]  # select all rows that have effect size value

#dev.off()

studies <- df2 %>% group_by(EMF_source) %>% summarise(n = n_distinct(study)) # sample sizes: number of studies
samp <- df2 %>% group_by(EMF_source) %>% summarise(k = n()) # sample sizes: number of experiments
samp <- full_join(studies, samp)
samp

# Ratio of means (ROM, of treatment vs control) boxplot with WiFi

r1 <- ggplot(aes(x = EMF_source, y = ROM_norm, col=EMF_source), data = df2) + 
  geom_boxplot() + xlab(NULL) + labs(col="EMF", title="A") +
  scale_y_continuous(limits = c(0.3, 4)) +  theme_classic2() +
  ylab(expression(paste("Effect size (ROM)")))+
  theme(axis.title.y = element_text(size = 11), plot.margin = unit(c(0,0,0,0), 'lines')) +
  scale_color_manual(values = color_table, drop = TRUE)
r1

ggstat1 <- ggplot_build(r1)$data # retrieving the data from the plot using ggplot_build
ggwhisk1 <- ggstat1[[1]]$ymax
ggwhisk1 <- data.frame(samp, whisk = ggwhisk1) # combining that with the sample size, and call that data in geom_text

r1 <- r1 + geom_text(data = ggwhisk1, size = 3, aes(x = EMF_source, y = whisk, label = paste0("n =", n,", k =", k), vjust = -0.4, hjust = -0.3)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_hline(yintercept = 1.5, linetype="dotted")

r1

# Percent change vs control boxplot (with WiFi)

r1b <- ggplot(aes(x = EMF_source, y = ROM_norm-1, col=EMF_source), data = df2) + 
  geom_boxplot() + xlab(NULL) +
  labs(col="EMF", title="A") +
  scale_y_continuous(limits = c(0.3-1, 4-1), labels = percent) +  theme_classic2() +
  ylab(expression(paste("Percent change vs. control, norm."))) +
  theme(axis.title.y = element_text(size = 11), plot.margin = unit(c(0.1,0,0.1,0), 'lines')) + 
  scale_color_manual(values = color_table, drop = TRUE)

r1b

ggplot_build(r1b)$data
ggstat1b <- ggplot_build(r1b)$data # retrieving the data from the plot using ggplot_build

ggwhisk1b <- ggstat1b[[1]]$ymax
ggwhisk1b <- data.frame(samp, whisk = ggwhisk1b) # combining that with the sample size, and call that data in geom_text

r1b <- r1b + geom_text(data = ggwhisk1b, size = 3, aes(x = EMF_source, y = whisk, label = paste0("n =", n,", k =", k), vjust = -0.4, hjust = -0.3)) +
           geom_hline(yintercept = 0, linetype = 'dashed') +
           geom_hline(yintercept = 0.5, linetype="dotted")
r1b

#dev.off()


ggarrange(r1, r1b, ncol = 1, nrow = 2, 
          labels = c("",""),  hjust=0, align="v",
          legend="none", common.legend = T)



df3 <- subset(df2, EMF_source  != "WiFi") # removing WiFi data, comment out if you want to plot WiFi data also
table(df3$EMF_source)
df3$EMF_source <- factor(df3$EMF_source, levels=c("Base station","DECT","Mobile phone","Signal generator","Coil system"), ordered=TRUE)
table(df3$EMF_source)

# ROM boxplot without WIFI

r2 <- ggplot(aes(x = EMF_source, y = ROM_norm, col=EMF_source), data = df3) + 
  geom_boxplot() + xlab(NULL) + labs(col="EMF", title="A") +
  scale_y_continuous(limits = c(0.3, 4)) +  theme_classic2() +
  ylab(expression(paste("Effect size (ROM)")))+
  theme(axis.title.y = element_text(size = 11), plot.margin = unit(c(0,0,0,0), 'lines')) +
  scale_color_manual(values = color_table, drop = TRUE)

r2


ggstat2 <- ggplot_build(r2)$data # retrieving the data from the plot using ggplot_build
ggwhisk2 <- ggstat2[[1]]$ymax
samp2 <- samp[samp$EMF_source != "WiFi",]
ggwhisk2 <- data.frame(samp2, whisk = ggwhisk2) # combining that with the sample size, and call that data in geom_text
ggwhisk2

r2 <- r2 + geom_text(data = ggwhisk2, size = 3, aes(x = EMF_source, y = whisk, label = paste0("n =", n,", k =", k), vjust = -0.4, hjust = -0.3)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_hline(yintercept = 1.5, linetype="dotted")

r2

#####################################################################
# adding effect sizes from meta-analysis (see part 09_meta_analysis)
#####################################################################

es_data <- read.xlsx("tables/all_effect_size_estimates_from_meta-analysis.xlsx", sheet = 1, na.strings = "NA") 

str(es_data)
es_data

es_data <- arrange(es_data, estimate, EMF_source)

# assigning confidence stars according to p value using the "meta" value; change to "metafor" value by changing estimate to 1
es_data$sig <- ""
es_data <- mutate(es_data, sig = if_else (estimate == 0 & pval < 0.05, "*", sig))
es_data <- mutate(es_data, sig = if_else (estimate == 0 & pval < 0.01, "**", sig))
es_data <- mutate(es_data, sig = if_else (estimate == 0 & pval < 0.001, "***", sig))
es_data

#es_data <- subset(es_data, EMF_source  != "WiFi") # removing the WiFi studies - so far too few
es_data <- subset(es_data, estimate  != 2) # removing the Bayesmeta non-clustered estimate (lvl1, not at study-level)
es_data$EMF_source <- factor(es_data$EMF_source, levels=c("Base station","Base station (tox)","Base station (behavior)","DECT","Mobile phone","Signal generator","WiFi","Coil system"), ordered=TRUE)


# Plot 1: plotting base station without subgroups, with average of 4 estimates for each EMF source

plain_estimates <- subset(es_data, EMF_source  != "Base station (tox)")
plain_estimates <- subset(plain_estimates, EMF_source  != "Base station (behavior)")
plain_estimates # base station subgroups removed
plain_estimates$EMF_source <- factor(plain_estimates$EMF_source, levels=c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system"), ordered=TRUE)
str(plain_estimates$EMF_source)

plain_estimates

mean_estimate <- plain_estimates %>%
  group_by(EMF_source) %>% 
  summarise(n = mean(n), k = mean(k), ROM = mean(ROM), LLCI = mean(LLCI), ULCI = mean(ULCI))
mean_estimate$estimate <- "mean"
mean_estimate

#dev.off()

p1 <- ggplot(plain_estimates, aes(x = EMF_source, y = ROM, col=EMF_source, group = estimate)) +
  geom_pointrange(aes(group = estimate, y = ROM, ymin = LLCI, ymax = ULCI), position = position_dodge(width = 0.35)) +
  xlab(NULL) + labs(col="EMF source", title="B") +
  scale_y_continuous(limits = c(0.65, 3.2)) +  theme_classic2() +
  ylab(expression(paste("Pooled effect size (ROM)")))+
  theme(axis.title.y = element_text(size = 11)) +
  scale_color_manual(values = color_table, drop = TRUE) +
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_hline(yintercept = 1.5, linetype="dotted")

p1

p1 <- p1 + geom_pointrange(data = mean_estimate, aes(group = estimate, y = ROM, ymin = LLCI, ymax = ULCI), linewidth = 1, position = position_nudge(x = 0.25))
p1 <- p1 + geom_text(data = mean_estimate, size = 3, aes(x = EMF_source, y = ULCI, label = paste0("n =", n,", k =", k), vjust = -2.5))

p1

ggarrange(r1, p1, ncol = 1, nrow = 2, align="v", legend="none", common.legend = T)

ggsave("figures/Estimated_Toxicity_ROM_with_mean.jpg", width = 12, height = 6)
ggsave("figures/Estimated_Toxicity_ROM_with_mean.png", width = 12, height = 6)


# Plot 2: plotting base station subgroups separately, with average of 4 estimates for each EMF source

# select only base station subgroups (altered behavior vs direct markers of toxicity)
split_base_station <- subset(es_data, EMF_source  != "Base station")
#split_base_station <- subset(split_base_station, EMF_source  != "WiFi")
split_base_station
str(split_base_station$EMF_source)

mean_estimate_split <- split_base_station %>%
  group_by(EMF_source) %>% 
  summarise(n = mean(n), k = mean(k), ROM = mean(ROM), LLCI = mean(LLCI), ULCI = mean(ULCI))
mean_estimate_split$estimate <- "mean"
mean_estimate_split

#dev.off()

p2 <- ggplot(split_base_station, aes(x = EMF_source, y = ROM, col=EMF_source, group = estimate)) +  
  geom_pointrange(aes(group = estimate, y = ROM, ymin = LLCI, ymax = ULCI), position = position_dodge(width = 0.35)) +
  xlab(NULL) + labs(col="Field source", title="C") +
  scale_y_continuous(limits = c(0.37, 5.6), labels = waiver()) +  theme_classic2() +
  ylab(expression(paste("Pooled effect size (ROM)"))) +
  theme(axis.title.y = element_text(size = 11)) +
  scale_color_manual(values = color_table, drop = TRUE) +
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_hline(yintercept = 1.5, linetype="dotted")

p2

p2 <- p2 + geom_text(data = mean_estimate_split, size = 3, aes(x = EMF_source, y = ROM, label = paste0("n =", n,", k =", k), vjust = -2, hjust = -1.1))
p2 <- p2 + geom_pointrange(data = mean_estimate_split, aes(group = estimate, y = ROM, ymin = LLCI, ymax = ULCI), linewidth = 1.3, position = position_nudge(x = 0.25))
p2


# two values are missing because they exceed plot limits: draw manually with arrows on end
# draw segment with arrow to shorten the very large CI Nr.4 (brms) of Base station (behavior)
p2 <- p2 + geom_segment(
  x = 2.13, y = 0.706,
  xend = 2.13, yend = 5.6,
  lineend = "round",
  linejoin = "round",
  linewidth = 0.5,
  size = 1, 
  arrow = arrow(length = unit(0.3, "cm")),
  colour = "red"
) 

# draw segment with arrow to shorten the very large CI Nr.5 (average) of Base station (behavior)
p2 <- p2 + geom_segment(
  x = 2.25, y = 0.372,
  xend = 2.25, yend = 5.6,
  lineend = "round",
  linejoin = "round",
  linewidth = 1.3,
  size = 1, 
  arrow = arrow(length = unit(0.3, "cm")),
  colour = "red"
) 

p2

#dev.off()
ggarrange(r1, p2, ncol = 1, nrow = 2, align="v",
          legend="none", common.legend = T)

ggsave("figures/Estimated_Toxicity_ROM_base_split.jpg", width = 12, height = 6)
ggsave("figures/Estimated_Toxicity_ROM_base_split.png", width = 12, height = 6)


# same plot but expressing values in percent change vs control instead of ROM
p2b <- ggplot(split_base_station, aes(x = EMF_source, y = ROM-1, col=EMF_source, group = estimate)) +  
  geom_pointrange(aes(group = estimate, y = ROM-1, ymin = LLCI-1, ymax = ULCI-1), position = position_dodge(width = 0.35)) +
  xlab(NULL) + labs(col="Field source", title="C") +
  scale_y_continuous(limits = c(0.37-1, 5.6-1), labels = percent) +  theme_classic2() +
  ylab(expression(paste("Pooled effect size (percent change)")))+
  theme(axis.title.y = element_text(size = 11)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_hline(yintercept = 0.5, linetype="dotted") +
  scale_color_manual(values = color_table, drop = TRUE) +
  geom_text(data = mean_estimate_split, size = 3, aes(x = EMF_source, y = ROM-1, label = paste0("n =", n,", k =", k), vjust = -2, hjust = -1.1))

# adding average of 4 estimates
p2b <- p2b + geom_pointrange(data = mean_estimate_split, aes(group = estimate, y = ROM-1, ymin = LLCI-1, ymax = ULCI-1), linewidth = 1.3, position = position_nudge(x = 0.25))
p2b

# draw segment with arrow to shorten the very large CI Nr.4 of Base station (behavior)
p2b <- p2b + geom_segment(
  x = 2.13, y = 0.706-1,
  xend = 2.13, yend = 5.6-1,
  lineend = "round", linejoin = "round", linewidth = 0.5, size = 1, 
  arrow = arrow(length = unit(0.3, "cm")),
  colour = "red"
) 

# draw segment with arrow to shorten the very large CI Nr.5 (average) of Base station (behavior)
p2b <- p2b + geom_segment(
  x = 2.25, y = 0.372-1,
  xend = 2.25, yend = 5.6-1,
  lineend = "round", linejoin = "round", linewidth = 1.3, size = 1, 
  arrow = arrow(length = unit(0.3, "cm")),
  colour = "red"
) 

p2b


ggarrange(r1b, p2b, ncol = 1, nrow = 2, align="v",
          legend="none", common.legend = T)

ggsave("figures/Estimated_Toxicity_percent_base_split.jpg", width = 12, height = 6)
ggsave("figures/Estimated_Toxicity_percent_base_split.png", width = 12, height = 6)


# same plot but excluding metafor and brms estimates (so only clustered meta and bayesmeta estimates are shown, 
# as it should have been in the original publication)
# coding for estimate: meta = 0, rma (metafor) = 1, bayesmeta = 2, bayesmeta studylevel = 3, brms = 4

split_base_station2 <- subset(split_base_station, EMF_source  != "WiFi")
split_base_station2 <- subset(split_base_station2, estimate != 1)
split_base_station2 <- subset(split_base_station2, estimate != 4)
split_base_station2

p3 <- ggplot(split_base_station2, aes(x = EMF_source, y = ROM, col=EMF_source, group = estimate)) +  
  geom_pointrange(aes(group = estimate, y = ROM, ymin = LLCI, ymax = ULCI), position = position_dodge(width = 0.4)) +
  xlab(NULL) + labs(col="Field source", title="C") +
  scale_y_continuous(limits = c(0.37, 3.5), labels = waiver()) +  theme_classic2() +
  ylab(expression(paste("Pooled effect size (ROM)")))+
  theme(axis.title.y = element_text(size = 11)) +
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_hline(yintercept = 1.5, linetype="dotted")

p3

p3 <- p3 + scale_color_manual(values = color_table, drop = TRUE)
p3 <- p3 + geom_text(data = split_base_station2[split_base_station2$estimate == 0,], size = 3, aes(x = EMF_source, y = ROM, label = paste0("n =", n,", k =", k, " ", sig), vjust = -1, hjust = -0.5))

p3

# draw segment with arrow to shorten the very large CI Nr.1 of Base station (behavior)
p3 <- p3 + geom_segment(
  x = 1.9, y = 0.761,
  xend = 1.9, yend = 3.5,
  position = position_nudge(y = "1"),
  lineend = "round",
  linejoin = "round",
  linewidth = 0.5,
  size = 1, 
  arrow = arrow(length = unit(0.3, "cm")),
  colour = "red"
)

# draw segment with arrow to shorten the very large CI Nr.5 (average) of Base station (behavior)
p3 <- p3 + geom_segment(
  x = 2.1, y = 0.876,
  xend = 2.1, yend = 3.5,
  lineend = "round",
  linejoin = "round",
  linewidth = 0.5,
  size = 1, 
  arrow = arrow(length = unit(0.3, "cm")),
  colour = "red"
) 


ggarrange(r2, p3, ncol = 1, nrow = 2, align="v",
          legend="none", common.legend = T)

ggsave("figures/Estimated_Toxicity_like_original_publication.jpg", width = 12, height = 6)
ggsave("figures/Estimated_Toxicity_like_original_publication.png", width = 12, height = 6)



