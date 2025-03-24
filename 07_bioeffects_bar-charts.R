
# run this first

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
getwd() # working directory set to folder of .R file

if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(dplyr, tidyr, ggplot2, ggpubr, scales, RColorBrewer, openxlsx, sqldf)


############################################################
# Standardizing Bioeffects 
############################################################

df <- openxlsx::read.xlsx("tables/data_table_HFLF_4.xlsx", sheet = 1, na.strings = "NA") 
df <- relocate(df, experiment_id, .after = "experiment")

names(df)

old_bioeffects <- df$Bioeffects # save for comparison

df$Bioeffect_cat <- ""
df <- relocate(df, Bioeffect_cat, .after = "Bioeffects") # place new column behind "Bioeffects" column
#View(df)

df$Bioeffect_cat[grep("reduced reproductive", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("Reduced number of F1", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("F1 pupae per maternal fly", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("Reduced ovicity", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("ovarian cell death", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("Ovarian Cell Death", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("Decreased Ovarian Size", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("Increased ovarian apoptosis", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("Reduced Reproductive Capacity", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("Reduced egg laying", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("developmental success reduced", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("winter survival rate", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("survival rate for eggs", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("decrease in number of eggs", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("reduced egg laying", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("Reduced egg laying", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("Increased mortality", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("decreased egg to adult viability", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("Reduced reproductive capacity", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("Reduced number of offspring", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("Reduced brood", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("increased mortality in damp", df$Bioeffects)] <- "Reduced reproductive capacity"
df$Bioeffect_cat[grep("Reduced number of bees compared to control group", df$Bioeffects)] <- "Reduced reproductive capacity"


df$Bioeffect_cat[grep("ovarian DNA fragmentation", df$Bioeffects)] <- "DNA damage"
df$Bioeffect_cat[grep("Fragmentation", df$Bioeffects)] <- "DNA damage"
df$Bioeffect_cat[grep("DNA Damage", df$Bioeffects)] <- "DNA damage"
df$Bioeffect_cat[grep("mutation", df$Bioeffects)] <- "DNA damage" ###
df$Bioeffect_cat[grep("mutant", df$Bioeffects)] <- "DNA damage" ###
df$Bioeffect_cat[grep("mutagenic effects", df$Bioeffects)] <- "DNA damage" ###

df$Bioeffect_cat[grep("pupation time", df$Bioeffects)] <- "Developmental effects"
df$Bioeffect_cat[grep("Developmental Effects", df$Bioeffects)] <- "Developmental effects"
df$Bioeffect_cat[grep("developmental time", df$Bioeffects)] <- "Developmental effects"
df$Bioeffect_cat[grep("defects", df$Bioeffects)] <- "Developmental effects"
df$Bioeffect_cat[grep("developmental abnormalities", df$Bioeffects)] <- "Developmental effects"
df$Bioeffect_cat[grep("Increase in hatching time", df$Bioeffects)] <- "Developmental effects"
df$Bioeffect_cat[grep("apoptotic-like cells", df$Bioeffects)] <- "Developmental effects"
df$Bioeffect_cat[grep("lower developmental stability", df$Bioeffects)] <- "Developmental effects"
df$Bioeffect_cat[grep("Altered brain", df$Bioeffects)] <- "Developmental effects"
df$Bioeffect_cat[grep("Reduced nymphal gut mass", df$Bioeffects)] <- "Developmental effects"

df$Bioeffect_cat[grep("upregulation", df$Bioeffects)] <- "Altered DNA / transcription"
df$Bioeffect_cat[grep("downregulation", df$Bioeffects)] <- "Altered DNA / transcription"
df$Bioeffect_cat[grep("gene transcription", df$Bioeffects)] <- "Altered DNA / transcription"
df$Bioeffect_cat[grep("Gene transcription", df$Bioeffects)] <- "Altered DNA / transcription"
df$Bioeffect_cat[grep("altered DNA", df$Bioeffects)] <- "Altered DNA / transcription"
df$Bioeffect_cat[grep("Altered DNA", df$Bioeffects)] <- "Altered DNA / transcription"
df$Bioeffect_cat[grep("altered polyteny degree", df$Bioeffects)] <- "Altered DNA / transcription"

df$Bioeffect_cat[grep("Increased Reactive Oxygen Species", df$Bioeffects)] <- "Oxidative stress"
df$Bioeffect_cat[grep("ROS accumulation", df$Bioeffects)] <- "Oxidative stress"
df$Bioeffect_cat[grep("Decreased GST activity", df$Bioeffects)] <- "Oxidative stress"
df$Bioeffect_cat[grep("Increased SOD", df$Bioeffects)] <- "Oxidative stress"
df$Bioeffect_cat[grep("Increased SOD and CAT activity", df$Bioeffects)] <- "Oxidative stress"
df$Bioeffect_cat[grep("Glutathione S-transferase", df$Bioeffects)] <- "Oxidative stress"
df$Bioeffect_cat[grep("Superoxide Dismutase", df$Bioeffects)] <- "Oxidative stress"
df$Bioeffect_cat[grep("Superoxide dismutase", df$Bioeffects)] <- "Oxidative stress"
df$Bioeffect_cat[grep("TBARS", df$Bioeffects)] <- "Oxidative stress"
df$Bioeffect_cat[grep("Catalase", df$Bioeffects)] <- "Oxidative stress"


df$Bioeffect_cat[grep("Increased AchE", df$Bioeffects)] <- "Altered enzyme activity / metabolism"
df$Bioeffect_cat[grep("AchE activity", df$Bioeffects)] <- "Altered enzyme activity / metabolism"
df$Bioeffect_cat[grep("ALP activity", df$Bioeffects)] <- "Altered enzyme activity / metabolism"
df$Bioeffect_cat[grep("enzyme activity", df$Bioeffects)] <- "Altered enzyme activity / metabolism"
df$Bioeffect_cat[grep("proteases activity", df$Bioeffects)] <- "Altered enzyme activity / metabolism"
df$Bioeffect_cat[grep("metabolite levels", df$Bioeffects)] <- "Altered enzyme activity / metabolism"

df$Bioeffect_cat[grep("memory", df$Bioeffects)] <- "Impaired memory"
df$Bioeffect_cat[grep("conditioned response", df$Bioeffects)] <- "Impaired memory"
df$Bioeffect_cat[grep("learning response", df$Bioeffects)] <- "Impaired memory"
df$Bioeffect_cat[grep("proboscis extension reflex", df$Bioeffects)] <- "Impaired memory"

#View(df)
df$Bioeffect_cat[grep("jerking", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("movement", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("piping", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("level of activity", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("linear speed", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("angular speed", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("reduced mobility", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("Reduced mobility", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("locomotion frequency", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("latency to escape", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("walking distance", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("leaving", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("aggressiveness", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("agitation", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("Agitation", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("wingbeat frequency", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("wing beat frequency", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("circadian rhythm", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("magnetic compass", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("Impaired orientation", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("disturbed magnetic orientation", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("sleep duration", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("circadian period", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("move toward", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("avoidance behavior", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("Mosquitoes flee", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("Mosquitoes stop", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("Preference for lower", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("Preference for irradiated", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("neuronal activity", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("at higher power densities", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("Reduced number of returning bees", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("Reduced abundance of insects", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("Reduced abundance", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("Reduced rate of climbing", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("More flying insects caught", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("Reduced homing rate", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("Reduced percentage of returning", df$Bioeffects)] <- "Altered behavior"
df$Bioeffect_cat[grep("Increased percentage of returning", df$Bioeffects)] <- "Altered behavior"

df$Bioeffect_cat[grep("Increased egg laying", df$Bioeffects)] <- "Other"
df$Bioeffect_cat[grep("Increased number of offspring", df$Bioeffects)] <- "Other"

df$Bioeffect_cat[grep("no effect", df$Bioeffects)] <- "No effect"
df$Bioeffect_cat[grep("Not conclusive", df$Bioeffects)] <- "No effect"
df$Bioeffect_cat[grep("No effect", df$Bioeffects)] <- "No effect"
df$Bioeffect_cat[grep("no significant effect", df$Bioeffects)] <- "No effect"

# better not to remove non-significant studies, since they can be used in meta-analysis
# omitting those studies causes increased bias

table(df$Bioeffect_cat)
#View(df)

# for Vilic2024: assigning direction of effect for oxidative stress, and reassigning reduced oxidative stress to bioeffect category "other"
df <- mutate(df, Direction_of_effect = if_else ( Ratio_of_means < 1 & study == "Vilic2024", "beneficial", Direction_of_effect))
df <- mutate(df, Direction_of_effect = if_else ( Ratio_of_means > 1 & study == "Vilic2024", "detrimental", Direction_of_effect))
df <- mutate(df, Bioeffect_cat = if_else ( Ratio_of_means < 1 & study == "Vilic2024", "Other", Bioeffect_cat))

# for experiments, where bioeffect category was described, but no actual change was observed, classify as "No effect"
df <- mutate(df, Bioeffect_cat = if_else ( Percent.change.vs.control == "0%" & is.na(Direction_of_effect), "No effect", Bioeffect_cat))

df <- mutate(df, Bioeffect_cat = if_else ( is.na(Bioeffect_cat), "No effect", Bioeffect_cat))
df <- mutate(df, Direction_of_effect = if_else ( is.na(Direction_of_effect), "none", Direction_of_effect))

# for the left-over experiments that do not fall into any of the 9 categories, classify as "Other"
df <- mutate(df, Bioeffect_cat = if_else ( Bioeffect_cat== '', "Other", Bioeffect_cat))

# reassigning base station studies that found reduced or changed abundance as "uncertain"
# those were mostly classified as detrimental previously (that is, for the effect size plot figure 5A)
# the base station experiments will be analysed twice by meta-analysis, 
# either including or excluding the experiments that found reduced or changed abundance
df$Direction_of_effect[grep("Altered behavior", df$Bioeffect_cat)] <- "uncertain"

df$Direction_of_effect[grep("positive effects", df$Bioeffects)] <- "beneficial"

table(df$Bioeffect_cat)
table(df$Direction_of_effect)

names(df)

# to check if assigning of bioeffect category was correctly done
check_bioeffects <- cbind(old_bioeffects, df$Bioeffect_cat)
View(check_bioeffects)

# save expanded table
openxlsx::write.xlsx(df, "tables/data_table_HFLF_5.xlsx", rowNames = F, colWidths = 15, firstRow = T, firstCol = T)


#######################################################
# Pie- and Bar-Charts of Bioeffects by category: LF
#######################################################


pacman::p_load(dplyr, tidyr, ggplot2, ggpubr, scales, RColorBrewer, openxlsx, sqldf)

df <- openxlsx::read.xlsx("tables/data_table_HFLF_5.xlsx", sheet = 1, na.strings = "NA") 

#View(df)
names(df)

# select only LF experiments
df <- df[grep("LF", df$EMF_type),]

df$Bioeffect_cat <- as.factor(df$Bioeffect_cat)
levels(df$Bioeffect_cat)
length(df$Bioeffect_cat) # 139 LF experiments
sort(table(df$Bioeffect_cat), decreasing = F)

table_bioeffects <- df %>%
  group_by(Bioeffect_cat) %>% 
#  filter(n() > 2) %>% 
  summarise(n = n())

table_bioeffects

Other <- length(df$Bioeffect_cat) - sum(table_bioeffects$n)
Other # control to see if numbers match (should be 0)


amount <- table_bioeffects$n
labels <- table_bioeffects$Bioeffect_cat
amount <- as.vector(amount)
pie(amount, labels)


count.data <- data.frame(
  Bioeffects = labels,
  n = amount,
  prop = amount*100/sum(amount))

count.data <- count.data %>%
  arrange(desc(Bioeffects)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
count.data <- count.data %>% mutate(prop = round(prop, digits = 1))

count.data
maximum <- max(count.data$n)

p1 <- ggplot(count.data, aes(x = "", y = prop, fill = Bioeffects)) +
  geom_bar(width = 1, stat = "identity", color = "white", alpha=1) +
  geom_text(aes(y = lab.ypos, label = percent(signif(prop/100), digits = 1)), size=3) +
  coord_polar("y", start = 0) + theme_void() +
  labs(title = "LF-EMF", ) + theme(legend.position = "bottom") +
  scale_fill_brewer(palette="Set3")
p1

p1b <- ggplot(count.data, aes(x = n, y = reorder(Bioeffects, +n), fill = Bioeffects)) +
  geom_bar(stat = "identity", alpha=1) +
  geom_text(aes(x = 3+n/2, label = percent(prop/100)), size=3) +
  labs(title = "A", subtitle = "LF-EMF", y = "", x = "Number of experiments") +
  scale_fill_brewer(palette="Set3") + theme_classic() + theme(legend.position = "none")+
  geom_vline(xintercept = c(maximum), linetype = 3)
p1b
#ggsave("figures/barchart_bioeffect_categories_LF.png", width = 4, height = 4)

leg <- get_legend(p1)


#######################################################
# Pie- and Bar-Charts of Bioeffects by category: HF
#######################################################

df <- openxlsx::read.xlsx("tables/data_table_HFLF_5.xlsx", sheet = 1, na.strings = "NA") 

# select only HF experiments
df <- df[grep("HF", df$EMF_type),]

names(df)

df$Bioeffect_cat <- as.factor(df$Bioeffect_cat)
levels(df$Bioeffect_cat)
length(df$Bioeffect_cat) # 348 HF experiments
sort(table(df$Bioeffect_cat), decreasing = F)


table_bioeffects <- df %>%
  group_by(Bioeffect_cat) %>% 
  #filter(n() > 2) %>% 
  summarise(n = n())
table_bioeffects

amount <- table_bioeffects$n
labels <- table_bioeffects$Bioeffect_cat
amount <- as.vector(amount)
pie(amount,labels)


count.data2 <- data.frame(
  Bioeffects = labels,
  n = amount,
  prop = amount*100/sum(amount))

count.data2 <- count.data2 %>%
  arrange(desc(Bioeffects)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
count.data2 <- count.data2 %>% mutate(prop = round(prop, digits = 1))

count.data2
maximum2 <- max(count.data2$n) 


p2 <- ggplot(count.data2, aes(x = "", y = prop, fill = Bioeffects)) +
  geom_bar(width = 1, stat = "identity", color = "white", alpha=1) +
  geom_text(aes(y = lab.ypos, label = percent(prop/100)), size=3) +
  coord_polar("y", start = 0) +
  labs(title = "HF-EMF", ) +
  scale_fill_brewer(palette="Set3") +
  theme_void()
p2

p2b <- ggplot(count.data2, aes(x = n, y = reorder(Bioeffects, +n), fill = Bioeffects)) +
  geom_bar(stat = "identity", alpha=1) + 
  geom_text(aes(x = 6+n/2, label = percent(prop/100)), size=3) +
  labs(title = "B", subtitle = "HF-EMF", y = "", x = "Number of experiments") +
  scale_fill_brewer(palette="Set3") + theme_classic() + theme(legend.position = "none")+
  geom_vline(xintercept = c(maximum2), linetype = 3)
p2b
#ggsave("figures/barchart_bioeffect_categories_HF.png", width = 4, height = 4)


ggarrange(p1, p2, ncol = 2, nrow = 1, align = "hv", legend.grob = leg, legend = "bottom")
#ggsave("HFLF_Bioeffect_categories_pie-charts.jpg", width = 8, height = 4)


ggarrange(p1b, p2b, ncol = 2, nrow = 1, align = "hv", legend = "none")
ggsave("figures/barcharts_HFLF_Bioeffect_categories.png", width = 8, height = 4)
ggsave("figures/barcharts_HFLF_Bioeffect_categories.jpg", width = 8, height = 4)

ggarrange(p1b, p2b, ncol = 1, nrow = 2, align = "v", legend = "none")
ggsave("figures/barcharts_HFLF_Bioeffect_categories(vertical).png", width = 4, height = 8)
ggsave("figures/barcharts_HFLF_Bioeffect_categories(vertical).jpg", width = 4, height = 8)



#######################################################
# Pie-charts of direction of effect
#######################################################

pacman::p_load(dplyr, tidyr, ggplot2, ggpubr, scales, RColorBrewer, openxlsx, sqldf)

df <- openxlsx::read.xlsx("tables/data_table_HFLF_5.xlsx", sheet = 1, na.strings = "NA") 

df$Direction_of_effect <- as.character(df$Direction_of_effect)
table(df$Direction_of_effect)

df$Direction_of_effect <- as.factor(df$Direction_of_effect)
HF <- df[grep("HF", df$EMF_type),]
LF <- df[grep("LF", df$EMF_type),]

table(df$Direction_of_effect)
table(LF$Direction_of_effect)
table(HF$Direction_of_effect)

amount <- summary(LF$Direction_of_effect)
labels <- names(amount)
labels
amount <- as.vector(amount)
amount
pie(amount,labels)

LF_data <- data.frame("labels" = labels, "prop" = amount*100/sum(amount))
LF_data$prop <- round(LF_data$prop, digits = 1)
LF_data

p1 <- ggplot(LF_data, aes(x="", y=prop, fill=labels)) +
  geom_bar(stat="identity", width=1, alpha=3/4) +
  coord_polar("y", start=0) +
  geom_text(aes(label = paste0(prop, "%")), position = position_stack(vjust=0.5)) +
  theme_void() +  scale_fill_brewer(palette="Dark2") +
  labs(title ="A", subtitle = "LF-EMF", fill='Effect') + theme(legend.position = "bottom")
p1

amount <- summary(HF$Direction_of_effect)
labels <- names(amount)
amount <- as.vector(amount)
pie(amount,labels)

HF_data <- data.frame("labels" = labels, "prop" = amount*100/sum(amount))
HF_data$prop <- round(HF_data$prop, digits = 1)

p2 <- ggplot(HF_data, aes(x="", y=prop, fill=labels)) +
  geom_bar(stat="identity", width=1, alpha=3/4) +
  coord_polar("y", start=0) +
  geom_text(aes(label = paste0(prop, "%")), position = position_stack(vjust=0.5)) +
  theme_void() +  scale_fill_brewer(palette="Dark2") +
  labs(title ="B", subtitle = "HF-EMF", fill='Effect') + theme(legend.position = "bottom")
p2

ggarrange(p1, p2, ncol = 2, nrow = 1, align = "hv", legend = "bottom", common.legend = T, hjust = 0)
ggsave("figures/piecharts_effects.jpg", width = 8, height = 4)
ggsave("figures/piecharts_effects.png", width = 8, height = 4)

ggarrange(p1, p2, ncol = 1, nrow = 2, align = "hv", legend = "bottom", common.legend = T, hjust = 0)
ggsave("figures/piecharts_effects(vertical).jpg", width = 4, height = 8)
ggsave("figures/piecharts_effects(vertical).png", width = 4, height = 8)


sort(table(LF$Direction_of_effect), decreasing = T)
sort(table(HF$Direction_of_effect), decreasing = T)
length(LF$Direction_of_effect)
length(HF$Direction_of_effect)

