
# if "pacman" doesn't automatically install properly, run this first!
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("survcomp")
#install.packages("remotes")
#remotes::install_github("MathiasHarrer/dmetar")


# run this first!
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

dir.create(file.path(getwd(), "suppl_figures")) # create subfolder for saving supplemental figures

#rm(list=ls())  # clear global environment of stuff from previous sessions


pacman::p_load(dplyr, tidyr, ggplot2, ggpubr, scales, RColorBrewer, openxlsx, sqldf, stringr,
               dmetar, metafor, meta, bayesmeta, RoBMA, brms, tidyverse, tidybayes, ggridges)


# function to make a study-level (bayesian) forest plot (credit to Matti Vuorre)

gg_forestplot <- function(brm_obj,
                          ci = 0.95,
                          ...) {
  # Study-specific effects are deviations + average
  out_r <- spread_draws(brm_obj, r_study[study,term], b_Intercept) %>% 
    mutate(b_Intercept = r_study + b_Intercept) 
  # Average effect
  out_f <- spread_draws(brm_obj, b_Intercept) %>% 
    mutate(study = "Average")
  # Combine average and study-specific effects' data frames
  out_all <- bind_rows(out_r, out_f) %>% 
    ungroup() %>%
    # Ensure that Average effect is on the bottom of the forest plot
    mutate(study = fct_relevel(study, "Average"))
  # Data frame of summary numbers
  out_all_sum <- group_by(out_all, study) %>% mean_qi(b_Intercept)
  #out_all_sum$b_Intercept <- exp(out_all_sum$b_Intercept)
  out_all_sum$.lower <- exp(out_all_sum$.lower)
  out_all_sum$.upper <- exp(out_all_sum$.upper)
  # Draw plot
  out_all %>%   
    ggplot(aes(exp(b_Intercept), study)) +
    geom_density_ridges(
      rel_min_height = 0.01, 
      col = NA,
      scale = 1
    ) +
    geom_pointinterval(data = out_all_sum, size = 1, xmin = out_all_sum$.lower, xmax = out_all_sum$.upper) +
    theme_classic2() +
    geom_vline(xintercept = 1, linetype="dashed") +
    geom_vline(xintercept = 1.5, linetype="dotted") +
    xlab("Effect size") +
    geom_text(
      data = mutate_if(out_all_sum, is.numeric, round, 2),
      # Use glue package to combine strings
      aes(label = glue::glue("{round(exp(b_Intercept), digits=2)} [{.lower}, {.upper}]"), x = Inf),
      hjust = "inward"
    )
}


####################################################################
# Meta-analysis of Drosophila reproductive toxicity studies
####################################################################

df <- read.xlsx("tables/HFLF_meta_table.xlsx", sheet = 1, na.strings = "NA") 
names(df)

sort(table(df$EMF_source), decreasing = F)
sort(table(df$Insect_type), decreasing = F)
table(df$Bioeffect_cat)

# Retain only experiments (within studies) that have effect size and standard error values
df <- as_tibble(df)
df <- df[complete.cases(df[,c("log_ROM","log_SE")]),]

# Retain only HF studies
df <- subset(df, EMF_type == "HF")

# Retain only Drosophila studies
df <- subset(df, Insect_type == "Drosophila")
table(df$SignalGen)
table(df$Bioeffect_cat)

# Retain only studies of reproductive effects
df <- subset(df, Bioeffect_cat == "Reduced reproductive capacity" |
               Bioeffect_cat == "Other" |
               Bioeffect_cat == "No effect")

# the filtering is correct, except for one study - due to the imprecise classifier "Other"
df <- subset(df, study != "Manta2013")

str(df)

# Subdivision of all data points into high, medium, and low E-field groups.
droso_high <- subset(df, E_Field > 7)
droso_mid <- subset(df, E_Field > 2.2 & E_Field < 7)
droso_low <- subset(df, E_Field < 2.2)


# Meta-analysis of all HF Drosophila studies over 7 V/m

# Frequentist models
# Package "meta", three-level, i.e. clustered model

madata1 <- metagen(TE = log_ROM,
                   seTE = log_SE,
                   data = droso_high,
                   cluster = study,
                   #subgroup = study,
                   comb.fixed = T,
                   comb.random = T,
                   hakn = FALSE,
                   prediction = TRUE,
                   sm = "ROM")
madata1
meta::forest(madata1)

# nicer
meta::forest(madata1, layout = "RevMan5",
       label.right = "detrimental", col.label.right = "red",
       label.left = "beneficial", col.label.left = "green",
       prediction = TRUE, test.overall.random = TRUE,
       file = "suppl_figures/forestplot_repro_tox_over_7Vm.png", width = 3500, rows.gr = 150, args.gr = list(res = 300), dev.off = T,
       leftcols = c("Author_Year","EMF_source","E_Field","effect","ci","w.random"),
       leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
       just = "center", just.addcols = "center", xlim = c(0.75,4))

madata1 <- update(madata1, subgroup = study) # show subgroup means

meta::forest(madata1, layout = "RevMan5",
            label.right = "detrimental", col.label.right = "red",
            label.left = "beneficial", col.label.left = "green",
            prediction = TRUE, test.overall.random = TRUE,
            print.subgroup.labels = F, print.subgroup.name = F,
            file = "suppl_figures/forestplot_repro_tox_over_7Vm_subgroups.png", width = 3500, rows.gr = 200, args.gr = list(res = 300), dev.off = T,
            leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
            leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
            just = "center", just.addcols = "center", xlim = c(0.75,4))


meta::funnel(madata1)
radial(madata1)

exp(madata1$TE.random) # ROM = "ratio of means" effect size
percent((1/exp(madata1$TE.random))-1) # Effect size as percentage change compared to control
madata1$pval.random

meta.est_high <- c(madata1$TE.random, madata1$lower.random, madata1$upper.random)
meta.est_high <- meta:::backtransf(meta.est_high, sm="ROM")
meta.est_high


# three-level model using the "Metafor" package (or 4-level model, depending on counting)

full.model <- rma.mv(yi, vi, random = ~ 1 | study/experiment, data = droso_high, slab = Author_Year)
#full.model <- rma.mv(yi, vi, random = list(~ 1 | experiment, ~ 1 | study), data = droso_high) # same
exp(c(full.model$b,full.model$ci.lb,full.model$ci.ub)) # About same effect size and CI as the clustered "meta" model above

meta::forest(full.model, atransf = exp, at = log(c(0.5, 1, 4)), xlim = c(-1, 2))

funnel.rma(full.model)

rma.est_high <- exp(c(full.model$b,full.model$ci.lb,full.model$ci.ub))
rma.est_high

full.model2 <- rma.mv(yi, vi, random = ~ 1 | study, data = droso_high)
exp(c(full.model2$b,full.model2$ci.lb,full.model2$ci.ub)) # nearly the same
full.model2 <- rma.mv(yi, vi, random = ~ 1 | experiment_id, data = droso_high) 
exp(c(full.model2$b,full.model2$ci.lb,full.model2$ci.ub)) # not the same


# Bayesian models

# "brms" model: correct form of LME4 formula for meta-analysis at study-level and within-study "experiment" subgroups
# ~ 1 + (1|study) + (1|experiment) is the same as ~ 1 + (1|study/experiment)

# takes a long time to run
brm1 <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study) + (1|experiment), data = droso_high,
            #tau.prior=function(t){dhalfnormal(t, scale=0.5)},
            iter = 40000, warmup = 1000,
            control = list(adapt_delta = 0.99),
            core = 6, chains = 6)

summary(brm1)
exp(fixef(brm1))
# old estimate 1.402960 1.233168 1.601654
# new estimate 1.342872 1.168114 1.552806

plot(brm1, ask = FALSE)
loo(brm1)
get_variables(brm1)

gg_forestplot(brm1) + xlim(0.5, 2.5)

ggsave("suppl_figures/brms_forestplot_repro_tox_over_7Vm.png", width = 12, height = 6)


# simple "bayesmeta" model
bma1 <- bayesmeta(y = droso_high$log_ROM, sigma = droso_high$log_SE,
                  labels = droso_high$study, mu.prior.mean = 0, mu.prior.sd = 4,
                  tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

bma1

forestplot(bma1)
exp(bma1$summary[,2])

png(file = "suppl_figures/bayesmeta_forestplot_repro_tox_over_7Vm.png", width = 3500, height = 2600, res = 300)
forestplot(bma1, expo=TRUE, clip=c(0.5,2.5))
dev.off()


# study-level "bayesmeta" model
X_h <- model.matrix( ~ -1 + study, data = droso_high)
X_h

bmr1 <- bmr(y = droso_high$log_ROM, sigma = droso_high$log_SE, X=X_h,
                  labels = droso_high$study, mu.prior.mean = 0, mu.prior.sd = 4,
                  tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

# transform meta results mean + CI into mean + SE
# Standard Error (SE) = Standard Deviation (SD) / √(Sample Size)
# For 95% CI: SE = (Upper Limit - Lower Limit) / 3.92
bmr1_means <- bmr1$summary[3,]
bmr1_means <- bmr1_means[-1]
bmr1_ses <- (bmr1$summary[6,]-bmr1$summary[5,])/3.92
bmr1_ses <- bmr1_ses[-1]

bmr1_studylevel <- bayesmeta(y = bmr1_means, sigma = bmr1_ses,
            labels = unique(droso_high$study), mu.prior.mean = 0, mu.prior.sd = 4,
            tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

exp(bmr1_studylevel$summary)

forestplot(bmr1_studylevel, expo=TRUE, clip=c(0.5,2.5))


png(file = "suppl_figures/bayesmeta_study-level_forestplot_repro_tox_over_7Vm.png", width = 3500, height = 2600, res = 300)
forestplot(bmr1_studylevel, expo=TRUE, clip=c(0.5,2.5))
dev.off()

bmr1_studylevel$likelihood(mu=0)
bmr1_studylevel$bayesfactor[1, "mu=0"]
bmr1_studylevel$bayesfactor

# extract estimates
brm.est_high <- c(fixef(brm1)[1],fixef(brm1)[3],fixef(brm1)[4])
brm.est_high <- meta:::backtransf(brm.est_high, sm="ROM")
bayes.est_high <- c(bma1$summary["mode",2],bma1$summary["95% lower",2],bma1$summary["95% upper",2])
bayes.est_high <- meta:::backtransf(bayes.est_high, sm="ROM")
bayes.est_high2 <- exp(bmr1_studylevel$summary[c(3,5,6),2])

# list all estimates
brm.est_high
bayes.est_high
bayes.est_high2
meta.est_high
rma.est_high


# Meta-analysis of all RF Drosophila studies between 2 and 7 V/m.
# "meta" model, clustered at study-level
madata2 <- metagen(TE = log_ROM,
                   seTE = log_SE,
                   data = droso_mid,
                   cluster = study,
                   #subgroup = study,
                   comb.fixed = TRUE,
                   comb.random = T,
                   #method.tau = "DL",
                   hakn = FALSE,
                   prediction = TRUE,
                   sm = "ROM")
madata2

exp(madata2$TE.random) # ROM = "ratio of means" Effect size
percent((1/exp(madata2$TE.random))-1) # Effect size as percentage change over control
madata2$pval.random

meta::forest(madata2)

meta::forest(madata2, layout = "RevMan5",
            label.right = "detrimental", col.label.right = "red",
            label.left = "beneficial", col.label.left = "green",
            prediction = TRUE, test.overall.random = TRUE,
            file = "suppl_figures/forestplot_repro_tox_over_2Vm_under_7Vm.png", width = 3500, rows.gr = 150, args.gr = list(res = 300), dev.off = T,
            leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
            leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
            just = "center", just.addcols = "center")


madata2 <- update(madata2, subgroup = study) # subgroups shown

meta::forest(madata2, layout = "RevMan5",
             label.right = "detrimental", col.label.right = "red",
             label.left = "beneficial", col.label.left = "green",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = F, print.subgroup.name = F,
             file = "suppl_figures/forestplot_repro_tox_over_2Vm_under_7Vm_subgroups.png", width = 3500, rows.gr = 200, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
             leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
             just = "center", just.addcols = "center")

meta::funnel(madata2)

meta.est_mid <- c(madata2$TE.random, madata2$lower.random, madata2$upper.random)
meta.est_mid <- meta:::backtransf(meta.est_mid, sm="ROM")
meta.est_mid


# "metafor" model
full.model2 <- rma.mv(yi, vi, random = ~ 1 | study/experiment, data = droso_mid, slab = Author_Year) 
exp(c(full.model2$b,full.model2$ci.lb,full.model2$ci.ub)) # basically the same as study-level "meta" model

droso_mid$experiment

meta::forest(full.model2, atransf = exp, at = log(c(0.5, 1, 4)), xlim = c(-1, 2))

funnel.rma(full.model2)

rma.est_mid <- exp(c(full.model2$b,full.model2$ci.lb,full.model2$ci.ub))
rma.est_mid



# Bayesian models

# "brms" model
brm2 <- brm(formula = yi | se(sei) ~ 1 + (1|study) + (1|experiment), data = droso_mid,
             iter = 40000, warmup = 1000,
             control = list(adapt_delta = 0.99),
             core = 6, chains = 6)

summary(brm2)
exp(fixef(brm2))

gg_forestplot(brm2) + xlim(0.5, 3)
ggsave("suppl_figures/brms_forestplot_repro_tox_over_2Vm_under_7Vm.png", width = 12, height = 6)


# "bayesmeta" model
bma2 <- bayesmeta(y = droso_mid$log_ROM, sigma = droso_mid$log_SE,
                 labels = droso_mid$study, mu.prior.mean = 0, mu.prior.sd = 4,
                 tau.prior = function(t) dhalfnormal(t, scale = 0.25))

exp(ma2$summary)
forestplot(ma2, expo=TRUE, clip=c(0.5,3))

png(file = "suppl_figures/bayesmeta_forestplot_repro_tox_over_2Vm_under_7Vm.png", width = 3500, height = 2600, res = 300)
forestplot(ma2, expo=TRUE, clip=c(0.5,3))
dev.off()


# "bayesmeta" study-level  model
X_m <- model.matrix(~ -1 + study, data = droso_mid)

bmr2 <- bmr(y = droso_mid$log_ROM, sigma = droso_mid$log_SE, X=X_m,
            labels = droso_mid$study, mu.prior.mean = 0, mu.prior.sd = 4,
            tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

bmr2_means <- bmr2$summary[3,]
bmr2_means <- bmr2_means[-1]
bmr2_ses <- (bmr2$summary[6,]-bmr2$summary[5,])/3.92
bmr2_ses <- bmr2_ses[-1]

bmr2_studylevel <- bayesmeta(y = bmr2_means, sigma = bmr2_ses,
                             labels = unique(droso_mid$study), mu.prior.mean = 0, mu.prior.sd = 4,
                             tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

exp(bmr2_studylevel$summary)
forestplot(bmr2_studylevel, expo=TRUE, clip=c(-1, 3))

png(file = "suppl_figures/bayesmeta_study-level_forestplot_repro_tox_over_2Vm_under_7Vm.png", width = 3300, height = 2800, res = 300)
forestplot(bmr2_studylevel, expo=TRUE, clip=c(-1, 3))
dev.off()

bmr2_studylevel$likelihood(mu=0)
bmr2_studylevel$bayesfactor[1, "mu=0"]


brm.est_mid <- c(fixef(brm2)[1],fixef(brm2)[3],fixef(brm2)[4])
brm.est_mid <- meta:::backtransf(brm.est_mid, sm="ROM")
bayes.est_mid <- c(ma2$summary["mode",2],ma2$summary["95% lower",2],ma2$summary["95% upper",2])
bayes.est_mid <- meta:::backtransf(bayes.est_mid, sm="ROM")
bayes.est_mid2 <- exp(bmr2_studylevel$summary[c(3,5,6),2])

# all estimates
brm.est_mid
bayes.est_mid
bayes.est_mid2
meta.est_mid
rma.est_mid


# Meta-analysis of all RF Drosophila studies below 2 V/m
# "meta" model
madata3 <- metagen(TE = log_ROM,
                   seTE = log_SE,
                   data = droso_low,
                   cluster = study,
                   comb.fixed = TRUE,
                   comb.random = T,
                   hakn = FALSE,
                   prediction = TRUE,
                   sm = "ROM")

#View(df3)
meta::forest(madata3)

exp(madata3$TE.random)
percent((1/exp(madata3$TE.random))-1) # Effect size as percentage change over control
madata3$pval.random


meta::forest(madata3, layout = "RevMan5",
            label.right = "detrimental", col.label.right = "red",
            label.left = "beneficial", col.label.left = "green",
            prediction = TRUE, test.overall.random = TRUE,
            file = "suppl_figures/forestplot_repro_tox_under_2Vm.png", width = 3300, rows.gr = 150, args.gr = list(res = 300), dev.off = T,
            leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
            leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
            just = "center", just.addcols = "center")

madata3 <- update(madata3, subgroup = study) # subgroups shown

meta::forest(madata3, layout = "RevMan5",
             label.right = "detrimental", col.label.right = "red",
             label.left = "beneficial", col.label.left = "green",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = F, print.subgroup.name = F,
             file = "suppl_figures/forestplot_repro_tox_under_2Vm_subgroups.png", width = 3300, rows.gr = 150, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
             leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
             just = "center", just.addcols = "center")

meta::funnel(madata3)
meta.est_low <- c(madata3$TE.random, madata3$lower.random, madata3$upper.random)
meta.est_low <- meta:::backtransf(meta.est_low, sm="ROM")


# "metafor" model
full.model3 <- rma.mv(yi, vi, random = ~ 1 | study/experiment, data = droso_low, slab = Author_Year) 
exp(c(full.model3$b,full.model3$ci.lb,full.model3$ci.ub)) # basically the same as study-level "meta" model

meta::forest(full.model3, atransf = exp, at = log(c(0.5, 1, 4)), xlim = c(-1, 2))

rma.est_low <- exp(c(full.model3$b,full.model3$ci.lb,full.model3$ci.ub))
rma.est_low


# Bayesian models

# "brms" model
brm3 <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study) + (1|experiment), data = droso_low,
             iter = 40000, warmup = 1000,
             control = list(adapt_delta = 0.99),
             core = 6, chains = 6)

exp(fixef(brm3))

gg_forestplot(brm3) + xlim(0.5, 2.5)
ggsave("suppl_figures/brms_forestplot_repro_tox_under_2Vm.png", width = 12, height = 6)


# "bayesmeta" model
bma3 <- bayesmeta(y = droso_low$log_ROM, sigma = droso_low$log_SE,
                 labels = droso_low$study, mu.prior.mean = 0, mu.prior.sd = 4,
                 tau.prior = function(t) dhalfnormal(t, scale = 0.5),
)

forestplot(ma3, expo=TRUE, clip=c(-1, 2))

png(file = "suppl_figures/bayesmeta_forestplot_repro_tox_under_2Vm.png", width = 3100, height = 2200, res = 300)
forestplot(ma3, expo=TRUE, clip=c(-1, 2))
dev.off()


# "bayesmeta" study-level model

X_l <- model.matrix(~ -1 + droso_low$study, data = droso_low)
X_l

bmr3 <- bmr(y = droso_low$log_ROM, sigma = droso_low$log_SE, X=X_l,
                 labels = droso_low$study, mu.prior.mean = 0, mu.prior.sd = 4,
                 tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

bmr3_means <- bmr3$summary[3,]
bmr3_means <- bmr3_means[-1]
bmr3_ses <- (bmr3$summary[6,]-bmr3$summary[5,])/3.92
bmr3_ses <- bmr3_ses[-1]

bmr3_studylevel <- bayesmeta(y = bmr3_means, sigma = bmr3_ses,
                             labels = unique(droso_low$study), mu.prior.mean = 0, mu.prior.sd = 4,
                             tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

exp(bmr3_studylevel$summary)
forestplot(bmr3_studylevel, expo=TRUE, clip=c(-1, 3))


png(file = "suppl_figures/bayesmeta_study-level_forestplot_repro_tox_under_2Vm.png", width = 3300, height = 2800, res = 300)
forestplot(bmr3_studylevel, expo=TRUE, clip=c(-1, 3))
dev.off()

bmr3_studylevel$likelihood(mu=0)
bmr3_studylevel$bayesfactor[1, "mu=0"]

brm.est_low <- c(fixef(brm3)[1],fixef(brm3)[3],fixef(brm3)[4])
brm.est_low <- meta:::backtransf(brm.est_low, sm="ROM")
bayes.est_low <- c(ma3$summary["mode",2],ma3$summary["95% lower",2],ma3$summary["95% upper",2])
bayes.est_low <- meta:::backtransf(bayes.est_low, sm="ROM")
bayes.est_low2 <- exp(bmr3_studylevel$summary[c(3,5,6),2])

# all estimates
brm.est_low
bayes.est_low
bayes.est_low2
meta.est_low
rma.est_low


####################################################################
#  Estimates of effect size by EMF source based on meta-analysis
####################################################################

df <- read.xlsx("tables/HFLF_meta_table.xlsx", sheet = 1, na.strings = "NA") 

table(df$EMF_source)
df <- as_tibble(df)
df <- df[complete.cases(df[,c("log_ROM","log_SE")]),]

table(df$EMF_source)
table(df$Bioeffect_cat)

df2 <- subset(df, EMF_source == "Base station" |
                EMF_source == "DECT" |
                EMF_source == "Mobile phone" |
                EMF_source == "RF signal generator" |
                EMF_source == "WiFi" |
                EMF_source == "Coil system 50/60 Hz"
              )

#View(df2)
table(df2$EMF_source)
str(df)


tower <- sqldf("SELECT * FROM df WHERE EMF_source LIKE '%Base station%'")
DECT <- sqldf("SELECT * FROM df WHERE EMF_source LIKE '%DECT%'")
mobile <- sqldf("SELECT * FROM df WHERE EMF_source LIKE '%Mobile phone%'")
sig_gen <- sqldf("SELECT * FROM df WHERE EMF_source LIKE '%Signal Generator%'")
coilsystem <- sqldf("SELECT * FROM df WHERE EMF_source LIKE '%Coil system%'")
wifi <- sqldf("SELECT * FROM df WHERE EMF_source LIKE '%WiFi%'")


# subgroup of all studies that found direct evidence of toxicity
df_basestation_tox <- subset(tower, Bioeffect_cat != "Altered behavior" & Bioeffect_cat != "No effect")
tower2 <- sqldf("SELECT * FROM df_basestation_tox WHERE EMF_source LIKE '%Base station%'")

# subgroup of all studies that found altered behavior
df_basestation_alt_behav <- subset(df, Bioeffect_cat == "Altered behavior" | Bioeffect_cat == "No effect")
tower3 <- sqldf("SELECT * FROM df_basestation_alt_behav WHERE EMF_source LIKE '%Base station%'")

#View(tower)
#View(tower2)
#View(tower3)
#View(DECT)
#View(coilsystem)
#View(wifi)


######################################################################
# Meta-analysis of all EMF groups together: faster but less accurate
######################################################################

madata_all <- metagen(TE = log_ROM,
                      seTE = log_SE,
                      data = df2,
                      cluster = study,
                      comb.fixed = TRUE,
                      comb.random = T,
                      #method.tau = "DL",
                      hakn = FALSE,
                      prediction = TRUE,
                      sm = "ROM")
#forest(madata)

meta::funnel(madata_all)
# the observed outcomes (experiment effect sizes) in two "rays" on the right side are those derived from "se.from.p" (p = 0.5 or p = 0.05)
# this somewhat reduces the explanatory power of this funnel plot, hence no judgement on study bias will be made

res0 <- rma(TE, sei=seTE, data=madata_all)
regtest(res0)
funnel.rma(res0, label = "out")

results <- update(madata_all, subgroup = EMF_source, tau.common = TRUE, common = FALSE)

summary(results)
# See "Results for subgroups"

# k    ROM           95%-CI  tau^2    tau       Q   I^2
# EMF_source = Base station          51 1.4379 [1.1212; 1.8443] 0.1604 0.4005  489.51 89.8%
# EMF_source = Mobile phone          79 1.5719 [1.3334; 1.8531] 0.1604 0.4005 3045.01 97.4%
# EMF_source = RF signal generator   38 1.2872 [1.0165; 1.6299] 0.1604 0.4005 1231.40 97.0%
# EMF_source = DECT                  20 1.4776 [1.0928; 1.9978] 0.1604 0.4005  478.00 96.0%
# EMF_source = Coil system 50/60 Hz  39 1.3533 [1.1349; 1.6137] 0.1604 0.4005  122.80 69.1%

# with new data added post publication of the review
# k    ROM           95%-CI  tau^2    tau        Q   I^2
# EMF_source = Base station         126 1.3973 [1.1287; 1.7298] 0.1724 0.4152 10594.99 98.8%
# EMF_source = RF signal generator   42 1.2318 [0.9807; 1.5473] 0.1724 0.4152  1235.88 96.7%
# EMF_source = Mobile phone          79 1.5367 [1.2979; 1.8195] 0.1724 0.4152  3046.87 97.4%
# EMF_source = DECT                  20 1.4414 [1.0608; 1.9585] 0.1724 0.4152   478.37 96.0%
# EMF_source = WiFi                   9 1.1784 [0.8807; 1.5767] 0.1724 0.4152   185.38 95.7%
# EMF_source = Coil system 50/60 Hz  39 1.1731 [0.9787; 1.4061] 0.1724 0.4152   127.06 70.1%

#################################################################
# Meta analysis of all EMF groups individually and sequentially
#################################################################
###########################################
# Base station (all experimental findings)
###########################################

m1 <- metagen(TE = log_ROM,
              seTE = log_SE,
              data = tower,
              cluster = study,
              comb.fixed = TRUE,
              comb.random = T,
              hakn = FALSE,
              prediction = TRUE,
              sm = "ROM")
#meta::forest(m1)

meta::forest(m1, layout = "RevMan5",
            label.right = "detrimental", col.label.right = "red",
            label.left = "beneficial", col.label.left = "green",
            prediction = TRUE, test.overall.random = TRUE,
            file = "suppl_figures/forestplot_base_station.png", width = 3300, rows.gr = 430, args.gr = list(res = 300), dev.off = T,
            leftcols = c("Author_Year","E_Field","CummHrs","Bioeffect_cat","effect","ci","w.random"),
            leftlabs = c("Author, Year","E-Field [V/m]","Duration [h]","Bioeffect category"),
            just = "center", just.addcols = "center", xlim = c(0.2,10))


m1 <- update(m1, subgroup = study)

meta::forest(m1, layout = "RevMan5",
             label.right = "detrimental", col.label.right = "red",
             label.left = "beneficial", col.label.left = "green",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = F, print.subgroup.name = F,
             #file = "suppl_figures/forestplot_base_station_subgroups.png", width = 3300, rows.gr = 520, args.gr = list(res = 300), dev.off = T,
             file = "suppl_figures/forestplot_base_station_subgroups.png", width = 3300, rows.gr = 250, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
             leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
             just = "center", just.addcols = "center", xlim = c(0.2,10))

summary(m1)
meta.est1 <- c(m1$TE.random, m1$lower.random, m1$upper.random)
meta.est1 <- meta:::backtransf(meta.est1, sm="ROM")
meta.est1
percent((1/exp(m1$TE.random))-1) # Effect size as percentage change over control
m1$pval.random

meta::funnel(m1, label = "out")

res1 <- rma(TE, sei=seTE, data=m1)
inf <- influence(res1)
plot(inf)
regtest(res1)
funnel.rma(res1, label = "out")
funnel.rma(res1, label = "out", refline=0) # centered on effect size = 0
res1


full.model1 <- rma.mv(yi, vi, random = ~ 1 | study/experiment, data = tower)
#full.model1 <- rma.mv(yi, vi, random = ~ 1 | study, data = tower) # gives basically the same mean + CI
exp(c(full.model1$b, full.model1$ci.lb, full.model1$ci.ub))
rma.est1 <- exp(c(full.model1$b, full.model1$ci.lb, full.model1$ci.ub))
rma.est1.pval <- full.model1$pval

percent((1/exp(full.model1$b))-1)
meta::forest(full.model1)


# Bayesian models

# "brms" model
brm_m1 <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study/experiment), data = tower,
             iter = 30000, warmup = 1000,
             control = list(adapt_delta = 0.98),
             core = 6, chains = 6)
exp(fixef(brm_m1))

gg_forestplot(brm_m1) + xlim(0.5, 6)
ggsave("suppl_figures/brms_forestplot_base_station.png", width = 12, height = 6)


# "bayesmeta" model
bm1 <- bayesmeta(y = tower$log_ROM, sigma = tower$log_SE,
                 labels = tower$study, mu.prior.mean = 0, mu.prior.sd = 4,
                 tau.prior = function(t) dhalfnormal(t, scale = 0.5))
forestplot(bm1, expo=TRUE, clip=c(-1, 3))

png(file = "suppl_figures/bayesmeta_forestplot_base_station.png", width = 3500, height = 3700, res = 300)
forestplot(bm1, expo=TRUE, clip=c(-1, 3))
dev.off()

exp(bm1$summary)
bm1$bayesfactor
bm1$bayesfactor[1, "mu=0"]


# "bayesmeta" model at study-level
X_tower <- model.matrix(~ -1 + tower$study, data = tower)
X_tower

bmr_1 <- bmr(y = tower$log_ROM, sigma = tower$log_SE, X=X_tower,
             labels = tower$study, #mu.prior.mean = 0, mu.prior.sd = 4,
             tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)


bmr_1_means <- bmr_1$summary[3,]
bmr_1_means <- bmr_1_means[-1]
bmr_1_ses <- (bmr_1$summary[6,]-bmr_1$summary[5,])/3.92
bmr_1_ses <- bmr_1_ses[-1]

bmr_1_studylevel <- bayesmeta(y = bmr_1_means, sigma = bmr_1_ses,
                              labels = sort(unique(tower$study)),
                              #labels = bmr1$study,
                              mu.prior.mean = 0, mu.prior.sd = 4,
                              tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

exp(bmr_1_studylevel$summary)
forestplot(bmr_1_studylevel, expo=T, clip=c(-1, 3))

png(file = "suppl_figures/bayesmeta_study-level_forestplot_tower.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_1_studylevel, expo=T, clip=c(-1, 3))
dev.off()

bmr_1_studylevel$likelihood(mu=0)
bmr_1_studylevel$bayesfactor[1, "mu=0"]


brm.est1 <- c(fixef(brm_m1)[1],fixef(brm_m1)[3],fixef(brm_m1)[4])
brm.est1 <- meta:::backtransf(brm.est1, sm="ROM")

bayes.est1 <- c(bm1$summary["mode",2], bm1$summary["95% lower",2], bm1$summary["95% upper",2])
bayes.est1 <- meta:::backtransf(bayes.est1, sm="ROM")
bayes.est1_studylvl <- exp(bmr_1_studylevel$summary[c(3,5,6),2])

brm.est1
bayes.est1
bayes.est1_studylvl
meta.est1
rma.est1

# Robust bayesian model
# rbm1 <- RoBMA(z = tower$z, se = tower$se_z, seed = 1, model = "PSMA", parallel = TRUE)
# 
# forest(rbm1, output_scale = "logOR")
# summary(rbm1)
# robust.bay.est1 <- summary(rbm1, output_scale = "logOR")$estimates[1,]
# robust.bay.est1 <- meta:::backtransf(robust.bay.est1, sm="ROM")
# RoBMA.est1 <- c(robust.bay.est1$Mean, robust.bay.est1$"0.025", robust.bay.est1$"0.975")
# RoBMA.est1


# Base station (without "reduced abundance" findings)
m1b <- metagen(TE = log_ROM,
              seTE = log_SE,
              data = tower2,
              cluster = study,
              comb.fixed = TRUE,
              comb.random = T,
              hakn = FALSE,
              prediction = TRUE,
              sm = "ROM")
meta::forest(m1b)

m1b <- update(m1b, subgroup = study)

meta::forest(m1b, layout = "RevMan5",
            label.right = "detrimental", col.label.right = "red",
            label.left = "beneficial", col.label.left = "green",
            prediction = TRUE, test.overall.random = TRUE,
            #file = "suppl_figures/forestplot_base_station2_subgroups.png", width = 3300, rows.gr = 420, args.gr = list(res = 300), dev.off = T,
            file = "suppl_figures/forestplot_base_station2_subgroups.png", width = 4200, rows.gr = 160, args.gr = list(res = 300), dev.off = T,
            #leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
            #leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
            leftcols = c("Author_Year","E_Field","CummHrs","Insect_type","Bioeffects","effect","w.random"),
            leftlabs = c("Author, Year","E-Field [V/m]","Duration [h]","Insect","Bioeffects"),
            just = "center", just.addcols = "center", xlim = c(0.2,5))

meta.est1b <- c(m1b$TE.random,m1b$lower.random,m1b$upper.random)
meta.est1b <- meta:::backtransf(meta.est1b, sm="ROM")
percent((1/exp(m1b$TE.random))-1) # Effect size as percentage change versus control
m1b$pval.random

res1 <- rma(TE, sei=seTE, data=m1b)
inf <- influence(res1)
plot(inf)
regtest(res1)
meta::funnel(m1b)

full.model1b <- rma.mv(yi, vi, random = ~ 1 | study/experiment, data = tower2) 
exp(c(full.model1b$b, full.model1b$ci.lb, full.model1b$ci.ub))
rma.est1b <- exp(c(full.model1b$b, full.model1b$ci.lb, full.model1b$ci.ub))
rma.est1b.pval <- full.model1b$pval
meta.est1b


# Bayesian models
# "brms" model
brm_m1b <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study/experiment), data = tower2,
              iter = 30000, warmup = 1000,
              control = list(adapt_delta = 0.98),
              core = 6, chains = 6)
exp(fixef(brm_m1b))


gg_forestplot(brm_m1b) + xlim(0.5, 2.5)
ggsave("suppl_figures/brms_forestplot_base_station2.png", width = 12, height = 6)


# "bayesmeta" model
bm1b <- bayesmeta(y = tower2$log_ROM, sigma = tower2$log_SE,
                 labels = tower2$study, mu.prior.mean = 0, mu.prior.sd = 4,
                 tau.prior = function(t) dhalfnormal(t, scale = 0.5))
forestplot(bm1b, expo=T)

png(file = "suppl_figures/bayesmeta_forestplot_base_station2.png", width = 3500, height = 3000, res = 300)
forestplot(bm1b, expo=T, clip=c(-1, 3))
dev.off()

bm1b$bayesfactor
bm1b$bayesfactor[1, "mu=0"]

# "bayesmeta" study-level model
X_tower2 <- model.matrix(~ -1 + tower2$study, data = tower2)
X_tower2

bmr_1b <- bmr(y = tower2$log_ROM, sigma = tower2$log_SE, X=X_tower2,
             labels = tower2$study, #mu.prior.mean = 0, mu.prior.sd = 4,
             tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

bmr_1b_means <- bmr_1b$summary[3,]
bmr_1b_means <- bmr_1b_means[-1]
bmr_1b_ses <- (bmr_1b$summary[6,]-bmr_1b$summary[5,])/3.92
bmr_1b_ses <- bmr_1b_ses[-1]

bmr_1b_studylevel <- bayesmeta(y = bmr_1b_means, sigma = bmr_1b_ses,
                              labels = sort(unique(tower2$study)),
                              mu.prior.mean = 0, mu.prior.sd = 4,
                              tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

exp(bmr_1b_studylevel$summary)
forestplot(bmr_1b_studylevel, expo=T,)

png(file = "suppl_figures/bayesmeta_study-level_forestplot_tower2.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_1b_studylevel, expo=T, clip=c(-1, 3))
dev.off()

bmr_1b_studylevel$likelihood(mu=0)
bmr_1b_studylevel$bayesfactor[1, "mu=0"]


# grabbing estimates of mean effect size + 95% CI
brm.est1b <- c(fixef(brm_m1b)[1],fixef(brm_m1b)[3],fixef(brm_m1b)[4])
brm.est1b <- meta:::backtransf(brm.est1b, sm="ROM")
bayes.est1b <- c(bm1b$summary["mode",2], bm1b$summary["95% lower",2], bm1b$summary["95% upper",2])
bayes.est1b <- meta:::backtransf(bayes.est1b, sm="ROM")
bayes.est1b_studylvl <- exp(bmr_1b_studylevel$summary[c(3,5,6),2])

# all estimates
brm.est1b
bayes.est1b
bayes.est1b_studylvl
meta.est1b
rma.est1b


# Base station (only "reduced abundance" findings)
m1c <- metagen(TE = log_ROM,
               seTE = log_SE,
               data = tower3,
               cluster = study,
               comb.fixed = TRUE,
               comb.random = T,
               hakn = FALSE,
               prediction = TRUE,
               sm = "ROM")
meta::forest(m1c)

m1c <- update(m1c, subgroup = study)

meta::forest(m1c, layout = "RevMan5",
             label.right = "detrimental", col.label.right = "red",
             label.left = "beneficial", col.label.left = "green",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = F, print.subgroup.name = T,
             #file = "suppl_figures/forestplot_base_station3_subgroups.png", width = 3600, rows.gr = 220, args.gr = list(res = 300), dev.off = T,
             file = "suppl_figures/forestplot_base_station3_subgroups.png", width = 4200, rows.gr = 160, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","E_Field","Insect_type","Bioeffects","effect","w.random"),
             leftlabs = c("Author, Year","E-Field [V/m]","Insect","Bioeffects"),
             just = "center", just.addcols = "center", xlim = c(0.2,12))

meta.est1c <- c(m1c$TE.random,m1c$lower.random,m1c$upper.random)
meta.est1c <- meta:::backtransf(meta.est1c, sm="ROM")
percent((1/exp(m1c$TE.random))-1) # Effect size as percentage change over control
m1c$pval.random

meta::funnel(m1c)

# "metafor" model
full.model1c <- rma.mv(yi, vi, random = ~ 1 | study/experiment, data = tower3) 
exp(c(full.model1c$b, full.model1c$ci.lb, full.model1c$ci.ub))
rma.est1c <- exp(c(full.model1c$b, full.model1c$ci.lb, full.model1c$ci.ub))
rma.est1c.pval <- full.model1c$pval

# Bayesian models
# "brms" model
brm_m1c <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study/experiment), data = tower3,
               iter = 30000, warmup = 1000,
               control = list(adapt_delta = 0.98),
               core = 6, chains = 6)
exp(fixef(brm_m1c))

gg_forestplot(brm_m1c) + xlim(0.2, 8)

ggsave("suppl_figures/brms_forestplot_base_station3.png", width = 12, height = 6)


# "bayesmeta" model
bm1c <- bayesmeta(y = tower3$log_ROM, sigma = tower3$log_SE,
                  labels = tower3$study, mu.prior.mean = 0, mu.prior.sd = 4,
                  tau.prior = function(t) dhalfnormal(t, scale = 0.5))
forestplot(bm1c, expo=T)

png(file = "suppl_figures/bayesmeta_forestplot_base_station3.png", width = 3500, height = 3000, res = 300)
forestplot(bm1c, expo=T, clip=c(-1, 3))
dev.off()

summary(bm1c)
bm1c$bayesfactor
bm1c$bayesfactor[1, "mu=0"]

# "bayesmeta" study-level model
X_tower3 <- model.matrix(~ -1 + tower3$study, data = tower3)
X_tower3

bmr_1c <- bmr(y = tower3$log_ROM, sigma = tower3$log_SE, X=X_tower3,
              labels = tower3$study, #mu.prior.mean = 0, mu.prior.sd = 4,
              tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

bmr_1c_means <- bmr_1c$summary[3,]
bmr_1c_means <- bmr_1c_means[-1]
bmr_1c_ses <- (bmr_1c$summary[6,]-bmr_1c$summary[5,])/3.92
bmr_1c_ses <- bmr_1c_ses[-1]

bmr_1c_studylevel <- bayesmeta(y = bmr_1c_means, sigma = bmr_1c_ses,
                               labels = sort(unique(tower3$study)),
                               mu.prior.mean = 0, mu.prior.sd = 4,
                               tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

exp(bmr_1c_studylevel$summary)
forestplot(bmr_1c_studylevel, expo=T)

png(file = "suppl_figures/bayesmeta_study-level_forestplot_tower3.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_1c_studylevel, expo=T, clip=c(-1, 3))
dev.off()

bmr_1c_studylevel$likelihood(mu=0)
bmr_1c_studylevel$bayesfactor[1, "mu=0"]

brm.est1c <- c(fixef(brm_m1c)[1],fixef(brm_m1c)[3],fixef(brm_m1c)[4])
brm.est1c <- meta:::backtransf(brm.est1c, sm="ROM")
bayes.est1c <- c(bm1c$summary["mode",2], bm1c$summary["95% lower",2], bm1c$summary["95% upper",2])
bayes.est1c <- meta:::backtransf(bayes.est1c, sm="ROM")
bayes.est1c_studylvl <- exp(bmr_1c_studylevel$summary[c(3,5,6),2])

# all estimates
brm.est1c
bayes.est1c
bayes.est1c_studylvl
meta.est1c
rma.est1c


##########
# DECT
##########
m2 <- metagen(TE = log_ROM,
              seTE = log_SE,
              data = DECT,
              cluster = study,
              comb.fixed = TRUE,
              comb.random = T,
              hakn = FALSE,
              prediction = TRUE,
              sm = "ROM")

summary(m2)
meta.est2 <- c(m2$TE.random, m2$lower.random, m2$upper.random)
meta.est2 <- meta:::backtransf(meta.est2, sm="ROM")
meta.est2

meta::forest(m2)

#png(file = "suppl_figures/forestplot_DECT.png", width = 3500, height = 2000, res = 300)
meta::forest(m2, layout = "RevMan5",
            label.right = "detrimental", col.label.right = "red",
            label.left = "beneficial", col.label.left = "green",
            prediction = TRUE, test.overall.random = TRUE,
            leftcols = c("Author_Year","E_Field","CummHrs","CellsAnimal","Bioeffects"),
            leftlabs = c("Author, Year","E-Field [V/m]","Duration [h]","Insect","Bioeffects"),
            just = "center", just.addcols = "center", xlim = c(0.25,4))
#dev.off()

m2 <- update(m2, subgroup = study)

meta::forest(m2, layout = "RevMan5",
             label.right = "detrimental", col.label.right = "red",
             label.left = "beneficial", col.label.left = "green",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = F, print.subgroup.name = T,
             file = "suppl_figures/forestplot_DECT_subgroups.png", width = 3300, rows.gr = 100, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
             leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
             just = "center", just.addcols = "center", xlim = c(0.25,4))

summary(m2)
meta.est2 <- c(m2$TE.random, m2$lower.random, m2$upper.random)
meta.est2 <- meta:::backtransf(meta.est2, sm="ROM")
percent((1/exp(m2$TE.random))-1) # Effect size as percentage change over control
m2$pval.random

meta::funnel(m2)
res2 <- rma(TE, sei=seTE, data=m2)
regtest(res2)


# three-level model using the "Metafor" package
full.model <- rma.mv(yi, vi, random = ~ 1 | study/experiment, data = DECT) # correct full model
exp(c(full.model$b, full.model$ci.lb, full.model$ci.ub))
rma.est2 <- exp(c(full.model$b, full.model$ci.lb, full.model$ci.ub))
rma.est2.pval <- full.model$pval

full.model2 <- rma.mv(yi, vi, random = ~ 1 | study, data = DECT) 
exp(c(full.model2$b, full.model2$ci.lb, full.model2$ci.ub))
full.model3 <- rma.mv(yi, vi, random = ~ 1 | experiment, data = DECT) 
exp(c(full.model3$b, full.model3$ci.lb, full.model3$ci.ub))
percent((1/exp(full.model$b))-1)
percent((1/exp(full.model2$b))-1)
percent((1/exp(full.model3$b))-1)


# Bayesian models
# "brms" model

brm_m2 <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|experiment) + (1|study), data = DECT,
               iter = 40000, warmup = 1000,
               control = list(adapt_delta = 0.99),
               core = 6, chains = 6)
exp(fixef(brm_m2))
#Intercept 1.449571  1.266467 0.9225131 2.421757


gg_forestplot(brm_m2) + xlim(0,3)
ggsave("suppl_figures/brms_forestplot_DECT.png", width = 12, height = 6)


# "bayesmeta" model
bm2 <- bayesmeta(y = DECT$log_ROM, sigma = DECT$log_SE,
                 labels = DECT$study, mu.prior.mean = 0, mu.prior.sd = 4,
                 tau.prior = function(t) dhalfnormal(t, scale = 0.5))

exp(bm2$summary)
forestplot(bm2, expo=T)

png(file = "suppl_figures/bayesmeta_forestplot_DECT.png", width = 3500, height = 2000, res = 300)
forestplot(bm2, expo=T, clip=c(-1, 3))
dev.off()


#bm2$likelihood(mu=0)
bm2$bayesfactor
bm2$bayesfactor[1, "mu=0"]


# "bayesmeta" model at study-level
X_DECT <- model.matrix(~ -1 + DECT$study, data = DECT)
X_DECT

bmr_2 <- bmr(y = DECT$log_ROM, sigma = DECT$log_SE, X=X_DECT,
            labels = DECT$study, mu.prior.mean = 0, mu.prior.sd = 4,
            tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

bmr_2_means <- bmr_2$summary[3,]
bmr_2_means <- bmr_2_means[-1]
# For 95% CI: SE = (Upper Limit - Lower Limit) / 3.92
bmr_2_ses <- (bmr_2$summary[6,]-bmr_2$summary[5,])/3.92
bmr_2_ses <- bmr_2_ses[-1]

bmr_2_studylevel <- bayesmeta(y = bmr_2_means, sigma = bmr_2_ses,
                             labels = unique(DECT$study), mu.prior.mean = 0, mu.prior.sd = 4,
                             tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

exp(bmr_2_studylevel$summary)
forestplot(bmr_2_studylevel, expo=T, clip=c(-1, 3))

png(file = "suppl_figures/bayesmeta_study-level_forestplot_DECT.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_2_studylevel, expo=T, clip=c(-1, 3))
dev.off()

bmr_2_studylevel$likelihood(mu=0)
bmr_2_studylevel$bayesfactor[1, "mu=0"]


brm.est2 <- c(fixef(brm_m2)[1], fixef(brm_m2)[3], fixef(brm_m2)[4])
brm.est2 <- meta:::backtransf(brm.est2, sm="ROM")
bayes.est2 <- c(bm2$summary["mode",2], bm2$summary["95% lower",2], bm2$summary["95% upper",2])
bayes.est2 <- meta:::backtransf(bayes.est2, sm="ROM")
bayes.est2_studylvl <- exp(bmr_2_studylevel$summary[c(3,5,6),2])

# all estimates
brm.est2
bayes.est2
bayes.est2_studylvl
meta.est2
rma.est2


###############
# Mobile phone
###############

m3 <- metagen(TE = log_ROM,
              seTE = log_SE,
              data = mobile,
              cluster = study,
              comb.fixed = TRUE,
              comb.random = T,
              hakn = FALSE,
              prediction = TRUE,
              sm = "ROM")
meta::forest(m3)

#png(file = "suppl_figures/forestplot_mobile_phone.png", width = 3500, height = 5500, res = 300)
meta::forest(m3, layout = "RevMan5",
            label.right = "detrimental", col.label.right = "red",
            label.left = "beneficial", col.label.left = "green",
            prediction = TRUE, test.overall.random = TRUE,
            leftcols = c("Author_Year","E_Field","CummHrs","Bioeffect_cat","effect","ci","w.random"),
            leftlabs = c("Author, Year","E-Field [V/m]","Duration [h]","Bioeffect category"),
            just = "center", just.addcols = "center")
#dev.off()

m3 <- update(m3, subgroup = study)

meta::forest(m3, layout = "RevMan5",
             label.right = "detrimental", col.label.right = "red",
             label.left = "beneficial", col.label.left = "green",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = F, print.subgroup.name = F,
             file = "suppl_figures/forestplot_mobile_phone_subgroups.png", width = 3300, rows.gr = 340, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
             leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
             just = "center", just.addcols = "center")


summary(m3)

meta.est3 <- c(m3$TE.random, m3$lower.random, m3$upper.random)
meta.est3 <- meta:::backtransf(meta.est3, sm="ROM")
m3$pval.random

meta::funnel(m3)
res3 <- rma(TE, sei=seTE, data=m3)
regtest(res3)

# three-level model using the "Metafor" package
full.model3 <- rma.mv(yi, vi, random = ~ 1 | study/experiment, data = mobile) # correct full model
exp(c(full.model3$b, full.model3$ci.lb, full.model3$ci.ub))
rma.est3 <- exp(c(full.model3$b, full.model3$ci.lb, full.model3$ci.ub))
rma.est3.pval <- full.model3$pval


# Bayesian models
# "brms" model
brm_m3 <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study/experiment), data = mobile,
              iter = 30000, warmup = 1000,
              control = list(adapt_delta = 0.95),
              core = 6, chains = 6)
exp(fixef(brm_m3))
#Intercept 1.436919  1.051347 1.303984 1.591406

gg_forestplot(brm_m3) + xlim(0.5, 2.5)

ggsave("suppl_figures/brms_forestplot_mobile_phone.png", width = 12, height = 6)


# "bayesmeta" models
bm3 <- bayesmeta(y = mobile$log_ROM, sigma = mobile$log_SE,
                 labels = mobile$study, mu.prior.mean = 0, mu.prior.sd = 4,
                 tau.prior = function(t) dhalfnormal(t, scale = 0.5))
forestplot(bm3, expo=T)

png(file = "suppl_figures/bayesmeta_forestplot_mobile_phone.png", width = 3500, height = 5500, res = 300)
forestplot(bm3, expo=T, clip=c(-1, 3))
dev.off()

bm3$bayesfactor
bm3$bayesfactor[1, "mu=0"]

# study-level model
X_mobile <- model.matrix(~ -1 + mobile$study, data = mobile)
X_mobile

bmr_3 <- bmr(y = mobile$log_ROM, sigma = mobile$log_SE, X=X_mobile,
             labels = mobile$study, #mu.prior.mean = 0, mu.prior.sd = 4,
             tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

bmr_3_means <- bmr_3$summary[3,]
bmr_3_means <- bmr_3_means[-1]
bmr_3_ses <- (bmr_3$summary[6,]-bmr_3$summary[5,])/3.92
bmr_3_ses <- bmr_3_ses[-1]

bmr_3_studylevel <- bayesmeta(y = bmr_3_means, sigma = bmr_3_ses,
                              labels = unique(mobile$study), mu.prior.mean = 0, mu.prior.sd = 4,
                              tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

exp(bmr_3_studylevel$summary)
bayes.est3_studylvl <- exp(bmr_3_studylevel$summary[c(3,5,6),2])
forestplot(bmr_3_studylevel, expo=T)

png(file = "suppl_figures/bayesmeta_study-level_forestplot_mobile.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_3_studylevel, expo=T, clip=c(-1, 3))
dev.off()

bmr_3_studylevel$likelihood(mu=0)
bmr_3_studylevel$bayesfactor[1, "mu=0"]

brm.est3 <- c(fixef(brm_m3)[1], fixef(brm_m3)[3], fixef(brm_m3)[4])
brm.est3 <- meta:::backtransf(brm.est3, sm="ROM")
bayes.est3 <- c(bm3$summary["mode",2], bm3$summary["95% lower",2], bm3$summary["95% upper",2])
bayes.est3 <- meta:::backtransf(bayes.est3, sm="ROM")

# all estimates
brm.est3
bayes.est3
bayes.est3_studylvl
meta.est3
rma.est3


###################
# Signal generator
###################
m4 <- metagen(TE = log_ROM,
              seTE = log_SE,
              data = sig_gen,
              #studlab = paste(key),
              cluster = study,
              comb.fixed = TRUE,
              comb.random = T,
              hakn = FALSE,
              prediction = TRUE,
              sm = "ROM")
meta::forest(m4)

m4 <- update(m4, subgroup = study)

meta::forest(m4, layout = "RevMan5",
             label.right = "detrimental", col.label.right = "red",
             label.left = "beneficial", col.label.left = "green",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = F, print.subgroup.name = F,
             file = "suppl_figures/forestplot_signal_generator_subgroups.png", width = 3300, rows.gr = 200, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
             leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
             just = "center", just.addcols = "center", xlim = c(0.25,5))


summary(m4)
meta.est4 <- c(m4$TE.random, m4$lower.random, m4$upper.random)
meta.est4 <- meta:::backtransf(meta.est4, sm="ROM")
meta.est4
meta::funnel(m4)
radial(m4)

# three-level model using the "Metafor" package
full.model4 <- rma.mv(yi, vi, random = ~ 1 | study/experiment, data = sig_gen) # correct full model
exp(c(full.model4$b, full.model4$ci.lb, full.model4$ci.ub))
rma.est4 <- exp(c(full.model4$b, full.model4$ci.lb, full.model4$ci.ub))
rma.est4.pval <- full.model4$pval


# Bayesian models
# "brms" model
brm_m4 <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study/experiment), data = sig_gen,
              iter = 30000, warmup = 1000,
              control = list(adapt_delta = 0.99),
              core = 6, chains = 6)
exp(fixef(brm_m4))
#Intercept 1.289902  1.128159 1.018425 1.649492

gg_forestplot(brm_m4) + xlim(0.5, 2.5)

ggsave("suppl_figures/brms_forestplot_signal_generator.png", width = 12, height = 6)


# "bayesmeta" models
bm4 <- bayesmeta(y = sig_gen$log_ROM, sigma = sig_gen$log_SE,
                 labels = sig_gen$study, mu.prior.mean = 0, mu.prior.sd = 4,
                 tau.prior = function(t) dhalfnormal(t, scale = 0.5))
forestplot(bm4, expo=t)

png(file = "suppl_figures/bayesmeta_forestplot_signal_generator.png", width = 3500, height = 2500, res = 300)
forestplot(bm4, expo=T, clip=c(-1, 3))
dev.off()

bm4$bayesfactor[1, "mu=0"]

# study-level model
X_sig_gen <- model.matrix(~ -1 + sig_gen$study, data = sig_gen)
X_sig_gen

bmr_4 <- bmr(y = sig_gen$log_ROM, sigma = sig_gen$log_SE, X=X_sig_gen,
             labels = sig_gen$study, #mu.prior.mean = 0, mu.prior.sd = 4,
             tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

bmr_4_means <- bmr_4$summary[3,]
bmr_4_means <- bmr_4_means[-1]
bmr_4_ses <- (bmr_4$summary[6,]-bmr_4$summary[5,])/3.92
bmr_4_ses <- bmr_4_ses[-1]

bmr_4_studylevel <- bayesmeta(y = bmr_4_means, sigma = bmr_4_ses,
                              labels = unique(sig_gen$study), mu.prior.mean = 0, mu.prior.sd = 4,
                              tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

exp(bmr_4_studylevel$summary)
bayes.est4_studylvl <- exp(bmr_4_studylevel$summary[c(3,5,6),2])
forestplot(bmr_4_studylevel, expo=T)

png(file = "suppl_figures/bayesmeta_study-level_forestplot_signal_generator.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_4_studylevel, expo=T, clip=c(-1, 3))
dev.off()

bmr_4_studylevel$likelihood(mu=0)
bmr_4_studylevel$bayesfactor[1, "mu=0"]

brm.est4 <- c(fixef(brm_m4)[1], fixef(brm_m4)[3], fixef(brm_m4)[4])
brm.est4 <- meta:::backtransf(brm.est4, sm="ROM")
bayes.est4 <- c(bm4$summary["mode",2], bm4$summary["95% lower",2], bm4$summary["95% upper",2])
bayes.est4 <- meta:::backtransf(bayes.est4, sm="ROM")

# all estimates
brm.est4
bayes.est4
bayes.est4_studylvl
meta.est4
rma.est4


###########################
# WiFi
###########################
m5 <- metagen(TE = log_ROM,
              seTE = log_SE,
              data = wifi,
              #studlab = paste(key),
              cluster = study,
              comb.fixed = TRUE,
              comb.random = T,
              hakn = FALSE,
              prediction = TRUE,
              sm = "ROM")
meta::forest(m5)

m5 <- update(m5, subgroup = study)

meta::forest(m5, layout = "RevMan5",
             label.right = "detrimental", col.label.right = "red",
             label.left = "beneficial", col.label.left = "green",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = F, print.subgroup.name = F,
             file = "suppl_figures/forestplot_wifi_subgroups.png", width = 3300, rows.gr = 100, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
             leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
             just = "center", just.addcols = "center")

summary(m5)

m5$I2
meta.est5 <- c(m5$TE.random, m5$lower.random, m5$upper.random)
meta.est5 <- meta:::backtransf(meta.est5, sm="ROM")
meta::funnel(m5)

# three-level model using the "Metafor" package
full.model5 <- rma.mv(yi, vi, random = ~ 1 | study/experiment, data = wifi) # correct full model
exp(c(full.model5$b, full.model5$ci.lb, full.model5$ci.ub))
rma.est5 <- exp(c(full.model5$b, full.model5$ci.lb, full.model5$ci.ub))
rma.est5.pval <- full.model5$pval


# Bayesian models
#"brms" model

brm_m5 <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study/experiment), data = wifi,
              iter = 30000, warmup = 1000,
              control = list(adapt_delta = 0.99),
              core = 6, chains = 6)
exp(fixef(brm_m5))
#Intercept  1.23627  1.246618 0.8031737 1.97265

gg_forestplot(brm_m5) + xlim(0.5, 2.5)

ggsave("suppl_figures/brms_forestplot_wifi.png", width = 12, height = 6)


# "bayesmeta" model
bm5 <- bayesmeta(y = wifi$log_ROM, sigma = wifi$log_SE,
                 labels = wifi$study, mu.prior.mean = 0, mu.prior.sd = 4,
                 tau.prior = function(t) dhalfnormal(t, scale = 0.5))
forestplot(bm5, expo=T)

png(file = "suppl_figures/bayesmeta_forestplot_wifi.png", width = 3500, height = 2500, res = 300)
forestplot(bm5, expo=T, clip=c(-1, 3))
dev.off()

bm5$bayesfactor[1, "mu=0"]

# "bayesmeta" study-level model
X_wifi <- model.matrix(~ -1 + wifi$study, data = wifi)
X_wifi

bmr_5 <- bmr(y = wifi$log_ROM, sigma = wifi$log_SE, X=X_wifi,
             labels = wifi$study, #mu.prior.mean = 0, mu.prior.sd = 4,
             tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

bmr_5_means <- bmr_5$summary[3,]
bmr_5_means <- bmr_5_means[-1]
bmr_5_ses <- (bmr_5$summary[6,]-bmr_5$summary[5,])/3.92
bmr_5_ses <- bmr_5_ses[-1]

bmr_5_studylevel <- bayesmeta(y = bmr_5_means, sigma = bmr_5_ses,
                              labels = unique(wifi$study), mu.prior.mean = 0, mu.prior.sd = 4,
                              tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

exp(bmr_5_studylevel$summary)
bayes.est5_studylvl <- exp(bmr_5_studylevel$summary[c(3,5,6),2])
forestplot(bmr_5_studylevel, expo=T)

png(file = "suppl_figures/bayesmeta_study-level_forestplot_wifi.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_5_studylevel, expo=T, clip=c(-1, 3))
dev.off()

bmr_5_studylevel$likelihood(mu=0)
bmr_5_studylevel$bayesfactor[1, "mu=0"]

brm.est5 <- c(fixef(brm_m5)[1], fixef(brm_m5)[3], fixef(brm_m5)[4])
brm.est5 <- meta:::backtransf(brm.est5, sm="ROM")
bayes.est5 <- c(bm5$summary["mode",2], bm5$summary["95% lower",2], bm5$summary["95% upper",2])
bayes.est5 <- meta:::backtransf(bayes.est5, sm="ROM")

# all estimates
brm.est5
bayes.est5
bayes.est5_studylvl
meta.est5
rma.est5


###########################
# (Helmholtz) coil systems
###########################
coilsystem <- coilsystem %>% arrange(experiment_id)

m6 <- metagen(TE = log_ROM,
              seTE = log_SE,
              data = coilsystem,
              #studlab = paste(key),
              cluster = study,
              comb.fixed = TRUE,
              comb.random = T,
              hakn = FALSE,
              prediction = TRUE,
              sm = "ROM")
meta::forest(m6)

m6 <- update(m6, subgroup = study)

meta::forest(m6, layout = "RevMan5",
             label.right = "detrimental", col.label.right = "red",
             label.left = "beneficial", col.label.left = "green",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = F, print.subgroup.name = F,
             file = "suppl_figures/forestplot_coil_system_subgroups.png", width = 3300, rows.gr = 200, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
             leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
             just = "center", just.addcols = "center", xlim = c(0.25,5))


summary(m6)

meta.est6 <- c(m6$TE.random, m6$lower.random, m6$upper.random)
meta.est6 <- meta:::backtransf(meta.est6, sm="ROM")
meta::funnel(m6)

# three-level model using the "Metafor" package
full.model <- rma.mv(yi, vi, random = ~ 1 | study/experiment, data = coilsystem)
exp(c(full.model$b, full.model$ci.lb, full.model$ci.ub))
rma.est6 <- exp(c(full.model$b, full.model$ci.lb, full.model$ci.ub))
rma.est6.pval <- full.model$pval


# Bayesian models
#"brms" model

brm_m6 <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study/experiment), data = coilsystem,
              iter = 30000, warmup = 1000,
              control = list(adapt_delta = 0.99),
              core = 6, chains = 6)
exp(fixef(brm_m6))
#Intercept 1.363545  1.095532 1.15304 1.655338

gg_forestplot(brm_m6) + xlim(0.5, 2.5)

ggsave("suppl_figures/brms_forestplot_coil_system.png", width = 12, height = 6)


# "bayesmeta" model
bm6 <- bayesmeta(y = coilsystem$log_ROM, sigma = coilsystem$log_SE,
                 labels = coilsystem$study, mu.prior.mean = 0, mu.prior.sd = 4,
                 tau.prior = function(t) dhalfnormal(t, scale = 0.5))
forestplot(bm6, expo=T)

png(file = "suppl_figures/bayesmeta_forestplot_coil_system.png", width = 3500, height = 2500, res = 300)
forestplot(bm6, expo=T, clip=c(-1, 3))
dev.off()

bm6$bayesfactor[1, "mu=0"]

# "bayesmeta" study-level model
X_coilsystem <- model.matrix(~ -1 + coilsystem$study, data = coilsystem)
X_coilsystem

bmr_6 <- bmr(y = coilsystem$log_ROM, sigma = coilsystem$log_SE, X=X_coilsystem,
             labels = coilsystem$study, #mu.prior.mean = 0, mu.prior.sd = 4,
             tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

bmr_6_means <- bmr_6$summary[3,]
bmr_6_means <- bmr_6_means[-1]
bmr_6_ses <- (bmr_6$summary[6,]-bmr_6$summary[5,])/3.92
bmr_6_ses <- bmr_6_ses[-1]

bmr_6_studylevel <- bayesmeta(y = bmr_6_means, sigma = bmr_6_ses,
                              labels = unique(coilsystem$study), mu.prior.mean = 0, mu.prior.sd = 4,
                              tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

exp(bmr_6_studylevel$summary)
forestplot(bmr_6_studylevel, expo=T, clip=c(-1, 5))

png(file = "suppl_figures/bayesmeta_study-level_forestplot_coil_system.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_6_studylevel, expo=T, clip=c(-1, 5))
dev.off()

bmr_6_studylevel$likelihood(mu=0)
bmr_6_studylevel$bayesfactor[1, "mu=0"]

brm.est6 <- c(fixef(brm_m6)[1], fixef(brm_m6)[3], fixef(brm_m6)[4])
brm.est6 <- meta:::backtransf(brm.est6, sm="ROM")
bayes.est6 <- c(bm6$summary["mode",2], bm6$summary["95% lower",2], bm6$summary["95% upper",2])
bayes.est6 <- meta:::backtransf(bayes.est6, sm="ROM")
bayes.est6_studylvl <- exp(bmr_6_studylevel$summary[c(3,5,6),2])

# all estimates
brm.est6
bayes.est6
bayes.est6_studylvl
meta.est6
rma.est6


########################################
# Gathering all estimates in a table
#######################################

meta_es_data <- data.frame(
  EMF_source = factor(c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system")),
  k = c(m1$k, m2$k, m3$k, m4$k, m5$k, m6$k),
  ROM = c(meta.est1[1],meta.est2[1],meta.est3[1],meta.est4[1],meta.est5[1],meta.est6[1]),
  LLCI = c(meta.est1[2],meta.est2[2],meta.est3[2],meta.est4[2],meta.est5[2],meta.est6[2]),
  ULCI = c(meta.est1[3],meta.est2[3],meta.est3[3],meta.est4[3],meta.est5[3],meta.est6[3]),
  estimate = c(0,0,0,0,0,0),
  pval = c(m1$pval.random, m2$pval.random, m3$pval.random, m4$pval.random, m5$pval.random, m6$pval.random))

bayes_es_data <- data.frame(
  EMF_source = factor(c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system")),
  k = c(m1$k, m2$k, m3$k, m4$k, m5$k, m6$k),
  ROM = c(bayes.est1[1],bayes.est2[1],bayes.est3[1],bayes.est4[1],bayes.est5[1],bayes.est6[1]),
  LLCI = c(bayes.est1[2],bayes.est2[2],bayes.est3[2],bayes.est4[2],bayes.est5[2],bayes.est6[2]),
  ULCI = c(bayes.est1[3],bayes.est2[3],bayes.est3[3],bayes.est4[3],bayes.est5[3],bayes.est6[3]),
  estimate = c(1,1,1,1,1,1),
  pval = c(bm1$bayesfactor[1, "mu=0"],bm2$bayesfactor[1, "mu=0"],bm3$bayesfactor[1, "mu=0"],bm4$bayesfactor[1, "mu=0"],bm5$bayesfactor[1, "mu=0"],bm6$bayesfactor[1, "mu=0"]))

bayes_es_studylvl_data <- data.frame(
  EMF_source = factor(c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system")),
  k = c(m1$k, m2$k, m3$k, m4$k, m5$k, m6$k),
  ROM = c(bayes.est1_studylvl[1],bayes.est2_studylvl[1],bayes.est3_studylvl[1],bayes.est4_studylvl[1],bayes.est5_studylvl[1],bayes.est6_studylvl[1]),
  LLCI = c(bayes.est1_studylvl[2],bayes.est2_studylvl[2],bayes.est3_studylvl[2],bayes.est4_studylvl[2],bayes.est5_studylvl[2],bayes.est6_studylvl[2]),
  ULCI = c(bayes.est1_studylvl[3],bayes.est2_studylvl[3],bayes.est3_studylvl[3],bayes.est4_studylvl[3],bayes.est5_studylvl[3],bayes.est6_studylvl[3]),
  estimate = c(2,2,2,2,2,2),
  pval = c(bmr_1_studylevel$bayesfactor[1, "mu=0"],bmr_2_studylevel$bayesfactor[1, "mu=0"],bmr_3_studylevel$bayesfactor[1, "mu=0"],bmr_4_studylevel$bayesfactor[1, "mu=0"],bmr_5_studylevel$bayesfactor[1, "mu=0"],bmr_6_studylevel$bayesfactor[1, "mu=0"]))

brm_es_data <- data.frame(
  EMF_source = factor(c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system")),
  k = c(m1$k, m2$k, m3$k, m4$k, m5$k, m6$k),
  ROM = c(brm.est1[1],brm.est2[1],brm.est3[1],brm.est4[1],brm.est5[1],brm.est6[1]),
  LLCI = c(brm.est1[2],brm.est2[2],brm.est3[2],brm.est4[2],brm.est5[2],brm.est6[2]),
  ULCI = c(brm.est1[3],brm.est2[3],brm.est3[3],brm.est4[3],brm.est5[3],brm.est6[3]),
  estimate = c(3,3,3,3,3,3))

rma_es_data <- data.frame(
  EMF_source = factor(c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system")),
  k = c(m1$k, m2$k, m3$k, m4$k, m5$k, m6$k),
  ROM = c(rma.est1[1],rma.est2[1],rma.est3[1],rma.est4[1],rma.est5[1],rma.est6[1]),
  LLCI = c(rma.est1[2],rma.est2[2],rma.est3[2],rma.est4[2],rma.est5[2],rma.est6[2]),
  ULCI = c(rma.est1[3],rma.est2[3],rma.est3[3],rma.est4[3],rma.est5[3],rma.est6[3]),
  estimate = c(4,4,4,4,4,4),
  pval = c(rma.est1.pval, rma.est2.pval, rma.est3.pval, rma.est4.pval, rma.est5.pval, rma.est6.pval))


base_data <- data.frame(
  EMF_source = factor(c("Base station","Base station","Base station","Base station","Base station","Base station","Base station","Base station","Base station","Base station")),
  k = c(m1b$k, m1b$k, m1b$k, m1b$k, m1b$k, m1c$k, m1c$k, m1c$k, m1c$k, m1c$k),
  ROM = c(meta.est1b[1],bayes.est1b[1],bayes.est1b_studylvl[1],brm.est1b[1],rma.est1b[1],meta.est1c[1],bayes.est1c[1],bayes.est1c_studylvl[1],brm.est1c[1],rma.est1c[1]),
  LLCI = c(meta.est1b[2],bayes.est1b[2],bayes.est1b_studylvl[2],brm.est1b[2],rma.est1b[2],meta.est1c[2],bayes.est1c[2],bayes.est1c_studylvl[2],brm.est1c[2],rma.est1c[2]),
  ULCI = c(meta.est1b[3],bayes.est1b[3],bayes.est1b_studylvl[3],brm.est1b[3],rma.est1b[3],meta.est1c[3],bayes.est1c[3],bayes.est1c_studylvl[3],brm.est1c[3],rma.est1c[3]),
  estimate = c(5,6,7,8,9,10,11,12,13,14))

es_data <- full_join(meta_es_data, base_data)
es_data <- full_join(es_data, bayes_es_data)
es_data <- full_join(es_data, bayes_es_studylvl_data)
es_data <- full_join(es_data, brm_es_data)

all_es_data <- full_join(es_data, rma_es_data)

all_es_data <- arrange(all_es_data, EMF_source, estimate)

all_es_data$ROM <- round(all_es_data$ROM, digits = 3)
all_es_data$LLCI <- round(all_es_data$LLCI, digits = 3)
all_es_data$ULCI <- round(all_es_data$ULCI, digits = 3)
all_es_data$pval <- round(all_es_data$pval, digits = 6)

all_es_data

# saving all estimates into spreadsheet
openxlsx::write.xlsx(all_es_data, "tables/all_es_data.xlsx", rowNames = F, colWidths = 15, overwrite = T, firstRow = T, firstCol = T)


############################################################################
# Summary table with results of three-level meta-analysis (package meta)
############################################################################


# Building a table
`EMF type` <- c("HF","HF","HF","Base station (all)","Base station (tox)","DECT","Mobile phone","Signal generator","Coil system")
`E-field strength` <- c("> 7 V/m", "2.2--7 V/m", "< 2.2 V/m",
                        "All","All","All","All","All","All")
`E-field (median)` <- c(median(madata1$data$E_Field),median(madata2$data$E_Field),median(madata3$data$E_Field),
                        median(m1$data$E_Field), median(m1b$data$E_Field), median(m2$data$E_Field), median(m3$data$E_Field), median(m4$data$E_Field),median(m6$data$E_Field))
`E-field (median)` <- round(`E-field (median)`, digits = 2)
`Insect` <- c("Drosophila","Drosophila","Drosophila","All","All","All","All","All","All")
`N studies` <- c(madata1$k.study,madata2$k.study,madata3$k.study, m1$k.study, m1b$k.study, m2$k.study,m3$k.study,m4$k.study,m6$k.study)
`N experiments` <- c(madata1$k.TE, madata2$k.TE, madata3$k.TE, m1$k.TE, m1b$k.TE, m2$k.TE, m3$k.TE, m4$k.TE, m6$k.TE)
`ROM` <- c(exp(madata1$TE.random), exp(madata2$TE.random), exp(madata3$TE.random),
                     exp(m1$TE.random), exp(m1b$TE.random), exp(m2$TE.random), exp(m3$TE.random), exp(m4$TE.random), exp(m6$TE.random))
`ROM` <- round(`ROM`, digits = 2)
LLCI <- c(meta.est_high[2],meta.est_mid[2],meta.est_low[2],meta.est1[2],meta.est1b[2],meta.est2[2],meta.est3[2],meta.est4[2],meta.est6[2])
ULCI <- c(meta.est_high[3],meta.est_mid[3],meta.est_low[3],meta.est1[3],meta.est1b[3], meta.est2[3],meta.est3[3],meta.est4[3],meta.est6[3])
LLCI <- round(LLCI, digits=2)
ULCI <- round(ULCI, digits=2)
`95% CI` <- paste0(LLCI,"--",ULCI)
pvalue <- c(madata1$pval.random, madata2$pval.random, madata3$pval.random, m1$pval.random, m1b$pval.random, m2$pval.random, m3$pval.random, m4$pval.random, m6$pval.random)
pvalue <- round(pvalue, digits=3)
pvalue2 <- as.character(pvalue)

meta_summary <- tibble(`EMF type`,`E-field strength`,`E-field (median)`,`Insect`,`N studies`,`N experiments`,`ROM`,`95% CI`,pvalue,pvalue2)

meta_summary <- mutate(meta_summary, pvalue2 = if_else ( pvalue < 0.0001, "< 0.0001", pvalue2))
meta_summary <-  meta_summary %>% select(!pvalue)

names(meta_summary) <- c("EMF type", "E-field strength", "E-field (median)", "Insect", "Studies (n= )", "Experiments (n= )", "ROM", "95% CI", "p-value") 

meta_summary
View(meta_summary)

openxlsx::write.xlsx(meta_summary, "tables/table_meta_summary.xlsx", rowNames = F, overwrite = T)


##############################
# table of probabilities
##############################

library(gridExtra)
#meta_summary <- openxlsx::read.xlsx("tables/table_meta_summary.xlsx", sheet = 1, na.strings = "NA") 

pdf("tables/table_meta_summary.pdf", height=3, width=11) # save to pdf
gridExtra::grid.table(meta_summary, rows=NULL)
dev.off()

# latex code for table
latexcode <- xtable::xtable(meta_summary, caption = 'Summary table of meta-analysis results', auto = T)
print(latexcode, include.rownames = F)


############################################################################
# Summary table with results of Bayesian meta-analysis (package bayesmeta)
############################################################################


# Building a table
`EMF type` <- c("HF","HF","HF","Base station (all)","Base station (tox)","DECT","Mobile phone","Signal generator","Coil system")
`E-field strength` <- c("> 7 V/m", "2.2--7 V/m", "< 2.2 V/m",
                        "All","All","All","All","All","All")
`E-field (median)` <- c(median(madata1$data$E_Field),median(madata2$data$E_Field),median(madata3$data$E_Field),
                        median(m1$data$E_Field), median(m1b$data$E_Field), median(m2$data$E_Field), median(m3$data$E_Field), median(m4$data$E_Field),median(m6$data$E_Field))
`E-field (median)` <- round(`E-field (median)`, digits = 2)
`Insect` <- c("Drosophila","Drosophila","Drosophila","All","All","All","All","All","All")
`N studies` <- c(madata1$k.study,madata2$k.study,madata3$k.study, m1$k.study, m1b$k.study, m2$k.study,m3$k.study,m4$k.study,m6$k.study)
`N experiments` <- c(madata1$k.TE, madata2$k.TE, madata3$k.TE, m1$k.TE, m1b$k.TE, m2$k.TE, m3$k.TE, m4$k.TE, m6$k.TE)
`ROM` <- c(bayes.est_high[1],bayes.est_mid[1],bayes.est_low[1],
           bayes.est1[1],bayes.est1b[1],bayes.est2[1],bayes.est3[1],bayes.est4[1],bayes.est5[1])
`ROM` <- round(`ROM`, digits = 2)
LLCI <- c(bayes.est_high[2],bayes.est_mid[2],bayes.est_low[2],bayes.est1[2],bayes.est1b[2],bayes.est2[2],bayes.est3[2],bayes.est4[2],bayes.est5[2])
ULCI <- c(bayes.est_high[3],bayes.est_mid[3],bayes.est_low[3],bayes.est1[3],bayes.est1b[3], bayes.est2[3],bayes.est3[3],bayes.est4[3],bayes.est5[3])
LLCI <- round(LLCI, digits=2)
ULCI <- round(ULCI, digits=2)
`95% CI` <- paste0(LLCI,"--",ULCI)
pvalue <- c(bma1$bayesfactor[1, "mu=0"],bma2$bayesfactor[1, "mu=0"],bma3$bayesfactor[1, "mu=0"], 
               bm1$bayesfactor[1, "mu=0"],bm1b$bayesfactor[1, "mu=0"],bm2$bayesfactor[1, "mu=0"],bm3$bayesfactor[1, "mu=0"],bm4$bayesfactor[1, "mu=0"],bm6$bayesfactor[1, "mu=0"])
pvalue2 <- round(pvalue, digits=4)
pvalue2 <- as.character(pvalue2)


bayesmeta_summary <- tibble(`EMF type`,`E-field strength`,`E-field (median)`,`Insect`,`N studies`,`N experiments`,`ROM`,`95% CI`,pvalue)

bayesmeta_summary <- mutate(bayesmeta_summary, pvalue2 = if_else ( pvalue < 0.0001, "< 0.0001", pvalue2))
bayesmeta_summary <-  bayesmeta_summary %>% select(!pvalue)

names(bayesmeta_summary) <- c("EMF type", "E-field strength", "E-field (median)", "Insect", "Studies (n= )", "Experiments (n= )", "ROM", "95% CI", "Bayes factor") 

bayesmeta_summary
View(bayesmeta_summary)

openxlsx::write.xlsx(bayesmeta_summary, "tables/table_bayesmeta_summary.xlsx", rowNames = F, overwrite = T)


##############################
# table of probabilities
##############################

library(gridExtra)
meta_summary <- openxlsx::read.xlsx("tables/table_bayesmeta_summary.xlsx", sheet = 1, na.strings = "NA") 

pdf("tables/table_bayesmeta_summary.pdf", height=3, width=11) # save to pdf
gridExtra::grid.table(bayesmeta_summary, rows=NULL)
dev.off()

# latex code for table
latexcode <- xtable::xtable(bayesmeta_summary, caption = 'Summary table of meta-analysis results', auto = T)
print(latexcode, include.rownames = F)




####################################################################
# Meta-analysis of bee, ant and wasp studies (optional)
####################################################################

df <- read.xlsx("tables/HFLF_meta_table.xlsx", sheet = 1, na.strings = "NA") 
names(df)

df <- as_tibble(df)
df <- df[complete.cases(df[,c("log_ROM","log_SE")]),]

sort(table(df$EMF_source), decreasing = F)
sort(table(df$Insect_type), decreasing = F)
table(df$Bioeffect_cat)

# Retain only HF studies
df <- subset(df, EMF_type == "HF")

# Retain only bee, ant and wasp studies
df <- subset(df, Insect_type == "Honeybee" |
               Insect_type == "Mason bee" |
               Insect_type == "Ant" |
               Insect_type == "Wasp" |
               Insect_type == "Wild bees"
)

table(df$Insect_type)
table(df$EMF_source)
table(df$Bioeffect_cat)
table(df$study)


