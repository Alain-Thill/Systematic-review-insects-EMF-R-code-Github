
# run this first! (up to line 135)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

dir.create(file.path(getwd(), "suppl_figures")) # create subfolder for saving supplemental figures

#rm(list=ls())  # clear global environment of stuff from previous sessions


if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(dplyr, tidyr, ggplot2, ggpubr, scales, RColorBrewer, openxlsx, sqldf, stringr,
               dmetar, metafor, meta, bayesmeta, RoBMA, brms, tidyverse, tidybayes, ggridges, gridExtra)
# if "pacman" doesn't automatically install "dmetar" properly, run this!
#if (!require("remotes")) {  install.packages("remotes")}
#remotes::install_github("MathiasHarrer/dmetar")



### function to make a study-level (bayesian) BRMS forest plot (credit to Matti Vuorre, code slightly adapted to work with newest R version)

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


### changing the global settings for package "meta" 

#settings.meta("reset")
settings.meta("RevMan5")
#settings.meta("print")
# some additional settings
settings.meta(test.effect.subgroup = FALSE,
              test.subgroup = FALSE,
              common = FALSE,
              random = TRUE,
              prediction = TRUE,
              col.subgroup = "darkgray",
              label.right = "detrimental", col.label.right = "red",
              label.left = "beneficial", col.label.left = "green")


### helper function to add I2 estimate (derived from tau2) for three-level "Meta" model:

meta_I2_multilevel <- function(metaobj) {
n <- length(metaobj$seTE)
list.inverse.variances <- 1 / (metaobj$seTE ^2)
sum.inverse.variances <- sum(list.inverse.variances)
squared.sum.inverse.variances <- (sum.inverse.variances) ^ 2
list.inverse.variances.square <- 1 / (metaobj$seTE ^4)
sum.inverse.variances.square <- sum(list.inverse.variances.square)
numerator <- (n - 1) * sum.inverse.variances
denominator <- squared.sum.inverse.variances - sum.inverse.variances.square
estimated.sampling.variance <- numerator / denominator
I2_1 <- (estimated.sampling.variance) / (metaobj$tau2[[1]] + metaobj$tau2[[2]] + estimated.sampling.variance)
I2_2 <- (metaobj$tau2[2]) / (metaobj$tau2[1] + metaobj$tau2[2] + estimated.sampling.variance)
I2_3 <- (metaobj$tau2[1]) / (metaobj$tau2[1] + metaobj$tau2[2] + estimated.sampling.variance)
names(I2_1) <- "sampling variance"
total_I2 <- percent(I2_2+I2_3)
I2_2 <- percent(I2_2)
I2_3 <- percent(I2_3)
paste_command <- list(bquote(paste('Total I² = ',.(total_I2), "; within-study variance = ", .(I2_2[[1]]), "; between-study variance = ", .(I2_3[[1]])
                                   )))
eval(parse(text = paste_command))
}


### helper function to add I2 estimate for three-level "Metafor" model:

mlabfun <- function(text, metaforobj) {
  list(bquote(paste(.(text),
                    I["Level 2"]^2," = ", .(var.comp(metaforobj)$results$I2[2]), "%, ",
                    I["Level 3"]^2," = ", .(var.comp(metaforobj)$results$I2[3]), "%")))
  }


### function to make a study-level "bayesmeta" model

bayesmeta_studylevel <- function(bmr_obj) {
X_studies <- model.matrix( ~ -1 + study, data = bmr_obj)

bmr_within_study <- bmr(y = bmr_obj$log_ROM, sigma = bmr_obj$log_SE, X=X_studies,
                labels = bmr_obj$study, mu.prior.mean = 0, mu.prior.sd = 4,
                tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

# transform individual study mean + 95% CI into mean + SE:  SE = (Upper Limit - Lower Limit) / 3.92
bmr_means <- bmr_within_study$summary[3,]
bmr_means <- bmr_means[-1]
bmr_ses <- (bmr_within_study$summary[6,]-bmr_within_study$summary[5,])/3.92
bmr_ses <- bmr_ses[-1]
# rerun bayesmeta with means and standard errors for each study
bmr_studylevel <- bayesmeta(y = bmr_means, sigma = bmr_ses,
                                 labels = unique(bmr_obj$study), mu.prior.mean = 0, mu.prior.sd = 4,
                                 tau.prior = function(t) dhalfnormal(t, scale = 0.5)
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
df <- df[complete.cases(df[,c("log_ROM","log_SE")]),]

# Retain only High-Frequency EMF studies
df <- subset(df, EMF_type == "HF")

# Retain only Drosophila studies
df <- subset(df, Insect_type == "Drosophila")
table(df$SignalGen)
table(df$Bioeffect_cat)

# Retain only studies of reproductive effects
df <- subset(df, Bioeffect_cat == "Reduced reproductive capacity" |
               Bioeffect_cat == "Other" |
               Bioeffect_cat == "No effect")

#View(df)
# the filtering is correct, except for one study - due to the imprecise classifier "Other"
df <- subset(df, study != "Manta2013")

str(df)
df <- as_tibble(df)


# visualizing all Drosophila experiments by EMF source, as function of effect size and electric field
ggplot(aes(y = exp(log_ROM), x = E_Field, col=EMF_source), data = df) + 
  geom_point() + scale_color_brewer(palette="Set1") + labs(col="Field source") +
  xlab("E-Field [V/m]") + labs(col="Field source") +
  scale_x_continuous(limits = c(0, 40)) +
  scale_y_continuous(limits = c(0.5, 3)) +  #theme_classic2() +
  ylab(expression(paste("Effect size (ROM)")))+
  theme(axis.title.y = element_text(size = 11)) +
  geom_vline(xintercept = 2, linetype="dashed") +
  geom_vline(xintercept = 7, linetype="dotted") +
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_hline(yintercept = 1.5, linetype="dotted")


# visualizing all Drosophila experiments by EMF source, as function of effect size and dose
ggplot(aes(y = log_ROM, x = log(E_Field * CummHrs), col=EMF_source), data = df) + 
  geom_point() + scale_color_brewer(palette="Set1") + labs(col="Field source")


# Subdivision of all data points into high, medium, and low E-field groups.
droso_high <- subset(df, E_Field > 7)
droso_mid <- subset(df, E_Field > 2.2 & E_Field < 7)
droso_low <- subset(df, E_Field < 2.2)

#View(droso_high)
#View(droso_mid)
#View(droso_low)


#####################################################################
# Meta-analysis of all HF Drosophila studies over 7 V/m
#####################################################################

# Frequentist models
# Package "meta", three-level, i.e. clustered model

madata1 <- metagen(TE = log_ROM,
                   seTE = log_SE,
                   data = droso_high,
                   cluster = study,
                   common = F,
                   random = T,
                   prediction = TRUE,
                   sm = "ROM"
                   )
madata1
meta::forest(madata1)

# nicer
meta::forest(madata1, layout = "Revman5",
       label.right = "detrimental", col.label.right = "red",
       label.left = "beneficial", col.label.left = "green",
       prediction = TRUE, test.overall.random = TRUE,
       #file = "suppl_figures/forestplot_repro_tox_over_7Vm.png", width = 3500, rows.gr = 150, args.gr = list(res = 300), dev.off = T,
       leftcols = c("Author_Year","EMF_source","E_Field","effect","ci","w.random"),
       leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
       just = "center", just.addcols = "center", xlim = c(0.75,4),
       text.addline1 = meta_I2_multilevel(madata1)
       )

madata1 <- update(madata1, subgroup = study) # show subgroup means

meta::forest(madata1, layout = "RevMan5")

# save customized forest plot to file
meta::forest(madata1, layout = "RevMan5",
            prediction = TRUE, test.overall.random = TRUE,
            print.subgroup.labels = T, print.subgroup.name = F,
            file = "suppl_figures/forestplot_repro_tox_over_7Vm_subgroups.png", width = 3500, rows.gr = 250, args.gr = list(res = 300), dev.off = T,
            leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
            leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
            just = "center", just.addcols = "center", xlim = c(0.75,4),
            text.addline1 = meta_I2_multilevel(madata1)
            )


summary(madata1)
meta::funnel(madata1)
radial(madata1)

exp(madata1$TE.random) # ROM = "ratio of means" effect size
percent((1/exp(madata1$TE.random))-1) # Effect size as percentage change compared to control
madata1$pval.random
meta_I2_multilevel(madata1)

meta.est_high <- c(madata1$TE.random, madata1$lower.random, madata1$upper.random)
meta.est_high <- meta:::backtransf(meta.est_high, sm="ROM")
meta.est_high


# four-level model using the "Metafor" package
full.model_high <- rma.mv(yi, vi, random = ~ 1 | study/experiment/experiment_id, tdist = T, data = droso_high) # correct full model
l3.model_high <- rma.mv(yi, vi, random = ~ 1 | study/experiment_id, data = droso_high) # 
l3.removed_high <- rma.mv(yi, vi, random = ~ 1 | study/experiment_id, sigma2=c(NA,0), data = droso_high) # gives similar but not equal mean + CI
l2.model_high <- rma.mv(yi, vi, random = ~ 1 | study, data = droso_high) 


exp(c(full.model_high$b, full.model_high$ci.lb, full.model_high$ci.ub))
exp(c(l3.model_high$b, l3.model_high$ci.lb, l3.model_high$ci.ub))  # same as clustered "meta" model above
exp(c(l3.removed_high$b, l3.removed_high$ci.lb, l3.removed_high$ci.ub))
exp(c(l2.model_high$b, l2.model_high$ci.lb, l2.model_high$ci.ub)) # same as l3.removed


var.comp(l3.model_high)

# comparing full model (4-level) with study-level model indicates full model is equivalent / not superior
anova(full.model_high, l3.model_high)

profile(full.model_high)
profile(l3.model_high)

rma.est_high <- exp(c(full.model_high$b, full.model_high$ci.lb, full.model_high$ci.ub))
rma.est_high.pval <- full.model_high$pval

percent((1/exp(full.model_high$b))-1) # Effect size as percentage change compared to control
rma.est_high.pval

metafor::forest(full.model_high)

# draw customized forest plot and save to file
png(file = "suppl_figures/rma_forestplot_repro_tox_over_7Vm.png", width = 3000, height = 2500, res = 300)
metafor::forest(full.model_high, annotate=T, addfit=T, addpred=T, showweights=T, atrans=exp, refline=0,
                at=log(c(0.25, 1, 4, 16)), xlim=c(-2.5,4.5), shade="zebra", slab = experiment, cex=0.75,
                mlab=mlabfun("Variance distrib.: ", l3.model_high),
                header=c("Author and Year","Weight        ROM      95% CI"))
dev.off()

# funnel plots
funnel.rma(full.model_high)

inf <- influence(res_high)
plot(inf)

# Regression Test for Funnel Plot Asymmetry
res_high <- rma(TE, sei=seTE, data=madata1)
regtest(res_high)

# contour-enhanced funnel plot (credit to Wolfgang Viechtbauer)
funnel(full.model_high, yaxis="seinv", refline = 0,
       xaxs="i", yaxs="i", las=1,
       level=c(.10, .05, .01), shade=c("white", "gray55", "gray75"),
       legend=TRUE, back="gray90", hlines=NULL, ylab="Precision (1/se)")


# Bayesian models

# "brms" model: correct form of LME4 formula for meta-analysis at study-level and within-study "experiment" subgroups
# ~ 1 + (1|study) + (1|experiment) + (1|experiment_id)  is the same as ~ 1 + (1|study/experiment/experiment_id)

# takes a long time to run
brm1 <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study) + (1|experiment)+ (1|experiment_id), data = droso_high,
            #tau.prior=function(t){dhalfnormal(t, scale=0.5)},
            iter = 40000, warmup = 1000,
            control = list(adapt_delta = 0.99),
            core = 6, chains = 6)

summary(brm1)
exp(fixef(brm1))
percent((1/exp(fixef(brm1)[1]))-1) # Effect size as percentage change compared to control


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

#bmr_high_studylevel <- bayesmeta_studylevel(droso_high) # one-line version using function defined above

X_h <- model.matrix( ~ -1 + study, data = droso_high)
X_h

bmr_high <- bmr(y = droso_high$log_ROM, sigma = droso_high$log_SE, X=X_h,
                  labels = droso_high$study, mu.prior.mean = 0, mu.prior.sd = 4,
                  tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

# transform meta results mean + CI into mean + SE
# For 95% CI: SE = (Upper Limit - Lower Limit) / 3.92
bmr_high_means <- bmr_high$summary[3,]
bmr_high_means <- bmr_high_means[-1]
bmr_high_ses <- (bmr_high$summary[6,]-bmr_high$summary[5,])/3.92
bmr_high_ses <- bmr_high_ses[-1]

bmr_high_studylevel <- bayesmeta(y = bmr_high_means, sigma = bmr_high_ses,
            labels = unique(droso_high$study), mu.prior.mean = 0, mu.prior.sd = 4,
            tau.prior = function(t) dhalfnormal(t, scale = 0.5)
)

# results:
exp(bmr_high_studylevel$summary)
percent((1/exp(bmr_high_studylevel$summary[3,2]))-1) # Effect size as percentage change in reproductive capacity compared to control

# this could also have been achieved using the function "bayesmeta_studylevel" previously defined (start of this file) in just one line:
#bmr_high_studylevel <- bayesmeta_studylevel(droso_high)

forestplot(bmr_high_studylevel, expo=TRUE, clip=c(0.5,2.5))

png(file = "suppl_figures/bayesmeta_study-level_forestplot_repro_tox_over_7Vm.png", width = 3500, height = 2600, res = 300)
forestplot(bmr_high_studylevel, expo=TRUE, clip=c(0.5,2.5))
dev.off()

# # compute p value for effect greater or less than zero
# prob_detrimental_effect <- pppvalue(bmr_high_studylevel, parameter="mu", value=0, alternative = "greater", n=100)
# prob_detrimental_effect
# # p-value = 0.01
# prob_beneficial_effect <- pppvalue(bmr_high_studylevel, parameter="mu", value=0, alternative = "less", n=100)
# prob_beneficial_effect
# #p-value = 0.43

# bmr_high_studylevel$pposterior(mu = 0)
# bmr_high_studylevel$likelihood(mu=0)

bmr_high_studylevel$bayesfactor
bmr_high_studylevel$bayesfactor[1, "mu=0"]

# grab estimates
brm.est_high <- c(fixef(brm1)[1],fixef(brm1)[3],fixef(brm1)[4])
brm.est_high <- meta:::backtransf(brm.est_high, sm="ROM")
bayes.est_high <- c(bma1$summary["mode",2],bma1$summary["95% lower",2],bma1$summary["95% upper",2])
bayes.est_high <- meta:::backtransf(bayes.est_high, sm="ROM")
bayes.est_high2 <- exp(bmr_high_studylevel$summary[c(3,5,6),2])

# list all effect size estimates with 95% CI
meta.est_high
rma.est_high
brm.est_high
bayes.est_high
bayes.est_high2


#####################################################################
# Meta-analysis of all RF Drosophila studies between 2 and 7 V/m.
#####################################################################

# "meta" model, clustered at study-level
madata2 <- metagen(TE = log_ROM,
                   seTE = log_SE,
                   data = droso_mid,
                   cluster = study,
                   common = F,
                   random = T,
                   prediction = TRUE,
                   sm = "ROM")
madata2

exp(madata2$TE.random) # ROM = "ratio of means" Effect size
percent((1/exp(madata2$TE.random))-1) # Effect size as percentage change in reproduction over control
madata2$pval.random

meta::forest(madata2)

meta::forest(madata2, layout = "RevMan5",
            prediction = TRUE, test.overall.random = TRUE,
            #file = "suppl_figures/forestplot_repro_tox_over_2Vm_under_7Vm.png", width = 3500, rows.gr = 150, args.gr = list(res = 300), dev.off = T,
            leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
            leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
            just = "center", just.addcols = "center",
            text.addline1 = meta_I2_multilevel(madata2))


madata2 <- update(madata2, subgroup = study) # subgroups shown

meta::forest(madata2, layout = "RevMan5",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = T, print.subgroup.name = F,
             file = "suppl_figures/forestplot_repro_tox_over_2Vm_under_7Vm_subgroups.png", width = 3500, rows.gr = 200, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
             leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
             just = "center", just.addcols = "center",
             text.addline1 = meta_I2_multilevel(madata2))

meta::funnel(madata2)

meta.est_mid <- c(madata2$TE.random, madata2$lower.random, madata2$upper.random)
meta.est_mid <- meta:::backtransf(meta.est_mid, sm="ROM")
meta.est_mid



# four-level model using the "Metafor" package
full.model_mid <- rma.mv(yi, vi, random = ~ 1 | study/experiment/experiment_id, tdist = T, data = droso_mid) # correct full model
l3.model_mid <- rma.mv(yi, vi, random = ~ 1 | study/experiment_id, data = droso_mid) # 

exp(c(full.model_mid$b, full.model_mid$ci.lb, full.model_mid$ci.ub))
exp(c(l3.model_mid$b, l3.model_mid$ci.lb, l3.model_mid$ci.ub))

var.comp(l3.model_mid)

# comparing full model 4-level model with study-level model indicates full model is equivalent or slightly worse
anova(full.model_mid, l3.model_mid)

rma.est_mid <- exp(c(full.model_mid$b,full.model_mid$ci.lb,full.model_mid$ci.ub))
percent((1/exp(full.model_mid$b))-1) # Effect size as percentage change compared to control
rma.est_mid.pval <- full.model_mid$pval
rma.est_mid.pval

png(file = "suppl_figures/rma_forestplot_repro_tox_over_2Vm_under_7Vm.png", width = 3000, height = 2500, res = 300)
metafor::forest(full.model_mid, annotate=T, addfit=T, addpred=T, showweights=T, atrans=exp, refline=0,
                at=log(c(0.25, 1, 4, 16)), xlim=c(-2.5,4.5), shade="zebra", slab = experiment, cex=0.75,
                mlab=mlabfun("Variance distrib.: ", l3.model_mid),
                header=c("Author and Year","Weight        ROM      95% CI"))
dev.off()

funnel.rma(full.model_mid)

res_mid <- rma(TE, sei=seTE, data=madata2)
inf <- influence(res_mid)
plot(res_mid)
regtest(res_mid)

# contour-enhanced funnel plot (credit to Wolfgang Viechtbauer)
funnel(full.model_mid, yaxis="seinv", refline = 0,
       xaxs="i", yaxs="i", las=1,
       level=c(.10, .05, .01), shade=c("white", "gray55", "gray75"),
       legend=TRUE, back="gray90", hlines=NULL, ylab="Precision (1/se)")


# Bayesian models

# "brms" model
brm2 <- brm(formula = yi | se(sei) ~ 1 + (1|study/experiment/experiment_id), data = droso_mid,
             iter = 40000, warmup = 1000,
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             core = 6, chains = 6)

summary(brm2)
exp(fixef(brm2))
percent((1/exp(fixef(brm2)[1]))-1) # Effect size as percentage change compared to control

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
bmr_mid_studylevel <- bayesmeta_studylevel(droso_mid)

exp(bmr_mid_studylevel$summary)
percent((1/exp(bmr_mid_studylevel$summary[3,2]))-1) # Effect size as percentage change compared to control
forestplot(bmr_mid_studylevel, expo=TRUE, clip=c(-1, 3))

png(file = "suppl_figures/bayesmeta_study-level_forestplot_repro_tox_over_2Vm_under_7Vm.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_mid_studylevel, expo=TRUE, clip=c(-1, 3))
dev.off()

#bmr_mid_studylevel$pposterior(mu = 0)
#bmr_mid_studylevel$likelihood(mu=0)
bmr_mid_studylevel$bayesfactor[1, "mu=0"]

# grab estimates
brm.est_mid <- c(fixef(brm2)[1],fixef(brm2)[3],fixef(brm2)[4])
brm.est_mid <- meta:::backtransf(brm.est_mid, sm="ROM")
bayes.est_mid <- c(ma2$summary["mode",2],ma2$summary["95% lower",2],ma2$summary["95% upper",2])
bayes.est_mid <- meta:::backtransf(bayes.est_mid, sm="ROM")
bayes.est_mid2 <- exp(bmr_mid_studylevel$summary[c(3,5,6),2])

# all estimates
meta.est_mid
rma.est_mid
brm.est_mid
bayes.est_mid
bayes.est_mid2


#####################################################################
# Meta-analysis of all RF Drosophila studies below 2 V/m
#####################################################################

# "meta" model
madata3 <- metagen(TE = log_ROM,
                   seTE = log_SE,
                   data = droso_low,
                   cluster = study,
                   comb.fixed = F,
                   comb.random = T,
                   prediction = TRUE,
                   sm = "ROM")

#View(droso_low)
meta::forest(madata3)

exp(madata3$TE.random)
percent((1/exp(madata3$TE.random))-1) # Effect size as percentage change over control
madata3$pval.random

meta::forest(madata3, layout = "RevMan5",
            prediction = TRUE, test.overall.random = TRUE,
            #file = "suppl_figures/forestplot_repro_tox_under_2Vm.png", width = 3300, rows.gr = 150, args.gr = list(res = 300), dev.off = T,
            leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
            leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
            just = "center", just.addcols = "center",
            text.addline1 = meta_I2_multilevel(madata3))

madata3 <- update(madata3, subgroup = study) # subgroups shown

meta::forest(madata3, layout = "RevMan5",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = T, print.subgroup.name = F,
             file = "suppl_figures/forestplot_repro_tox_under_2Vm_subgroups.png", width = 3300, rows.gr = 175, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
             leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
             just = "center", just.addcols = "center",
             text.addline1 = meta_I2_multilevel(madata3))

meta::funnel(madata3)
meta.est_low <- c(madata3$TE.random, madata3$lower.random, madata3$upper.random)
meta.est_low <- meta:::backtransf(meta.est_low, sm="ROM")



# four-level model using the "Metafor" package
full.model_low <- rma.mv(yi, vi, random = ~ 1 | study/experiment/experiment_id, tdist = T, data = droso_low) # correct full model
l3.model_low <- rma.mv(yi, vi, random = ~ 1 | study/experiment_id, data = droso_low) # 

exp(c(full.model_low$b, full.model_low$ci.lb, full.model_low$ci.ub))
exp(c(l3.model_low$b, l3.model_low$ci.lb, l3.model_low$ci.ub))

var.comp(l3.model_low)

# comparing full model 3-level model with study-level model indicates full model is slightly worse
anova(full.model_low, l3.model_low)

rma.est_low <- exp(c(full.model_low$b,full.model_low$ci.lb,full.model_low$ci.ub))
percent((1/exp(full.model_low$b))-1) # Effect size as percentage change compared to control
rma.est_low.pval <- full.model_low$pval
rma.est_low.pval

png(file = "suppl_figures/rma_forestplot_repro_tox_under_2Vm.png", width = 3100, height = 2200, res = 300)
metafor::forest(full.model_low, annotate=T, addfit=T, addpred=T, showweights=T, atrans=exp, refline=0,
                at=log(c(0.25, 1, 4, 16)), xlim=c(-2.5,4.5), shade="zebra", slab = experiment, cex=0.75,
                mlab=mlabfun("Variance distrib.: ", l3.model_low),
                header=c("Author and Year","Weight        ROM      95% CI"))
dev.off()

funnel.rma(full.model_low)

res_low <- rma(TE, sei=seTE, data=madata3)
inf <- influence(res_low)
plot(res_low)
regtest(res_low)

# contour-enhanced funnel plot (credit to Wolfgang Viechtbauer)
funnel(full.model_low, yaxis="seinv", refline = 0,
       xaxs="i", yaxs="i", las=1,
       level=c(.10, .05, .01), shade=c("white", "gray55", "gray75"),
       legend=TRUE, back="gray90", hlines=NULL, ylab="Precision (1/se)")



# Bayesian models

# "brms" model
brm3 <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study/experiment/experiment_id), data = droso_low,
             iter = 40000, warmup = 1000,
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             core = 6, chains = 6)

exp(fixef(brm3))
percent((1/exp(fixef(brm3)[1]))-1) # Effect size as percentage change compared to control

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
bmr_low_studylevel <- bayesmeta_studylevel(droso_low)

exp(bmr_low_studylevel$summary)
percent((1/exp(bmr_low_studylevel$summary[3,2]))-1) # Effect size as percentage change compared to control
forestplot(bmr_low_studylevel, expo=TRUE, clip=c(-1, 3))

png(file = "suppl_figures/bayesmeta_study-level_forestplot_repro_tox_under_2Vm.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_low_studylevel, expo=TRUE, clip=c(-1, 3))
dev.off()

# bmr_low_studylevel$pposterior(mu = 0)
# bmr_low_studylevel$likelihood(mu=0)
bmr_low_studylevel$bayesfactor[1, "mu=0"]

# grab estimates
brm.est_low <- c(fixef(brm3)[1],fixef(brm3)[3],fixef(brm3)[4])
brm.est_low <- meta:::backtransf(brm.est_low, sm="ROM")
bayes.est_low <- c(ma3$summary["mode",2],ma3$summary["95% lower",2],ma3$summary["95% upper",2])
bayes.est_low <- meta:::backtransf(bayes.est_low, sm="ROM")
bayes.est_low2 <- exp(bmr_low_studylevel$summary[c(3,5,6),2])

# all estimates
meta.est_low
rma.est_low
brm.est_low
bayes.est_low
bayes.est_low2



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

studies <- df2 %>% group_by(EMF_source) %>% summarise(n = n_distinct(study)) # sample sizes: number of studies
studies

tower <- sqldf("SELECT * FROM df WHERE EMF_source LIKE '%Base station%'")
# remove Thill2018 study
tower <- subset(tower, study != "Thill2018")
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

View(tower)
#View(tower2)
#View(tower3)
View(DECT)
View(mobile)
View(sig_gen)
View(wifi) # too few observations so far, even with newly added data
View(coilsystem)



######################################################################
# Meta-analysis of all EMF groups together: faster but less accurate
######################################################################

madata_all <- metagen(TE = log_ROM,
                      seTE = log_SE,
                      data = df2,
                      cluster = study,
                      common = F,
                      random = T,
                      prediction = TRUE,
                      sm = "ROM")


# various funnel plots
meta::funnel(madata_all)
# the observed outcomes (effect sizes) in two "rays" on the right side are those derived from "se.from.p" (for p = 0.5 or p = 0.05)
# this somewhat reduces the explanatory power of this funnel plot, hence no judgement on study bias will be made

res0 <- rma(TE, sei=seTE, data=madata_all)
regtest(res0)
funnel.rma(res0)

full.model_alldata <- rma.mv(yi, vi, random = ~ 1 | study/experiment/experiment_id, tdist = T, data = df2) # correct full model

# contour-enhanced funnel plot (credit to Wolfgang Viechtbauer)
funnel(full.model_alldata, yaxis="seinv", refline = 0,
       xaxs="i", yaxs="i", las=1,
       level=c(.10, .05, .01), shade=c("white", "gray55", "gray75"),
       legend=TRUE, back="gray90", hlines=NULL, ylab="Precision (1/se)", 
       ylim = c(0.5,100))



results <- update(madata_all, subgroup = EMF_source, tau.common = F, common = FALSE)

summary(results)
# Results for subgroups (random effects model):
#   k  ROM       95%-CI tau^2    tau       Q I^2
# EMF_source = Base station          50 1.60 [0.99, 2.60] 0.966 0.9827  639.24 92%
# EMF_source = RF signal generator   42 1.37 [1.14, 1.64] 0.275 0.5241 1231.44 97%
# EMF_source = Mobile phone          79 1.48 [1.35, 1.62] 0.091 0.3012 3046.87 97%
# EMF_source = WiFi                  11 1.49 [1.03, 2.16] 0.256 0.5056  175.43 94%
# EMF_source = DECT                  20 1.55 [1.31, 1.82] 0.128 0.3573  478.37 96%
# EMF_source = Coil system 50/60 Hz  39 1.35 [1.16, 1.56] 0.079 0.2803  123.07 69%


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
              comb.fixed = F,
              comb.random = T,
              prediction = TRUE,
              sm = "ROM")

meta::forest(m1, layout = "RevMan5",
            prediction = TRUE, test.overall.random = TRUE,
            #file = "suppl_figures/forestplot_base_station.png", width = 3300, rows.gr = 250, args.gr = list(res = 300), dev.off = T,
            leftcols = c("Author_Year","E_Field","CummHrs","Bioeffect_cat","effect","ci","w.random"),
            leftlabs = c("Author, Year","E-Field [V/m]","Duration [h]","Bioeffect category"),
            just = "center", just.addcols = "center", xlim = c(0.2,10),
            text.addline1 = meta_I2_multilevel(m1))


m1 <- update(m1, subgroup = study)

meta::forest(m1, layout = "RevMan5")

meta::forest(m1, layout = "RevMan5",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = T, print.subgroup.name = F,
             file = "suppl_figures/forestplot_base_station_subgroups.png", width = 4000, rows.gr = 275, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","Insect_type","EMF_source","E_Field","CummHrs","Bioeffect_cat","effect","ci","w.random"),
             leftlabs = c("Author, Year","Insect","EMF","[V/m]","Duration [h]","Bioeffect"),
             just = "center", just.addcols = "center", xlim = c(0.2,10),
             text.addline1 = meta_I2_multilevel(m1))

summary(m1)
meta_I2_multilevel(m1)

meta.est1 <- c(m1$TE.random, m1$lower.random, m1$upper.random)
meta.est1 <- meta:::backtransf(meta.est1, sm="ROM")
meta.est1
percent((1/exp(m1$TE.random))-1) # Effect size as percentage change over control
m1$pval.random

res1 <- rma(TE, sei=seTE, data=m1)
inf <- influence(res1)
plot(inf)
regtest(res1)

funnel.rma(res1)
funnel.rma(res1, label = "out", refline=0) # centered on effect size = 0


# four-level model using the "Metafor" package
full.model1 <- rma.mv(yi, vi, random = ~ 1 | study/experiment/experiment_id, tdist = T, data = tower) # correct full model
l3.model1 <- rma.mv(yi, vi, random = ~ 1 | study/experiment_id, data = tower) # 

exp(c(full.model1$b, full.model1$ci.lb, full.model1$ci.ub))
exp(c(l3.model1$b, l3.model1$ci.lb, l3.model1$ci.ub))

var.comp(l3.model1)

# comparing full model with study-level (3-level) model indicates full model is equivalent
anova(full.model1, l3.model1)

rma.est1 <- exp(c(full.model1$b, full.model1$ci.lb, full.model1$ci.ub))
percent((1/exp(full.model1$b))-1) # Effect size as percentage change compared to control
rma.est1.pval <- full.model1$pval


png(file = "suppl_figures/rma_forestplot_base_station.png", width = 3500, height = 3700, res = 300)
metafor::forest(full.model1, annotate=T, addfit=T, addpred=F, showweights=T, atrans=exp, refline=0,
                at=log(c(0.25, 1, 4, 16)), xlim=c(-2.5,4.5), shade="zebra", slab = tower$experiment, cex=0.75,
                mlab=mlabfun("Variance distrib.: ", l3.model1),
                header=c("Author and Year","Weight        ROM      95% CI"))
dev.off()


# Bayesian models

# "brms" model
brm_m1 <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study/experiment/experiment_id), data = tower,
             iter = 30000, warmup = 1000,
             control = list(adapt_delta = 0.98, max_treedepth = 15),
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


# "bayesmeta" model at study-level
bmr_1_studylevel <- bayesmeta_studylevel(tower)

exp(bmr_1_studylevel$summary)
forestplot(bmr_1_studylevel, expo=T, clip=c(-1, 3))

png(file = "suppl_figures/bayesmeta_study-level_forestplot_base_station.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_1_studylevel, expo=T, clip=c(-1, 3))
dev.off()

# bmr_1_studylevel$pposterior(mu = 0)
# bmr_1_studylevel$likelihood(mu=0)
bmr_1_studylevel$bayesfactor[1, "mu=0"]

# grabbing estimates of mean effect size + 95% CI
brm.est1 <- c(fixef(brm_m1)[1],fixef(brm_m1)[3],fixef(brm_m1)[4])
brm.est1 <- meta:::backtransf(brm.est1, sm="ROM")
bayes.est1 <- c(bm1$summary["mode",2], bm1$summary["95% lower",2], bm1$summary["95% upper",2])
bayes.est1 <- meta:::backtransf(bayes.est1, sm="ROM")
bayes.est1_studylvl <- exp(bmr_1_studylevel$summary[c(3,5,6),2])

meta.est1
rma.est1
brm.est1
bayes.est1
bayes.est1_studylvl


# Base station (without "reduced abundance" findings)
m1b <- metagen(TE = log_ROM,
              seTE = log_SE,
              data = tower2,
              cluster = study,
              comb.fixed = F,
              comb.random = T,
              prediction = TRUE,
              sm = "ROM")
meta::forest(m1b)

m1b <- update(m1b, subgroup = study)

meta::forest(m1b, layout = "RevMan5",
            prediction = TRUE, test.overall.random = TRUE,
            file = "suppl_figures/forestplot_base_station2_subgroups.png", width = 4200, rows.gr = 160, args.gr = list(res = 300), dev.off = T,
            #leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
            #leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
            leftcols = c("Author_Year","Insect_type","EMF_source","E_Field","CummHrs","Bioeffect_cat","effect","ci","w.random"),
            leftlabs = c("Author, Year","Insect","EMF","[V/m]","Duration [h]","Bioeffect"),
            just = "center", just.addcols = "center", xlim = c(0.2,5),
            text.addline1 = meta_I2_multilevel(m1b))

meta.est1b <- c(m1b$TE.random,m1b$lower.random,m1b$upper.random)
meta.est1b <- meta:::backtransf(meta.est1b, sm="ROM")
percent((1/exp(m1b$TE.random))-1) # Effect size as percentage change versus control
m1b$pval.random

res1b <- rma(TE, sei=seTE, data=m1b)
regtest(res1b)

meta::funnel(m1b)

# four-level model using the "Metafor" package
full.model1b <- rma.mv(yi, vi, random = ~ 1 | study/experiment/experiment_id, tdist = T, data = tower2) # correct full model
l3.model1b <- rma.mv(yi, vi, random = ~ 1 | study/experiment_id, data = tower2) # 

exp(c(full.model1b$b, full.model1b$ci.lb, full.model1b$ci.ub))
exp(c(l3.model1b$b, l3.model1b$ci.lb, l3.model1b$ci.ub))

var.comp(l3.model1b)

# comparing full model with study-level (3-level) model indicates full model is equivalent
anova(full.model1b, l3.model1b)

rma.est1b <- exp(c(full.model1b$b, full.model1b$ci.lb, full.model1b$ci.ub))
rma.est1b.pval <- full.model1b$pval


png(file = "suppl_figures/rma_forestplot_base_station2.png", width = 3000, height = 3500, res = 300)
metafor::forest(full.model1b, annotate=T, addfit=T, addpred=T, showweights=T, atrans=exp, refline=0,
                at=log(c(0.25, 1, 4, 16)), xlim=c(-2.5,4.5), shade="zebra", slab = tower2$experiment, cex=0.75,
                mlab=mlabfun("Variance distrib.: ", l3.model1b),
                header=c("Author and Year","Weight        ROM      95% CI"))
dev.off()


# Bayesian models
# "brms" model
brm_m1b <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study/experiment/experiment_id), data = tower2,
              iter = 60000, warmup = 1000,
              control = list(adapt_delta = 0.999, max_treedepth = 15),
              core = 6, chains = 6)

exp(fixef(brm_m1b))
#pairs(brm_m1b)

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


# "bayesmeta" study-level model
bmr_1b_studylevel <- bayesmeta_studylevel(tower2)

exp(bmr_1b_studylevel$summary)
forestplot(bmr_1b_studylevel, expo=T,)

png(file = "suppl_figures/bayesmeta_study-level_forestplot_base_station2.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_1b_studylevel, expo=T, clip=c(-1, 3))
dev.off()

#bmr_1b_studylevel$likelihood(mu=0)
bmr_1b_studylevel$bayesfactor[1, "mu=0"]


# grabbing estimates of mean effect size + 95% CI
brm.est1b <- c(fixef(brm_m1b)[1],fixef(brm_m1b)[3],fixef(brm_m1b)[4])
brm.est1b <- meta:::backtransf(brm.est1b, sm="ROM")
bayes.est1b <- c(bm1b$summary["mode",2], bm1b$summary["95% lower",2], bm1b$summary["95% upper",2])
bayes.est1b <- meta:::backtransf(bayes.est1b, sm="ROM")
bayes.est1b_studylvl <- exp(bmr_1b_studylevel$summary[c(3,5,6),2])

# all estimates
meta.est1b
rma.est1b
bayes.est1b
bayes.est1b_studylvl
brm.est1b


# Base station (only "reduced abundance" findings)
m1c <- metagen(TE = log_ROM,
               seTE = log_SE,
               data = tower3,
               cluster = study,
               comb.fixed = F,
               comb.random = T,
               prediction = TRUE,
               sm = "ROM")
meta::forest(m1c)

m1c <- update(m1c, subgroup = study)

meta::forest(m1c, layout = "RevMan5",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = T, print.subgroup.name = T,
             file = "suppl_figures/forestplot_base_station3_subgroups.png", width = 4200, rows.gr = 160, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","Insect_type","EMF_source","E_Field","CummHrs","Bioeffect_cat","effect","ci","w.random"),
             leftlabs = c("Author, Year","Insect","EMF","[V/m]","Duration [h]","Bioeffect"),
             just = "center", just.addcols = "center", xlim = c(0.2,12),
             text.addline1 = meta_I2_multilevel(m1c))

meta.est1c <- c(m1c$TE.random,m1c$lower.random,m1c$upper.random)
meta.est1c <- meta:::backtransf(meta.est1c, sm="ROM")
meta_I2_multilevel(m1c)
percent((1/exp(m1c$TE.random))-1) # Effect size as percentage change over control
m1c$pval.random

meta::funnel(m1c)
res1c <- rma(TE, sei=seTE, data=m1c)
regtest(res1c)


# four-level model using the "Metafor" package
full.model1c <- rma.mv(yi, vi, random = ~ 1 | study/experiment/experiment_id, tdist = T, data = tower3) # correct full model
l3.model1c <- rma.mv(yi, vi, random = ~ 1 | study/experiment_id, data = tower3) # 

exp(c(full.model1c$b, full.model1c$ci.lb, full.model1c$ci.ub))
exp(c(l3.model1c$b, l3.model1c$ci.lb, l3.model1c$ci.ub))

var.comp(l3.model1c)

# comparing full model with study-level (3-level) model indicates full model is equivalent
anova(full.model1c, l3.model1c)

rma.est1c <- exp(c(full.model1c$b, full.model1c$ci.lb, full.model1c$ci.ub))
rma.est1c.pval <- full.model1c$pval

png(file = "suppl_figures/rma_forestplot_base_station3.png", width = 3500, height = 3000, res = 300)
metafor::forest(full.model1c, annotate=T, addfit=T, addpred=T, showweights=T, atrans=exp, refline=0,
                at=log(c(0.25, 1, 4, 16)), xlim=c(-2.5,4.5), shade="zebra", slab = tower3$experiment, cex=0.75,
                mlab=mlabfun("Variance distrib.: ", l3.model1c),
                header=c("Author and Year","Weight        ROM      95% CI"))
dev.off()


# Bayesian models
# "brms" model
brm_m1c <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study/experiment/experiment_id), data = tower3,
               iter = 55000, warmup = 1000,
               control = list(adapt_delta = 0.99, max_treedepth = 15),
               core = 6, chains = 6)

exp(fixef(brm_m1c))
#pairs(brm_m1c)

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


# "bayesmeta" study-level model
bmr_1c_studylevel <- bayesmeta_studylevel(tower3)

exp(bmr_1c_studylevel$summary)
forestplot(bmr_1c_studylevel, expo=T)

png(file = "suppl_figures/bayesmeta_study-level_forestplot_base_station3.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_1c_studylevel, expo=T, clip=c(-1, 3))
dev.off()


bmr_1c_studylevel$bayesfactor[1, "mu=0"]

# grab estimates
brm.est1c <- c(fixef(brm_m1c)[1],fixef(brm_m1c)[3],fixef(brm_m1c)[4])
brm.est1c <- meta:::backtransf(brm.est1c, sm="ROM")
bayes.est1c <- c(bm1c$summary["mode",2], bm1c$summary["95% lower",2], bm1c$summary["95% upper",2])
bayes.est1c <- meta:::backtransf(bayes.est1c, sm="ROM")
bayes.est1c_studylvl <- exp(bmr_1c_studylevel$summary[c(3,5,6),2])

# all estimates
meta.est1c
rma.est1c
brm.est1c
bayes.est1c
bayes.est1c_studylvl


##########
# DECT
##########
m2 <- metagen(TE = log_ROM,
              seTE = log_SE,
              data = DECT,
              cluster = study,
              #cluster = experiment,
              sm = "ROM")

summary(m2)
#View(DECT)
meta.est2 <- c(m2$TE.random, m2$lower.random, m2$upper.random)
meta.est2 <- meta:::backtransf(meta.est2, sm="ROM")
meta.est2

meta::forest(m2)

#png(file = "suppl_figures/forestplot_DECT.png", width = 4000, height = 2000, res = 300)
meta::forest(m2, layout = "RevMan5",
            prediction = TRUE, test.overall.random = TRUE,
            leftcols = c("Author_Year","E_Field","CummHrs","Insect_type","Bioeffects"),
            leftlabs = c("Author, Year","E-Field [V/m]","Duration [h]","Insect","Bioeffects"),
            just = "center", just.addcols = "center", xlim = c(0.25,4),
            text.addline1 = meta_I2_multilevel(m2) )
#dev.off()

m2 <- update(m2, subgroup = study) # optionally: subgroup = experiment

meta::forest(m2, layout = "RevMan5",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = T, print.subgroup.name = F,
             file = "suppl_figures/forestplot_DECT_subgroups.png", width = 4000, rows.gr = 130, args.gr = list(res = 300), dev.off = T,
             #leftcols = c("Author_Year","EMF_source","E_Field","CummHrs","effect","ci","w.random"),
             #leftlabs = c("Author, Year","EMF","E-Field [V/m]","Duration [h]"),
             leftcols = c("Author_Year","Insect_type","X","EMF_source","E_Field","CummHrs","Bioeffect_cat","effect","ci","w.random"),
             leftlabs = c("Author, Year","Insect","Details","EMF","[V/m]","Duration [h]","Bioeffect"),
             just = "center", just.addcols = "center", xlim = c(0.25,4),
             text.addline1 = meta_I2_multilevel(m2))

summary(m2)
meta.est2 <- c(m2$TE.random, m2$lower.random, m2$upper.random)
meta.est2 <- meta:::backtransf(meta.est2, sm="ROM")
percent((1/exp(m2$TE.random))-1) # Effect size as percentage change over control
m2$pval.random

meta::funnel(m2)
res2 <- rma(TE, sei=seTE, data=m2)
regtest(res2)


# four-level model using the "Metafor" package
full.model2 <- rma.mv(yi, vi, random = ~ 1 | study/experiment/experiment_id, tdist = T, data = DECT) # correct full model
l3.model2 <- rma.mv(yi, vi, random = ~ 1 | study/experiment_id, data = DECT) # 

exp(c(full.model2$b, full.model2$ci.lb, full.model2$ci.ub))
exp(c(l3.model2$b, l3.model2$ci.lb, l3.model2$ci.ub))

var.comp(l3.model2)

# comparing full model (4-level) with study-level (3-level) model indicates full model is equivalent
anova(full.model2, l3.model2)

profile(full.model2)
profile(l3.model2)


rma.est2 <- exp(c(full.model2$b, full.model2$ci.lb, full.model2$ci.ub))
rma.est2
rma.est2.pval <- full.model2$pval
full.model2$pval

png(file = "suppl_figures/rma_forestplot_DECT.png", width = 3000, height = 2500, res = 300)
metafor::forest(full.model2, annotate=T, addfit=T, addpred=T, showweights=T, atrans=exp, refline=0,
                at=log(c(0.25, 1, 4, 16)), xlim=c(-2.5,4.5), shade="zebra", slab = DECT$experiment, cex=0.75,
                mlab=mlabfun("Variance distrib.: ", l3.model2),
                header=c("Author and Year","Weight        ROM      95% CI"))
dev.off()


# contour-enhanced funnel plot (credit to Wolfgang Viechtbauer)
funnel(full.model2, yaxis="seinv", refline = 0,
       xaxs="i", yaxs="i", las=1,
       level=c(.10, .05, .01), shade=c("white", "gray55", "gray75"),
       legend=TRUE, back="gray90", hlines=NULL, ylab="Precision (1/se)")


# Bayesian models
# "brms" model

brm_m2 <- brm(formula = log_ROM | se(log_SE) ~ 1+ (1|study/experiment/experiment_id), data = DECT,
               iter = 60000, warmup = 1000,
               control = list(adapt_delta = 0.995, max_treedepth = 12),
               core = 6, chains = 6,
              )

exp(fixef(brm_m2))

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


# "bayesmeta" model at study-level
bmr_2_studylevel <- bayesmeta_studylevel(DECT)

exp(bmr_2_studylevel$summary)
forestplot(bmr_2_studylevel, expo=T, clip=c(-1, 3))

png(file = "suppl_figures/bayesmeta_study-level_forestplot_DECT.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_2_studylevel, expo=T, clip=c(-1, 3))
dev.off()

# bmr_2_studylevel$pposterior(mu = 0)
# bmr_2_studylevel$likelihood(mu = 0)
bmr_2_studylevel$bayesfactor[1, "mu=0"]

# # compute p values for effect greater or less than zero
p_detrimental_effect_DECT <- pppvalue(bmr_2_studylevel, parameter="mu", value=0, alternative = "greater", n=100)
p_detrimental_effect_DECT


# grab estimates
brm.est2 <- c(fixef(brm_m2)[1], fixef(brm_m2)[3], fixef(brm_m2)[4])
brm.est2 <- meta:::backtransf(brm.est2, sm="ROM")
bayes.est2 <- c(bm2$summary["mode",2], bm2$summary["95% lower",2], bm2$summary["95% upper",2])
bayes.est2 <- meta:::backtransf(bayes.est2, sm="ROM")
bayes.est2_studylvl <- exp(bmr_2_studylevel$summary[c(3,5,6),2])

# all estimates
meta.est2
rma.est2
bayes.est2
bayes.est2_studylvl
brm.est2


###############
# Mobile phone
###############

m3 <- metagen(TE = log_ROM,
              seTE = log_SE,
              data = mobile,
              cluster = study,
              sm = "ROM")

summary(m3)

#png(file = "suppl_figures/forestplot_mobile_phone.png", width = 3500, height = 5500, res = 300)
meta::forest(m3, layout = "RevMan5",
            prediction = TRUE, test.overall.random = TRUE,
            leftcols = c("Author_Year","E_Field","CummHrs","Bioeffect_cat","effect","ci","w.random"),
            leftlabs = c("Author, Year","E-Field [V/m]","Duration [h]","Bioeffect category"),
            just = "center", just.addcols = "center",
            text.addline1 = meta_I2_multilevel(m3))
#dev.off()

m3 <- update(m3, subgroup = study)

meta::forest(m3, layout = "RevMan5",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = FALSE, print.subgroup.name = FALSE,
             file = "suppl_figures/forestplot_mobile_phone_subgroups.png", width = 4000, rows.gr = 340, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","Insect_type","EMF_source","E_Field","CummHrs","Bioeffect_cat","effect","ci","w.random"),
             leftlabs = c("Author, Year","Insect","EMF","[V/m]","Duration [h]","Bioeffect"),
             just = "center", just.addcols = "center",
             text.addline1 = meta_I2_multilevel(m3))


summary(m3)
meta_I2_multilevel(m3)

m3$I2
m3$tau2
m3$I2 * m3$tau2[1]/(m3$tau2[1]+m3$tau2[2])
m3$I2 * m3$tau2[2]/(m3$tau2[1]+m3$tau2[2])
# approximately the same as meta_I2_multilevel


meta.est3 <- c(m3$TE.random, m3$lower.random, m3$upper.random)
meta.est3 <- meta:::backtransf(meta.est3, sm="ROM")
m3$pval.random

meta::funnel(m3)
res3 <- rma(TE, sei=seTE, data=m3)
regtest(res3)


# four-level model using the "Metafor" package
full.model3 <- rma.mv(yi, vi, random = ~ 1 | study/experiment/experiment_id, data = mobile) # correct full model
l3.model3 <- rma.mv(yi, vi, random = ~ 1 | study/experiment_id, data = mobile) # 

exp(c(full.model3$b, full.model3$ci.lb, full.model3$ci.ub))
exp(c(l3.model3$b, l3.model3$ci.lb, l3.model3$ci.ub))

var.comp(l3.model3)

# comparing full  model (4-level) with study-level model indicates full model is equivalent
anova(full.model3, l3.model3)


rma.est3 <- exp(c(full.model3$b, full.model3$ci.lb, full.model3$ci.ub))
rma.est3
rma.est3.pval <- full.model3$pval
rma.est3.pval

png(file = "suppl_figures/rma_forestplot_mobile_phone.png", width = 3500, height = 5500, res = 300)
metafor::forest(full.model3, annotate=T, addfit=T, addpred=T, showweights=T, atrans=exp, refline=0,
                at=log(c(0.25, 1, 4)), xlim=c(-1,2), shade="zebra", slab = mobile$experiment, cex=0.75,
                mlab=mlabfun("Variance distrib.: ", l3.model3),
                header=c("Author and Year","Weight        ROM      95% CI"))
dev.off()

# contour-enhanced funnel plot (credit to Wolfgang Viechtbauer)
funnel(full.model3, yaxis="seinv", refline = 0,
       xaxs="i", yaxs="i", las=1,
       level=c(.10, .05, .01), shade=c("white", "gray55", "gray75"),
       legend=TRUE, back="gray90", hlines=NULL, ylab="Precision (1/se)")


# Bayesian models
# "brms" model
brm_m3 <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study/experiment/experiment_id), data = mobile,
              iter = 30000, warmup = 1000,
              control = list(adapt_delta = 0.95),
              core = 6, chains = 6)

exp(fixef(brm_m3))

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


# study-level model
bmr_3_studylevel <- bayesmeta_studylevel(mobile)

exp(bmr_3_studylevel$summary)
forestplot(bmr_3_studylevel, expo=T)

png(file = "suppl_figures/bayesmeta_study-level_forestplot_mobile_phone.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_3_studylevel, expo=T, clip=c(-1, 3))
dev.off()

# bmr_3_studylevel$likelihood(mu = 0)
# bmr_3_studylevel$pposterior(mu = 0)
bmr_3_studylevel$bayesfactor[1, "mu=0"]

# # compute p values for effect greater or less than zero
p_detrimental_effect_mobile <- pppvalue(bmr_3_studylevel, parameter="mu", value=0, alternative = "greater", n=100)
p_detrimental_effect_mobile

# grab estimates
brm.est3 <- c(fixef(brm_m3)[1], fixef(brm_m3)[3], fixef(brm_m3)[4])
brm.est3 <- meta:::backtransf(brm.est3, sm="ROM")
bayes.est3 <- c(bm3$summary["mode",2], bm3$summary["95% lower",2], bm3$summary["95% upper",2])
bayes.est3 <- meta:::backtransf(bayes.est3, sm="ROM")
bayes.est3_studylvl <- exp(bmr_3_studylevel$summary[c(3,5,6),2])

# all estimates
meta.est3
rma.est3
brm.est3
bayes.est3
bayes.est3_studylvl



###################
# Signal generator
###################
m4 <- metagen(TE = log_ROM,
              seTE = log_SE,
              data = sig_gen,
              #studlab = paste(key),
              cluster = study,
              sm = "ROM")

m4
meta::forest(m4)

m4 <- update(m4, subgroup = study)

meta::forest(m4, layout = "RevMan5",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = TRUE, print.subgroup.name = F,
             file = "suppl_figures/forestplot_signal_generator_subgroups.png", width = 4200, rows.gr = 270, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","Insect_type","EMF_source","E_Field","CummHrs","Bioeffect_cat","effect","ci","w.random"),
             leftlabs = c("Author, Year","Insect","EMF","[V/m]","Duration [h]","Bioeffect"),
             just = "center", just.addcols = "center", xlim = c(0.25,5),
             text.addline1 = meta_I2_multilevel(m4))


summary(m4)
#View(sig_gen)
meta_I2_multilevel(m4)

meta.est4 <- c(m4$TE.random, m4$lower.random, m4$upper.random)
meta.est4 <- meta:::backtransf(meta.est4, sm="ROM")
meta.est4

meta::funnel(m4)
radial(m4)
res4 <- rma(TE, sei=seTE, data=m4)
regtest(res4)


# four-level model using the "Metafor" package
full.model4 <- rma.mv(yi, vi, random = ~ 1 | study/experiment/experiment_id, data = sig_gen) # correct full model
l3.model4 <- rma.mv(yi, vi, random = ~ 1 | study/experiment_id, data = sig_gen) # 

exp(c(full.model4$b, full.model4$ci.lb, full.model4$ci.ub))
exp(c(l3.model4$b, l3.model4$ci.lb, l3.model4$ci.ub))

var.comp(l3.model4)

# comparing full model 4-level model with 3-level model indicates full model is equivalent
anova(full.model4, l3.model4)

rma.est4 <- exp(c(full.model4$b, full.model4$ci.lb, full.model4$ci.ub))
rma.est4
rma.est4.pval <- full.model4$pval
rma.est4.pval

png(file = "suppl_figures/rma_forestplot_signal_generator.png", width = 3500, height = 2500, res = 300)
metafor::forest(full.model4, annotate=T, addfit=T, addpred=T, showweights=T, atrans=exp, refline=0,
                at=log(c(0.25, 1, 4, 16)), xlim=c(-2,4), shade="zebra", slab = sig_gen$experiment, cex=0.75,
                mlab=mlabfun("Variance distrib.: ", l3.model4),
                header=c("Author and Year","Weight        ROM      95% CI"))
dev.off()


# contour-enhanced funnel plot (credit to Wolfgang Viechtbauer)
funnel(full.model4, yaxis="seinv", refline = 0,
       xaxs="i", yaxs="i", las=1,
       level=c(.10, .05, .01), shade=c("white", "gray55", "gray75"),
       legend=TRUE, back="gray90", hlines=NULL, ylab="Precision (1/se)")


# Bayesian models
# "brms" model
brm_m4 <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study/experiment/experiment_id), data = sig_gen,
              iter = 30000, warmup = 1000,
              control = list(adapt_delta = 0.99),
              core = 6, chains = 6)

exp(fixef(brm_m4))
#Intercept 1.384759  1.125584 1.08561 1.741763

gg_forestplot(brm_m4) + xlim(0.5, 2.5)

ggsave("suppl_figures/brms_forestplot_signal_generator.png", width = 12, height = 6)


# "bayesmeta" models
bm4 <- bayesmeta(y = sig_gen$log_ROM, sigma = sig_gen$log_SE,
                 labels = sig_gen$study, mu.prior.mean = 0, mu.prior.sd = 4,
                 tau.prior = function(t) dhalfnormal(t, scale = 0.5))
forestplot(bm4, expo=T)

png(file = "suppl_figures/bayesmeta_forestplot_signal_generator.png", width = 3500, height = 2500, res = 300)
forestplot(bm4, expo=T, clip=c(-1, 3))
dev.off()


# study-level model
bmr_4_studylevel <- bayesmeta_studylevel(sig_gen)

exp(bmr_4_studylevel$summary)
forestplot(bmr_4_studylevel, expo=T)

png(file = "suppl_figures/bayesmeta_study-level_forestplot_signal_generator.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_4_studylevel, expo=T, clip=c(-1, 3))
dev.off()

# bmr_4_studylevel$likelihood(mu = 0)
# bmr_4_studylevel$pposterior(mu = 0)
bmr_4_studylevel$bayesfactor[1, "mu=0"]

# grab estimates
brm.est4 <- c(fixef(brm_m4)[1], fixef(brm_m4)[3], fixef(brm_m4)[4])
brm.est4 <- meta:::backtransf(brm.est4, sm="ROM")
bayes.est4 <- c(bm4$summary["mode",2], bm4$summary["95% lower",2], bm4$summary["95% upper",2])
bayes.est4 <- meta:::backtransf(bayes.est4, sm="ROM")
bayes.est4_studylvl <- exp(bmr_4_studylevel$summary[c(3,5,6),2])

# all estimates
meta.est4
rma.est4
brm.est4
bayes.est4
bayes.est4_studylvl



###########################
# WiFi
###########################
m5 <- metagen(TE = log_ROM,
              seTE = log_SE,
              data = wifi,
              #studlab = paste(key),
              cluster = study,
              sm = "ROM")

m5
meta::forest(m5)

m5 <- update(m5, subgroup = study)

meta::forest(m5, layout = "RevMan5",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = F, print.subgroup.name = F,
             file = "suppl_figures/forestplot_wifi_subgroups.png", width = 4000, rows.gr = 100, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","Insect_type","EMF_source","E_Field","CummHrs","Bioeffect_cat","effect","ci","w.random"),
             leftlabs = c("Author, Year","Insect","EMF","[V/m]","Duration [h]","Bioeffect"),
             just = "center", just.addcols = "center")

summary(m5)


meta.est5 <- c(m5$TE.random, m5$lower.random, m5$upper.random)
meta.est5 <- meta:::backtransf(meta.est5, sm="ROM")
meta::funnel(m5)
res5 <- rma(TE, sei=seTE, data=m5)
regtest(res5)


# four-level model using the "Metafor" package
full.model5 <- rma.mv(yi, vi, random = ~ 1 | study/experiment/experiment_id, data = wifi) # correct full model
l3.model5 <- rma.mv(yi, vi, random = ~ 1 | study/experiment_id, data = wifi) # 

exp(c(full.model5$b, full.model5$ci.lb, full.model5$ci.ub))
exp(c(l3.model5$b, l3.model5$ci.lb, l3.model5$ci.ub))

var.comp(l3.model5)

# comparing full model 4-level model with study-level model indicates full model is not superior: higher AIC and BIC
anova(full.model5, l3.model5)


rma.est5 <- exp(c(full.model5$b, full.model5$ci.lb, full.model5$ci.ub))
rma.est5
rma.est5.pval <- full.model5$pval
rma.est5.pval

# Bayesian models
#"brms" model

brm_m5 <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study/experiment/experiment_id), data = wifi,
              iter = 40000, warmup = 1000,
              control = list(adapt_delta = 0.999, max_treedepth = 15),
              core = 6, chains = 6)

exp(fixef(brm_m5))

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


# "bayesmeta" study-level model
bmr_5_studylevel <- bayesmeta_studylevel(wifi)

exp(bmr_5_studylevel$summary)
forestplot(bmr_5_studylevel, expo=T)

png(file = "suppl_figures/bayesmeta_study-level_forestplot_wifi.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_5_studylevel, expo=T, clip=c(-1, 3))
dev.off()

# bmr_5_studylevel$likelihood(mu = 0)
# bmr_5_studylevel$pposterior(mu = 0)
bmr_5_studylevel$bayesfactor[1, "mu=0"]

# grab estimates
brm.est5 <- c(fixef(brm_m5)[1], fixef(brm_m5)[3], fixef(brm_m5)[4])
brm.est5 <- meta:::backtransf(brm.est5, sm="ROM")
bayes.est5 <- c(bm5$summary["mode",2], bm5$summary["95% lower",2], bm5$summary["95% upper",2])
bayes.est5 <- meta:::backtransf(bayes.est5, sm="ROM")
bayes.est5_studylvl <- exp(bmr_5_studylevel$summary[c(3,5,6),2])

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
              comb.fixed = F,
              comb.random = T,
              prediction = TRUE,
              sm = "ROM")

meta::forest(m6)

m6 <- update(m6, subgroup = study)

meta::forest(m6, layout = "RevMan5",
             prediction = TRUE, test.overall.random = TRUE,
             print.subgroup.labels = FALSE, print.subgroup.name = FALSE,
             file = "suppl_figures/forestplot_coil_system_subgroups.png", width = 4300, rows.gr = 230, args.gr = list(res = 300), dev.off = T,
             leftcols = c("Author_Year","Insect_type","EMF_source","B_Field","CummHrs","Bioeffect_cat","effect","ci","w.random"),
             leftlabs = c("Author, Year","Insect","EMF","[uT]","Duration [h]","Bioeffect"),
             just = "center", just.addcols = "center", xlim = c(0.25,5),
             text.addline1 = meta_I2_multilevel(m6))

summary(m6)
#View(coilsystem)
meta_I2_multilevel(m6)

meta.est6 <- c(m6$TE.random, m6$lower.random, m6$upper.random)
meta.est6 <- meta:::backtransf(meta.est6, sm="ROM")

meta::funnel(m6)
res6 <- rma(TE, sei=seTE, data=m6)
regtest(res6)

# four-level model using the "Metafor" package
full.model6 <- rma.mv(yi, vi, random = ~ 1 | study/experiment/experiment_id, data = coilsystem) # correct full model
l3.model6 <- rma.mv(yi, vi, random = ~ 1 | study/experiment_id, data = coilsystem) # 

exp(c(full.model6$b, full.model6$ci.lb, full.model6$ci.ub))
exp(c(l3.model6$b, l3.model6$ci.lb, l3.model6$ci.ub))

var.comp(l3.model6)

# comparing full model (4-level) with study-level (3-level) model indicates full model is equivalent
anova(full.model6, l3.model6)


rma.est6 <- exp(c(full.model6$b, full.model6$ci.lb, full.model6$ci.ub))
rma.est6
rma.est6.pval <- full.model6$pval
rma.est6.pval

png(file = "suppl_figures/rma_forestplot_coil_system.png", width = 3500, height = 2500, res = 300)
metafor::forest(full.model6, annotate=T, addfit=T, addpred=T, showweights=T, atrans=exp, refline=0,
                at=log(c(0.25, 1, 4, 16)), xlim=c(-2,4), shade="zebra", slab = coilsystem$experiment, cex=0.75,
                mlab=mlabfun("Variance distrib.: ", l3.model6),
                header=c("Author and Year","Weight        ROM      95% CI"))
dev.off()


# contour-enhanced funnel plot
funnel(full.model6, yaxis="seinv", refline = 0,
       xaxs="i", yaxs="i", las=1,
       level=c(.10, .05, .01), shade=c("white", "gray55", "gray75"),
       legend=TRUE, back="gray90", hlines=NULL, ylab="Precision (1/se)")


# Bayesian models
#"brms" model

brm_m6 <- brm(formula = log_ROM | se(log_SE) ~ 1 + (1|study/experiment/experiment_id), data = coilsystem,
              iter = 30000, warmup = 1000,
              control = list(adapt_delta = 0.99, max_treedepth = 15),
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


# "bayesmeta" study-level model
bmr_6_studylevel <- bayesmeta_studylevel(coilsystem)

exp(bmr_6_studylevel$summary)
forestplot(bmr_6_studylevel, expo=T, clip=c(-1, 5))

png(file = "suppl_figures/bayesmeta_study-level_forestplot_coil_system.png", width = 3300, height = 2800, res = 300)
forestplot(bmr_6_studylevel, expo=T, clip=c(-1, 5))
dev.off()

# bmr_6_studylevel$likelihood(mu = 0)
# bmr_6_studylevel$pposterior(mu = 0)
bmr_6_studylevel$bayesfactor[1, "mu=0"]

# grab estimates
brm.est6 <- c(fixef(brm_m6)[1], fixef(brm_m6)[3], fixef(brm_m6)[4])
brm.est6 <- meta:::backtransf(brm.est6, sm="ROM")
bayes.est6 <- c(bm6$summary["mode",2], bm6$summary["95% lower",2], bm6$summary["95% upper",2])
bayes.est6 <- meta:::backtransf(bayes.est6, sm="ROM")
bayes.est6_studylvl <- exp(bmr_6_studylevel$summary[c(3,5,6),2])

# all estimates
meta.est6
rma.est6
brm.est6
bayes.est6
bayes.est6_studylvl



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

rma_es_data <- data.frame(
  EMF_source = factor(c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system")),
  k = c(m1$k, m2$k, m3$k, m4$k, m5$k, m6$k),
  ROM = c(rma.est1[1],rma.est2[1],rma.est3[1],rma.est4[1],rma.est5[1],rma.est6[1]),
  LLCI = c(rma.est1[2],rma.est2[2],rma.est3[2],rma.est4[2],rma.est5[2],rma.est6[2]),
  ULCI = c(rma.est1[3],rma.est2[3],rma.est3[3],rma.est4[3],rma.est5[3],rma.est6[3]),
  estimate = c(1,1,1,1,1,1),
  pval = c(rma.est1.pval, rma.est2.pval, rma.est3.pval, rma.est4.pval, rma.est5.pval, rma.est6.pval))

bayes_es_data <- data.frame(
  EMF_source = factor(c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system")),
  k = c(m1$k, m2$k, m3$k, m4$k, m5$k, m6$k),
  ROM = c(bayes.est1[1],bayes.est2[1],bayes.est3[1],bayes.est4[1],bayes.est5[1],bayes.est6[1]),
  LLCI = c(bayes.est1[2],bayes.est2[2],bayes.est3[2],bayes.est4[2],bayes.est5[2],bayes.est6[2]),
  ULCI = c(bayes.est1[3],bayes.est2[3],bayes.est3[3],bayes.est4[3],bayes.est5[3],bayes.est6[3]),
  estimate = c(2,2,2,2,2,2),
  pval = c(bm1$bayesfactor[1, "mu=0"],bm2$bayesfactor[1, "mu=0"],bm3$bayesfactor[1, "mu=0"],bm4$bayesfactor[1, "mu=0"],bm5$bayesfactor[1, "mu=0"],bm6$bayesfactor[1, "mu=0"]))

bayes_es_studylvl_data <- data.frame(
  EMF_source = factor(c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system")),
  k = c(m1$k, m2$k, m3$k, m4$k, m5$k, m6$k),
  ROM = c(bayes.est1_studylvl[1],bayes.est2_studylvl[1],bayes.est3_studylvl[1],bayes.est4_studylvl[1],bayes.est5_studylvl[1],bayes.est6_studylvl[1]),
  LLCI = c(bayes.est1_studylvl[2],bayes.est2_studylvl[2],bayes.est3_studylvl[2],bayes.est4_studylvl[2],bayes.est5_studylvl[2],bayes.est6_studylvl[2]),
  ULCI = c(bayes.est1_studylvl[3],bayes.est2_studylvl[3],bayes.est3_studylvl[3],bayes.est4_studylvl[3],bayes.est5_studylvl[3],bayes.est6_studylvl[3]),
  estimate = c(3,3,3,3,3,3),
  pval = c(bmr_1_studylevel$bayesfactor[1, "mu=0"],bmr_2_studylevel$bayesfactor[1, "mu=0"],bmr_3_studylevel$bayesfactor[1, "mu=0"],bmr_4_studylevel$bayesfactor[1, "mu=0"],bmr_5_studylevel$bayesfactor[1, "mu=0"],bmr_6_studylevel$bayesfactor[1, "mu=0"]))

brm_es_data <- data.frame(
  EMF_source = factor(c("Base station","DECT","Mobile phone","Signal generator","WiFi","Coil system")),
  k = c(m1$k, m2$k, m3$k, m4$k, m5$k, m6$k),
  ROM = c(brm.est1[1],brm.est2[1],brm.est3[1],brm.est4[1],brm.est5[1],brm.est6[1]),
  LLCI = c(brm.est1[2],brm.est2[2],brm.est3[2],brm.est4[2],brm.est5[2],brm.est6[2]),
  ULCI = c(brm.est1[3],brm.est2[3],brm.est3[3],brm.est4[3],brm.est5[3],brm.est6[3]),
  estimate = c(4,4,4,4,4,4))

base_data <- data.frame(
  EMF_source = factor(c("Base station (tox)","Base station (tox)","Base station (tox)","Base station (tox)","Base station (tox)","Base station (behavior)","Base station (behavior)","Base station (behavior)","Base station (behavior)","Base station (behavior)")),
  k = c(m1b$k, m1b$k, m1b$k, m1b$k, m1b$k, m1c$k, m1c$k, m1c$k, m1c$k, m1c$k),
  ROM = c(meta.est1b[1],rma.est1b[1],bayes.est1b[1],bayes.est1b_studylvl[1],brm.est1b[1],meta.est1c[1],rma.est1c[1],bayes.est1c[1],bayes.est1c_studylvl[1],brm.est1c[1]),
  LLCI = c(meta.est1b[2],rma.est1b[2],bayes.est1b[2],bayes.est1b_studylvl[2],brm.est1b[2],meta.est1c[2],rma.est1c[2],bayes.est1c[2],bayes.est1c_studylvl[2],brm.est1c[2]),
  ULCI = c(meta.est1b[3],rma.est1b[3],bayes.est1b[3],bayes.est1b_studylvl[3],brm.est1b[3],meta.est1c[3],rma.est1c[3],bayes.est1c[3],bayes.est1c_studylvl[3],brm.est1c[3]),
  estimate = c(0,1,2,3,4, 0,1,2,3,4),
  pval = c(m1b$pval.random, rma.est1b.pval, NA, bmr_1b_studylevel$bayesfactor[1, "mu=0"], NA, m1c$pval.random, rma.est1c.pval, NA, bmr_1c_studylevel$bayesfactor[1, "mu=0"], NA))


es_data <- full_join(meta_es_data, base_data)
es_data <- full_join(es_data, rma_es_data)
es_data <- full_join(es_data, bayes_es_data)
es_data <- full_join(es_data, bayes_es_studylvl_data)
all_es_data <- full_join(es_data, brm_es_data)
all_es_data

# counting number of studies by EMF source and inputing in table, to merge with all_es_data
studies <- df2 %>% group_by(EMF_source) %>% summarise(n = n_distinct(study)) # sample sizes: number of studies
studies$EMF_source <- as.character(studies$EMF_source)
studies$EMF_source <- gsub("RF signal generator", "Signal generator", studies$EMF_source)
studies$EMF_source <- gsub("Coil system 50/60 Hz", "Coil system", studies$EMF_source)
new_row <- c(EMF_source = "Base station (tox)", tower2 %>% summarise(n = n_distinct(study)))
studies <- rbind(studies, new_row)

new_row <- c(EMF_source = "Base station (behavior)", tower3 %>% summarise(n = n_distinct(study)))
studies <- rbind(studies, new_row)

studies <- arrange(studies, EMF_source)
studies


all_es_data <- full_join(all_es_data, studies)
all_es_data

all_es_data <- arrange(all_es_data, EMF_source, estimate)

# rounding up digits to sensible values (this isn't nuclear physics)
all_es_data$ROM <- round(all_es_data$ROM, digits = 3)
all_es_data$LLCI <- round(all_es_data$LLCI, digits = 3)
all_es_data$ULCI <- round(all_es_data$ULCI, digits = 3)
all_es_data$pval <- round(all_es_data$pval, digits = 6)


all_es_data


# saving all estimates into spreadsheet
openxlsx::write.xlsx(all_es_data, "tables/all_effect_size_estimates_from_meta-analysis.xlsx", rowNames = F, colWidths = 15, overwrite = T, firstRow = T, firstCol = T)


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

meta_summary <- openxlsx::read.xlsx("tables/table_meta_summary.xlsx", sheet = 1, na.strings = "NA") 

pdf("tables/table_meta_summary.pdf", height=3, width=11) # save to pdf
gridExtra::grid.table(meta_summary, rows=NULL)
dev.off()

# latex code for table
latexcode <- xtable::xtable(meta_summary, caption = 'Summary table of meta-analysis results', auto = T)
print(latexcode, include.rownames = F)


#######################################################################################
# Summary table with results of three-level Bayesian meta-analysis (package bayesmeta)
#######################################################################################


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
`ROM` <- c(bayes.est_high2[1],bayes.est_mid2[1],bayes.est_low2[1],
           bayes.est1_studylvl[1],bayes.est1b_studylvl[1],bayes.est2_studylvl[1],bayes.est3_studylvl[1],bayes.est4_studylvl[1],bayes.est6_studylvl[1])
`ROM` <- round(`ROM`, digits = 2)
LLCI <- c(bayes.est_high2[2],bayes.est_mid2[2],bayes.est_low2[2],bayes.est1_studylvl[2],bayes.est1b_studylvl[2],bayes.est2_studylvl[2],bayes.est3_studylvl[2],bayes.est4_studylvl[2],bayes.est6_studylvl[2])
ULCI <- c(bayes.est_high2[3],bayes.est_mid2[3],bayes.est_low2[3],bayes.est1_studylvl[3],bayes.est1b_studylvl[3],bayes.est2_studylvl[3],bayes.est3_studylvl[3],bayes.est4_studylvl[3],bayes.est6_studylvl[3])
LLCI <- round(LLCI, digits=2)
ULCI <- round(ULCI, digits=2)
`95% CI` <- paste0(LLCI,"--",ULCI)
pvalue <- c(bmr_high_studylevel$bayesfactor[1, "mu=0"],bmr_mid_studylevel$bayesfactor[1, "mu=0"],bmr_low_studylevel$bayesfactor[1, "mu=0"], 
            bmr_1b_studylevel$bayesfactor[1, "mu=0"],bmr_1c_studylevel$bayesfactor[1, "mu=0"],bmr_2_studylevel$bayesfactor[1, "mu=0"],bmr_3_studylevel$bayesfactor[1, "mu=0"],bmr_4_studylevel$bayesfactor[1, "mu=0"],bmr_6_studylevel$bayesfactor[1, "mu=0"])
pvalue2 <- round(pvalue, digits=4)
pvalue2 <- as.character(pvalue2)

# chances for no effect given as bayes factor value (1 = 50/50 chance of no effect)
bmr_1_studylevel$bayesfactor[1, "mu=0"]
bmr_1b_studylevel$bayesfactor[1, "mu=0"]
bmr_1c_studylevel$bayesfactor[1, "mu=0"]
bmr_2_studylevel$bayesfactor[1, "mu=0"]
bmr_3_studylevel$bayesfactor[1, "mu=0"]
bmr_4_studylevel$bayesfactor[1, "mu=0"]
bmr_6_studylevel$bayesfactor[1, "mu=0"]


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

meta_summary <- openxlsx::read.xlsx("tables/table_bayesmeta_summary.xlsx", sheet = 1, na.strings = "NA") 

pdf("tables/table_bayesmeta_summary.pdf", height=3, width=11) # save to pdf
gridExtra::grid.table(bayesmeta_summary, rows=NULL)
dev.off()

# latex code for table
latexcode <- xtable::xtable(bayesmeta_summary, caption = 'Summary table of meta-analysis results', auto = T)
print(latexcode, include.rownames = F)



################################################################################
# Summary textfile of multilevel heterogeneity (for easy adding to latex file)
################################################################################

meta_I2_multilevel(madata1)
var.comp(l3.model_high)$results$I2[2]

paste_command_high <- list(bquote(paste("Meta-analysis of all HF Drosophila studies over 7 V/m.\n Variance distribution (meta):",
                        .(meta_I2_multilevel(madata1)), ".\n", 
                        "Variance distribution (metafor):", 
                        "Level 2 (within-study) I2 = ", .(var.comp(l3.model_high)$results$I2[2]), "%,",
                        "Level 3 (between-study) I2 = ", .(var.comp(l3.model_high)$results$I2[3]), "%."
                        )))

paste_command_mid <- list(bquote(paste("Meta-analysis of all HF Drosophila studies between 2 and 7 V/m.\n Variance distribution (meta):",
                                       .(meta_I2_multilevel(madata2)), ".\n", 
                                       "Variance distribution (metafor):", 
                                       "Level 2 (within-study) I2 = ", .(var.comp(l3.model_mid)$results$I2[2]), "%,",
                                       "Level 3 (between-study) I2 = ", .(var.comp(l3.model_mid)$results$I2[3]), "%."
                                       )))

paste_command_low <- list(bquote(paste("Meta-analysis of all HF Drosophila studies under 2 V/m.\n Variance distribution (meta):",
                                       .(meta_I2_multilevel(madata3)), ".\n", 
                                       "Variance distribution (metafor):", 
                                       "Level 2 (within-study) I2 = ", .(var.comp(l3.model_low)$results$I2[2]), "%,",
                                       "Level 3 (between-study) I2 = ", .(var.comp(l3.model_low)$results$I2[3]), "%."
                                       )))

paste_command1 <- list(bquote(paste("Base station (all studies).\n Variance distribution (meta):",
                                     .(meta_I2_multilevel(m1)), ".\n", 
                                     "Variance distribution (metafor):", 
                                     "Level 2 (within-study) I2 = ", .(var.comp(l3.model1)$results$I2[2]), "%,",
                                     "Level 3 (between-study) I2 = ", .(var.comp(l3.model1)$results$I2[3]), "%."
                                     )))

paste_command1b <- list(bquote(paste("Base station (only repro. tox. findings).\n Variance distribution (meta):",
                                    .(meta_I2_multilevel(m1b)), ".\n", 
                                    "Variance distribution (metafor):", 
                                    "Level 2 (within-study) I2 = ", .(var.comp(l3.model1b)$results$I2[2]), "%,",
                                    "Level 3 (between-study) I2 = ", .(var.comp(l3.model1b)$results$I2[3]), "%."
                                    )))

paste_command1c <- list(bquote(paste("Base station (only reduced abundance findings).\n Variance distribution (meta):",
                                    .(meta_I2_multilevel(m1c)), ".\n", 
                                    "Variance distribution (metafor):", 
                                    "Level 2 (within-study) I2 = ", .(var.comp(l3.model1c)$results$I2[2]), "%,",
                                    "Level 3 (between-study) I2 = ", .(var.comp(l3.model1c)$results$I2[3]), "%."
                                    )))

paste_command2 <- list(bquote(paste("DECT.\n Variance distribution (meta):",
                                    .(meta_I2_multilevel(m2)), ".\n", 
                                    "Variance distribution (metafor):", 
                                    "Level 2 (within-study) I2 = ", .(var.comp(l3.model2)$results$I2[2]), "%,",
                                    "Level 3 (between-study) I2 = ", .(var.comp(l3.model2)$results$I2[3]), "%."
                                    )))

paste_command3 <- list(bquote(paste("Mobile phone.\n Variance distribution (meta):",
                                    .(meta_I2_multilevel(m3)), ".\n", 
                                    "Variance distribution (metafor):", 
                                    "Level 2 (within-study) I2 = ", .(var.comp(l3.model3)$results$I2[2]), "%,",
                                    "Level 3 (between-study) I2 = ", .(var.comp(l3.model3)$results$I2[3]), "%."
                                    )))

paste_command4 <- list(bquote(paste("Signal generator.\n Variance distribution (meta):",
                                    .(meta_I2_multilevel(m4)), ".\n", 
                                    "Variance distribution (metafor):", 
                                    "Level 2 (within-study) I2 = ", .(var.comp(l3.model4)$results$I2[2]), "%,",
                                    "Level 3 (between-study) I2 = ", .(var.comp(l3.model4)$results$I2[3]), "%."
                                    )))

paste_command5 <- list(bquote(paste("WiFi.\n Variance distribution (meta):",
                                    .(meta_I2_multilevel(m5)), ".\n", 
                                    "Variance distribution (metafor):", 
                                    "Level 2 (within-study) I2 = ", .(var.comp(l3.model5)$results$I2[2]), "%,",
                                    "Level 3 (between-study) I2 = ", .(var.comp(l3.model5)$results$I2[3]), "%."
                                    )))

paste_command6 <- list(bquote(paste("Coil system.\n Variance distribution (meta):",
                                    .(meta_I2_multilevel(m6)), ".\n", 
                                    "Variance distribution (metafor):", 
                                    "Level 2 (within-study) I2 = ", .(var.comp(l3.model6)$results$I2[2]), "%,",
                                    "Level 3 (between-study) I2 = ", .(var.comp(l3.model6)$results$I2[3]), "%."
                                    )))


eval(parse(text = paste_command_high))
eval(parse(text = paste_command3))


heterogeneity_text <- c(
  eval(parse(text = paste_command_high)), 
  eval(parse(text = paste_command_mid)), 
  eval(parse(text = paste_command_low)),
  eval(parse(text = paste_command1)),   
  eval(parse(text = paste_command1b)), 
  eval(parse(text = paste_command1c)), 
  eval(parse(text = paste_command2)),  
  eval(parse(text = paste_command3)),
  eval(parse(text = paste_command4)),
  eval(parse(text = paste_command5)), 
  eval(parse(text = paste_command6))
       )

heterogeneity_text

# writing all the multilevel heterogeneity values into a text file, preformatted for easily copying into Latex
cat(heterogeneity_text, file="tables/multilevel_heterogeneity.txt", sep="\n\n")


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


