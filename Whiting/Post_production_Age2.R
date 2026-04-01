#############
## whiting fits
## diagnostics, plots, index
#############

library(sdmTMB)
library(sdmTMBextra)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(DHARMa)
library(patchwork)
library(purrr)

path <- "C:/Users/astroh/Desktop/Chapter 2/sdmTMB/WHG/"
plot_path <- "C:/Users/astroh/Desktop/Chapter 2/sdmTMB/WHG/all plots"

# load helper functions
source("funs_whg.R")

### Check AIC ------

load(paste0(path, "fits/2/2_aic.RData"))
aic <- aic[order(aic$AIC), ]
min(aic$AIC) # m_dlogn6

sink(file = paste0(path, "/fits/2/", "aic.txt"))
aic
sink(file = NULL)

### Read data ------

mods_2 <- readRDS(paste0(path, "fits/2/2_model_object.rds"))
mods_2$m_dlogn6[[1]]$formula

sanity(mods_2$m_dlogn_spattemp[[1]])
sanity(mods_2$m_dlogn6[[1]]) 

mod_spattemp <- mods_2$m_dlogn_spattemp[[1]]
mod_depth_subst <- mods_2$m_dlogn6[[1]]

sink(file = paste0(path, "/fits/2/", "base_model.txt"))
mod_spattemp

print("SDREPORT")
mod_spattemp$sd_report
sink(file = NULL)

sink(file = paste0(path, "/fits/2/", "cov_model.txt"))
mod_depth_subst

print("SDREPORT")
mod_depth_subst$sd_report
sink(file = NULL)


########## -------------------------
### GET OUTPUTS AND PRODUCTS
########## -------------------------

### Make predictions ------
load("prediction_grid.RData")
head(pred_grid)
pred_grid <- add_utm_columns(pred_grid, c("lon", "lat"), units = "km")
grid_yrs <- data.frame(X = pred_grid$X,
                       Y = pred_grid$Y,
                       middepth_scaled = pred_grid$middepth_scaled,
                       substrate2 = pred_grid$substrate2,
                       fyear = as.factor(pred_grid$fyear),
                       year = pred_grid$year)
grid_yrs <- distinct(grid_yrs)

# With spatial random effects and spatiotemporal rf in both predictors
predictions_spattemp <- predict(mod_spattemp, newdata = grid_yrs[, c("X", "Y", "fyear", "year")], return_tmb_object = TRUE)
# Including depth
predictions_d_s <- predict(mod_depth_subst, newdata = grid_yrs, return_tmb_object = TRUE)

preds <- list(predictions_spattemp, 
              predictions_d_s)
modelIndices <- purrr::map(preds, ~ list(
  get_index(.x, bias_correct = TRUE)
))

names(modelIndices) <- c("mod_spattemp", "mod_depth_subst")
modelIndices$mod_spattemp[[1]]$model <- "base"
modelIndices$mod_depth_subst[[1]]$model <- "base+depth+substrate"

all_indices <- as.data.frame(do.call(rbind, lapply(modelIndices, as.data.frame)))
rownames(all_indices) <- NULL

### Visualise indices
index_2 <- ggplot(all_indices, aes(year, est, colour = model), size=2) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model), alpha = 0.2) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  labs(fill = "Model", colour = "Model") +
  xlab('Year') + ylab('Biomass estimate (kg)') +
  theme_bw() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15)) +
  ggtitle("Abundance index by model approach: Whiting age 2")
ggsave(filename = paste0(path, "/fits/2/","logn_indices.jpg"), plot = index_2)

index_2_2 <- ggplot(all_indices, aes(year, est, colour = model, linetype = model), size=2) + 
  geom_line(linewidth = 2) +
  scale_color_manual(values = c("grey", "black")) +
  labs(fill = "Model", colour = "Model", linetype = "Model") +
  xlab('Year') + ylab('Relative Index of Abundance') +
  theme_classic() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15),
        legend.position = c(.3,.9),
        legend.background = element_rect(linewidth = 1, colour = "black"), 
        plot.margin = margin(5, 20, 5, 5))
ggsave(filename = paste0(path, "/fits/2/","indices2.jpg"), plot = index_2_2, 
       width = 2100, height = 2100, units = "px")
saveRDS(index_2_2, file = paste0(plot_path, "/2_index2.rds"))

all_indices$CI_width <- with(all_indices, upr - lwr)
index_2_error <- ggplot(all_indices, aes(year, CI_width, colour = model, linetype = model), size=2) + 
  geom_line(linewidth = 2) +
  #scale_color_manual(values = c("red", "blue", "green")) +
  scale_color_manual(values = c("grey", "black")) +
  #scale_color_manual(values = c("red", "blue")) +
  labs(fill = "Model", colour = "Model", linetype = "Model") +
  xlab('Year') + ylab('CI width of Index') +
  theme_bw() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15),
        legend.position = c(.7,.9),
        legend.background = element_rect(linewidth = 1, colour = "black"), 
        plot.margin = margin(5, 20, 5, 5))
ggsave(filename = paste0(path, "/fits/2/","2_index_CI.jpg"), plot = index_2_error)
saveRDS(index_2_error, file = paste0(plot_path, "/2_index_CI.rds"))

index_2_se <- ggplot(all_indices, aes(year, log(se), colour = model, linetype = model), size=2) + 
  geom_line(linewidth = 2) +
  scale_color_manual(values = c("grey", "black")) +
  labs(fill = "Model", colour = "Model", linetype = "Model") +
  xlab('Year') + ylab('Log-SE of the Index') +
  theme_classic() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15),
        legend.position = c(.3,.9),
        legend.background = element_rect(linewidth = 1, colour = "black"), 
        plot.margin = margin(5, 20, 5, 5))
ggsave(filename = paste0(path, "/fits/2/","indices_se.jpg"), plot = index_2_se, 
       width = 2100, height = 2100, units = "px")
saveRDS(index_2_se, file = paste0(plot_path, "/2_index_se.rds"))

#save(all_indices, file = paste0(path, "/fits/2/", "allindices.RData"))
load(paste0(path, "/fits/2/", "allindices.RData"))

### Internal consistency
abundance <- read.delim("C:/Users/astroh/Desktop/Chapter 2/sdmTMB/WHG/VAST whg 2025.txt", 
                        sep = "")
abundance <- abundance[-nrow(abundance),]
abundance$year <- 2003:2023 # my index data only estimates until 2023
abundance <- abundance[, -c(1,4)] # also remove VAST-based age 2 index
long <- tidyr::pivot_longer(abundance, names_to = "age", cols = 1:5)
colnames(long)[3] <- "index"
long$age <- as.numeric(gsub("X", "", long$age))
long <- long[order(long$age, long$year), ]

#base
base_sub <- all_indices[all_indices$model == "base", c("year", "est")]
base_index <- data.frame(year = base_sub$year, 
                         age = 2, 
                         index = base_sub$est)
df_base <- rbind(long[c(1:42),], base_index, long[c(43:nrow(long)),])
df_base$cohort <- with(df_base, year - age)
df_base_wide <- reshape::cast(df_base, cohort ~ age, value = "index")

png(filename = paste0(path, "/fits/2/", "intconsistency_base2.png"))
PerformanceAnalytics::chart.Correlation(df_base_wide, histogram=TRUE, pch=19)
dev.off()

# cov
cov_sub <- all_indices[all_indices$model != "base", c("year", "est")]
cov_index <- data.frame(year = cov_sub$year, 
                        age = 2, 
                        index = cov_sub$est)
df_cov <- rbind(long[c(1:42),], cov_index, long[c(43:nrow(long)),])
df_cov$cohort <- with(df_cov, year - age)
df_cov_wide <- reshape::cast(df_cov, cohort ~ age, value = "index")

png(filename = paste0(path, "/fits/2/", "intconsistency_cov2.png"))
PerformanceAnalytics::chart.Correlation(df_cov_wide, histogram=TRUE, pch=19)
dev.off()


### Make plot of error
base_rf_enc <- tidy(mod_spattemp, model = 1, effects = "ran_pars")
base_rf_pos <- tidy(mod_spattemp, model = 2, effects = "ran_pars")
t_base <- rbind(base_rf_enc, base_rf_pos)
colnames(t_base)[1] <- "model_comp"
t_base$model_comp[t_base$model_comp == 1] <- "Binomial component"
t_base$model_comp[t_base$model_comp == 2] <- "Lognormal component"
t_base$model <- "base"
head(t_base)

cov_rf_enc <- tidy(mod_depth_subst, model = 1, effects = "ran_pars")
cov_rf_pos <- tidy(mod_depth_subst, model = 2, effects = "ran_pars")
t_cov <- rbind(cov_rf_enc, cov_rf_pos)
colnames(t_cov)[1] <- "model_comp"
t_cov$model_comp[t_cov$model_comp == 1] <- "Binomial component"
t_cov$model_comp[t_cov$model_comp == 2] <- "Lognormal component"
t_cov$model <- "base+depth+substrate"
head(t_cov)

all_rfs <- rbind(t_base, t_cov)

all_rfs_sub <- subset(all_rfs, term %in% c("sigma_O", "sigma_E"))
all_rfs_sub$term[all_rfs_sub$term == "sigma_O"] <- "spatial"
all_rfs_sub$term[all_rfs_sub$term == "sigma_E"] <- "spatiotemporal"
save(all_rfs_sub, file = paste0(path, "/fits/2/", "rf_data.RData"))

rfs_plot2D <- ggplot(all_rfs_sub, aes(term, estimate, fill = model)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width=.2, position=position_dodge(.9), 
                colour = "grey40", linewidth = 2) +
  scale_fill_manual(values = c("grey", "black")) +
  labs(x = "Random effect", y = "SD Estimate", fill = "Model") +
  facet_wrap(~ model_comp) +
  theme_bw() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15),
        legend.position = "top",
        #legend.position = c(.8,.9),
        legend.background = element_rect(linewidth = 1, colour = "black"), 
        plot.margin = margin(5, 20, 5, 5))
ggsave(filename = paste0(path, "/fits/2/","rf_2D.jpg"), plot = rfs_plot2D, 
       width = 2100, height = 2100, units = "px")
saveRDS(rfs_plot2D, file = paste0(path, "fits/2",  "/whg_2_rf.rds"))

### MARGINAL model effects ----
### Base model
# year
y_base2a <- visreg_delta(mod_spattemp, xvar = "fyear", 
                         model = 1, scale = "response", gg = TRUE, 
                         partial = FALSE, rug = FALSE)
y_base2b <- visreg_delta(mod_spattemp, xvar = "fyear", 
                         model = 2, trans = exp, gg = TRUE, 
                         partial = FALSE, rug = FALSE)
y_base2adat <- visreg_delta(mod_spattemp, xvar = "fyear", 
                            model = 1, scale = "response", gg = TRUE, 
                            partial = FALSE, rug = FALSE, plot = FALSE)
y_base2bdat <- visreg_delta(mod_spattemp, xvar = "fyear", 
                            model = 2, trans = exp, gg = TRUE, 
                            partial = FALSE, rug = FALSE, plot = FALSE)
y_p_base <- format_mareff_plots(mareff1 = y_base2a, mareff2 = y_base2b)
ggsave(filename = paste0(path, "/fits/2/" ,"year_effect_base.jpg"), plot = y_p_base)

# Covariate model
# year
y_cova <- visreg_delta(mod_depth_subst, xvar = "fyear", model = 1, 
                       gg = TRUE)
y_covb <- visreg_delta(mod_depth_subst, xvar = "fyear", model = 2, 
                       gg = TRUE)
y_cov2a <- visreg_delta(mod_depth_subst, xvar = "fyear", 
                        model = 1, scale = "response", gg = TRUE, 
                        partial = FALSE, rug = FALSE)
y_cov2b <- visreg_delta(mod_depth_subst, xvar = "fyear", 
                        model = 2, trans = exp, gg = TRUE, 
                        partial = FALSE, rug = FALSE)
y_cov2adat <- visreg_delta(mod_depth_subst, xvar = "fyear", 
                           model = 1, scale = "response", gg = TRUE, 
                           partial = FALSE, rug = FALSE, plot = FALSE)
y_cov2bdat <- visreg_delta(mod_depth_subst, xvar = "fyear", 
                           model = 2, trans = exp, gg = TRUE, 
                           partial = FALSE, rug = FALSE, plot = FALSE)
y_p <- format_mareff_plots(mareff1 = y_cov2a, mareff2 = y_cov2b)
ggsave(filename = paste0(path, "/fits/2/", "year_effect_cov.jpg"), plot = y_p)

# depth
d_cova <- visreg_delta(mod_depth_subst, xvar = "middepth_scaled", model = 1,
                       gg = TRUE)
d_covb <- visreg_delta(mod_depth_subst, xvar = "middepth_scaled", model = 2,
                       gg = TRUE)
d_cov2a <- visreg_delta(mod_depth_subst, xvar = "middepth_scaled", 
                        model = 1, scale = "response", gg = TRUE, 
                        partial = FALSE, rug = FALSE)
d_cov2b <- visreg_delta(mod_depth_subst, xvar = "middepth_scaled", 
                        model = 2, trans = exp, gg = TRUE, 
                        partial = FALSE, rug = FALSE)
d_p <- format_mareff_plots(mareff1 = d_cov2a, mareff2 = d_cov2b)
ggsave(filename = paste0(path, "/fits/2/", "depth_effect.jpg"), plot = d_p, 
       width = 2100, height = 2100, units = "px")

# substrate
substr_cova <- visreg_delta(mod_depth_subst, xvar = "substrate2", model = 1,
                            gg = TRUE)
substr_covb <- visreg_delta(mod_depth_subst, xvar = "substrate2", model = 2,
                            gg = TRUE)
substr_cov2a <- visreg_delta(mod_depth_subst, xvar = "substrate2", 
                             model = 1, scale = "response", gg = TRUE, 
                             partial = FALSE, rug = FALSE)
substr_cov2b <- visreg_delta(mod_depth_subst, xvar = "substrate2", 
                             model = 2, trans = exp, gg = TRUE, 
                             partial = FALSE, rug = FALSE)
substr_p <- format_mareff_plots(mareff1 = substr_cov2a, mareff2 = substr_cov2b)
ggsave(filename = paste0(path, "/fits/2/", "substr_effect.jpg"), plot = substr_p)

# save all
mar_eff_list <- list("fyear2_base_m1" = y_base2a, "fyear2_base_m2" = y_base2b,
                     "fyear_m1" = y_cova, "fyear_m2" = y_covb,
                     "fyear2_m1" = y_cov2a, "fyear2_m2" = y_cov2b,
                     "depth_m1" = d_cova, "depth_m2" = d_covb,
                     "depth2_m1" = d_cov2a, "depth2_m2" = d_cov2b,
                     "substr_m1" = substr_cova, "substr_m2" = substr_covb,
                     "substr2_m1" = substr_cov2a, "substr2_m2" = substr_cov2b)
saveRDS(mar_eff_list, paste0(path, "/fits/2/", "allmarginaleffects.rds"))

year_eff_list <- list("fyear2_base_m1_dat" = y_base2adat, 
                      "fyear2_base_m2_dat" = y_base2bdat,
                      "fyear2_cov_m1_dat" = y_cov2adat, 
                      "fyear2_cov_m2_dat" = y_cov2bdat)
saveRDS(year_eff_list, paste0(path, "/fits/2/", "year_marginaleffects.rds"))

### Partial deviance explained (PDE) -----
mod_depth <- mods_2$m_dlogn3[[1]]
mod_substr <- mods_2$m_dlogn5[[1]]

null_model <- sdmTMB(
  data = mod_spattemp$data,
  list(biomass ~ 1,
       biomass ~ 1),
  mesh = mod_spattemp$spde,
  family = delta_lognormal(),
  offset = log(mod_spattemp$data$areakmsqadj)
)

sink(file = paste0(path, "/fits/2/", "pde_report.txt"))
deviance_report(full_cov_model = mod_depth_subst, null_model = null_model, 
                base_model = mod_spattemp, 
                reduced_models = list(mod_depth, mod_substr))
sink(file = NULL)

# Check significance of fixed covariates -----
# https://github.com/pbs-assess/sdmTMB/discussions/360
# https://github.com/pbs-assess/sdmTMB/discussions/228
# base - year effects
SDrepFixed_base <- as.data.frame(summary(mod_spattemp$sd_report, 
                                         p.value=TRUE, select = "fixed"))
SDrepFixed_base$Variable <- row.names(SDrepFixed_base)
SDrepFixed_base$Variable <- gsub("\\.[0-9]+", "", SDrepFixed_base$Variable)
row.names(SDrepFixed_base) <- NULL
table(SDrepFixed_base$Variable)

terms_base <- rbind(tidy(mod_spattemp, model = 1), 
                    tidy(mod_spattemp, model = 2))
print(terms_base, n = nrow(terms_base))
SDrepFixedFac_base <- SDrepFixed_base[SDrepFixed_base$Variable %in% c("b_j", "b_j2"),]
SDrepFixedFac_base$ModelComp <- ifelse(SDrepFixedFac_base$Variable == "b_j",
                                       1, 2)

nrow(SDrepFixedFac_base) == nrow(terms_base)
all_terms_base <- cbind(terms_base[, "term"], SDrepFixedFac_base[, c(1:4, 6)], 
                        terms_base[, c("conf.low", "conf.high")])
all_terms_base$sig001 <- ifelse(all_terms_base$`Pr(>|z^2|)` < 0.001,
                                "sig", "insig")
all_terms_base$sig05 <- ifelse(all_terms_base$`Pr(>|z^2|)` < 0.05,
                               "sig", "insig")
write.csv(all_terms_base, 
          file = paste0(path, "/fits/2/", "base_factor_level_pvalues.csv"))

# cov - year + covariates
# factor effects
SDrepFixed_cov <- as.data.frame(summary(mod_depth_subst$sd_report, 
                                        p.value=TRUE, select = "fixed"))
SDrepFixed_cov$Variable <- row.names(SDrepFixed_cov)
SDrepFixed_cov$Variable <- gsub("\\.[0-9]+", "", SDrepFixed_cov$Variable)
row.names(SDrepFixed_cov) <- NULL
table(SDrepFixed_cov$Variable)
SDrepFixedFac_cov <- SDrepFixed_cov[SDrepFixed_cov$Variable %in% c("b_j", "b_j2"),]
SDrepFixedFac_cov$ModelComp <- ifelse(SDrepFixedFac_cov$Variable == "b_j", 
                                      1, 2)

terms_cov <- rbind(tidy(mod_depth_subst, model = 1), 
                   tidy(mod_depth_subst, model = 2))
print(terms_cov, n = nrow(terms_cov))

nrow(SDrepFixedFac_cov) == nrow(terms_cov)
all_terms_cov <- cbind(terms_cov[, "term"], SDrepFixedFac_cov[, c(1:4, 6)], 
                       terms_cov[, c("conf.low", "conf.high")])
all_terms_cov$sig001 <- ifelse(all_terms_cov$`Pr(>|z^2|)` < 0.001,
                               "sig", "insig")
all_terms_cov$sig05 <- ifelse(all_terms_cov$`Pr(>|z^2|)` < 0.05,
                              "sig", "insig")
write.csv(all_terms_cov, 
          file = paste0(path, "/fits/2/", "cov_factor_level_pvalues.csv"))


### Internal consistency -----
abundance <- read.delim("C:/Users/astroh/Desktop/Chapter 2/sdmTMB/WHG/VAST whg 2025.txt", 
                        sep = "")
abundance <- abundance[-nrow(abundance),]
abundance$year <- 2003:2023 # my index data only estimates until 2023
abundance <- abundance[, -c(1,4)] # also remove VAST-based age 0 index
long <- tidyr::pivot_longer(abundance, names_to = "age", cols = 1:5)
colnames(long)[3] <- "index"
long$age <- as.numeric(gsub("X", "", long$age))
long <- long[order(long$age, long$year), ]

#base
base_sub <- all_indices[all_indices$model == "base", c("year", "est")]
base_index <- data.frame(year = base_sub$year, 
                         age = 1, 
                         index = base_sub$est)
df_base <- rbind(base_index, long)
df_base$cohort <- with(df_base, year - age)
df_base_wide <- reshape::cast(df_base, cohort ~ age, value = "index")

png(filename = paste0(path, "/fits/2/", "intconsistency_base2.png"))
PerformanceAnalytics::chart.Correlation(df_base_wide, histogram=TRUE, pch=19)
dev.off()

# cov
cov_sub <- all_indices[all_indices$model != "base", c("year", "est")]
cov_index <- data.frame(year = cov_sub$year, 
                        age = 1, 
                        index = cov_sub$est)
df_cov <- rbind(cov_index, long)
df_cov$cohort <- with(df_cov, year - age)
df_cov_wide <- reshape::cast(df_cov, cohort ~ age, value = "index")

png(filename = paste0(path, "/fits/2/", "intconsistency_cov2.png"))
PerformanceAnalytics::chart.Correlation(df_cov_wide, histogram=TRUE, pch=19)
dev.off()

consistency_dat <- list("base" = df_base_wide, "cov" = df_cov_wide)
saveRDS(consistency_dat, file = paste0(path, "/fits/2/", "consistency_data.rds"))


########## -------------------------
### TESTING MODEL DESIGN AND FIT
########## -------------------------

### Crossvalidation ------

# separate data in test and training folds
set.seed(5)
seeds <- sample(1:1500, 20, replace = TRUE)
partitions <- holdout(model = mod_spattemp, 
                      training_cutoff = 0.7, seeds = seeds)

# run CVs and calculate evaluation scores
baseCVs <- map_dfr(partitions, 
                   rep_cv, formula = "base", 
                   .id = "randomPartition_nr")
head(baseCVs)
colnames(baseCVs)[2:3] <- c("RMSE_base", "MAE_base")

covCVs <- map_dfr(partitions, 
                  rep_cv, formula = "depth_substrate", 
                  .id = "randomPartition_nr")
head(covCVs)
colnames(covCVs)[2:3] <- c("RMSE_cov", "MAE_cov")

CV_results <- merge(baseCVs, covCVs)
CV_results$RMSE_dif <- CV_results$RMSE_cov - CV_results$RMSE_base
CV_results$MAE_dif <- CV_results$MAE_cov - CV_results$MAE_base
CV_results$RMSE_dif_mean <- mean(CV_results$RMSE_dif)
CV_results$MAE_dif_mean <- mean(CV_results$MAE_dif)

write.table(CV_results, 
            file = paste0(path, "fits/2/", "crossvalidation_results.txt"), 
            append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = TRUE)

### Check residuals ------
#Looking at tails
# randomised quantile residuals - MCMC-MLE
mod_base_res1 <- residuals(mod_spattemp, type = "mle-mvn", model = 1)
mod_base_res2 <- residuals(mod_spattemp, type = "mle-mvn", model = 2)
mod_cov_res1 <- residuals(mod_depth_subst, type = "mle-mvn", model = 1)
mod_cov_res2 <- residuals(mod_depth_subst, type = "mle-mvn", model = 2)

png(paste0(path, "fits/2/", "randomised_quantile_base_bestfit.png"))
par(mfrow = c(2,2))
qqnorm(mod_base_res1, main = "Base model - binomial: Normal Q-Q Plot");abline(0, 1)
qqnorm(mod_base_res2, main = "Base model - lognormal: Normal Q-Q Plot");abline(0, 1)
qqnorm(mod_cov_res1, main = "Covariate model - binomial: Normal Q-Q Plot");abline(0, 1)
qqnorm(mod_cov_res2, main = "Covariate model - lognormal: Normal Q-Q Plot");abline(0, 1)
dev.off()
# simulation-based - uniform
s_mod_base <- make_dharma(model = mod_spattemp, model_comp = c(1, 2))
s_mod_cov <- make_dharma(model = mod_depth_subst, model_comp = c(1, 2))

# Uniformity, dispersion, outliers
png(paste0(path, "fits/2/", "dharma_covbestfit_test_m1.png"))
sink(file = paste0(path, "fits/2/", "dharma_standardtests_covmod_m1.txt"))
DHARMa::testResiduals(s_mod_cov$m1)
sink(file = NULL)
dev.off()
png(paste0(path, "fits/2/", "dharma_covbestfit_test_m2.png"))
sink(file = paste0(path, "fits/2/", "dharma_standardtests_covmod_m2.txt"))
DHARMa::testResiduals(s_mod_cov$m2)
sink(file = NULL)
dev.off()
png(paste0(path, "fits/2/", "dharma_base_test_m1.png"))
sink(file = paste0(path, "fits/2/", "dharma_standardtests_basemod_m1.txt"))
DHARMa::testResiduals(s_mod_base$m1)
sink(file = NULL)
dev.off()
png(paste0(path, "fits/2/", "dharma_base_test_m2.png"))
sink(file = paste0(path, "fits/2/", "dharma_standardtests_basemod_m2.txt"))
DHARMa::testResiduals(s_mod_base$m2)
sink(file = NULL)
dev.off()

# Res vs preds
png(paste0(path, "fits/2/", "dharma_mods_vs_preds.png"))
par(mfrow = c(2,2))
DHARMa::plotResiduals(s_mod_base$m1)
DHARMa::plotResiduals(s_mod_base$m2)
DHARMa::plotResiduals(s_mod_cov$m1)
DHARMa::plotResiduals(s_mod_cov$m2)
title(main = list("Base model (upper), Covariate model (lower)", cex = 1.4), 
      line = -1, outer = TRUE)
dev.off()

# Residuals vs covariates
mod_plot_data <- subset(mod_depth_subst$data, biomass > 0)
png(paste0(path, "fits/2/", "dharma_bestfit_res_vs_covs.png"))
par(mfrow = c(3,2))
DHARMa::plotResiduals(s_mod_cov$m1, mod_depth_subst$data$fyear)
DHARMa::plotResiduals(s_mod_cov$m2, mod_plot_data$fyear)
DHARMa::plotResiduals(s_mod_cov$m1, mod_depth_subst$data$middepth_scaled)
DHARMa::plotResiduals(s_mod_cov$m2, mod_plot_data$middepth_scaled)
DHARMa::plotResiduals(s_mod_cov$m1, mod_depth_subst$data$substrate2)
DHARMa::plotResiduals(s_mod_cov$m2, mod_plot_data$substrate2)
title(main = list("Covariate model residuals vs covariates", cex = 1.4), 
      line = -1, outer = TRUE)
dev.off()

mod_plot_data2 <- subset(mod_spattemp$data, biomass > 0)
png(paste0(path, "fits/2/", "dharma_base_res_vs_years.png"))
par(mfrow = c(1,2))
DHARMa::plotResiduals(s_mod_base$m1, mod_spattemp$data$fyear)
DHARMa::plotResiduals(s_mod_base$m2, mod_plot_data2$fyear)
title(main = list("Base model residuals vs year", cex = 1.4), 
      line = -1, outer = TRUE)
dev.off()

# Zeroinflation
png(paste0(path, "fits/2/", "zeroinflation_mods.png"))
par(mfrow = c(2,2))
DHARMa::testZeroInflation(s_mod_base$m1)
DHARMa::testZeroInflation(s_mod_cov$m1)
DHARMa::testZeroInflation(s_mod_base$m2)
DHARMa::testZeroInflation(s_mod_cov$m2)
title(main = list("Base model (left), Covariate model (right)", cex = 1.4), 
      line = -3.3, outer = TRUE)
dev.off()
sink(file = paste0(path, "fits/2/", "dharma_zeroinf_basemod_m1.txt"))
DHARMa::testZeroInflation(s_mod_base$m1)
sink(file = NULL)
sink(file = paste0(path, "fits/2/", "dharma_zeroinf_basemod_m2.txt"))
DHARMa::testZeroInflation(s_mod_base$m2)
sink(file = NULL)
sink(file = paste0(path, "fits/2/", "dharma_zeroinf_covmod_m1.txt"))
DHARMa::testZeroInflation(s_mod_cov$m1)
sink(file = NULL)
sink(file = paste0(path, "fits/2/", "dharma_zeroinf_covmod_m2.txt"))
DHARMa::testZeroInflation(s_mod_cov$m2)
sink(file = NULL)

# Uniformity
png(paste0(path, "fits/2/", "uniformity_mods.png"))
par(mfrow = c(2,2))
DHARMa::testUniformity(s_mod_base$m1)
DHARMa::testUniformity(s_mod_cov$m1)
DHARMa::testUniformity(s_mod_base$m2)
DHARMa::testUniformity(s_mod_cov$m2)
title(main = list("Base model (left), Covariate model (right)", cex = 1.4), 
      line = -3.3, outer = TRUE)
dev.off()
sink(file = paste0(path, "fits/2/", "dharma_uniformity_basemod_m1.txt"))
DHARMa::testUniformity(s_mod_base$m1)
sink(file = NULL)
sink(file = paste0(path, "fits/2/", "dharma_uniformity_basemod_m2.txt"))
DHARMa::testUniformity(s_mod_base$m2)
sink(file = NULL)
sink(file = paste0(path, "fits/2/", "dharma_uniformity_covmod.txt"))
DHARMa::testUniformity(s_mod_cov$m1)
sink(file = NULL)
sink(file = paste0(path, "fits/2/", "dharma_uniformity_covmod.txt"))
DHARMa::testUniformity(s_mod_cov$m2)
sink(file = NULL)


# Dispersion
png(paste0(path, "fits/2/", "dispersion_mods.png"))
par(mfrow = c(2,2))
DHARMa::testDispersion(s_mod_base$m1)
DHARMa::testDispersion(s_mod_cov$m1)
DHARMa::testDispersion(s_mod_base$m2)
DHARMa::testDispersion(s_mod_cov$m2)
title(main = list("Base model (left), Covariate model (right)", cex = 1.4), 
      line = -3.3, outer = TRUE)
dev.off()
sink(file = paste0(path, "fits/2/", "dharma_dispersion_basemod_m1.txt"))
DHARMa::testDispersion(s_mod_base$m1)
sink(file = NULL)
sink(file = paste0(path, "fits/2/", "dharma_dispersion_basemod_m2.txt"))
DHARMa::testDispersion(s_mod_base$m2)
sink(file = NULL)
sink(file = paste0(path, "fits/2/", "dharma_dispersion_covmod.txt"))
DHARMa::testDispersion(s_mod_cov$m1)
sink(file = NULL)
sink(file = paste0(path, "fits/2/", "dharma_dispersion_covmod.txt"))
DHARMa::testDispersion(s_mod_cov$m2)
sink(file = NULL)

# Temporal autocorrelation
grouped_basemod1 <- recalculateResiduals(s_mod_base$m1, group = unique(mod_depth_subst$data$fyear))
grouped_basemod2 <- recalculateResiduals(s_mod_base$m2, group = unique(mod_plot_data$fyear))

grouped_covmod1<- recalculateResiduals(s_mod_cov$m1, group = unique(mod_depth_subst$data$fyear))
grouped_covmod2 <- recalculateResiduals(s_mod_cov$m2, group = unique(mod_plot_data$fyear))

png(paste0(path, "fits/2/", "tempautocorr_basemod_m1.png"))
DHARMa::testTemporalAutocorrelation(grouped_basemod1, time = unique(mod_depth_subst$data$fyear))
title(main = list("Base model", cex = 1.4), 
      line = -1, outer = TRUE)
dev.off()
png(paste0(path, "fits/2/", "tempautocorr_basemod_m2.png"))
DHARMa::testTemporalAutocorrelation(grouped_basemod2, time = unique(mod_plot_data$fyear))
title(main = list("Base model", cex = 1.4), 
      line = -1, outer = TRUE)
dev.off()
png(paste0(path, "fits/2/", "tempautocorr_covmod_m1.png"))
DHARMa::testTemporalAutocorrelation(grouped_covmod1, time = unique(mod_depth_subst$data$fyear))
title(main = list("Covariate model", cex = 1.4), 
      line = -1, outer = TRUE)
dev.off()
png(paste0(path, "fits/2/", "tempautocorr_covmod_m2.png"))
DHARMa::testTemporalAutocorrelation(grouped_covmod2, time = unique(mod_plot_data$fyear))
title(main = list("Covariate model", cex = 1.4), 
      line = -1, outer = TRUE)
dev.off()
sink(file = paste0(path, "fits/2/", "dharma_dispersion_basemod_m1.txt"))
testTemporalAutocorrelation(grouped_basemod1, time = unique(mod_depth_subst$data$fyear))
sink(file = NULL)
sink(file = paste0(path, "fits/2/", "dharma_dispersion_basemod_m2.txt"))
testTemporalAutocorrelation(grouped_basemod2, time = unique(mod_plot_data$fyear))
sink(file = NULL)
sink(file = paste0(path, "fits/2/", "dharma_dispersion_covmod_m1.txt"))
testTemporalAutocorrelation(grouped_covmod1, time = unique(mod_depth_subst$data$fyear))
sink(file = NULL)
sink(file = paste0(path, "fits/2/", "dharma_dispersion_covmod_m2.txt"))
testTemporalAutocorrelation(grouped_covmod2, time = unique(mod_plot_data$fyear))
sink(file = NULL)

## Spatial autocorrelation
mod_spattemp$data$resids <- mod_base_res1
mod_spattemp$data$resids2 <- mod_base_res2
mod_depth_subst$data$resids <- mod_cov_res1
mod_depth_subst$data$resids2 <- mod_cov_res2

# base
ggplot(data=mod_spattemp$data, aes(x=resids, group=fyear, fill=fyear)) +
  geom_density(adjust=1.5, alpha=.4) # binomial
ggplot(data=mod_spattemp$data, aes(x=resids2, group=fyear, fill=fyear)) +
  geom_density(adjust=1.5, alpha=.4) # lognormal

# cov
ggplot(data=mod_depth_subst$data, aes(x=resids, group=fyear, fill=fyear)) +
  geom_density(adjust=1.5, alpha=.4) # binomial
ggplot(data=mod_depth_subst$data, aes(x=resids2, group=fyear, fill=fyear)) +
  geom_density(adjust=1.5, alpha=.4) # lognormal

spatial_basemod1 <- ggplot(mod_spattemp$data, aes(X, Y, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year) + coord_fixed()
ggsave(filename = paste0(path, "fits/2/", "spatial_res_basemod_m1.png"), plot = spatial_basemod1)
spatial_basemod2 <- ggplot(mod_spattemp$data, aes(X, Y, col = resids2)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year) + coord_fixed()
ggsave(filename = paste0(path, "fits/2/", "spatial_res_basemod_m2.png"), plot = spatial_basemod2)
spatial_covmod1 <- ggplot(mod_depth_subst$data, aes(X, Y, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year) + coord_fixed()
ggsave(filename = paste0(path, "fits/2/", "spatial_res_covmod_m1.png"), plot = spatial_covmod1)
spatial_covmod2 <- ggplot(mod_depth_subst$data, aes(X, Y, col = resids2)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year) + coord_fixed()
ggsave(filename = paste0(path, "fits/2/", "spatial_res_covmod_m2.png"), plot = spatial_covmod2)

## Global Moran's I
# base
knw_spattemp_binom <- make_kn_dist_obj_binom(coords = mod_spattemp$data)
knw_spattemp_pos <- make_kn_dist_obj_pos(coords = mod_spattemp$data)

morans_spattemp_m1 <- purrr::map(knw_spattemp_binom, ~ list(
  moran.mc(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw, nsim = 999)
  #moran.test(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw)
))
morans_spattemp_m2 <- purrr::map(knw_spattemp_pos, ~ list(
  moran.mc(.x[[1]]$coordsdat$resids2, .x[[1]]$nblistw, nsim = 999)
  #moran.test(.x[[1]]$coordsdat$resids2, .x[[1]]$nblistw)
))

# we also need to record the sample size for  
n_spattemp_m1 <- imap_dfr(knw_spattemp_binom, ~ {
  data.frame(year = .y, N = length(.x[[1]]$nb))
})
n_spattemp_m2 <- imap_dfr(knw_spattemp_pos, ~ {
  data.frame(year = .y, N = length(.x[[1]]$nb))
})

sink(file = paste0(path, "fits/2/","global_morans_results_base_binom.txt"))
cat("Sample size in each year \n")
n_spattemp_m1
cat("\n")
cat("Moran's I in each year \n")
morans_spattemp_m1
sink(file = NULL)
sink(file = paste0(path, "fits/2/","global_morans_results_base_logn.txt"))
cat("Sample size in each year \n")
n_spattemp_m2
cat("\n")
cat("Moran's I in each year \n")
morans_spattemp_m2
sink(file = NULL)

# finally we extract plots
morans_spattemp_m1_plots <- purrr::map(knw_spattemp_binom, ~ 
                                         moran.plot(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw, return_df = TRUE)
)
morans_spattemp_m1_gg <- purrr::map(morans_spattemp_m1_plots, ~ globm_plot_fun(
  moran_data = .x,
  model_comp = 1
))
globmorans_spattemp_m1_grid <- cowplot::plot_grid(plotlist = morans_spattemp_m1_gg,
                                                  align = "hv", axis = "tb")
buffered_spattemp_m1_grid <- cowplot::ggdraw() +
  cowplot::draw_plot(globmorans_spattemp_m1_grid, 
                     x = 0.05, y = 0.05, width = 0.9, height = 0.9)
ggsave(filename = paste0(path, "fits/2/", "globm_plots_base_binom.jpg"), 
       buffered_spattemp_m1_grid)

morans_spattemp_m2_plots <- purrr::map(knw_spattemp_pos, ~ 
                                         moran.plot(.x[[1]]$coordsdat$resids2, .x[[1]]$nblistw, return_df = TRUE)
)
morans_spattemp_m2_gg <- purrr::map(morans_spattemp_m2_plots, ~ globm_plot_fun(
  moran_data = .x,
  model_comp = 2
))
globmorans_spattemp_m2_grid <- cowplot::plot_grid(plotlist = morans_spattemp_m2_gg,
                                                  align = "hv", axis = "tb")
buffered_spattemp_m2_grid <- cowplot::ggdraw() +
  cowplot::draw_plot(globmorans_spattemp_m2_grid, 
                     x = 0.05, y = 0.05, width = 0.9, height = 0.9)
ggsave(filename = paste0(path, "fits/2/", "globm_plots_base_logn.jpg"), 
       buffered_spattemp_m2_grid)


# covariate
knw_cov_binom <- make_kn_dist_obj_binom(coords = mod_depth_subst$data)
knw_cov_pos <- make_kn_dist_obj_pos(coords = mod_depth_subst$data)

morans_cov_m1 <- purrr::map(knw_cov_binom, ~ list(
  moran.mc(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw, nsim = 999)
  #moran.test(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw)
))
morans_cov_m2 <- purrr::map(knw_cov_pos, ~ list(
  moran.mc(.x[[1]]$coordsdat$resids2, .x[[1]]$nblistw, nsim = 999)
  #moran.test(.x[[1]]$coordsdat$resids2, .x[[1]]$nblistw)
))

# we also need to record the sample size for  
n_cov_m1 <- imap_dfr(knw_cov_binom, ~ {
  data.frame(year = .y, N = length(.x[[1]]$nb))
})
n_cov_m2 <- imap_dfr(knw_cov_pos, ~ {
  data.frame(year = .y, N = length(.x[[1]]$nb))
})

sink(file = paste0(path, "fits/2/","global_morans_results_cov_binom.txt"))
cat("Sample size in each year \n")
n_cov_m1
cat("\n")
cat("Moran's I in each year \n")
morans_cov_m1
sink(file = NULL)
sink(file = paste0(path, "fits/2/","global_morans_results_cov_logn.txt"))
cat("Sample size in each year \n")
n_cov_m2
cat("\n")
cat("Moran's I in each year \n")
morans_cov_m2
sink(file = NULL)

# finally we extract plots
morans_cov_m1_plots <- purrr::map(knw_cov_binom, ~ 
                                    moran.plot(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw, return_df = TRUE)
)
morans_cov_m1_gg <- purrr::map(morans_cov_m1_plots, ~ globm_plot_fun(
  moran_data = .x,
  model_comp = 1
))
globmorans_cov_m1_grid <- cowplot::plot_grid(plotlist = morans_cov_m1_gg,
                                             align = "hv", axis = "tb")
buffered_cov_m1_grid <- cowplot::ggdraw() +
  cowplot::draw_plot(globmorans_cov_m1_grid, 
                     x = 0.05, y = 0.05, width = 0.9, height = 0.9)
ggsave(filename = paste0(path, "fits/2/", "globm_plots_cov_binom.jpg"), 
       buffered_cov_m1_grid)

morans_cov_m2_plots <- purrr::map(knw_cov_pos, ~ 
                                    moran.plot(.x[[1]]$coordsdat$resids2, .x[[1]]$nblistw, return_df = TRUE)
)
morans_cov_m2_gg <- purrr::map(morans_cov_m2_plots, ~ globm_plot_fun(
  moran_data = .x,
  model_comp = 2
))
globmorans_cov_m2_grid <- cowplot::plot_grid(plotlist = morans_cov_m2_gg,
                                             align = "hv", axis = "tb")
buffered_cov_m2_grid <- cowplot::ggdraw() +
  cowplot::draw_plot(globmorans_cov_m2_grid, 
                     x = 0.05, y = 0.05, width = 0.9, height = 0.9)
ggsave(filename = paste0(path, "fits/2/", "globm_plots_cov_logn.jpg"), 
       buffered_cov_m2_grid)


## Local Moran's I
# apply functions to base and covariate models
library(parallel)
invisible(spdep::set.coresOption(max(detectCores()-1L, 1L)))
# base
locm_spattemp_binom <- purrr::map(knw_spattemp_binom, make_local_morans_obj_binom)
locm_spattemp_binom_df <- do.call(rbind, locm_spattemp_binom)
save(locm_spattemp_binom_df, 
     file = paste0(path, "fits/2/","local_morans_psignificance_base_binom.RData"))
locm_spattemp_binom_sig <- purrr::map(locm_spattemp_binom, ~ check_psig_locm(.x))
sink(file = paste0(path, "fits/2/","local_morans_psignificance_base_binom.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(locm_spattemp_binom_sig)
sink(file = NULL)

locm_spattemp_pos <- purrr::map(knw_spattemp_pos, make_local_morans_obj_pos)
locm_spattemp_pos_df <- do.call(rbind, locm_spattemp_pos)
save(locm_spattemp_pos_df, 
     file = paste0(path, "fits/2/","local_morans_psignificance_base_pos.RData"))
locm_spattemp_pos_sig <- purrr::map(locm_spattemp_pos, ~ check_psig_locm(.x))
sink(file = paste0(path, "fits/2/","local_morans_psignificance_base_pos.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(locm_spattemp_pos_sig)
sink(file = NULL)


# covariate
locm_cov_binom <- purrr::map(knw_cov_binom, make_local_morans_obj_binom)
locm_cov_binom_df <- do.call(rbind, locm_cov_binom)
save(locm_cov_binom_df, 
     file = paste0(path, "fits/2/","local_morans_psignificance_cov_binom.RData"))
locm_cov_binom_sig <- purrr::map(locm_cov_binom, ~ check_psig_locm(.x))
sink(file = paste0(path, "fits/2/","local_morans_psignificance_cov_binom.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(locm_cov_binom_sig)
sink(file = NULL)

locm_cov_pos <- purrr::map(knw_cov_pos, make_local_morans_obj_pos)
locm_cov_pos_df <- do.call(rbind, locm_cov_pos)
save(locm_cov_pos_df, 
     file = paste0(path, "fits/2/","local_morans_psignificance_cov_pos.RData"))
locm_cov_pos_sig <- purrr::map(locm_cov_pos, ~ check_psig_locm(.x))
sink(file = paste0(path, "fits/2/","local_morans_psignificance_cov_pos.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(locm_cov_pos_sig)
sink(file = NULL)

## Getis-Ord G
# base
getisord_spattemp_binom <- purrr::map(knw_spattemp_binom, make_getis_ord_obj_binom)
getisord_spattemp_binom_df <- do.call(rbind, getisord_spattemp_binom)
save(getisord_spattemp_binom_df, 
     file = paste0(path, "fits/2/","getis_ord_psignificance_spattemp_binom.RData"))
getisord_spattemp_binom_sig <- purrr::map(getisord_spattemp_binom, ~ check_psig_locg(.x))
sink(file = paste0(path, "fits/2/","getis_ord_psignificance_base_binom.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(getisord_spattemp_binom_sig)
sink(file = NULL)

getisord_spattemp_pos <- purrr::map(knw_spattemp_pos, make_getis_ord_obj_pos)
getisord_spattemp_pos_df <- do.call(rbind, getisord_spattemp_pos)
save(getisord_spattemp_pos_df, 
     file = paste0(path, "fits/2/","getis_ord_psignificance_spattemp_pos.RData"))
getisord_spattemp_pos_sig <- purrr::map(getisord_spattemp_pos, ~ check_psig_locg(.x))
sink(file = paste0(path, "fits/2/","getis_ord_psignificance_base_pos.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(getisord_spattemp_pos_sig)
sink(file = NULL)

# covariate
getisord_cov_binom <- purrr::map(knw_cov_binom, make_getis_ord_obj_binom)
getisord_cov_binom_df <- do.call(rbind, getisord_cov_binom)
save(getisord_cov_binom_df, 
     file = paste0(path, "fits/2/","getis_ord_psignificance_cov_binom.RData"))
getisord_cov_binom_sig <- purrr::map(getisord_cov_binom, ~ check_psig_locg(.x))
sink(file = paste0(path, "fits/2/","getis_ord_psignificance_cov_binom.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(getisord_cov_binom_sig)
sink(file = NULL)

getisord_cov_pos <- purrr::map(knw_cov_pos, make_getis_ord_obj_pos)
getisord_cov_pos_df <- do.call(rbind, getisord_cov_pos)
save(getisord_cov_pos_df, 
     file = paste0(path, "fits/2/","getis_ord_psignificance_cov_pos.RData"))
getisord_cov_pos_sig <- purrr::map(getisord_cov_pos, ~ check_psig_locg(.x))
sink(file = paste0(path, "fits/2/","getis_ord_psignificance_cov_pos.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(getisord_cov_pos_sig)
sink(file = NULL)

# Make the Voronoi tesselation and hull
ire <- st_read("C:/Users/astroh/OneDrive - Marine Institute/Chapter 1/Plotting canvases/ie survey canvas.shp")
ire_utm <- st_transform(ire, st_crs(32629))

# Make a final dataframe encompassing all local autocorrelation measures
# Base model
# binomial
locm_spattemp_df_binom <- locm_spattemp_binom_df
rownames(locm_spattemp_df_binom) <- NULL
locg_spattemp_df_binom <- getisord_spattemp_binom_df
rownames(locg_spattemp_df_binom) <- NULL
nrow(locm_spattemp_df_binom) == nrow(locg_spattemp_df_binom)
local_autocor_spattemp_binom <- cbind(locm_spattemp_df_binom, 
                                      locg_spattemp_df_binom[, c("Gi", "E.Gi",
                                                                 "Var.Gi", "StdDev.Gi",
                                                                 "Pr(z != E(Gi))",
                                                                 "Pr(z != E(Gi)) Sim",
                                                                 "locg_p_pv", 
                                                                 "locg_significance",
                                                                 "locg_p_clust")])

# lognormal
locm_spattemp_df_logn <- locm_spattemp_pos_df
rownames(locm_spattemp_df_logn) <- NULL
locg_spattemp_df_logn <- getisord_spattemp_pos_df
rownames(locg_spattemp_df_logn) <- NULL
nrow(locm_spattemp_df_logn) == nrow(locg_spattemp_df_logn)
local_autocor_spattemp_logn <- cbind(locm_spattemp_df_logn, 
                                     locg_spattemp_df_logn[, c("Gi", "E.Gi",
                                                               "Var.Gi", "StdDev.Gi",
                                                               "Pr(z != E(Gi))",
                                                               "Pr(z != E(Gi)) Sim",
                                                               "locg_p_pv", 
                                                               "locg_significance",
                                                               "locg_p_clust")])


# Covariate model
# binomial
locm_cov_df_binom <- locm_cov_binom_df
rownames(locm_cov_df_binom) <- NULL
locg_cov_df_binom <- getisord_cov_binom_df
rownames(locg_cov_df_binom) <- NULL
nrow(locm_cov_df_binom) == nrow(locg_cov_df_binom)
local_autocor_cov_binom <- cbind(locm_cov_df_binom, 
                                 locg_cov_df_binom[, c("Gi", "E.Gi", "Var.Gi", 
                                                       "StdDev.Gi", "Pr(z != E(Gi))",
                                                       "Pr(z != E(Gi)) Sim",
                                                       "locg_p_pv", 
                                                       "locg_significance",
                                                       "locg_p_clust")])

# lognormal
locm_cov_df_logn <- locm_cov_pos_df
rownames(locm_cov_df_logn) <- NULL
locg_cov_df_logn <- getisord_cov_pos_df
rownames(locg_cov_df_logn) <- NULL
nrow(locm_cov_df_logn) == nrow(locg_cov_df_logn)
local_autocor_cov_logn <- cbind(locm_cov_df_logn, 
                                locg_cov_df_logn[, c("Gi", "E.Gi", "Var.Gi", 
                                                     "StdDev.Gi", "Pr(z != E(Gi))",
                                                     "Pr(z != E(Gi)) Sim",
                                                     "locg_p_pv", 
                                                     "locg_significance",
                                                     "locg_p_clust")])

h_v_spattemp_binom <- make_hull_and_vortess(local_autocor_spattemp_binom, ire)
h_v_spattemp_logn <- make_hull_and_vortess(local_autocor_spattemp_logn, ire)
h_v_cov_binom <- make_hull_and_vortess(local_autocor_cov_binom, ire)
h_v_cov_logn <- make_hull_and_vortess(local_autocor_cov_logn, ire)

# Save all outputs 
spat_autocor_res <- list("base-binom" = local_autocor_spattemp_binom,
                         "base-lognormal" = local_autocor_spattemp_logn,
                         "cov-binom" = local_autocor_cov_binom,
                         "cov-lognormal" = local_autocor_cov_logn)


saveRDS(spat_autocor_res, file = paste0(path, "fits/2/", "spatialAutocorResults.rds"))


# Plot out
# plotting local Moran's
# Base
# binomial
all_locm_plots_base_binom <- purrr::map(h_v_spattemp_binom, ~ locm_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Ii = "Ii"
))
locm_plots_base_binom <- cowplot::plot_grid(plotlist = all_locm_plots_base_binom)
cowplot::save_plot(paste0(path, "fits/2/","locm_plots_base_binom.jpg"), 
                   locm_plots_base_binom, base_asp = 1.6)

# lognormal
all_locm_plots_base_logn <- purrr::map(h_v_spattemp_logn, ~ locm_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Ii = "Ii"
))
locm_plots_base_logn <- cowplot::plot_grid(plotlist = all_locm_plots_base_logn)
cowplot::save_plot(paste0(path, "fits/2/","locm_plots_base_logn.jpg"), 
                   locm_plots_base_logn, base_asp = 1.6)


# Covariate
# binomial
all_locm_plots_cov_binom <- purrr::map(h_v_cov_binom, ~ locm_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Ii = "Ii"
))
locm_plots_cov_binom <- cowplot::plot_grid(plotlist = all_locm_plots_cov_binom)
cowplot::save_plot(paste0(path, "fits/2/", "locm_plots_cov_binom.jpg"), locm_plots_cov_binom)

# lognormal
all_locm_plots_cov_logn <- purrr::map(h_v_cov_logn, ~ locm_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Ii = "Ii"
))
locm_plots_cov_logn <- cowplot::plot_grid(plotlist = all_locm_plots_cov_logn)
cowplot::save_plot(paste0(path, "fits/2/", "locm_plots_cov_logn.jpg"), locm_plots_cov_logn)


# plotting Getis Ord
# Base
# binomial
all_locg_plots_base_binom <- purrr::map(h_v_spattemp_binom, ~ locg_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Gi = "Gi"
))
gord_plots_base_binom <- cowplot::plot_grid(plotlist = all_locg_plots_base_binom)
cowplot::save_plot(paste0(path, "fits/2/", "getis_ord_plots_base_binom.jpg"), gord_plots_base_binom)

# lognormal
all_locg_plots_base_logn <- purrr::map(h_v_spattemp_logn, ~ locg_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Gi = "Gi"
))
gord_plots_base_logn <- cowplot::plot_grid(plotlist = all_locg_plots_base_logn)
cowplot::save_plot(paste0(path, "fits/2/", "getis_ord_plots_base_logn.jpg"), gord_plots_base_logn)

# Covariate
# binomial
all_locg_plots_cov_binom <- purrr::map(h_v_cov_binom, ~ locg_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Gi = "Gi"
))
gord_plots_cov_binom <- cowplot::plot_grid(plotlist = all_locg_plots_cov_binom)
cowplot::save_plot(paste0(path, "fits/2/", "getis_ord_plots_cov_binom.jpg"), gord_plots_cov_binom)

# lognormal
all_locg_plots_cov_logn <- purrr::map(h_v_cov_logn, ~ locg_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Gi = "Gi"
))
gord_plots_cov_logn <- cowplot::plot_grid(plotlist = all_locg_plots_cov_logn)
cowplot::save_plot(paste0(path, "fits/2/", "getis_ord_plots_cov_logn.jpg"), gord_plots_cov_logn)

# plotting Getis Ord hotspots
# Base
# binomial
all_locg_hotspot_plots_base_binom <- purrr::map(h_v_spattemp_binom, ~ locg_hotspot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  locg_p_clust = "locg_p_clust"
))
hotspot_plots_base_binom <- cowplot::plot_grid(plotlist = all_locg_hotspot_plots_base_binom)
cowplot::save_plot(paste0(path, "fits/2/","getis_ord_hotspot_plots_base_binom.jpg"), hotspot_plots_base_binom)

# lognormal
all_locg_hotspot_plots_base_logn <- purrr::map(h_v_spattemp_logn, ~ locg_hotspot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  locg_p_clust = "locg_p_clust"
))
hotspot_plots_base_logn <- cowplot::plot_grid(plotlist = all_locg_hotspot_plots_base_logn)
cowplot::save_plot(paste0(path, "fits/2/","getis_ord_hotspot_plots_base_logn.jpg"), 
                   hotspot_plots_base_logn)

# Covariate
# binomial
all_locg_hotspot_plots_cov_binom <- purrr::map(h_v_cov_binom, ~ locg_hotspot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  locg_p_clust = "locg_p_clust"
))
hotspot_plots_cov_binom <- cowplot::plot_grid(plotlist = all_locg_hotspot_plots_cov_binom)
cowplot::save_plot(paste0(path, "fits/2/","getis_ord_hotspot_plots_cov_binom.jpg"), 
                   hotspot_plots_cov_binom, base_asp = 1.6)

# lognormal
all_locg_hotspot_plots_cov_logn <- purrr::map(h_v_cov_logn, ~ locg_hotspot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  locg_p_clust = "locg_p_clust"
))
hotspot_plots_cov_logn <- cowplot::plot_grid(plotlist = all_locg_hotspot_plots_cov_logn)
cowplot::save_plot(paste0(path, "fits/2/","getis_ord_hotspot_plots_cov_logn.jpg"), 
                   hotspot_plots_cov_logn, base_asp = 1.6)

### Map spatial and spatiotemporal variability -----
### Map outputs ------
library(ggthemes)

# Save predictions
predictions_spattemp$data$density <- plogis(predictions_spattemp$data$est1) * exp(predictions_spattemp$data$est2)
saveRDS(predictions_spattemp, file = paste0(path, "fits/2/", "fullobj_base_model_predictions.rds"))
write.csv(predictions_spattemp$data, file = paste0(path, "fits/2/", "base_model_predictions.csv"), row.names = FALSE)

predictions_d_s$data$density <- plogis(predictions_d_s$data$est1) * exp(predictions_d_s$data$est2)
saveRDS(predictions_d_s, file = paste0(path, "fits/2/", "fullobj_cov_model_predictions.rds"))
write.csv(predictions_d_s$data, file = paste0(path, "fits/2/", "cov_model_predictions.csv"), row.names = FALSE)

# DENSITY
base_density <- plot_map(predictions_spattemp$data, density) + 
  ggtitle("Base model: Density") +
  scale_fill_viridis_c(option = "G", direction = -1) +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/", "density_base.jpg"), plot = base_density, dpi = 300)

cov_density <- plot_map(predictions_d_s$data, density) + 
  ggtitle("Covariate model: Density") +
  scale_fill_viridis_c(option = "G", direction = -1) +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/", "density_cov.jpg"), plot = cov_density, dpi = 300)

# Uncertainty around model estimates
predictions_base_est1 <- predict(mod_spattemp, 
                                 newdata = grid_yrs[, c("X", "Y", "fyear", "year")],
                                 model = 1, nsim = 500, 
                                 sims_var = "est")
grid_yrs$base_est1_se <- apply(predictions_base_est1, 1, sd) 
rm(predictions_base_est1)

predictions_base_est2 <- predict(mod_spattemp, 
                                 newdata = grid_yrs[, c("X", "Y", "fyear", "year")],
                                 model = 2, nsim = 500, 
                                 sims_var = "est")
grid_yrs$base_est2_se <- apply(predictions_base_est2, 1, sd)
rm(predictions_base_est2)

est_uncertainty_base1 <- plot_map(grid_yrs, base_est1_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, binomial component: Model estimate uncertainty", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","est_se_base1.jpg"), plot = est_uncertainty_base1)

est_uncertainty_base2 <- plot_map(grid_yrs, base_est2_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, lognormal component: Model estimate uncertainty", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","est_se_base2.jpg"), plot = est_uncertainty_base2)


# covariate model
predictions_cov_est1 <- predict(mod_depth_subst, 
                                newdata = grid_yrs,
                                model = 1, nsim = 500, 
                                sims_var = "est")
grid_yrs$cov_est1_se <- apply(predictions_cov_est1, 1, sd) # uncertainty of spatial effect
rm(predictions_cov_est1)

predictions_cov_est2 <- predict(mod_depth_subst, 
                                newdata = grid_yrs, 
                                model = 2, nsim = 500, 
                                sims_var = "est")
grid_yrs$cov_est2_se <- apply(predictions_cov_est2, 1, sd)
rm(predictions_cov_est2)

est_uncertainty_cov1 <- plot_map(grid_yrs, cov_est1_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, binomial component: Model estimate uncertainty", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","est_se_cov1.jpg"), plot = est_uncertainty_cov1)

est_uncertainty_cov2 <- plot_map(grid_yrs, cov_est2_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, lognormal component: Model estimate uncertainty", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","est_se_cov2.jpg"), plot = est_uncertainty_cov2)

# FIXED AN RANDOM EFFECTS
base_est1 <- plot_map(predictions_spattemp$data, plogis(est1)) + 
  ggtitle("Base model: Binomial model estimate") +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(fill = "Encounter\nprobability") +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","enc_prob_base.jpg"), plot = base_est1, dpi = 300)

base_est2 <- plot_map(predictions_spattemp$data, exp(est2)) + 
  ggtitle("Base model: Lognormal model estimate") +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(fill = "Biomass\n(kg)") +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","pos_catch_base.jpg"), plot = base_est2, dpi = 300)


cov_est1 <- plot_map(predictions_d_s$data, plogis(est1)) + 
  ggtitle("Covariate model: Binomial model estimate") +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(fill = "Encounter\nprobability") +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","enc_prob_cov.jpg"), plot = cov_est1)

cov_est2 <- plot_map(predictions_d_s$data, exp(est2)) + 
  ggtitle("Covariate model: Lognormal model estimate") +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(fill = "Biomass\n(kg)") +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","pos_catch_cov.jpg"), plot = cov_est2, dpi = 300)


# FIXED EFFECTS
base_ff1 <- plot_map2(predictions_spattemp$data, est_non_rf1) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, binomial component: Fixed effects only", 
       subtitle = "Fixed effects: fyear") +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","ff_effects_base1.jpg"), plot = base_ff1)

base_ff2 <- plot_map2(predictions_spattemp$data, est_non_rf2) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, lognormal component: Fixed effects only", 
       subtitle = "Fixed effects: fyear") +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","ff_effects_base2.jpg"), plot = base_ff2)

cov_ff1 <- plot_map2(predictions_d_s$data, est_non_rf1) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, binomial component: Fixed effects only", 
       subtitle = "Fixed effects: fyear + depth + substrate") +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","ff_effects_cov1.jpg"), plot = cov_ff1)

cov_ff2 <- plot_map2(predictions_d_s$data, est_non_rf2) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, lognormal component: Fixed effects only", 
       subtitle = "Fixed effects: fyear + depth + substrate") +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","ff_effects_cov2.jpg"), plot = cov_ff2)

# SPATIAL EFFECTS
base_omega1 <- plot_map2(predictions_spattemp$data, omega_s1) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, binomial component", fill = "Estimate") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        legend.position = c(.85,.65))

base_omega2 <- plot_map2(predictions_spattemp$data, omega_s2) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, lognormal component", fill = "Estimate") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        legend.position = c(.85,.65))

cov_omega1 <- plot_map2(predictions_d_s$data, omega_s1) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, binomial component", fill = "Estimate") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        legend.position = c(.85,.65)) 

cov_omega2 <- plot_map2(predictions_d_s$data, omega_s2) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, lognormal component", fill = "Estimate") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        legend.position = c(.85,.65))

omega <- (base_omega1 + base_omega2) / (cov_omega1 + cov_omega2) + 
  plot_annotation(title = "Spatial random effects only",
                  theme = theme(plot.title = element_text(hjust = 0.45, size = 16, face = "bold")))
ggsave(filename = paste0(path, "fits/2/","all_omega_base_cov.jpg"), plot = omega)

# Uncertainty around spatial random effects
# base model
predictions_base_om1 <- predict(mod_spattemp, 
                                newdata = grid_yrs[, c("X", "Y", "fyear", "year")],
                                model = 1, nsim = 500, 
                                sims_var = "omega_s")
grid_yrs$base_om1_se <- apply(predictions_base_om1, 1, sd) # uncertainty of spatial effect
rm(predictions_base_om1)

predictions_base_om2 <- predict(mod_spattemp, 
                                newdata = grid_yrs[, c("X", "Y", "fyear", "year")],
                                model = 2, nsim = 500, 
                                sims_var = "omega_s")
grid_yrs$base_om2_se <- apply(predictions_base_om2, 1, sd)
rm(predictions_base_om2)

omega_uncertainty_base1 <- plot_map2(grid_yrs, base_om1_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, binomial component", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        legend.position = c(.85,.65))

omega_uncertainty_base2 <- plot_map2(grid_yrs, base_om2_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, lognormal component", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        legend.position = c(.85,.65))

# covariate model
predictions_cov_om1 <- predict(mod_depth_subst, 
                               newdata = grid_yrs,
                               model = 1, nsim = 500, 
                               sims_var = "omega_s")
grid_yrs$cov_om1_se <- apply(predictions_cov_om1, 1, sd) # uncertainty of spatial effect
rm(predictions_cov_om1)

predictions_cov_om2 <- predict(mod_depth_subst, 
                               newdata = grid_yrs, 
                               model = 2, nsim = 500, 
                               sims_var = "omega_s")
grid_yrs$cov_om2_se <- apply(predictions_cov_om2, 1, sd)
rm(predictions_cov_om2)

omega_uncertainty_cov1 <- plot_map2(grid_yrs, cov_om1_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, binomial component", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        legend.position = c(.85,.65))

omega_uncertainty_cov2 <- plot_map2(grid_yrs, cov_om2_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, lognormal component", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
        legend.position = c(.85,.65))

omega_se <- (omega_uncertainty_base1 + omega_uncertainty_base2) / 
  (omega_uncertainty_cov1 + omega_uncertainty_cov2) + 
  plot_annotation(title = "Point-wise uncertainty in spatial predictions",
                  theme = theme(plot.title = element_text(hjust = 0.45, size = 16, face = "bold")))
ggsave(filename = paste0(path, "fits/2/","all_omega_se_base_cov.jpg"), plot = omega_se)

# SPATIOTEMPORAL EFFECTS
base_preds <- read.csv(file = paste0(path, "fits/2/","base_model_predictions.csv"))
base_eps1 <- plot_map(predictions_spattemp$data, epsilon_st1) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, binomial component: Spatiotemporal random effects only", 
       fill = "Estimate") +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","eps_effects_base1.jpg"), plot = base_eps1, dpi = 300)

base_eps2 <- plot_map(predictions_spattemp$data, epsilon_st2) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, lognormal component: Spatiotemporal random effects only", 
       fill = "Estimate") +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","eps_effects_base2.jpg"), plot = base_eps2, dpi = 300)

cov_preds <- read.csv(file = paste0(path, "fits/2/","cov_model_predictions.csv"))
cov_eps1 <- plot_map(predictions_d_s$data, epsilon_st1) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, binomial component: Spatiotemporal random effects only", 
       fill = "Estimate") +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","eps_effects_cov1.jpg"), plot = cov_eps1, dpi = 300)

cov_eps2 <- plot_map(predictions_d_s$data, epsilon_st2) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, lognormal component: Spatiotemporal random effects only", 
       fill = "Estimate") +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","eps_effects_cov2.jpg"), plot = cov_eps2, dpi = 300)

# Uncertainty of spatiotemporal effects
# base model
predictions_base_eps1 <- predict(mod_spattemp, 
                                 newdata = grid_yrs[, c("X", "Y", "fyear", "year")], 
                                 model = 1, nsim = 500, 
                                 sims_var = "epsilon_st")
grid_yrs$base_eps1_se <- apply(predictions_base_eps1, 1, sd) # uncertainty of spatial effect
rm(predictions_base_eps1)

predictions_base_eps2 <- predict(mod_spattemp, 
                                 newdata = grid_yrs[, c("X", "Y", "fyear", "year")], 
                                 model = 1, nsim = 500, 
                                 sims_var = "epsilon_st")
grid_yrs$base_eps2_se <- apply(predictions_base_eps2, 1, sd)
rm(predictions_base_eps2)

eps_uncertainty_base1 <- plot_map(grid_yrs, base_eps1_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, binom comp: Point-wise uncertainty in spatiotemporal predictions", 
       fill = "SE") +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","eps_se_base1.jpg"), plot = eps_uncertainty_base1, dpi = 300)

eps_uncertainty_base2 <- plot_map(grid_yrs, base_eps2_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model, logn comp: Point-wise uncertainty in spatiotemporal predictions", 
       fill = "SE") +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","eps_se_base2.jpg"), plot = eps_uncertainty_base2, dpi = 300)


# covariate model
predictions_cov_eps1 <- predict(mod_depth_subst, 
                                newdata = grid_yrs,
                                model = 1, nsim = 500, 
                                sims_var = "epsilon_st")
grid_yrs$cov_eps1_se <- apply(predictions_cov_eps1, 1, sd) # uncertainty of spatial effect
rm(predictions_cov_eps1)

predictions_cov_eps2 <- predict(mod_depth_subst, 
                                newdata = grid_yrs,
                                model = 2, nsim = 500, 
                                sims_var = "epsilon_st")
grid_yrs$cov_eps2_se <- apply(predictions_cov_eps2, 1, sd)
rm(predictions_cov_eps2)

eps_uncertainty_cov1 <- plot_map(grid_yrs, cov_eps1_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, binom comp: Point-wise uncertainty in spatiotemporal predictions", 
       fill = "SE") +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","eps_se_cov1.jpg"), plot = eps_uncertainty_cov1)

eps_uncertainty_cov2 <- plot_map(grid_yrs, cov_eps2_se) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model, logn comp: Point-wise uncertainty in spatiotemporal predictions", 
       fill = "SE") +
  theme_map() + 
  theme(legend.position = c(.85,.05))
ggsave(filename = paste0(path, "fits/2/","eps_se_cov2.jpg"), plot = eps_uncertainty_cov2)

# Save all predictions
#save(grid_yrs, file = paste0(path, "fits/2/", "uncertainty_mapping_data.RData"))
#load(paste0(path, "fits/2/", "uncertainty_mapping_data.RData"))


