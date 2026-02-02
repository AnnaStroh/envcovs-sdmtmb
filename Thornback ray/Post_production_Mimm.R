#############
## Mimm fits
## diagnostics, plots, index
#############

library(sdmTMB)
library(sdmTMBextra)
library(ggeffects)
library(visreg)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(DHARMa)
library(patchwork)
library(cowplot)
library(purrr)

path <- "C:/Users/astroh/Desktop/Chapter 2/sdmTMB/RAJ/"
plot_path <- "C:/Users/astroh/Desktop/Chapter 2/sdmTMB/RAJ/all plots"

source("funs_raj.R")

### Check AIC ------
mods_mimm <- readRDS(paste0(path, "fits/mimm/mimm_model_object4.rds"))
names(mods_mimm)

all_aic <- map_dfr(names(mods_mimm), function(name) {
  model <- mods_mimm[[name]]
  tidyr::tibble(
    Model = name,
    AIC = model$aic
    #cAIC = model$caic
  )
})
print(all_aic, n=30)

sink(file = paste0(path, "/fits/Mimm/","aic.txt"))
all_aic
sink(file = NULL)

### Read data ------
sanity(mods_mimm$m_d_2nb3[[1]])
sanity(mods_mimm$m_d_2nb9a[[1]])

mod_spattemp <- mods_mimm$m_d_2nb3[[1]]
mod_dss <- mods_mimm$m_d_2nb9a[[1]] # best

sink(file = paste0(path, "/fits/Mimm/","base_model.txt"))
mod_spattemp

print("SDREPORT")
mod_spattemp$sd_report
sink(file = NULL)

sink(file = paste0(path, "/fits/Mimm/","cov_model.txt"))
mod_dss

print("SDREPORT")
mod_dss$sd_report
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
                       substrate = pred_grid$substrate,
                       bottomT_scaled = pred_grid$bottomT_scaled,
                       fyear = as.factor(pred_grid$fyear),
                       fgear = as.factor(pred_grid$gear),
                       year = pred_grid$year)
grid_yrs <- distinct(grid_yrs)

# With spatial random effects and spatiotemporal rf in both predictors
predictions_spattemp <- predict(mod_spattemp, newdata = grid_yrs[, c("X", "Y", "fyear", "fgear", "year")], 
                                offset = rep(0, nrow(grid_yrs)), 
                                return_tmb_object = TRUE)

# Including all covariates
predictions_d_s_substr <- predict(mod_dss, newdata = grid_yrs, 
                                  offset = rep(0, nrow(grid_yrs)),
                                  return_tmb_object = TRUE)

preds <- list(predictions_spattemp, #predictions_depth, 
              predictions_d_s_substr)
modelIndices <- purrr::map(preds, ~ list(
  get_index(.x, bias_correct = TRUE)
))

names(modelIndices) <- c("mod_spattemp", #"mod_depth", 
                         "mod_dss")
modelIndices$mod_spattemp[[1]]$model <- "base"
#modelIndices$mod_depth[[1]]$model <- "base+depth"
modelIndices$mod_dss[[1]]$model <- "base+substrate+depth*gear+bottomT"

all_indices <- as.data.frame(do.call(rbind, lapply(modelIndices, as.data.frame)))
rownames(all_indices) <- NULL


### Visualise indices
index_mimm <- ggplot(all_indices, aes(year, est, colour = model), size=2) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model), alpha = 0.2) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  labs(fill = "Model", colour = "Model") +
  xlab('Year') + ylab('Count') +
  theme_bw() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.position = c(.25,.95),
        legend.background = element_rect(linewidth = 1, colour = "black"), 
        plot.margin = margin(5, 20, 5, 5)) +
  ggtitle("Abundance index by model approach: Male immature Thornback ray")
ggsave(filename = paste0(path, "/fits/Mimm/","nbinom_indices.jpg"), plot = index_mimm)
saveRDS(index_mimm, file = paste0(plot_path, "/Mimm_index_ribbon.rds"))

index_mimm2 <- ggplot(all_indices, aes(year, est, colour = model, linetype = model), size=2) + 
  geom_line() +
  scale_color_manual(values = c("grey", "black")) +
  labs(fill = "Model", colour = "Model", linetype = "Model") +
  xlab('Year') + ylab('Relative Index of Abundance') +
  theme_classic() +
  theme(axis.text = element_text(size = 13), 
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = c(.28,.9),
        legend.background = element_rect(linewidth = 1, colour = "black"), 
        plot.margin = margin(5, 20, 5, 5)) +
  ggtitle("Abundance index by model approach: Male immature Thornback ray")
ggsave(filename = paste0(path, "/fits/Mimm/","nbinom_indices2.jpg"), plot = index_mimm2)
saveRDS(index_mimm2, file = paste0(plot_path, "/Mimm_index.rds"))

all_indices$CI_width <- with(all_indices, upr - lwr)
head(all_indices)

index_mimm_error <- ggplot(all_indices, aes(year, CI_width, colour = model, linetype = model), size=2) + 
  geom_line() +
  scale_color_manual(values = c("grey", "black")) +
  labs(fill = "Model", colour = "Model", linetype = "Model") +
  xlab('Year') + ylab('log SE of index') +
  theme_bw() +
  theme(axis.text = element_text(size = 13), 
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = c(.7,.9),
        legend.background = element_rect(linewidth = 1, colour = "black"), 
        plot.margin = margin(5, 20, 5, 5)) +
  ggtitle("Index error by model approach: Male immature Thornback ray")
ggsave(filename = paste0(path, "/fits/Mimm/","nbinom_indices_error.jpg"), plot = index_mimm_error)
saveRDS(index_mimm_error, file = paste0(plot_path, "Mimm_index_CI.rds"))

index_mimm_se <- ggplot(all_indices, aes(year, log(se), colour = model, linetype = model), size=2) + 
  geom_line() +
  #geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
  #              width=.2, position=position_dodge(.9), 
  #              colour = "grey40", linewidth = 2) +
  scale_fill_manual(values = c("grey", "black")) +
  scale_colour_manual(values = c("grey", "black")) +
  labs(fill = "Model", colour = "Model", linetype = "Model") +
  xlab('Year') + ylab('log SE of index') +
  theme_bw() +
  theme(axis.text = element_text(size = 13), 
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = c(.7,.9),
        legend.background = element_rect(linewidth = 1, colour = "black"), 
        plot.margin = margin(5, 20, 5, 5)) +
  ggtitle("Index error by model approach: Male immature Thornback ray")
ggsave(filename = paste0(path, "/fits/Mimm/","nbinom_indices_se.jpg"), plot = index_mimm_se)
saveRDS(index_mimm_se , file = paste0(plot_path, "/Mimm_index_se.rds"))

save(all_indices, file = paste0(path, "/fits/Mimm/", "allindices.RData"))

### Make plot of error

mods_mimm$m_d_2nb3$rf$model <- "base"
t <- mods_mimm$m_d_2nb3$rf
mods_mimm$m_d_2nb9$rf$model <- "base+substrate+depth*+bottomT"
t2 <- mods_mimm$m_d_2nb9$rf
all_rfs <- rbind(t, t2)

all_rfs_sub <- subset(all_rfs, term %in% c("sigma_O", "sigma_E"))
all_rfs_sub$term[all_rfs_sub$term == "sigma_O"] <- "spatial"
all_rfs_sub$term[all_rfs_sub$term == "sigma_E"] <- "spatiotemporal"

rfs_plot2D <- ggplot(all_rfs_sub, aes(term, estimate, fill = model)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width=.2, position=position_dodge(.9), 
                colour = "grey40", linewidth = 2) +
  scale_fill_manual(values = c("grey", "black")) +
  #labs(title = "Spatial and spatiotemporal random effects by model: whiting age 1", 
  #     subtitle = "Each model contains a fixed effect of year", 
  #     x = "Random effect", y = "Estimate", fill = "Model") +
  labs(x = "Random effect", y = "Estimate", fill = "Model", 
       title = "Random effects: Male immature Thornback ray") +
  theme_bw() +
  theme(axis.text = element_text(size = 13), 
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = c(.7,.9),
        legend.background = element_rect(linewidth = 1, colour = "black"), 
        plot.margin = margin(5, 20, 5, 5))
ggsave(filename = paste0(path, "/fits/Mimm/","rf_2D.jpg"), plot = rfs_plot2D)
saveRDS(rfs_plot2D, file = paste0(plot_path, "/Mimm_rf_error.rds"))

### MARGINAL model effects ----
### Base model
# year
y_cov <- visreg(mod_spattemp, xvar = "fyear", 
                data = mod_spattemp$data, 
                trans = exp, gg = TRUE)
y_cov2 <- visreg(mod_spattemp, xvar = "fyear", 
                 data = mod_spattemp$data, 
                 trans = exp, gg = TRUE, partial = FALSE, rug = FALSE)

g_cov <- visreg(mod_spattemp, xvar = "fgear", 
                data = mod_spattemp$data, gg = TRUE)
g_cov2 <- visreg(mod_spattemp, xvar = "fgear", 
                 data = mod_spattemp$data, 
                 trans = exp, gg = TRUE, 
                 partial = FALSE, rug = FALSE)

mar_eff_list_base <- list("fyear" = y_cov, "fyear2" = y_cov2,
                          "fgear" = g_cov, "fgear2" = g_cov2)
saveRDS(mar_eff_list_base, paste0(path, "/fits/Mimm/", "allmarginaleffects_base.rds"))

mar_effs_base <- y_cov / g_cov 
ggsave(filename = paste0(path, "/fits/Mimm/","mareffplot_base.jpg"), plot = mar_effs_base)
mar_effs2_base <- y_cov2 / g_cov2
ggsave(filename = paste0(path, "/fits/Mimm/","mareffplot2_base.jpg"), plot = mar_effs2_base)

# Covariate model
y_cov <- visreg(mod_dss, xvar = "fyear", 
                data = mod_dss$data, gg = TRUE)
y_cov2 <- visreg(mod_dss, xvar = "fyear", 
                 data = mod_dss$data, 
                 trans = exp, gg = TRUE, partial = FALSE, rug = FALSE)

g_cov <- visreg(mod_dss, xvar = "fgear", 
                data = mod_dss$data, gg = TRUE)
g_cov2 <- visreg(mod_dss, xvar = "fgear", 
                 data = mod_dss$data, 
                 trans = exp, gg = TRUE, 
                 partial = FALSE, rug = FALSE)

d_cov <- visreg(mod_dss, xvar = "middepth_scaled",
                by = "fgear", data = mod_dss$data, gg = TRUE)
d_cov2 <- visreg(mod_dss, xvar = "middepth_scaled", 
                 data = mod_dss$data, by = "fgear", 
                 trans = exp, gg = TRUE, 
                 partial = FALSE, rug = FALSE)

substr_cov <- visreg(mod_dss, xvar = "substrate", 
                     data = mod_dss$data, gg = TRUE)
substr_cov2 <- visreg(mod_dss, xvar = "substrate", 
                      data = mod_dss$data, 
                      trans = exp, gg = TRUE, 
                      partial = FALSE, rug = FALSE)

sbt_cov <- visreg(mod_dss, "bottomT_scaled", 
                  data = mod_dss$data, gg = TRUE)
sbt_cov2 <- visreg(mod_dss, "bottomT_scaled", 
                   data = mod_dss$data, 
                   trans = exp, gg = TRUE, 
                   partial = FALSE, rug = FALSE)

mar_eff_list <- list("fyear" = y_cov, "fyear2" = y_cov2,
                     "fgear" = g_cov, "fgear2" = g_cov2,
                     "depth" = d_cov, "depth2" = d_cov2,
                     "substrate" = substr_cov, "substrate2" = substr_cov2,
                     "bottomT" = sbt_cov, "bottomT2" = sbt_cov2)
saveRDS(mar_eff_list, paste0(path, "/fits/Mimm/", "allmarginaleffects_cov.rds"))


mar_effs <- (y_cov + g_cov) / d_cov / (substr_cov + sbt_cov)
ggsave(filename = paste0(path, "/fits/Mimm/","mareffplot_cov.jpg"), plot = mar_effs)
mar_effs2 <- (y_cov2 + g_cov2) / d_cov2 / (substr_cov2 + sbt_cov2)
ggsave(filename = paste0(path, "/fits/Mimm/","mareffplot2_cov.jpg"), plot = mar_effs2)

### Partial deviance explained (PDE) -----
mod_depth <- mods_mimm$m_d_2nb4a[[1]]
mod_substr <- mods_mimm$m_d_2nb6[[1]]
mod_sbt <- mods_mimm$m_d_2nb5[[1]]

null_model <- sdmTMB(
  data = mod_spattemp$data,
  formula = count ~ 1,
  mesh = mod_spattemp$spde,
  family = nbinom2(),
  offset = mod_spattemp$data$log_sweptareakmsqadj,
)

sink(file = paste0(path, "/fits/Mimm/", "pde_report.txt"))
deviance_report(full_cov_model = mod_dss, 
                null_model = null_model, 
                base_model = mod_spattemp, 
                reduced_models = list(mod_depth, mod_substr, mod_sbt))
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

terms_base <- rbind(tidy(mod_spattemp))
print(terms_base, n = nrow(terms_base))
SDrepFixedFac_base <- SDrepFixed_base[SDrepFixed_base$Variable == "b_j",]

nrow(SDrepFixedFac_base) == nrow(terms_base)
all_terms_base <- cbind(terms_base[, "term"], SDrepFixedFac_base, 
                        terms_base[, c("conf.low", "conf.high")])
all_terms_base$sig001 <- ifelse(all_terms_base$`Pr(>|z^2|)` < 0.001,
                                "sig", "insig")
all_terms_base$sig05 <- ifelse(all_terms_base$`Pr(>|z^2|)` < 0.05,
                               "sig", "insig")
write.csv(all_terms_base, 
          file = paste0(path, "/fits/Mimm/", "base_factor_level_pvalues.csv"))

# cov - year + covariates
# factor effects
SDrepFixed_cov <- as.data.frame(summary(mod_dss$sd_report, 
                                        p.value=TRUE, select = "fixed"))
SDrepFixed_cov$Variable <- row.names(SDrepFixed_cov)
SDrepFixed_cov$Variable <- gsub("\\.[0-9]+", "", SDrepFixed_cov$Variable)
row.names(SDrepFixed_cov) <- NULL
table(SDrepFixed_cov$Variable)
SDrepFixedFac_cov <- SDrepFixed_cov[SDrepFixed_cov$Variable == "b_j",]

terms_cov <- rbind(tidy(mod_dss))
print(terms_cov, n = nrow(terms_cov))

nrow(SDrepFixedFac_cov) == nrow(terms_cov)
all_terms_cov <- cbind(terms_cov[, "term"], SDrepFixedFac_cov, 
                       terms_cov[, c("conf.low", "conf.high")])
all_terms_cov$sig001 <- ifelse(all_terms_cov$`Pr(>|z^2|)` < 0.001,
                               "sig", "insig")
all_terms_cov$sig05 <- ifelse(all_terms_cov$`Pr(>|z^2|)` < 0.05,
                              "sig", "insig")
write.csv(all_terms_cov, 
          file = paste0(path, "/fits/Mimm/", "cov_factor_level_pvalues.csv"))

########## -------------------------
### TESTING MODEL DESIGN AND FIT
########## -------------------------

### Crossvalidation ------

# separate data in test and training folds
set.seed(343) # required for random sampling
base_datasets <- holdout(model = mod_spattemp, training_cutoff = 0.7)
cov_datasets <- holdout(model = mod_dss, training_cutoff = 0.7)

# run cross validation on training datasets
mesh_base <- make_mesh(base_datasets$training, c("X", "Y"), # original mesh settings
                       fmesher_func = fmesher::fm_mesh_2d_inla,
                       cutoff = 10, max.edge = c(75, 100), offset = 10) 
base_cv <- sdmTMB(
  data = base_datasets$training,
  formula = count ~ fgear + fyear,
  mesh = mesh_base,
  family = nbinom2(),
  offset = base_datasets$training$log_sweptareakmsqadj,
  spatial = "on",
  time = "year", 
  spatiotemporal = "IID",
  anisotropy = TRUE,
  share_range = TRUE,
  silent = FALSE
)

mesh_cov <- make_mesh(cov_datasets$training, c("X", "Y"), # original mesh settings
                      fmesher_func = fmesher::fm_mesh_2d_inla,
                      cutoff = 10, max.edge = c(75, 100), offset = 10)
cov_cv <- sdmTMB( 
  data = cov_datasets$training,
  formula = count ~ fgear + fyear + s(bottomT_scaled, m = 1) + 
    s(middepth_scaled, m = 1, by = fgear) + substrate,
  mesh = mesh_cov,
  family = nbinom2(),
  offset = cov_datasets$training$log_sweptareakmsqadj,
  spatial = "on",
  time = "year", 
  spatiotemporal = "IID",
  anisotropy = TRUE,
  share_range = TRUE,
  silent = FALSE
)

# predict model estimates onto test data
base_cv_preds <- predict(base_cv, 
                         newdata = base_datasets$test, type = "response")
cov_cv_preds <- predict(cov_cv, 
                        newdata = cov_datasets$test, type = "response")


# save crossvalidation products
crossvalfits <- list("basemodel_cv" = base_cv, "covmodel_cv" = cov_cv)
saveRDS(file = paste0(path, "/fits/Mimm/", "cv_fits.rds"), crossvalfits)

sink(file = paste0(path, "/fits/Mimm/", "crossvalidation_results.txt"))
crossval_report(base_cv_preds, cov_cv_preds)
sink(file = NULL)

### Check residuals ------
#Looking at tails
# randomised quantile residuals - normal
mod_base_res <- residuals(mod_spattemp, type = "mle-mvn")
mod_cov_res <- residuals(mod_dss, type = "mle-mvn")

png(paste0(path, "/fits/Mimm/", "randomised_quantile_base_bestfit.png"))
par(mfrow = c(1,2))
qqnorm(mod_base_res, main = "Base model: Normal Q-Q Plot");abline(0, 1)
qqnorm(mod_cov_res, main = "Covariate model: Normal Q-Q Plot");abline(0, 1)
#title(main = list("Best fit: Randomised quantile residuals", cex = 1.4), 
#      line = -1, outer = TRUE)
dev.off()

# simulation-based - uniform
s_mod_base <- simulate(mod_spattemp, nsim = 500, type = "mle-mvn") |>
  dharma_residuals(mod_spattemp, return_DHARMa = TRUE)
s_mod_cov <- simulate(mod_dss, nsim = 500, type = "mle-mvn") |>
  dharma_residuals(mod_dss, return_DHARMa = TRUE)

# Uniformity, dispersion, outliers
png(paste0(path, "/fits/Mimm/", "dharma_covbestfit_test.png"))
sink(file = paste0(path, "/fits/Mimm/", "dharma_standardtests_covmod.txt"))
DHARMa::testResiduals(s_mod_cov)
sink(file = NULL)
dev.off()
png(paste0(path, "/fits/Mimm/", "dharma_base_test.png"))
sink(file = paste0(path, "/fits/Mimm/", "dharma_standardtests_basemod.txt"))
DHARMa::testResiduals(s_mod_base)
sink(file = NULL)
dev.off()

# Res vs preds
png(paste0(path, "/fits/Mimm/", "dharma_mods_vs_preds.png"))
par(mfrow = c(1,2))
DHARMa::plotResiduals(s_mod_base)
DHARMa::plotResiduals(s_mod_cov)
title(main = list("Base model (left), Covariate model (right)", cex = 1.4), 
      line = -1, outer = TRUE)
dev.off()

# Residuals vs covariates
png(paste0(path, "/fits/Mimm/", "dharma_bestfit_res_vs_covs.png"))
par(mfrow = c(3,1))
DHARMa::plotResiduals(s_mod_cov, mod_dss$data$middepth_scaled)
DHARMa::plotResiduals(s_mod_cov, mod_dss$data$bottomT_scaled)
DHARMa::plotResiduals(s_mod_cov, mod_dss$data$substrate)
dev.off()

# Zeroinflation
png(paste0(path, "/fits/Mimm/", "zeroinflation_mods.png"))
par(mfrow = c(1,2))
DHARMa::testZeroInflation(s_mod_base)
DHARMa::testZeroInflation(s_mod_cov)
title(main = list("Base model (left), Covariate model (right)", cex = 1.4), 
      line = -3.3, outer = TRUE)
dev.off()
sink(file = paste0(path, "/fits/Mimm/", "dharma_zeroinf_basemod.txt"))
DHARMa::testZeroInflation(s_mod_base)
sink(file = NULL)
sink(file = paste0(path, "/fits/Mimm/", "dharma_zeroinf_covmod.txt"))
DHARMa::testZeroInflation(s_mod_cov)
sink(file = NULL)

# Dispersion
png(paste0(path, "/fits/Mimm/", "dispersion_mods.png"))
par(mfrow = c(1,2))
DHARMa::testDispersion(s_mod_base)
DHARMa::testDispersion(s_mod_cov)
title(main = list("Base model (left), Covariate model (right)", cex = 1.4), 
      line = -3.3, outer = TRUE)
dev.off()
sink(file = paste0(path, "/fits/Mimm/", "dharma_dispersion_basemod.txt"))
DHARMa::testDispersion(s_mod_base)
sink(file = NULL)
sink(file = paste0(path, "/fits/Mimm/", "dharma_dispersion_covmod.txt"))
DHARMa::testDispersion(s_mod_cov)
sink(file = NULL)

# Temporal autocorrelation
grouped_basemod <- recalculateResiduals(s_mod_base, group = unique(mod_spattemp$data$fyear))
grouped_covmod <- recalculateResiduals(s_mod_cov, group = unique(mod_dss$data$fyear))

png(paste0(path, "/fits/Mimm/", "tempautocorr_basemod.png"))
DHARMa::testTemporalAutocorrelation(grouped_basemod, time = unique(mod_spattemp$data$fyear))
title(main = list("Base model", cex = 1.4), 
      line = -1, outer = TRUE)
dev.off()
png(paste0(path, "/fits/Mimm/", "tempautocorr_covmod.png"))
DHARMa::testTemporalAutocorrelation(grouped_covmod, time = unique(mod_dss$data$fyear))
title(main = list("Covariate model", cex = 1.4), 
      line = -1, outer = TRUE)
dev.off()
sink(file = paste0(path, "/fits/Mimm/", "dharma_dispersion_basemod.txt"))
testTemporalAutocorrelation(grouped_basemod, time = unique(mod_spattemp$data$fyear))
sink(file = NULL)
sink(file = paste0(path, "/fits/Mimm/", "dharma_dispersion_covmod.txt"))
testTemporalAutocorrelation(grouped_covmod, time = unique(mod_dss$data$fyear))
sink(file = NULL)

## Spatial autocorrelation
set.seed(123)
mod_spattemp$data$resids <- residuals(mod_spattemp, type = "mle-mvn")
mod_spattemp$data$cov_resids <- residuals(mod_dss, type = "mle-mvn")
mod_spattemp$data$dif_resids <- with(mod_spattemp$data,  resids - cov_resids)
mod_dss$data$resids <- residuals(mod_dss, type = "mle-mvn")

spatial_basemod <- ggplot(mod_spattemp$data, aes(X, Y, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year) + coord_fixed()
ggsave(filename = paste0(path, "/fits/Mimm/","spatial_res_basemod.png"), plot = spatial_basemod)
spatial_covmod <- ggplot(mod_dss$data, aes(X, Y, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~year) + coord_fixed()
ggsave(filename = paste0(path, "/fits/Mimm/","spatial_res_covmod.png"), plot = spatial_covmod)

## Global Moran's I
knw_spattemp <- make_kn_dist_obj(coords = mod_spattemp$data)
morans_spattemp <- purrr::map(knw_spattemp, ~ list(
  moran.mc(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw, nsim = 999)
))

knw_cov <- make_kn_dist_obj(coords = mod_dss$data)
morans_cov <- purrr::map(knw_cov, ~ list(
  moran.mc(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw, nsim = 999)
))

# we also need to record the sample size for  
n_base <- imap_dfr(knw_spattemp, ~ {
  data.frame(year = .y, N = length(.x[[1]]$nb))
})
n_cov <- imap_dfr(knw_cov, ~ {
  data.frame(year = .y, N = length(.x[[1]]$nb))
})

sink(file = paste0(path, "/fits/Mimm/","global_morans_results_base.txt"))
cat("Sample size in each year \n")
n_base
cat("\n")
cat("Moran's I in each year \n")
morans_spattemp
sink(file = NULL)
sink(file = paste0(path, "/fits/Mimm/","global_morans_results_cov.txt"))
cat("Sample size in each year \n")
n_cov
cat("\n")
cat("Moran's I in each year \n")
morans_cov
sink(file = NULL)

# finally we extract plots
morans_base_plots <- purrr::map(knw_spattemp, ~ 
                                  moran.plot(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw, return_df = TRUE)
)
morans_base_gg <- purrr::map(morans_base_plots, ~ globm_plot_fun(
  moran_data = .x,
  model_comp = 1
))
globmorans_base_grid <- cowplot::plot_grid(plotlist = morans_base_gg,
                                           align = "hv", axis = "tb")
buffered_base_grid <- cowplot::ggdraw() +
  cowplot::draw_plot(globmorans_base_grid, 
                     x = 0.05, y = 0.05, width = 0.9, height = 0.9)
ggsave(filename = paste0(path, "/fits/Mimm/", "globm_plots_base.jpg"), 
       buffered_base_grid)

morans_cov_plots <- purrr::map(knw_cov, ~ 
                                 moran.plot(.x[[1]]$coordsdat$resids, .x[[1]]$nblistw, return_df = TRUE)
)
morans_cov_gg <- purrr::map(morans_cov_plots, ~ globm_plot_fun(
  moran_data = .x,
  model_comp = 1
))
globmorans_cov_grid <- cowplot::plot_grid(plotlist = morans_cov_gg,
                                          align = "hv", axis = "tb")
buffered_cov_grid <- cowplot::ggdraw() +
  cowplot::draw_plot(globmorans_cov_grid, 
                     x = 0.05, y = 0.05, width = 0.9, height = 0.9)
ggsave(filename = paste0(path, "/fits/Mimm/", "globm_plots_cov.jpg"), 
       buffered_cov_grid)



# apply functions to base and covariate models
library(parallel)
invisible(spdep::set.coresOption(max(detectCores()-1L, 1L)))

locm_spattemp <- purrr::map(knw_spattemp, make_local_morans_obj)
locm_spattemp_sig <- purrr::map(locm_spattemp, ~ check_psig_locm(.x))
sink(file = paste0(path, "/fits/Mimm/","local_morans_psignificance_base.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(locm_spattemp_sig)
sink(file = NULL)

locm_cov <- purrr::map(knw_cov, make_local_morans_obj)
locm_cov_sig <- purrr::map(locm_cov, ~ check_psig_locm(.x))
sink(file = paste0(path, "/fits/Mimm/","local_morans_psignificance_cov.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(locm_cov_sig)
sink(file = NULL)

## Getis-Ord G
getisord_spattemp <- purrr::map(knw_spattemp, make_getis_ord_obj)
getisord_spattemp_sig <- purrr::map(getisord_spattemp, ~ check_psig_locg(.x))
sink(file = paste0(path, "/fits/Mimm/","getis_ord_psignificance_base.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(getisord_spattemp_sig)
sink(file = NULL)

getisord_cov <- purrr::map(knw_cov, make_getis_ord_obj)
getisord_cov_sig <- purrr::map(getisord_cov, ~ check_psig_locg(.x))
sink(file = paste0(path, "/fits/Mimm/","getis_ord_psignificance_cov.txt"))
print("Sum of coordinates with p-value significance below alpha level threshold 0.01")
print(getisord_cov_sig)
sink(file = NULL)

# Make the Voronoi tesselation and hull
ire <- st_read("C:/Users/astroh/OneDrive - Marine Institute/Chapter 1/Plotting canvases/ie survey canvas.shp")
ire_utm <- st_transform(ire, st_crs(32629))

# Make a final dataframe encompassing all local autocorrelation measures
# Base model
locm_spattemp_df <- do.call(rbind, locm_spattemp)
rownames(locm_spattemp_df) <- NULL
locg_spattemp_df <- do.call(rbind, getisord_spattemp)
rownames(locg_spattemp_df) <- NULL
nrow(locm_spattemp_df) == nrow(locg_spattemp_df)
local_autocor_spattemp <- cbind(locm_spattemp_df, locg_spattemp_df[, c("Gi", "E.Gi",
                                                                       "Var.Gi", "StdDev.Gi",
                                                                       "Pr(z != E(Gi))",
                                                                       "Pr(z != E(Gi)) Sim",
                                                                       "locg_p_pv", 
                                                                       "locg_significance",
                                                                       "locg_p_clust")])
# Covariate model
locm_cov_df <- do.call(rbind, locm_cov)
rownames(locm_cov_df) <- NULL
locg_cov_df <- do.call(rbind, getisord_cov)
rownames(locg_cov_df) <- NULL
nrow(locg_cov_df) == nrow(locg_cov_df)
local_autocor_cov <- cbind(locm_cov_df, locg_cov_df[, c("Gi", "E.Gi", "Var.Gi", 
                                                        "StdDev.Gi", "Pr(z != E(Gi))",
                                                        "Pr(z != E(Gi)) Sim",
                                                        "locg_p_pv", 
                                                        "locg_significance",
                                                        "locg_p_clust")])

h_v_spattemp <- make_hull_and_vortess(local_autocor_spattemp, ire)
h_v_cov <- make_hull_and_vortess(local_autocor_cov, ire)

# Save all outputs 
spat_autocor_res <- list("base" = local_autocor_spattemp,
                         "cov" = local_autocor_cov)
saveRDS(spat_autocor_res, file = paste0(path, "/fits/Mimm/", "spatialAutocorResults.rds"))


# Plot out
# Local Moran's
# plotting local Moran's
all_locm_plots_base <- purrr::map(h_v_spattemp, ~ locm_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Ii = "Ii"
))
locm_plots_base <- plot_grid(plotlist = all_locm_plots_base)
save_plot(paste0(path, "/fits/Mimm/","locm_plots_base.png"), 
          locm_plots_base, base_asp = 1.6)

all_locm_plots_cov <- purrr::map(h_v_cov, ~ locm_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Ii = "Ii"
))
locm_plots_cov <- plot_grid(plotlist = all_locm_plots_cov)
save_plot(paste0(path, "/fits/Mimm/","locm_plots_cov.png"), locm_plots_cov)


# plotting Getis Ord
all_locg_plots_base <- purrr::map(h_v_spattemp, ~ locg_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Gi = "Gi"
))
gord_plots_base <- plot_grid(plotlist = all_locg_plots_base)
save_plot(paste0(path, "/fits/Mimm/","getis_ord_plots_base.png"), gord_plots_base)


all_locg_plots_cov <- purrr::map(h_v_cov, ~ locg_plot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  Gi = "Gi"
))
gord_plots_cov <- plot_grid(plotlist = all_locg_plots_cov)
save_plot(paste0(path, "/fits/Mimm/","getis_ord_plots_cov.png"), gord_plots_cov)


# plotting Getis Ord hotspots
all_locg_hotspot_plots_base <- purrr::map(h_v_spattemp, ~ locg_hotspot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  locg_p_clust = "locg_p_clust"
))
hotspot_plots_base <- plot_grid(plotlist = all_locg_hotspot_plots_base)
save_plot(paste0(path, "/fits/Mimm/","getis_ord_hotspot_plots_base.png"), hotspot_plots_base)


all_locg_hotspot_plots_cov <- purrr::map(h_v_cov, ~ locg_hotspot_fun(
  full_obj = .x$full_obj,
  ire      = .x$ire,
  coords   = .x$coords,
  hull     = .x$hull,
  locg_p_clust = "locg_p_clust"
))
hotspot_plots_cov <- plot_grid(plotlist = all_locg_hotspot_plots_cov)
save_plot(paste0(path, "/fits/Mimm/","getis_ord_hotspot_plots_cov.png"), 
          hotspot_plots_cov, base_asp = 1.6)

### Map spatial and spatiotemporal variability 
# DENSITY
base_density <- plot_map(predictions_spattemp$data, exp(est)) + 
  ggtitle("Base model: Density") +
  scale_fill_viridis_c(option = "G", direction = -1) +
  theme_map()
ggsave(filename = paste0(path, "/fits/Mimm/", "density_base.jpg"), plot = base_density, dpi = 300)

cov_density <- plot_map(predictions_d_s_substr$data, exp(est)) + 
  ggtitle("Covariate model: Density") +
  scale_fill_viridis_c(option = "G", direction = -1) +
  theme_map()
ggsave(filename = paste0(path, "/fits/Mimm/", "density_cov.jpg"), plot = cov_density, dpi = 300)

# Uncertainty of model estimate
predictions_base_d <- predict(mod_spattemp, 
                              newdata = grid_yrs[, c("X", "Y", "fyear", "fgear", "year")], 
                              nsim = 500,
                              sims_var = "est")
grid_yrs$est_se_base <- apply(predictions_base_d, 1, sd) 
rm(predictions_base_d)

predictions_cov_d <- predict(mod_dss, 
                              newdata = grid_yrs, 
                              nsim = 500)
grid_yrs$est_se_cov <- apply(predictions_cov_d, 1, sd) 
rm(predictions_cov_d)

est_uncertainty_base <- plot_map(grid_yrs, exp(est_se_base)) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model", fill = "Exp(SE)") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
ggsave(filename = paste0(path, "/fits/Mimm/","density_se_base.jpg"), plot = est_uncertainty_base)

est_uncertainty_cov <- plot_map(grid_yrs, exp(est_se_cov)) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model", fill = "Exp(SE)") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
ggsave(filename = paste0(path, "/fits/Mimm/","density_se_cov.jpg"), plot = est_uncertainty_cov)

# FIXED EFFECTS
base_ff <- plot_map2(predictions_spattemp$data, est_non_rf) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model", fill = "Estimate") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
ggsave(filename = paste0(path, "/fits/Mimm/","ff_effects_base.jpg"), plot = base_ff)

cov_ff <- plot_map2(predictions_d_s_substr$data, est_non_rf) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model", fill = "Estimate") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
ggsave(filename = paste0(path, "/fits/Mimm/","ff_effects_cov.jpg"), plot = cov_ff)

ff <- base_ff + cov_ff + 
  plot_annotation(title = "Fixed effects only",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")))
ggsave(filename = paste0(path, "/fits/Mimm/","ff_effects_base_cov.jpg"), plot = ff)


# SPATIAL EFFECTS
base_omega <- plot_map2(predictions_spattemp$data, omega_s) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model", fill = "Estimate") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
#ggsave(filename = paste0(path, "/fits/Mimm/","omega_effects_base.png"), plot = base_omega)

cov_omega <- plot_map2(predictions_d_s_substr$data, omega_s) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model", fill = "Estimate") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
#ggsave(filename = paste0(path, "/fits/Mimm/","omega_effects_cov.png"), plot = cov_omega)

# Uncertainty around spatial random effects
predictions_base_om <- predict(mod_spattemp, 
                               newdata = grid_yrs[, c("X", "Y", "fyear", "fgear", "year")], 
                               nsim = 500, 
                               sims_var = "omega_s")
grid_yrs$omega_se_base <- apply(predictions_base_om, 1, sd) 
rm(predictions_base_om)

predictions_cov_om <- predict(mod_dss, 
                              newdata = grid_yrs, 
                              nsim = 500, 
                              sims_var = "omega_s")
grid_yrs$omega_se_cov <- apply(predictions_cov_om, 1, sd) 
rm(predictions_cov_om)

omega_uncertainty_base <- plot_map2(grid_yrs, omega_se_base) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
#ggsave(filename = paste0(path, "/fits/Mimm/","omega_se_base.png"), plot = omega_uncertainty_base)

omega_uncertainty_cov <- plot_map2(grid_yrs, omega_se_cov) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
#ggsave(filename = paste0(path, "/fits/Mimm/","omega_se_cov.png"), plot = omega_uncertainty_cov)

spat_pred_plots <- base_omega + cov_omega + 
  plot_annotation(title = "Spatial predictions",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")))
spat_pred_se_plots <- omega_uncertainty_base + omega_uncertainty_cov + 
  plot_annotation(title = "Point-wise uncertainty in spatial predictions",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")))

omega <- wrap_elements(spat_pred_plots) / wrap_elements(spat_pred_se_plots)
ggsave(filename = paste0(path, "/fits/Mimm/","all_omega_and_se_base_cov.jpg"), plot = omega)


# SPATIOTEMPORAL EFFECTS
base_eps <- plot_map(predictions_spattemp$data, epsilon_st) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model: Spatiotemporal random effects only", fill = "Estimate") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
ggsave(filename = paste0(path, "/fits/Mimm/","eps_effects_base.jpg"), plot = base_eps)

cov_eps <- plot_map(predictions_d_s_substr$data, epsilon_st) + 
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model: Spatiotemporal random effects only", fill = "Estimate") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
ggsave(filename = paste0(path, "/fits/Mimm/","eps_effects_cov.jpg"), plot = cov_eps)


# Uncertainty of spatiotemporal effects
predictions_base_eps <- predict(mod_spattemp, 
                                newdata = grid_yrs[, c("X", "Y", "fyear", "fgear", "year")], 
                                nsim = 500, 
                                sims_var = "epsilon_st")
grid_yrs$eps_se_base <- apply(predictions_base_eps, 1, sd) 
rm(predictions_base_eps)

predictions_cov_eps <- predict(mod_dss, 
                               newdata = grid_yrs, 
                               nsim = 500, 
                               sims_var = "epsilon_st")
grid_yrs$eps_se_cov <- apply(predictions_cov_eps, 1, sd) 
rm(predictions_cov_eps)

eps_uncertainty_base <- plot_map(grid_yrs, eps_se_base) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Base model: Point-wise uncertainty in spatiotemporal predictions", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
ggsave(filename = paste0(path, "/fits/Mimm/","eps_se_base.jpg"), plot = eps_uncertainty_base)

eps_uncertainty_cov <- plot_map(grid_yrs, eps_se_cov) +
  scale_fill_viridis_c(option = "G", direction = -1) +
  labs(title = "Covariate model: Point-wise uncertainty in spatiotemporal predictions", fill = "SE") +
  theme_map() +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))
ggsave(filename = paste0(path, "/fits/Mimm/","eps_se_cov.jpg"), plot = eps_uncertainty_cov)

# Save all predictions
save(grid_yrs, file = paste0(path, "/fits/Mimm/","uncertainty_mapping_data.RData"))
preds2 <- list("base" = predictions_spattemp, 
              "cov_full" = predictions_d_s_substr)
saveRDS(preds2, file = paste0(path, "/fits/Mimm/","prediction_mapping_data.rds"))


