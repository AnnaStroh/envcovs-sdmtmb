#######
## Internal consistency
##
#######


library(ggplot2); theme_set(theme_bw())
library(smplot2)
library(patchwork)

path <- "C:/Users/astroh/Desktop/Chapter 2/sdmTMB/WHG/"

# Load abundance indices for each age class
load(paste0(path, "/fits/0/", "allindices.RData"))
index_0 <- all_indices
index_0$age <- 0
load(paste0(path, "/fits/1/", "allindices.RData"))
index_1 <- all_indices
index_1$age <- 1
load(paste0(path, "/fits/2/", "allindices.RData"))
index_2 <- all_indices
index_2$age <- 2
rm(all_indices)

# Combine all indices into one dataframe
indices <- rbind(index_0, index_1, index_2)

# Subset dataframe into base and covariate model-bsed indices
indices_base <- subset(indices, model == "base")
indices_cov <- subset(indices, model != "base")

# Track cohort

# for base model
base_sub <- indices_base[, names(indices_base) %in% c("year", "age", "est")]
base_sub$cohort <- with(base_sub, year - age)
df_base_wide <- reshape::cast(base_sub, cohort ~ age, value = "est")
colnames(df_base_wide)[2:4] <- c("age0", "age1", "age2")

with(df_base_wide, plot(age0, age1))
with(df_base_wide, plot(age1, age2))


# for covariate model
cov_sub <- indices_cov[, names(indices_cov) %in% c("year", "age", "est")]
cov_sub$cohort <- with(cov_sub, year - age)
df_cov_wide <- reshape::cast(cov_sub, cohort ~ age, value = "est")
colnames(df_cov_wide)[2:4] <- c("age0", "age1", "age2")

with(df_cov_wide, plot(age0, age1))
with(df_cov_wide, plot(age1, age2))


# Calculate within-cohort correlation scores

base_int01 <- cor.test(df_base_wide$age0, df_base_wide$age1)
base_int12 <- cor.test(df_base_wide$age1, df_base_wide$age2)

cov_int01 <- cor.test(df_cov_wide$age0, df_cov_wide$age1)
cov_int12 <- cor.test(df_cov_wide$age1, df_cov_wide$age2)

# Visualise internal consistency

int01_base <- ggplot(df_base_wide, aes(age0, age1)) +
  geom_point() +
  sm_statCorr(text_size = 3) +
  #labs(title = "Base model: Age-0 to age-1 whiting") +
  scale_y_continuous(labels = scales::label_scientific(digits = 1)) +
  xlab("Age-0 Whiting Abundance") + ylab("Age-1 Whiting Abundance") +
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 6))

int12_base <- ggplot(df_base_wide, aes(age1, age2)) +
  geom_point() +
  sm_statCorr(text_size = 3, label_y = 2.7e+06) +
  #labs(title = "Base model: Age-1 to age-2 whiting") +
  scale_x_continuous(labels = scales::label_scientific(digits = 1)) +
  xlab("Age-1 Whiting Abundance") + ylab("Age-2 Whiting Abundance") +
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 6))

int01_cov <- ggplot(df_cov_wide, aes(age0, age1)) +
  geom_point() +
  sm_statCorr(text_size = 3) +
  #labs(title = "Covariate model: Age-0 to age-1 whiting") +
  xlab("Age-0 Whiting Abundance") + ylab("Age-1 Whiting Abundance") +
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 6))

int12_cov <- ggplot(df_cov_wide, aes(age1, age2)) +
  geom_point() +
  sm_statCorr(text_size = 3, label_y = 4.5e+06) +
  #labs(title = "Covariate model: Age-1 to age-2 whiting") +
  xlab("Age-1 Whiting Abundance") + ylab("Age-2 Whiting Abundance") +
  theme(axis.title = element_text(size = 7), 
        axis.text = element_text(size = 6))

# Create final plot

int_plot <- (int01_base + int12_base) / (int01_cov + int12_cov) +
  plot_annotation(tag_levels = 'a')
ggsave(filename = paste0(path, "/", "internal_consistency.png"), 
       plot = int_plot, dpi = 300, 
       height = 170, width = 159, units = "mm") 






















