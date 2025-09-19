# G4_new_Rhizo_DBH_2025.R
# ----
# DBH Modeling for Rhizophora (French Guiana)
# Author: [Insert Name]
# Date: 2025-09-19
# Usage: Place 'GD44_new.csv' in a 'data/' folder in repo root. Run in R or RStudio.

# ---- LOAD LIBARIES ----
library(ggplot2)
library(minpack.lm)
library(gridExtra)

# ---- DATA LOADING ----
DATA_PATH <- "C:/Users/agyekum/Documents/FrenchGuyiana/Script/GD44_new.csv"
GDS <- read.csv(DATA_PATH)

# ---- DATA FILTERING ----
remove_values <- c("LA1-2013", "LA1-2022", "TC6-2022", "KA17-2004")
GDSR <- subset(GDS, Species == "R" & age < 70 & !plot %in% remove_values)
d1 <- data.frame(age = GDSR$age, DBH = GDSR$DBH, plot = GDSR$plot)
d1 <- d1[d1$age != 0, ]
dage <- 0:100

# ---- MODEL FUNCTION DEFINITIONS ----
f1 <- function(t1, t2, t3, t) t1 - t2 * exp(-t3 * t)
f2 <- function(t1, t2, t) t1 * (1 - exp(-t2 * t))
f3 <- function(a, b, t) a * t^b
f5 <- function(t1, t2, t3, t) t1 * exp(-t2 * exp(-t3 * t))
f6 <- function(t1, t2, t3, t) t1 / (1 + t2 * exp(-t3 * t))
f7 <- function(t1, t2, t3, t) t1 * (1 - t2 * exp(-t3 * t))

# ---- MODEL FITTING ----
m1 <- nls(DBH ~ f1(t1, t2, t3, age), data = d1,
          start = c(t1 = 35, t2 = 15, t3 = 0.05),
          algorithm = "port", control = nls.control(maxiter = 2000, minFactor = 1/1024))
m2 <- nls(DBH ~ f2(t1, t2, age), data = d1,
          start = c(t1 = 35, t2 = 0.07),
          algorithm = "port", control = nls.control(maxiter = 1000, minFactor = 1/1024))
m3 <- nls(DBH ~ f3(a, b, age), data = d1,
          start = c(a = 2, b = 0.7),
          algorithm = "port", control = nls.control(maxiter = 1000, minFactor = 1/1024))
m5 <- nls(DBH ~ f5(t1, t2, t3, age), data = d1,
          start = c(t1 = 35, t2 = 0.17, t3 = 0.012),
          algorithm = "port", control = nls.control(maxiter = 1000, minFactor = 1/1024))
m6 <- nls(DBH ~ f6(t1, t2, t3, age), data = d1,
          start = c(t1 = 35, t2 = 0.11, t3 = 0.011),
          algorithm = "port", control = nls.control(maxiter = 1000, minFactor = 1/1024))
m7 <- nls(DBH ~ f7(t1, t2, t3, age), data = d1,
          start = c(t1 = 32, t2 = 0.9, t3 = 0.29),
          algorithm = "port", control = nls.control(maxiter = 1000, minFactor = 1/1024))

fit_list <- list(
  Weibull = m1,
  ChapmanRich = m2,
  Power = m3,
  Gompertz = m5,
  Logistic = m6,
  Monomolecular = m7
)
formula_list <- list(
  Weibull = DBH ~ f1(t1, t2, t3, age),
  ChapmanRich = DBH ~ f2(t1, t2, age),
  Power = DBH ~ f3(a, b, age),
  Gompertz = DBH ~ f5(t1, t2, t3, age),
  Logistic = DBH ~ f6(t1, t2, t3, age),
  Monomolecular = DBH ~ f7(t1, t2, t3, age)
)

# ---- METRICS ----
get_r2 <- function(obs, fit) cor(obs, fit)^2
get_rmse <- function(obs, fit) sqrt(mean((obs-fit)^2))
pred_train_list <- lapply(fit_list, function(fit) predict(fit, newdata = list(age = d1$age)))
metric_table <- data.frame(
  Model = names(fit_list),
  R2 = sapply(pred_train_list, function(p) get_r2(d1$DBH, p)),
  AIC = sapply(fit_list, AIC),
  RMSE = sapply(pred_train_list, function(p) get_rmse(d1$DBH, p)),
  row.names = NULL
)
print(metric_table)

# ---- BOOTSTRAP CI & PLOTTING ----
n_boot <- 200
set.seed(42)
plot_list <- list()
for (name in names(fit_list)) {
  fit <- fit_list[[name]]
  fmla <- formula_list[[name]]
  coef_vec <- as.numeric(coef(fit))
  names(coef_vec) <- names(coef(fit))
  boot_mat <- matrix(NA, n_boot, length(dage))
  for(i in 1:n_boot){
    idx <- sample(nrow(d1), replace=TRUE)
    boot_dat <- d1[idx,]
    pert <- coef_vec * runif(length(coef_vec), 0.95, 1.05)
    fitb <- tryCatch(
      nls(fmla, data=boot_dat, start=as.list(pert),
          algorithm="port", control=nls.control(maxiter=2000, minFactor=1e-6, warnOnly=TRUE)),
      error=function(e) NULL)
    boot_mat[i,] <- if(!is.null(fitb)) predict(fitb, newdata=list(age=dage)) else rep(NA, length(dage))
  }
  ci_mean <- apply(boot_mat, 2, mean, na.rm=TRUE)
  ci_lower <- apply(boot_mat, 2, quantile, probs=0.025, na.rm=TRUE)
  ci_upper <- apply(boot_mat, 2, quantile, probs=0.975, na.rm=TRUE)
  ci_df <- data.frame(age=dage, fit=ci_mean, lower=ci_lower, upper=ci_upper)
  plt <- ggplot() +
    geom_point(data=d1, aes(x=age, y=DBH), color='blue', size=2) +
    geom_text(data=d1, aes(x=age, y=DBH, label=plot), color='black', vjust=-1, size=2) +
    geom_ribbon(data=ci_df, aes(x=age, ymin=lower, ymax=upper), fill="grey60", alpha=0.33) +
    geom_line(data=ci_df, aes(x=age, y=fit), color="red", size=1.2) +
    labs(title=name, x="Stand Age (years)", y="Average DBH (cm)") +
    xlim(0,60) + ylim(0,30) +
    theme_minimal(base_size=13)
  plot_list[[name]] <- plt
}

# ---- MULTI-PANEL GRID ----
grid.arrange(
  plot_list$Weibull,
  plot_list$ChapmanRich,
  plot_list$Power,
  plot_list$Gompertz,
  plot_list$Logistic,
  plot_list$Monomolecular,
  nrow=2
)

# ---- END OF SCRIPT ----
