
# ----
# Stand DBH Modeling for Avicennia Germinans (French Guiana)
# Author: Michael Kyei Agyekum
# Date: 2025-09-19
# Description: Fits multiple growth models to DBH data; computes R2, AIC, RMSE, MAE; produces plots.
# Usage: Place 'GD44_new.csv' in a 'data' folder in repo root. Run script in R (>=4.x).

# ---- PACKAGE MANAGEMENT ----
required_packages <- c("ggplot2", "minpack.lm", "gridExtra")
new_packages <- setdiff(required_packages, rownames(installed.packages()))
if(length(new_packages) > 0) install.packages(new_packages)
lapply(required_packages, library, character.only = TRUE)

# ---- DATA LOADING ----
DATA_PATH <- "data/GD44_new.csv"
GDS <- read.csv(DATA_PATH)

# ---- DATA FILTERING ----
GDSA <- subset(GDS, Species == "A" & age < 70)
remove_values <- c("LA1-2013", "LA1-2022", "TC6-2022", "KA17-2004")
GDSA <- GDSA[!GDSA$plot %in% remove_values,]

d1 <- data.frame(age = GDSA$age, DBH = GDSA$DBH, plot = GDSA$plot)
d1 <- d1[d1$age != 0, ]
dage <- 0:100

# ---- MODEL FUNCTION DEFINITIONS ----
f1 <- function(t1, t2, t3, t) { t1 - t2 * exp(-t3 * t) }
f2 <- function(t1, t2, t) { t1 * (1 - exp(-t2 * t)) }
f3 <- function(a, b, t) { a * t^b }
f5 <- function(t1, t2, t3, t) { t1 * exp(-t2 * exp(-t3 * t)) }
f6 <- function(t1, t2, t3, t) { t1 / (1 + t2 * exp(-t3 * t)) }
f7 <- function(t1, t2, t3, t) { t1 * (1 - t2 * exp(-t3 * t)) }

# ---- METRIC HELPERS ----
get_r2 <- function(obs, pred) cor(obs, pred)^2
get_rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))
get_mae <- function(obs, pred) mean(abs(obs - pred))

# ---- MODEL FITTING ----
model_list <- list(
  Weibull = list(
    fit = nls(DBH ~ f1(t1, t2, t3, age), data=d1, start=c(t1=125, t2=50, t3=0.1),
              control=nls.control(maxiter=1000, minFactor=1/1024)),
    predict_fun = function(fit, ages) predict(fit, newdata=list(age=ages))
  ),
  ChapmanRich = list(
    fit = nls(DBH ~ f2(t1, t2, age), data=d1, start=c(t1=125, t2=0.1),
              control=nls.control(maxiter=1000, minFactor=1/1024)),
    predict_fun = function(fit, ages) predict(fit, newdata=list(age=ages))
  ),
  Power = list(
    fit = nls(DBH ~ f3(a, b, age), data=d1, start=c(a=125, b=0.1),
              control=nls.control(maxiter=1000, minFactor=1/1024)),
    predict_fun = function(fit, ages) predict(fit, newdata=list(age=ages))
  ),
  Gompertz = list(
    fit = nls(DBH ~ f5(t1, t2, t3, age), data=d1, start=c(t1=125, t2=0.1, t3=0.01),
              control=nls.control(maxiter=1000, minFactor=1/1024)),
    predict_fun = function(fit, ages) predict(fit, newdata=list(age=ages))
  ),
  Logistic = list(
    fit = nls(DBH ~ f6(t1, t2, t3, age), data=d1, start=c(t1=125, t2=0.1, t3=0.01),
              control=nls.control(maxiter=1000, minFactor=1/1024)),
    predict_fun = function(fit, ages) predict(fit, newdata=list(age=ages))
  ),
  Monomolecular = list(
    fit = nls(DBH ~ f7(t1, t2, t3, age), data=d1, start=c(t1=125, t2=1, t3=0.3),
              control=nls.control(maxiter=1000, minFactor=1/1024)),
    predict_fun = function(fit, ages) predict(fit, newdata=list(age=ages))
  )
)

# ---- METRIC CALCULATION & BOOTSTRAP CI ----
metric_table <- data.frame(Model = character(), R2 = numeric(), AIC = numeric(), RMSE = numeric(), MAE = numeric(), stringsAsFactors = FALSE)
plot_list <- list()
n_boot <- 200
set.seed(42)

for(model_name in names(model_list)) {
  fit <- model_list[[model_name]]$fit
  predict_fun <- model_list[[model_name]]$predict_fun
  
  pred_train <- predict_fun(fit, d1$age)
  obs <- d1$DBH
  
  r2 <- get_r2(obs, pred_train)
  rmse <- get_rmse(obs, pred_train)
  mae <- get_mae(obs, pred_train)
  aic <- AIC(lm(obs ~ pred_train - 1))
  
  metric_table <- rbind(metric_table, data.frame(Model = model_name, R2 = round(r2, 3), AIC = round(aic, 2), RMSE = round(rmse, 3), MAE = round(mae, 3)))
  
  # Bootstrap for CIs
  boot_preds <- replicate(n_boot, {
    ix <- sample(nrow(d1), replace = TRUE)
    d1b <- d1[ix, ]
    fitb <- tryCatch(
      suppressWarnings(nls(formula(fit), data = d1b, start = as.list(coef(fit)), control = nls.control(maxiter=1000, minFactor=1/1024))),
      error = function(e) NULL
    )
    if(!is.null(fitb)) predict_fun(fitb, dage) else rep(NA, length(dage))
  })
  
  ci_mean <- apply(boot_preds, 1, mean, na.rm=TRUE)
  ci_lower <- apply(boot_preds, 1, quantile, probs=0.025, na.rm=TRUE)
  ci_upper <- apply(boot_preds, 1, quantile, probs=0.975, na.rm=TRUE)
  ci_df <- data.frame(age = dage, fit = ci_mean, lower = ci_lower, upper = ci_upper)
  
  plt <- ggplot() +
    geom_point(data = d1, aes(x = age, y = DBH), color = "red", size = 2) +
    geom_text(data = d1, aes(x = age, y = DBH, label = plot), vjust = -1, size = 2) +
    geom_line(data = ci_df, aes(x = age, y = fit), color = "black", size = 1.1) +
    geom_ribbon(data = ci_df, aes(x = age, ymin = lower, ymax = upper), fill = "grey40", alpha = 0.30) +
    labs(title = model_name, x = "Stand Age (years)", y = "Average DBH (cm)") +
    xlim(0, 60) + ylim(0, 80) +
    theme_bw(base_size = 13)
  plot_list[[model_name]] <- plt
}

# ---- PRINT AND PLOT RESULTS ----
print(metric_table)

# Arrange all plots (2 rows x 3 columns)
grid.arrange(
  plot_list$Weibull,
  plot_list$ChapmanRich,
  plot_list$Power,
  plot_list$Gompertz,
  plot_list$Logistic,
  plot_list$Monomolecular,
  nrow = 2
)

# ---- END OF SCRIPT ----


