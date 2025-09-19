

# Stand data on Above-Ground Biomass (AGB) Modeling for Avicennia Germinans (French Guiana)
# Author: Michael Kyei Ageykum
# Date: 2025-09-19


# ---- LOAD LIBARIES ----
library(ggplot2)
library(minpack.lm)
library(gridExtra)

# ---- DATA LOADING ----
DATA_PATH <- "C:/Users/agyekum/Documents/FrenchGuyiana/Script/GD44_new.csv"
GDS <- read.csv(DATA_PATH, stringsAsFactors = FALSE)

# ---- DATA FILTERING ----
remove_values <- c("LA1-2013", "LA1-2022", "TC6-2022", "KA17-2004")
GDSA <- subset(GDS, Species == "A" & age < 70 & AGB > 0 & !plot %in% remove_values)
d1 <- data.frame(age = GDSA$age, AGB = GDSA$AGB, plot = GDSA$plot)
d1 <- subset(d1, age != 0)
dage <- 0:60

# ---- MODEL FUNCTION DEFINITIONS ----
f1 <- function(t1, t2, t3, t, m) t1 - t2 * exp(-t3 * t^m)
f2 <- function(t1, t2, t) t1 * (1 - exp(-t2 * t))
f3 <- function(a, b, t) a * t^b
f5 <- function(t1, t2, t3, t) t1 * exp(-t2 * exp(-t3 * t))
f6 <- function(t1, t2, t3, t) t1 / (1 + t2 * exp(-t3 * t))
f7 <- function(t1, t2, t3, t) t1 * (1 - t2 * exp(-t3 * t))

# ---- MODEL FITTING ----
m1 <- nls(AGB ~ f1(t1, t2, t3, age, m), data = d1,
          start = c(t1 = 100, t2 = 10, t3 = 0.01, m = 0.3),
          algorithm = "port", lower = c(0,0,0,0), upper = c(1000,1000,10,1),
          control = nls.control(maxiter = 2000, minFactor = 1e-6, warnOnly = TRUE))
m2 <- nls(AGB ~ f2(t1, t2, age), data = d1,
          start = c(t1 = 300, t2 = 0.1),
          algorithm = "port", control = nls.control(maxiter = 1000, minFactor = 1/1024))
m3 <- nls(AGB ~ f3(a, b, age), data = d1,
          start = c(a = 300, b = 0.1),
          algorithm = "port", control = nls.control(maxiter = 1000, minFactor = 1/1024))
m5 <- nls(AGB ~ f5(t1, t2, t3, age), data = d1,
          start = c(t1 = 300, t2 = 0.1, t3 = 0.01),
          algorithm = "port", control = nls.control(maxiter = 1000, minFactor = 1/1024))
m6 <- nls(AGB ~ f6(t1, t2, t3, age), data = d1,
          start = c(t1 = 300, t2 = 0.1, t3 = 0.01),
          algorithm = "port", control = nls.control(maxiter = 1000, minFactor = 1/1024))
m7 <- nls(AGB ~ f7(t1, t2, t3, age), data = d1,
          start = c(t1 = 300, t2 = 1, t3 = 0.3),
          algorithm = "port", control = nls.control(maxiter = 1000, minFactor = 1/1024))

# ---- METRICS ----
get_r2 <- function(obs, pred) cor(obs, pred)^2
get_rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))
get_mae <- function(obs, pred) mean(abs(obs - pred))

pred1 <- predict(m1, newdata = list(age = d1$age, m = coef(m1)["m"]))
pred2 <- predict(m2, newdata = list(age = d1$age))
pred3 <- predict(m3, newdata = list(age = d1$age))
pred5 <- predict(m5, newdata = list(age = d1$age))
pred6 <- predict(m6, newdata = list(age = d1$age))
pred7 <- predict(m7, newdata = list(age = d1$age))
obs <- d1$AGB

r2_vals <- c(
  Weibull = get_r2(obs, pred1),
  ChapmanRich = get_r2(obs, pred2),
  Power = get_r2(obs, pred3),
  Gompertz = get_r2(obs, pred5),
  Logistic = get_r2(obs, pred6),
  Monomolecular = get_r2(obs, pred7)
)
aic_vals <- setNames(
  AIC(
    lm(obs ~ pred1 - 1),
    lm(obs ~ pred2 - 1),
    lm(obs ~ pred3 - 1),
    lm(obs ~ pred5 - 1),
    lm(obs ~ pred6 - 1),
    lm(obs ~ pred7 - 1)
  )[,"AIC"],
  names(r2_vals)
)
rmse_vals <- c(
  Weibull = get_rmse(obs, pred1),
  ChapmanRich = get_rmse(obs, pred2),
  Power = get_rmse(obs, pred3),
  Gompertz = get_rmse(obs, pred5),
  Logistic = get_rmse(obs, pred6),
  Monomolecular = get_rmse(obs, pred7)
)


metrics <- data.frame(
  Model = names(r2_vals),
  R2 = round(as.numeric(r2_vals), 3),
  AIC = round(as.numeric(aic_vals), 2),
  RMSE = round(as.numeric(rmse_vals), 3),
  row.names = NULL,
  stringsAsFactors = FALSE
)
print(metrics)

# ---- BOOTSTRAP FUNCTION ----
robust_boot <- function(orig_fit, boot_fun, dage, d1, lower, upper) {
  n_boot <- 200
  pnames <- names(coef(orig_fit))
  base_start <- as.numeric(coef(orig_fit))
  boot_mat <- matrix(NA, n_boot, length(dage))
  for (i in 1:n_boot) {
    idx <- sample(nrow(d1), replace=TRUE)
    boot_data <- d1[idx, ]
    perturbed <- base_start * runif(length(base_start), 0.95, 1.05)
    names(perturbed) <- pnames
    fitb <- tryCatch(
      boot_fun(boot_data, as.list(perturbed)),
      error=function(e) NULL
    )
    boot_mat[i, ] <- if(!is.null(fitb)) tryCatch({
      if ("m" %in% names(coef(fitb)))
        predict(fitb, newdata = list(age = dage, m = coef(fitb)["m"]))
      else
        predict(fitb, newdata = list(age = dage))
    }, error=function(e) rep(NA, length(dage))) else rep(NA, length(dage))
  }
  boot_mat
}

# ---- MODEL DEFINITIONS FOR PLOTTING ----
model_defs <- list(
  Weibull = list(
    fit = m1,
    boot_fun = function(data, start) nls(AGB ~ f1(t1,t2,t3,age,m), data=data, start=start, 
                                         algorithm="port", lower=c(0,0,0,0), upper=c(1000,1000,10,1),
                                         control=nls.control(maxiter=3000,minFactor=1e-10,warnOnly=TRUE)),
    predict_plot = function() predict(m1, newdata = list(age = dage, m=coef(m1)["m"]))
  ),
  ChapmanRich = list(
    fit = m2,
    boot_fun = function(data, start) nls(AGB ~ f2(t1,t2,age), data=data, start=start, 
                                         algorithm="port", control=nls.control(maxiter=3000,minFactor=1e-10,warnOnly=TRUE)),
    predict_plot = function() predict(m2, newdata = list(age = dage))
  ),
  Power = list(
    fit = m3,
    boot_fun = function(data, start) nls(AGB ~ f3(a,b,age), data=data, start=start, 
                                         algorithm="port", control=nls.control(maxiter=3000,minFactor=1e-10,warnOnly=TRUE)),
    predict_plot = function() predict(m3, newdata = list(age = dage))
  ),
  Gompertz = list(
    fit = m5,
    boot_fun = function(data, start) nls(AGB ~ f5(t1,t2,t3,age), data=data, start=start, 
                                         algorithm="port", control=nls.control(maxiter=3000,minFactor=1e-10,warnOnly=TRUE)),
    predict_plot = function() predict(m5, newdata = list(age = dage))
  ),
  Logistic = list(
    fit = m6,
    boot_fun = function(data, start) nls(AGB ~ f6(t1,t2,t3,age), data=data, start=start, 
                                         algorithm="port", control=nls.control(maxiter=3000,minFactor=1e-10,warnOnly=TRUE)),
    predict_plot = function() predict(m6, newdata = list(age = dage))
  ),
  Monomolecular = list(
    fit = m7,
    boot_fun = function(data, start) nls(AGB ~ f7(t1,t2,t3,age), data=data, start=start, 
                                         algorithm="port", control=nls.control(maxiter=3000,minFactor=1e-10,warnOnly=TRUE)),
    predict_plot = function() predict(m7, newdata = list(age = dage))
  )
)

# ---- PLOTTING LOOP ----
plot_list <- list()
for(model_name in names(model_defs)) {
  def <- model_defs[[model_name]]
  orig_fit <- def$fit
  boot_mat <- robust_boot(orig_fit, def$boot_fun, dage, d1, c(0,0,0,0), c(1000,1000,10,1))
  ci_mean <- apply(boot_mat, 2, function(x) mean(x, na.rm=TRUE))
  ci_lower <- apply(boot_mat, 2, function(x) quantile(x, .025, na.rm=TRUE))
  ci_upper <- apply(boot_mat, 2, function(x) quantile(x, .975, na.rm=TRUE))
  ci_df <- data.frame(age = dage, fit = ci_mean, lower = ci_lower, upper = ci_upper)
  mean_curve <- data.frame(age = dage, mean = def$predict_plot())
  plt <- ggplot() +
    geom_ribbon(data = ci_df, aes(x = age, ymin = lower, ymax = upper), fill = "grey40", alpha = 0.3) +
    geom_line(data = mean_curve, aes(x = age, y = mean), color = 'black', size = 1.1) +
    geom_point(data = d1, aes(x = age, y = AGB), color = 'red', size = 2) +
    geom_text(data = d1, aes(x = age, y = AGB, label = plot), vjust = -1, size = 2) +
    labs(title = model_name, x = "Stand Age (years)", y = "AG Biomass (t dm/ha)") +
    xlim(0, 60) + ylim(0, 300) +
    theme_bw(base_size = 13)
  plot_list[[model_name]] <- plt
}

# ---- PRINT AND PLOT RESULTS ----
print(metrics)
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
