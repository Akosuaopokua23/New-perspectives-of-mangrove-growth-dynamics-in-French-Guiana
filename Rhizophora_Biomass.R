# scripts/biomass_models_bootstrap.R
# ----
# Prof. Dr. Marcelo Protazio & M.Sc. Michael Experiment
# Purpose: Fit six stand biomass models (FORET excluded), compute R²/AIC/RMSE,
# and draw bootstrap (B=200) prediction CIs in grey. Points are blue,
# curves red, each point labeled by dataset column "plot".
# ----

# ---- LOAD LIBARIES ----
library(ggplot2)
library(minpack.lm)
library(gridExtra)

# ---- DATA LOADING ----
DATA_PATH <- "C:/Users/agyekum/Documents/FrenchGuyiana/Script/GD44_new.csv"
amax <- 70
B <- 200
set.seed(123)
GDS <- read.csv(DATA_PATH, stringsAsFactors = FALSE)
remove_values <- c("LA1-2013", "LA1-2022", "TC6-2022", "KA17-2004", "KA8-2002")
GDS <- GDS[!(GDS$plot %in% remove_values), ]
sp <- unique(GDS$Species)
i3 <- which((GDS$Species == sp[3]) & (GDS$age < amax))
d3 <- data.frame(
  plot = GDS[i3, "plot"],
  age  = GDS[i3, "age"],
  AGB  = GDS[i3, "AGB"],    # expects AGB in column named "AGB"
  stringsAsFactors = FALSE
)
d3 <- d3[d3$age != 0 & is.finite(d3$AGB) & is.finite(d3$age), ]
row.names(d3) <- NULL
dage <- 0:100

# ---- [2] Model functions ----
f1 <- function(t1, t2, t3, t) t1 - t2 * exp(-t3 * t)
f2 <- function(t1, t2, t)     t1 * (1 - exp(-t2 * t))
f3 <- function(a, b, t)       a * t^b
f5 <- function(t1, t2, t3, t) t1 * exp(-t2 * exp(-t3 * t))
f6 <- function(t1, t2, t3, t) t1 / (1 + t2 * exp(-t3 * t))
f7 <- function(t1, t2, t3, t) t1 * (1 - t2 * exp(-t3 * t))

# ---- [3] Fit models ----
ctrl <- nls.control(maxiter = 1000, minFactor = 1/1024)
m1 <- nls(AGB ~ f1(t1, t2, t3, age),       start = c(t1=40, t2=5,   t3=0.1),  data=d3, algorithm="port", control=ctrl)
m2 <- nls(AGB ~ f2(t1, t2, age),           start = c(t1=40, t2=0.1),           data=d3, algorithm="port", control=ctrl)
m3 <- nls(AGB ~ f3(a, b, age),             start = c(a=40, b=0.1),             data=d3, algorithm="port", control=ctrl)
m5 <- nls(AGB ~ f5(t1, t2, t3, age),       start = c(t1=40, t2=0.2, t3=0.1),   data=d3, algorithm="port", control=ctrl)
m6 <- nls(AGB ~ f6(t1, t2, t3, age),       start = c(t1=40, t2=0.1, t3=0.01),  data=d3, algorithm="port", control=ctrl)
m7 <- nls(AGB ~ f7(t1, t2, t3, age),       start = c(t1=40, t2=1,   t3=0.3),   data=d3, algorithm="port", control=ctrl)

# ---- [4] Predictions (point estimates) ----
e1 <- predict(m1, newdata = list(age = dage))
e2 <- predict(m2, newdata = list(age = dage))
e3 <- predict(m3, newdata = list(age = dage))
e5 <- predict(m5, newdata = list(age = dage))
e6 <- predict(m6, newdata = list(age = dage))
e7 <- predict(m7, newdata = list(age = dage))

# ---- [5] Bootstrap CIs for prediction (B=200) ----
fit_once <- function(dat, model) {
  tryCatch({
    if (model == "weibull")       predict(nls(AGB ~ f1(t1, t2, t3, age), start=c(t1=40, t2=5, t3=0.1), data=dat, algorithm="port", control=ctrl), newdata=list(age=dage))
    else if (model == "chapman")  predict(nls(AGB ~ f2(t1, t2, age), start=c(t1=40, t2=0.1), data=dat, algorithm="port", control=ctrl), newdata=list(age=dage))
    else if (model == "power")    predict(nls(AGB ~ f3(a, b, age), start=c(a=40, b=0.1), data=dat, algorithm="port", control=ctrl), newdata=list(age=dage))
    else if (model == "gompertz") predict(nls(AGB ~ f5(t1, t2, t3, age), start=c(t1=40, t2=0.2, t3=0.1), data=dat, algorithm="port", control=ctrl), newdata=list(age=dage))
    else if (model == "logistic") predict(nls(AGB ~ f6(t1, t2, t3, age), start=c(t1=40, t2=0.1, t3=0.01), data=dat, algorithm="port", control=ctrl), newdata=list(age=dage))
    else if (model == "monomol")  predict(nls(AGB ~ f7(t1, t2, t3, age), start=c(t1=40, t2=1, t3=0.3), data=dat, algorithm="port", control=ctrl), newdata=list(age=dage))
    else stop("Unknown model key")
  }, error = function(e) rep(NA_real_, length(dage)))
}
boot_model <- function(model_key) {
  preds <- matrix(NA_real_, nrow = B, ncol = length(dage))
  n <- nrow(d3)
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    preds[b, ] <- fit_once(d3[idx, , drop = FALSE], model_key)
  }
  list(lower = apply(preds, 2, quantile, probs = 0.025, na.rm = TRUE),
       upper = apply(preds, 2, quantile, probs = 0.975, na.rm = TRUE),
       draws = preds)
}
ci1 <- boot_model("weibull")
ci2 <- boot_model("chapman")
ci3 <- boot_model("power")
ci5 <- boot_model("gompertz")
ci6 <- boot_model("logistic")
ci7 <- boot_model("monomol")

# ---- [6] Metrics table (R², AIC, RMSE) ----
metrics_of <- function(model, y, x_age) {
  yhat <- predict(model, newdata = list(age = x_age))
  sse <- sum((y - yhat)^2)
  sst <- sum((y - mean(y))^2)
  r2 <- 1 - sse / sst
  rmse <- sqrt(mean((y - yhat)^2))
  aic <- AIC(model)
  c(R2 = r2, AIC = aic, RMSE = rmse)
}
obs_y <- d3$AGB
obs_x <- d3$age
metrics_df <- as.data.frame(rbind(
  Weibull         = metrics_of(m1, obs_y, obs_x),
  ChapmanRichards = metrics_of(m2, obs_y, obs_x),
  Power           = metrics_of(m3, obs_y, obs_x),
  Gompertz        = metrics_of(m5, obs_y, obs_x),
  Logistic        = metrics_of(m6, obs_y, obs_x),
  Monomolecular   = metrics_of(m7, obs_y, obs_x)
))
metrics_df$model <- row.names(metrics_df)
metrics_df <- metrics_df[, c("model", "R2", "AIC", "RMSE")]
row.names(metrics_df) <- NULL
print(metrics_df)
write.csv(metrics_df, file = "model_metrics.csv", row.names = FALSE)

# ---- [7] Per-model panels with shaded CIs (points blue, curve red, labeled by 'plot') ----
shade_col <- grDevices::adjustcolor("grey70", alpha.f = 0.4)
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
plot_model <- function(x, y, labels, e_fit, ci, main_lbl) {
  plot(x, y, pch = 19, col = "blue", xlab = "Stand Age (years)", ylab = "AG Biomass (t dm/ha)")
  if (all(is.finite(ci$lower)) && all(is.finite(ci$upper))) {
    polygon(c(dage, rev(dage)), c(ci$lower, rev(ci$upper)), border = NA, col = shade_col)
  }
  lines(dage, e_fit, lwd = 2, col = "red")
  text(x, y, labels = labels, pos = 1, cex = 0.8)
  grid(nx = NULL, ny = NULL, col = 'gray', lty = 'dotted', lwd = 0.5)
  title(main = main_lbl)
}
plot_model(d3$age, d3$AGB, d3$plot, e1, ci1, "Weibull")
plot_model(d3$age, d3$AGB, d3$plot, e2, ci2, "Chapman-Richards")
plot_model(d3$age, d3$AGB, d3$plot, e3, ci3, "Power")
plot_model(d3$age, d3$AGB, d3$plot, e5, ci5, "Gompertz")
plot_model(d3$age, d3$AGB, d3$plot, e6, ci6, "Logistic")
plot_model(d3$age, d3$AGB, d3$plot, e7, ci7, "Monomolecular")
par(mfrow = c(1, 1))

# ---- [8] (Optional) Quick combined plot (no shading) ----
make_combined <- FALSE
if (make_combined) {
  pdf("dbh-Rhi.pdf", width = 8, height = 6)
  plot(d3$age, d3$AGB, pch = 19, col = "blue", xlab = "Age (years)", ylab = "AG Biomass (t dm/ha)")
  lines(dage, e1, lty = 1)
  lines(dage, e2, lty = 2)
  lines(dage, e3, lty = 3)
  lines(dage, e5, lty = 4)
  lines(dage, e6, lty = 5)
  lines(dage, e7, lty = 6)
  legend("topleft",
         legend = c("Observed", "Weibull", "Chapman-Richards", "Power", "Gompertz", "Logistic", "Monomolecular"),
         bty = "n", pch = c(19, NA, NA, NA, NA, NA, NA), lty = c(NA, 1, 2, 3, 4, 5, 6))
  grid()
  dev.off()
}

# ---- [END] ----

