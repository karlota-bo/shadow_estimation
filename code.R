library(readxl)
library(urca)
library(seasonal)
library(tempdisagg)
library(tidyr)
library(lavaan)
library(forecast)
library(zoo)
library(dynr)
library(tseries) 
library(forecast) 

# annual data

dataset <- read_excel("dataset.xlsx")

# quarterly data

ketv <- read_excel("dataset.xlsx", sheet = "Ketvirtiniai")

# variables that were only annual 2005-2006
VAT2005_2006 <- c(4841689, 6152214)/3.4528    #conversion from LT to EUR
CIT2005_2006 <- c(1507704, 1924521)/3.4528
PIT2005_2006 <- c(3566302, 4059245)/3.4528
excise2005_2006 <- c(2040088, 2374429)/3.4528
wage2005_2006 <- c(369.6, 433.2)*1.289          # wage adjustment
real_wage2005_2006 <- wage2005_2006/c(77.89, 74.59)*100 #adjusting for inflation

#interpolation
VAT <- predict(td(VAT2005_2006 ~ 1, to = 4, method = "denton-cholette"))
CIT <- predict(td(CIT2005_2006 ~ 1, to = 4, method = "denton-cholette"))
PIT <- predict(td(PIT2005_2006 ~ 1, to = 4, method = "denton-cholette"))
excise <- predict(td(excise2005_2006 ~ 1, to = 4, method = "denton-cholette"))
wage <- predict(td(wage2005_2006 ~ 1, to = 4, method = "denton-cholette",
                   conversion = "mean"))
real_wage <- predict(td(real_wage2005_2006 ~ 1, to = 4, method = 
                          "denton-cholette", conversion = "mean"))

# taxes as %GDP
ketv$VAT <- c(VAT, ketv$VAT[9:80])/1000/ketv$nom_GDP*100
ketv$CIT <- c(CIT, ketv$CIT[9:80])/1000/ketv$nom_GDP*100
ketv$PIT <- c(PIT, ketv$PIT[9:80])/1000/ketv$nom_GDP*100
ketv$excise <- c(excise, ketv$excise[9:80])/1000/ketv$nom_GDP*100

ketv$contr_of_corrupt <- c(predict(td(dataset$contr_of_corrupt ~ 1,
                                    to = 4, method = "denton-cholette",
                                    conversion = "mean")), rep(NA, 4))
ketv$gov_effect <- c(predict(td(dataset$gov_effect ~ 1,
                                to = 4, method = "denton-cholette",
                                conversion = "mean")), rep(NA, 4))
ketv$polit_stability <- c(predict(td(dataset$polit_stability ~ 1,
                                     to = 4, method = "denton-cholette",
                                     conversion = "mean")), rep(NA, 4))
ketv$regul_quality <- c(predict(td(dataset$regul_quality ~ 1,
                                     to = 4, method = "denton-cholette",
                                     conversion = "mean")), rep(NA, 4))
ketv$rule_of_law <- c(predict(td(dataset$rule_of_law ~ 1,
                                 to = 4, method = "denton-cholette",
                                 conversion = "mean")), rep(NA, 4))
ketv$kaitz <- predict(td(dataset$kaitz ~ 1, to = 4, method = "denton-cholette",
                           conversion = "mean"))

# converting to time series
ts_list <- list()
for (col in names(ketv)) {
  ts_list[[col]] <- ts(ketv[[col]], start = c(2005, 1), frequency = 4)
}

par(mfrow=c(2,3))
# seasonal adjustments
m <- seas(ts_list$participation_rate)
ketv$participation_rate <- final(m) 
plot(m, main = "Participation rate")
m <- seas(ts_list$CIT)
ketv$CIT <- final(m)
plot(m, main = "CIT")
m <- seas(ts_list$VAT)
ketv$VAT <- final(m)
plot(m, main = "VAT")
m <- seas(ts_list$excise)
ketv$excise <- final(m)
plot(m, main = "Excise")
m <- seas(ts_list$real_wage)
ketv$real_wage[9:80] <- final(m)
plot(m, main = "Real wage")
ketv$real_wage[1:8] <- real_wage

df_quart <- list()
for (col in names(ketv)) {
  df_quart[[col]] <- ts(ketv[[col]], start = c(2005, 1), frequency = 4)
}

df.quart <- data.frame(df_quart)
df.quart <- df.quart[1:76, 2:21] # dropping Period column and rows with NA

# Testing for unit roots
results <- data.frame(
  Variable = character(),
  Test_Statistic = numeric(),
  Critical_5pct = numeric(),
  Stationary = character(),
  stringsAsFactors = FALSE
)


for(col in names(df.quart)) {
  if(is.numeric(df.quart[[col]])) {
    test <- ur.df(df.quart[[col]], type = "drift")
    test_stat <- test@teststat[1]
    crit_5 <- test@cval[1, 2]
    results <- rbind(results, data.frame(
      Variable = col,
      Test_Statistic = test_stat,
      Critical_5pct = crit_5,
      Stationary = ifelse(test_stat < crit_5, "Yes", "No")
    ))
  }
}

results  #political stability, kaitz and VAT I(0)

# data transformation for variables with unit root
# skipping variables that are in percentages when taking logarithm
data <- cbind(df.quart[, c(1:8, 14:20)], data.frame(lapply(df.quart[,c(9:13)], log)))
data_transformed <- data.frame(lapply(data[,-c(6, 8, 13)], diff)) # - I(0) variables
data_transformed$kaitz <- data[2:76,8]
data_transformed$VAT <- data[2:76, 6]
data_transformed$polit_stability <- data[2:76, 13]

# testing tranformed data for unit roots
results <- data.frame(
  Variable = character(),
  Test_Statistic = numeric(),
  Critical_5pct = numeric(),
  Stationary = character(),
  stringsAsFactors = FALSE
)

df_test <- list()
for (col in names(data_transformed)) {
  df_test[[col]] <- ts(data_transformed[[col]], start = c(2005, 2), frequency = 4)
}
df_test <- data.frame(df_test)

for(col in names(df_test)) {
  if(is.numeric(df_test[[col]])) {
    test <- ur.df(df_test[[col]])
    test_stat <- test@teststat[1]
    crit_5 <- test@cval[1, 2]
    results <- rbind(results, data.frame(
      Variable = col,
      Test_Statistic = test_stat,
      Critical_5pct = crit_5,
      Stationary = ifelse(test_stat < crit_5, "Yes", "No")
    ))
  }
}
results

######################################################## MIMIC model
data_transformed$kaitz <- scale(data_transformed$kaitz)
model1 <- '
  shadow =~ -1 * real_GDP + real_GFCF + cash.GDP
  shadow ~ real_wage + kaitz + regul_quality
'

fit_mimic <- sem(model1, data = data_transformed, std.lv = F)

summary(fit_mimic, fit.measures = T)

plot(residuals(fit_mimic, type = "casewise")[,1], type = 'l')
checkresiduals(residuals(fit_mimic, type = "casewise")[,1])

estimates_mimic <- lavPredict(fit_mimic, type = "lv")  

# anchoring estimates
SE_MIMIC <- exp(cumsum(estimates_mimic))
base <- 27.5 * SE_MIMIC[1]   # external benchmark
SE_MIMIC <- SE_MIMIC * (base / SE_MIMIC[1])

Time_vec <- as.yearqtr(gsub("K", " Q", ketv$Period), format = "%Y Q%q")
plot(Time_vec[2:76], SE_MIMIC, type = "l",
     main = "Shadow Economy Estimate")


######################################################## State-space model
data_transformed$time <- c(1:75)  # needed for specification
data_transformed$id <- 1

dynr_data <- dynr.data(data = data_transformed, time = "time",
                       observed = c("real_GDP", "cash.GDP", "real_GFCF"),
                       covariates = c("real_wage", "excise", "kaitz",
                                      "regul_quality"))


measurement <- prep.measurement(
  values.load = matrix(c(-1, NA, NA), nrow = 3, byrow = TRUE), #real_GDP fixed to -1
  params.load = matrix(c("fixed", "free1", "free2"), nrow = 3, byrow = TRUE),
  state.names = "shadow",
  obs.names = c("real_GDP", "cash.GDP", "real_GFCF")
)


dynamics <- prep.formulaDynamics(
  formula = list(
    # lag is specified writing latent variable on right side of equation
    shadow ~ a * shadow + b1 * real_wage + b2 * kaitz +   
      b3 * regul_quality
  ),
  startval = c( a = 0, b1 = 0, 
               b2 = 0, b3 = 0),
  isContinuousTime = FALSE
)


noise <- prep.noise(
  values.latent = matrix(0.1, nrow = 1, ncol = 1),
  params.latent = matrix("dnoise", nrow = 1, ncol = 1),
  values.observed = diag(0.1, 3),
  params.observed = diag(c("eps_gdp", "eps_cash", "eps_gfcf"), 3)
)

initial <- prep.initial(
  values.inistate = c(0),           
  params.inistate = c("inipos"),   
  values.inicov = matrix(1),     
  params.inicov = matrix("fixed")
)


model2 <- dynr.model(
  dynamics = dynamics,
  measurement = measurement,
  noise = noise,
  initial = initial,
  data = dynr_data,
  outfile = "state_space_model"
)


fit_SS <- dynr.cook(model2, debug_flag = T)

summary(fit_SS)


# anchoring estimates
SE_SS <- fit_SS$eta_smooth_final[1,]
SE_SS <- exp(cumsum(SE_SS))
base <- 27.5 * SE_SS[1]   # external benchmark
SE_SS <- SE_SS * (base / SE_SS[1])

Time_vec <- as.yearqtr(gsub("K", " Q", ketv$Period), format = "%Y Q%q")
plot(SE_SS, type = "l", main = "Shadow Economy Estimate")

# checking residuals
innovations_matrix <- t(fit_SS@innov_vec) 
resid_series <- innovations_matrix[, 1]      #change column no for other variables
plot(resid_series, type = 'l')

par(mfrow=c(2,1))
acf(resid_series, main="ACF of Innovations")
pacf(resid_series, main="PACF of Innovations")

ljung_box_test <- Box.test(resid_series, lag=20, type="Ljung-Box")
print(ljung_box_test)

par(mfrow=c(1,2))
hist(resid_series, breaks="FD", prob=TRUE, main="Histogram of Innovations")
qqnorm(resid_series)
qqline(resid_series, col="red")
par(mfrow=c(1,1))

jarque.bera.test(resid_series)
shapiro.test(resid_series)

######################################################## regime switching
dynr_data <- dynr.data(
  data = data_transformed, 
  time = "time",
  observed = c("real_GDP", "real_GFCF", "cash.GDP"),
  covariates = c("real_wage", "excise", "kaitz", "regul_quality")
)

measurement <- prep.measurement(
  values.load = matrix(c(-1, NA, NA), nrow = 3, byrow = TRUE),
  params.load = matrix(c("fixed", "free1", "free2"), nrow = 3, byrow = TRUE),
  state.names = "shadow",
  obs.names = c("real_GDP", "real_GFCF", "cash.GDP")
)

# initial values similar to those obtained in SS model for faster convergence
dynamics <- prep.matrixDynamics(
  values.dyn = list(
    matrix(0.7, nrow = 1, ncol = 1), 
    matrix(0.7, nrow = 1, ncol = 1)
  ),
  params.dyn = list(
    matrix("a_const", nrow = 1, ncol = 1),  
    matrix("a_const", nrow = 1, ncol = 1)  
  ),
   
  values.exo = list(
    matrix(c(-0.1, 0, -0.1), nrow = 1, ncol = 3),
    matrix(c(-0.1, 0.1, -0.1), nrow = 1, ncol = 3)  
  ),
  params.exo = list(
    #excise fixed to 0 in regime 1
    matrix(c("b1_const", "fixed", "b3_const"), nrow = 1, ncol = 3),
    matrix(c("b1_const", "b2_r2", "b3_const"), nrow = 1, ncol = 3)
  ),
  
  covariates = c("real_wage", "excise", "kaitz"),
  isContinuousTime = FALSE
)

noise <- prep.noise(
  values.latent = matrix(0.01, nrow = 1, ncol = 1), 
  params.latent = matrix("dnoise", nrow = 1, ncol = 1),
  values.observed = diag(0.1, 3),
  params.observed = diag(c("eps_gdp", "eps_cash", "eps_gfcf"), 3)
)

initial <- prep.initial(
  values.inistate = c(0),
  params.inistate = c("fixed"),
  values.inicov = matrix(1), 
  params.inicov = matrix("fixed"),
  values.regimep = c(0.7, -1),
  params.regimep = c("fixed", "fixed")
)

regimes <- prep.regimes(
  values = matrix(c(0.7, -1, 0, 0), 2, 2), 
  params = matrix(c("p11", "p21", "fixed", "fixed"), 2, 2)
)


model3 <- dynr.model(
  dynamics = dynamics, measurement = measurement, noise = noise,
  initial = initial, regimes = regimes, data = dynr_data,
  outfile = "regime_switching_model"
)

fit_regimes <- dynr.cook(model3, debug_flag = T)
summary(fit_regimes)

# checking residuals
innovations_matrix <- t(fit_regimes@innov_vec) 
resid_series <- innovations_matrix[, 1]      #change column no for other variables
plot(resid_series, type = 'l')

par(mfrow=c(2,1))
acf(resid_series, main="ACF of Innovations")
pacf(resid_series, main="PACF of Innovations")

ljung_box_test <- Box.test(resid_series, lag=20, type="Ljung-Box")
print(ljung_box_test)

par(mfrow=c(1,2))
hist(resid_series, breaks="FD", prob=TRUE, main="Histogram of Innovations")
qqnorm(resid_series)
qqline(resid_series, col="red")
par(mfrow=c(1,1))

jarque.bera.test(resid_series)
shapiro.test(resid_series)

# anchoring estimates
se_regime <- fit_regimes$eta_smooth_final[1,]
SE_regimes <- exp(cumsum(se_regime))
base <- 27.5 * SE_regimes[1]   # external benchmark
SE_regimes <- SE_regimes * (base / SE_regimes[1])


#### plots

df <- data.frame(
  Time = Time_vec[2:76],
  SE_MIMIC = SE_MIMIC,
  SE_SS = SE_SS,
  SE_regimes = SE_regimes
)

df_long <- df |>
  pivot_longer(
    cols = -Time,         
    names_to = "Method",  
    values_to = "SE"       
  )

# Convert to long format for ggplot
ggplot(df_long, aes(x = Time, y = SE, color = Method)) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = c(
      "SE_MIMIC"   = "black",
      "SE_SS"      = "blue",
      "SE_regimes" = "green"
    ),
    breaks = c("SE_MIMIC", "SE_SS", "SE_regimes"),
    labels = c("MIMIC model", "State-space model", "Regime-switching SS model")
  ) +
  labs(x = "", y = "Shadow Economy, % GDP") +
  scale_x_continuous(breaks = seq(2005, 2025, by = 5)) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85, 0.85)
  )

# Probabilities plot
plot_df <- data.frame(Time = Time_vec[2:76], Prob = fit_regimes$pr_t_given_T[1,])
ggplot(plot_df, aes(x = Time, y = Prob)) +
  geom_line(color = "#2c3e50", size = 1) +  
  geom_area(fill = "#3498db", alpha = 0.2) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    y = "Probability of being in Regime 1"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
