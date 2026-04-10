# ------------------------------------------------------------------------------
# Linear & Generalized Linear Models for the LUCAS Dataset
# ------------------------------------------------------------------------------

rm(list=ls())

library(lme4)
library(lmtest)

# ------------------------------------------------------------------------------
# load data
source("load_lucas_data.R")

# ------------------------------------------------------------------------------
# non-linearity test
resettest(y~X+U)

# ------------------------------------------------------------------------------
# Linear model
# ------------------------------------------------------------------------------

# Fit model
linfit <- lm(y~X+U)

# Validation metrics
res_lm <- linfit$residuals
RSS_lm <- sum(res_lm^2)
df_lm <- linfit$df.residual
AIC_lm <- 2*n*log(RSS_lm)+2*df_lm

cat(
  "LM Validation | RSS:", RSS_lm, 
  "DF:", df_lm, 
  "AIC:", AIC_lm, 
  "Min fitted:", min(linfit$fitted.values), 
  "\n"
)

# Diagnostic plots
qqnorm(res_lm, main="y = x+u")
qqline(res_lm)
plot(
  linfit$fitted.values,
  res_lm,
  xlab="Fitted Values",
  ylab="Residuals",
  main="y = x+u"
)

# ------------------------------------------------------------------------------
# Generalized linear model
# ------------------------------------------------------------------------------

# Fit model
y_pos <- y
y_pos[y_pos==0] <- 1e-10 
expfit <- glm(y_pos~X+U, family=gaussian(link="log"))

# Validation metrics
res_glm <- expfit$residuals
RSS_glm <- sum(res_glm^2)
df_glm <- expfit$df.residual
AIC_glm <- 2*n*log(RSS_glm)+2*df_glm

cat(
  "GLM Validation | RSS:", RSS_glm, 
  "DF:", df_glm, 
  "AIC:", AIC_glm, 
  "Min fitted:", min(expfit$fitted.values), 
  "\n"
)

# Diagnostic plots
qqnorm(res_glm, main="y = exp(x+u)")
qqline(res_glm)
plot(
  expfit$fitted.values,
  res_glm,
  xlab="Fitted Values",
  ylab="Residuals",
  main="y = exp(x+u)"
)
