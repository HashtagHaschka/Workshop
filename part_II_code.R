# =============================================
# WORKSHOP: Copula-Based Endogeneity Correction
# =============================================
# https://github.com/HashtagHaschka/Copula-based-endogeneity-corrections
# =============================================


source("https://raw.githubusercontent.com/HashtagHaschka/Workshop/functions/CopReg2sCOPE-np_workshop.R")
source("https://raw.githubusercontent.com/HashtagHaschka/Workshop/functions/CopReg2sCOPE_workshop.R")
source("https://raw.githubusercontent.com/HashtagHaschka/Workshop/functions/CopRegBMW_workshop.R")
source("https://raw.githubusercontent.com/HashtagHaschka/Workshop/functions/CopRegIMA_workshop.R")
source("https://raw.githubusercontent.com/HashtagHaschka/Workshop/functions/CopRegJAMS_workshop.R")
source("https://raw.githubusercontent.com/HashtagHaschka/Workshop/functions/CopRegPG_workshop.R")
source("https://raw.githubusercontent.com/HashtagHaschka/Workshop/functions/ICAreg_workshop.R")


# ------------------- #
#  OMITTED VARIABLES  #
# ------------------- #

n <- 1000
W <- rnorm(n) # omitted variable
P <- W + (rgamma(n, shape = 1, rate = 1) - 1) # observed regressor
xi <- rnorm(n) # normal error

Y <- 3 - 1 * P + 1 * W + xi
data1 <- as.data.frame(cbind(Y, P, W))

mod1 <- lm(Y ~ P, data1) # variable W is omitted
summary(mod1) # bias because of omitted variable


# Estimate with Park-Gupta
pg_model <- CopRegPG(formula = Y ~ P,
                     data = data1, cdf = "ecdf")
pg_model[[1]] # still biased


### OUR PROPOSED ESTIMATOR
ica_model <- ica_reg (formula = Y ~ P, data = data1)
ica_model[[1]] # unbiased



# ----------------------------------- #
#  A TRUE OMITTED VARIABLES SCENARIO  #
# ----------------------------------- #

n <- 1000
W <- scale(rt(n, df = 4)) + scale(runif(n)) + 
  scale(rbinom(n, size = 1, prob = .5)) + scale(rpois(n, lambda = 1)) +
  scale(rexp(n)) + scale(rchisq(n, d = 2)) + scale(rlnorm(n)) + scale(rgamma(n, shape = 2))
W <- W/8
P <- W + rgamma(n, shape = 1, rate = 1) # observed regressor

Y <- 3 + 1 * P + 1 * W + xi
data1 <- as.data.frame(cbind(Y, P, W))
colnames(data1) <- c("Y", "P", "W")

mod1 <- lm(Y ~ P, data1) # variable W is omitted
summary(mod1) # bias because of omitted variable


# Estimate with Park-Gupta
pg_model <- CopRegPG(formula = Y ~ P,
                     data = data1, cdf = "ecdf")
pg_model[[1]] # still biased


mod1 <- ica_reg(formula = Y ~ P, data = data1)
mod1[[1]] # works



R <- 1000 
results <- matrix(data = NA, nrow = R, ncol = 3)
colnames(results) <- c("OLS", "Copula", "ICA")

for (r in 1:R) {
  
  n <- 1000
  W <- scale(rt(n, df = 4)) + scale(runif(n)) + 
    scale(rbinom(n, size = 1, prob = .5)) + scale(rpois(n, lambda = 1)) +
    scale(rexp(n)) + scale(rchisq(n, d = 2)) + scale(rlnorm(n)) + scale(rgamma(n, shape = 2))
  W <- W/8
  P <- W + rgamma(n, shape = 1, rate = 1) # observed regressor
  
  Y <- 3 - 1 * P + 1 * W + xi
  data1 <- as.data.frame(cbind(Y, P, W))
  colnames(data1) <- c("Y", "P", "W")
  
  mod1 <- lm(Y ~ P, data1) # variable W is omitted
  results[r, 1] <- coefficients(mod1)[2]
  
  pg_model <- CopRegPG(formula = Y ~ P, nboots = 1,
                       data = data1, cdf = "ecdf")
  results[r, 2] <- pg_model[[1]][2, 1]
  
  mod1 <- ica_reg(formula = Y ~ P, data = data1, nboots = 1)
  results[r, 3] <- mod1[[1]][2, 1]
  
}

mean(results[, 1]) # OLS is biased
mean(results[, 2]) # Copula overcorrects
mean(results[, 3]) # Proposed is unbiaseds



# ------------------------- #
#  REAL DATA DEMONSTRATION  #
# ------------------------- #

library(bayesm)
library(dplyr)

data("orangeJuice")
dat1 <- orangeJuice[[1]]
help(orangeJuice)

dat1_FloridaNatural <- dat1 %>% filter(brand == 3)
dat1_FloridaNatural <- dat1_FloridaNatural %>%
  mutate(across(starts_with("price"), log))


# Fit least squares
mod1 <- lm(formula = logmove ~ price3 + feat + as.factor(deal) +
             price1 + price2 + price4 + price5 + price6 + price7 +
             price8 + price9 + price10 + price11,
           data = dat1_FloridaNatural)
summary(mod1)


# 2sCOPE
Tscope_model <- CopRegPG(formula = logmove ~ price3 | feat + as.factor(deal) +
                           price1 + price2 + price4 + price5 + price6 + price7 +
                           price8 + price9 + price10 + price11,
                     data = dat1_FloridaNatural, cdf = "ecdf")
Tscope_model[[1]]


# Proposed ICA
ica_model <- ica_reg(formula = logmove ~ price3 | feat + as.factor(deal) +
                       price1 + price2 + price4 + price5 + price6 + price7 +
                       price8 + price9 + price10 + price11,
                         data = dat1_FloridaNatural)
ica_model[[1]]




