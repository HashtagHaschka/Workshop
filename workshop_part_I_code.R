## ================================================================================== ##
##  Copula-based endogeneity corrections & ICA-based correction for omitted variables ##
## ================================================================================== ##
##
## -----------------------------------------------------------------------------
##  Copula-based endogeneity corrections in R
##  https://github.com/HashtagHaschka/Copula-based-endogeneity-corrections
##
##  ICA-based endogeneity correction in R
##  https://github.com/HashtagHaschka/ICA-based-endogeneity-correction
##
##  Copyright (C) 2026 Rouven E. Haschka
##  ORCID: https://orcid.org/0000-0002-2916-9745
##
##  If this code contributes to work you publish, please cite the software
##
##    Haschka, R. E. (2026). Copula-based endogeneity corrections in R.
##    https://github.com/HashtagHaschka/Copula-based-endogeneity-corrections
##
##    Haschka, R. E. (2026). ICA-based endogeneity correction in R.
##    https://github.com/HashtagHaschka/ICA-based-endogeneity-correction
##
##  and, for the estimators it implements, the papers named in each file.
##
##  This program is free software: you can redistribute it and/or modify it
##  under the terms of the GNU General Public License as published by the Free
##  Software Foundation, either version 3 of the License, or (at your option)
##  any later version.
##
##  This program is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
##  more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program.  If not, see <https://www.gnu.org/licenses/>.
## -----------------------------------------------------------------------------


## Load the Gaussian copula-based endogeneity corrections

gh <- paste0("https://raw.githubusercontent.com/HashtagHaschka/",
             "Copula-based-endogeneity-corrections/functions/")

for (f in c("Copreg_core.R", # the core has to come first
            "Copreg_pg.R", "Copreg_2scope.R", "Copreg_ima.R",
            "Copreg_jams.R", "Copreg_bmw.R", "Copreg_2scope_np.R",
            "Copreg_panel.R", "Copreg_bayes.R"))
  source(paste0(gh, f))

need <- c("Formula", "np", "bayesm")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0L) install.packages(miss)


################################################################################


## load the ICA-based endogeneity correction 

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

source(paste0("https://raw.githubusercontent.com/HashtagHaschka/",
              "ICA-based-endogeneity-correction/main/IcaReg.R"))

pacman::p_load(dplyr, AER, ISLR, mlbench)


################################################################################


# ------------------------------------ #
#  A SIMPLE ENDOGENEITY DEMONSTRATION  #
# ------------------------------------ #

n <- 1000
rho <- 0.8  # Endogeneity strength

# Generate endogenous data (DGP from Slide 4)
latent <- mvtnorm::rmvnorm(n, mean = c(0, 0), 
                           sigma = matrix(c(1, rho, rho, 1), ncol = 2))
P <- qlnorm(pnorm(latent[, 1]))  # Endogenous regressor
xi <- latent[, 2]                # Error term
X <- rnorm(n)                    # Exogenous control
Y <- 2 + 1.5*X + 3*P + xi        # True model
data1 <- as.data.frame(cbind(Y, P, X))

# Biased OLS estimation
ols_model <- lm(Y ~ X + P, data1)
summary(ols_model)  
# Coefficient of P is biased (true is 3)
# Intercept is biased (true is 2)
# Coefficient of X is not biased (true is 1.5)



# ------------------------------ #
#  PARK-GUPTA (2012) CORRECTION  #
# ------------------------------ #


# Estimate with Park-Gupta
pg_model <- CopRegPG(formula = Y ~ P | X, data = data1)
summary(pg_model)
# not biased

# validity and diagnostic check
validity(pg_model)



# -------------------------------------------------- #
#  HASCHKA (2024) and YANG ET AL. (2025) CORRECTION  #
# -------------------------------------------------- #


# Simulate data with correlated regressors
n <- 1000
rho <- 0.8  # Endogeneity strength
r <- 0.5  # P-X correlation
Sigma <- matrix(c(1, r, rho,
                  r, 1, 0,
                  rho, 0, 1), ncol=3)
latent <- mvtnorm::rmvnorm(n, sigma = Sigma)
P <- qexp(pnorm(latent[,1]), rate = 0.5)
X <- latent[, 2]
xi <- latent[, 3]
Y <- 2 + 1.5*X + 3*P + xi
data1 <- as.data.frame(cbind(Y, P, X))

# Biased OLS estimation
ols_model <- lm(Y ~ X + P, data1)
summary(ols_model)  # Coefficients of both X and P are biased

# Estimate with Haschka (2024) method
ima_model <- CopRegIMA(formula = Y ~ P | X, data = data1)
summary(ima_model)

# Estimate with 2sCOPE method
twoscope_model <- CopReg2sCOPE(formula = Y ~ P | X, data = data1)
summary(twoscope_model)

# validity checks and diagnostics
validity(ima_model)
validity(twoscope_model)



# ------------------------------------------ #
#  HASCHKA (2022) CORRECTION FOR PANEL DATA  #
# ------------------------------------------ #


# Panel dimensions
N <- 200   # Number of cross-sectional units
Ti <- 5    # Number of time periods
n <- N * Ti

# Simulate data with correlated regressors
rho <- 0.8  # Endogeneity strength
r <- 0.5    # P-X correlation

Sigma <- matrix(c(1, r, rho,
                  r, 1, 0,
                  rho, 0, 1), ncol = 3)

latent <- mvtnorm::rmvnorm(n, sigma = Sigma)

P <- qexp(pnorm(latent[, 1]), rate = 0.5)
X <- latent[, 2]
xi <- latent[, 3]

# Panel identifiers
id <- rep(1:N, each = Ti)
time <- rep(1:Ti, times = N)

# Individual fixed effects
alpha_i <- rnorm(N)
FE <- alpha_i[id]

# Outcome variable
Y <- 2 + FE + 1.5 * X + 3 * P + xi

# Panel dataset
data1 <- data.frame(
  id = id,
  time = time,
  Y = Y,
  P = P,
  X = X
)

# Biased OLS estimation with fixed effects
ols_model <- lm(Y ~ X + P + as.factor(id), data = data1)
summary(ols_model)

# Estimate with Haschka (2022) panel method
panel_model <- CopRegPANEL(Y ~ P | X, data = data1, index = c("id", "time"))
summary(panel_model)

# Validity checks and diagnostics
validity(panel_model)



# ------------------------------------------------- #
#  LIENGAARD ET AL. (2025) - CATEGORICAL MODERATOR  #
# ------------------------------------------------- #


# Simulate data with endogeneity varying by group
n <- 1000
rho1 <- 0.4  # Endogeneity strength
rho2 <- 0.8  # Endogeneity strength
X <- sample(0:1, n, replace = TRUE)
P <- numeric(n)
xi <- numeric(n)

latent1 <- mvtnorm::rmvnorm(n/2, mean = c(0, 0), 
                            sigma = matrix(c(1, rho1, rho1, 1), ncol = 2))
latent2 <- mvtnorm::rmvnorm(n/2, mean = c(0, 0), 
                            sigma = matrix(c(1, rho2, rho2, 1), ncol = 2))

P <- c(qexp(pnorm(latent1[, 1]), rate = 0.5), 
       qunif(pnorm(latent2[, 1])))
X <- c(rep(0, n/2),
       rep(1, n/2))
xi <- c(latent1[, 2], 
        latent2[, 2])
Y <- 2 + 1.5*X + 3*P + xi
data1 <- as.data.frame(cbind(Y, P, X))

# Biased OLS estimation
ols_model <- lm(Y ~ X + P, data1)
summary(ols_model)  # Coefficients of both X and P are biased

# Estimate with Park-Gupta
pg_model <- CopRegPG(formula = Y ~ P | X, data = data1)
summary(pg_model) # Coefficients are still biased

# Estimate with Liengaard et al.
jams_model <- CopRegJAMS(formula = Y ~ P | as.factor(X), data = data1)
summary(jams_model) # unbiased



# ------------------------------------------- #
#  REAL DATA DEMONSTRATION - TRY IT YOURSELF  #
# ------------------------------------------- #


library(bayesm)
library(dplyr)

data("orangeJuice")
dat1 <- orangeJuice[[1]]
help(orangeJuice)

dat1_FloridaNatural <- dat1 %>% filter(brand == 3)
dat1_FloridaNatural <- dat1_FloridaNatural %>%
  mutate(across(starts_with("price"), log))


# Fit least squares
mod1 <- lm(formula = logmove ~ price3 + feat + as.factor(deal),
           data = dat1_FloridaNatural)
summary(mod1)

# Fit Park & Gupta (2012)
mod_PG <- CopRegPG(...)
summary(mod_PG)

# Fit Haschka (2024)
mod_IMA <- CopRegIMA(...)
summary(mod_IMA)

# Fit 2sCOPE
mod_2sCOPE <- CopReg2sCOPE(...)
summary(mod_2sCOPE)

# Fit Liengaard et al. (2025)
mod_JAMS <- CopRegJAMS(...)
summary(mod_JAMS)


# ----------------------------------- #
#  OMITTED VARIABLES I - TRY YOURSELF #
# ----------------------------------- #


n <- 10000
W <- rnorm(n) # omitted variable
Pstar<- rnorm(n)#(rgamma(n, shape = 1, rate = 1) - 1) #exogeneous part of the variance / latent instrument
P <- qgamma(pnorm(W + Pstar),shape = 1, rate = 1) # observed regressor (omitted variable plus latent instrument), both transformed to non-normal
xi <- rnorm(n) # normal error

Y <- 3 + 1 * P + 1 * W + xi
data1 <- as.data.frame(cbind(Y, P, W))

mod1 <- lm(Y ~ P, data1) # variable W is omitted
summary(mod1) # bias because of omitted variable

# Estimate with Park-Gupta
pg_model <- CopRegPG(...)
summary(pg_model) # Corrected?





# ------------------------------------ #
#  OMITTED VARIABLES II - TRY YOURSELF #
# ------------------------------------ #


n <- 10000
W <- rnorm(n) # omitted variable
Pstar<- (rgamma(n, shape = 1, rate = 1) - 1) #exogeneous part of the variance / latent instrument
P <- W + Pstar # observed regressor (omitted variable plus latent instrument)
xi <- rnorm(n) # normal error

Y <- 3 + 1 * P + 1 * W + xi
data1 <- as.data.frame(cbind(Y, P, W))

mod1 <- lm(Y ~ P, data1) # variable W is omitted
summary(mod1) # bias because of omitted variable

# Estimate with Park-Gupta
pg_model <- CopRegPG(...)
summary(pg_model) # Corrected ?


