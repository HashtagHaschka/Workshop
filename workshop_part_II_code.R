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


# ------------------- #
#  OMITTED VARIABLES I#
# ------------------- #


n <- 10000
W <- rnorm(n) # omitted variable
Pstar<- rnorm(n)#(rgamma(n, shape = 1, rate = 1) - 1) #exogeneous part of the variance / latent instrument
P <- qgamma(pnorm(W + Pstar),shape = 1, rate = 1) # observed regressor (omitted variable plus latent instrument)
xi <- rnorm(n) # normal error

Y <- 3 + 1 * P + 1 * W + xi
data1 <- as.data.frame(cbind(Y, P, W))

mod1 <- lm(Y ~ P, data1) # variable W is omitted
summary(mod1) # bias because of omitted variable

# Estimate with Park-Gupta
pg_model <- CopRegPG(formula = Y ~ P, data = data1,)
summary(pg_model) # still biased





# --------------------- #
#  OMITTED VARIABLES II #
# --------------------- #


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
pg_model <- CopRegPG(formula = Y ~ P, data = data1,)
summary(pg_model) # still biased


### AN EXAMPLE FOR 2sCOPE

n <- 1000
W <- rnorm(n) # omitted variable
X <- (rchisq(n, df = 1) - 1) # exogenous regressor
P <- X + W + (rgamma(n, shape = 1, rate = 1) - 1) # observed regressor
xi <- rnorm(n) # normal error

Y <- 3 + 1 * P + 1 * W + 1 * X + xi
data1 <- as.data.frame(cbind(Y, P, W, X))

mod1 <- lm(Y ~ P + X, data1) # variable W is omitted
summary(mod1) # bias because of omitted variable


# Estimate with 2sCOPE
scope_model <- CopReg2sCOPE(formula = Y ~ P | X, data = data1)
summary(scope_model) # still biased
validity(scope_model) # validity checks are passed! Crazy!

# Take home message: The Gaussian copula model cannot correct for endogeneity
# due to typical omitted variables



### OUR PROPOSED ESTIMATOR
ica_model <- IcaReg(formula = Y ~ P | X, data = data1)
ica_model[1] # unbiased. WOW!



# ----------------------------------- #
#  A TRUE OMITTED VARIABLES SCENARIO  #
# ----------------------------------- #


n <- 1000
W <- scale(rt(n, df = 4)) + scale(runif(n)) + 
  scale(rbinom(n, size = 1, prob = .5)) + scale(rpois(n, lambda = 1)) +
  scale(rexp(n)) + scale(rchisq(n, d = 2)) + scale(rlnorm(n)) + scale(rgamma(n, shape = 2))
W <- W/8 # multiple non-normal omitted variables
P <- W + rgamma(n, shape = 1, rate = 1) # observed regressor

Y <- 3 + 1 * P + 1 * W
data1 <- as.data.frame(cbind(Y, P, W))
colnames(data1) <- c("Y", "P", "W")

mod1 <- lm(Y ~ P, data1) # variable W is omitted
summary(mod1) # bias because of omitted variable


# Estimate with Park-Gupta
pg_model <- CopRegPG(formula = Y ~ P, data = data1)
summary(pg_model) # still biased

# Estimate with our proposed method
mod1 <- IcaReg(formula = Y ~ P, data = data1)
mod1[[1]] # unbiased. WOW!


### Let's scale it up!

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
  
  pg_model <- CopRegPG(formula = Y ~ P, nboots = 2, data = data1)
  results[r, 2] <- coefficients(pg_model)[2]
  
  mod1 <- IcaReg(formula = Y ~ P, data = data1, se = FALSE)
  results[r, 3] <- mod1[[1]][2, 1]
  
}

# TRUE VALUE IS -1
mean(results[, 1]) # OLS is biased
mean(results[, 2]) # Copula overcorrects
mean(results[, 3]) # Proposed is unbiaseds



################################################################################


# Comparison with IV in a scenario with a valid and strong instrument
library(data.table)

# download the dataset from GitHub
dat1 <- fread("https://raw.githubusercontent.com/HashtagHaschka/Workshop/main/dataset.csv")

# NO BOOTSTRAPPING OR JACKKNIFE TO SPEED UP EVERYTHING

# OLS - NO FE
ols_pooled <- lm(depvar ~ price + price_premium +
                   price_national + price_dominicks +
                   display + display_premium + 
                   display_national + display_dominicks + as.factor(quarter), 
                 dat1)

# OLS - FE
ols_fe <- lm(depvar ~ price + price_premium +
               price_national + price_dominicks +
               display + display_premium + 
               display_national + display_dominicks + as.factor(quarter) + 
               as.factor(store), dat1)

# IV - NO FE
preds1 <- lm(price ~ iv +
               price_dominicks + price_premium + 
               price_national + display_dominicks + display_premium + 
               display_national + display  + as.factor(quarter), 
             dat1)$fitted.values
iv_pooled <- lm(depvar ~ preds1 + 
                  price_premium + price_national + price_dominicks +
                  display + display_premium + display_national + 
                  display_dominicks + as.factor(quarter), dat1)

# IV - FE
preds2 <- lm(price ~ iv +
               price_dominicks + price_premium + price_national +
               display_dominicks + display_premium + display_national + 
               display+ as.factor(quarter) + as.factor(store), 
             dat1)$fitted.values
iv_fe <- lm(depvar ~ preds2 + 
              price_premium + price_national + price_dominicks +
              display + display_premium + display_national + 
              display_dominicks + as.factor(quarter) + as.factor(store), 
            dat1)

# 2sCOPE - NO FE
scope_pooled <- CopReg2sCOPE(formula = depvar ~ price | 
                               price_premium + price_national + 
                               price_dominicks + display + 
                               display_premium + display_national +
                               display_dominicks + as.factor(quarter), 
                             data = as.data.frame(dat1), 
                             nboots = 2)

# 2sCOPE - FE
scope_fe <- CopReg2sCOPE(formula = depvar ~ price | 
                           price_premium + price_national + 
                           price_dominicks + display + 
                           display_premium + display_national +
                           display_dominicks + as.factor(quarter) + 
                           as.factor(store), data = as.data.frame(dat1), 
                         nboots = 2)

# Panel-GC
mod_panel <- CopRegPANEL(formula = depvar ~ price | 
                           price_premium + price_national + 
                           price_dominicks + display + 
                           display_premium + display_national + 
                           display_dominicks + as.factor(quarter), 
                         data = as.data.frame(dat1), 
                         index = c("store", "week"), nboots = 2, 
                         intercept = TRUE)

ica_pooled <- IcaReg(formula = depvar ~ price | 
                       price_premium + price_national + 
                       price_dominicks + display + 
                       display_premium + display_national + 
                       display_dominicks + as.factor(quarter), 
                     as.data.frame(dat1), se = FALSE)

ica_fe <- IcaReg(formula = depvar ~ price | 
                   price_premium + price_national + 
                   price_dominicks + display + 
                   display_premium + display_national + 
                   display_dominicks + as.factor(quarter) + as.factor(store), 
                 as.data.frame(dat1), se = FALSE)


################################################################################

### Documenting the results in a Table

row_names <- c("Intercept", "log price", "log price_premium",
               "log price_national", "log price_Dominicks", "display",
               "display_premium", "display_national", "display_dominicks")

vars_std <- c("(Intercept)", "price", "price_premium",
              "price_national", "price_dominicks", "display",
              "display_premium", "display_national", "display_dominicks")

vars_iv1 <- vars_std
vars_iv1[2] <- "preds1"
vars_iv2 <- vars_std
vars_iv2[2] <- "preds2"

coef_list <- list(summary(ols_pooled)$coefficients,
                  summary(ols_fe)$coefficients,
                  summary(iv_pooled)$coefficients,
                  summary(iv_fe)$coefficients,
                  summary(scope_pooled)$coefficients,
                  summary(scope_fe)$coefficients,
                  summary(mod_panel)$coefficients,
                  ica_pooled[[1]],
                  ica_fe[[1]])

vars_list <- list(vars_std, vars_std, vars_iv1, vars_iv2, vars_std,
                  vars_std, vars_std, vars_std, vars_std)

col_names <- c("OLS_no", "OLS_yes", "IV_no", "IV_yes", "2sCOPE_no",
               "2sCOPE_yes", "PanelGC_yes", "ICA_no", "ICA_yes")

fmt <- function(x) {
  if (is.na(x)) return("")
  if (abs(x) >= 1) {
    out <- formatC(x, format = "f", digits = 3)
  } else {
    out <- formatC(x, format = "f", digits = 4)
    out <- sub("0\\.", ".", out)
  }
  out
}

tab <- matrix("", nrow = 2 * length(row_names), ncol = length(col_names))
rownames(tab) <- as.vector(rbind(row_names, paste0(row_names, "_se")))
colnames(tab) <- col_names

for (j in 1:length(coef_list)) {
  cf <- coef_list[[j]]
  vn <- vars_list[[j]]
  for (i in 1:length(row_names)) {
    if (vn[i] %in% rownames(cf)) {
      tab[2 * i - 1, j] <- fmt(cf[vn[i], "Estimate"])
      tab[2 * i,     j] <- paste0("(", fmt(cf[vn[i], "Std. Error"]), ")")
    }
  }
}

print(tab, quote = FALSE)


################################################################################