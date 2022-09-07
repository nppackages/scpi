################################################################################
## SCPI R Package
## R-file for Empirical Illustration - Multiple Treated Units
## Authors: Matias D. Cattaneo, Yingjie Feng, Filippo Palomba and Rocio Titiunik  
################################################################################
### Clear R environment
rm(list=ls(all=TRUE))

### Install R library
#install.packages('scpi')

### Load SCPI package
library(scpi)

###############################################################################
# MULTIPLE TREATED UNITS
###############################################################################

### Load data
data <- scpi_germany

# Create a second placebo treated unit
data$treatment <- 0
data[(data$country == "West Germany" & data$year >= 1991), 'treatment'] <- 1
data[(data$country == "Italy" & data$year >= 1992), 'treatment'] <- 1


###############################################
# unit-time treatment effect
###############################################

df <- scdataMulti(data, id.var = "country", outcome.var = "gdp", 
                  treatment.var = "treatment", time.var = "year", constant = TRUE, 
                  cointegrated.data = T, features = list(c("gdp","trade")),
                  cov.adj = list(c("constant", "trend")))

res <- scest(df, w.constr = list("name" = "simplex"))
scplotMulti(res)

res.pi <- scpi(df, w.constr = list("name" = "simplex"), cores = 1, sims = 50,
               e.method = "gaussian")

# plot series
scplotMulti(res.pi, type = "series")

# plot treatment
scplotMulti(res.pi, type = "series", joint = TRUE)


###############################################
# average unit post-treatment effect
###############################################

df <- scdataMulti(data, id.var = "country", outcome.var = "gdp", 
                  treatment.var = "treatment", time.var = "year", constant = T, 
                  cointegrated.data = T, features = list(c("gdp","trade")),
                  cov.adj = list(c("constant", "trend")), effect = "unit")

res <- scest(df, w.constr = list("name" = "simplex"))
scplotMulti(res)

res.pi <- scpi(df, w.constr = list("name" = "simplex"), cores = 1, sims = 50,
               e.method = "gaussian")

# plot series
scplotMulti(res.pi, type = "series")

# plot treatment
scplotMulti(res.pi, type = "series", joint = TRUE)
