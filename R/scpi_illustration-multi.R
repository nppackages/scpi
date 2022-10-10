################################################################################
## SCPI R Package
## R-file for Empirical Illustration - Multiple Treated Units
## Authors: Matias D. Cattaneo, Yingjie Feng, Filippo Palomba and Rocio Titiunik
################################################################################
### Clear R environment
rm(list = ls(all = TRUE))

### Install R library
#install.packages("scpi")

### Load SCPI package
library(scpi)

###############################################################################
# MULTIPLE TREATED UNITS
###############################################################################

### Load data
data <- scpi_germany

# Create a second placebo treated unit
data$treatment <- 0
data[(data$country == "West Germany" & data$year >= 1991), "treatment"] <- 1
data[(data$country == "Italy" & data$year >= 1992), "treatment"] <- 1


######################################################
# unit-time treatment effect (\tau_{ik})
######################################################

df <- scdataMulti(data, id.var = "country", outcome.var = "gdp",
                  treatment.var = "treatment", time.var = "year", constant = TRUE,
                  cointegrated.data = TRUE, features = list(c("gdp", "trade")),
                  cov.adj = list(c("constant", "trend")))

res <- scest(df, w.constr = list("name" = "simplex"))
scplotMulti(res)

respi <- scpi(df, w.constr = list("name" = "simplex"), cores = 1, sims = 50,
               e.method = "gaussian")

# plot series
scplotMulti(respi, type = "series")

# plot treatment
scplotMulti(respi, type = "series", joint = TRUE)


######################################################
# average unit treatment effect (\tau_{i.})
######################################################

df <- scdataMulti(data, id.var = "country", outcome.var = "gdp",
                  treatment.var = "treatment", time.var = "year", constant = TRUE,
                  cointegrated.data = TRUE, features = list(c("gdp", "trade")),
                  cov.adj = list(c("constant", "trend")), effect = "unit")

res <- scest(df, w.constr = list("name" = "simplex"))
scplotMulti(res)

respi <- scpi(df, w.constr = list("name" = "simplex"), cores = 1, sims = 50,
               e.method = "gaussian")

# plot series
scplotMulti(respi, type = "series")

# plot treatment
scplotMulti(respi, type = "series", joint = TRUE)

######################################################
# average treatment effect on the treated (\tau_{.k})
######################################################

df <- scdataMulti(data, id.var = "country", outcome.var = "gdp",
                  treatment.var = "treatment", time.var = "year", constant = TRUE,
                  cointegrated.data = TRUE, features = list(c("gdp", "trade")),
                  cov.adj = list(c("constant", "trend")), effect = "time")

res <- scest(df, w.constr = list("name" = "simplex"))
scplotMulti(res)

respi <- scpi(df, w.constr = list("name" = "simplex"), cores = 1, sims = 50,
               e.method = "gaussian")

# plot series
scplotMulti(respi, type = "series")

# plot treatment
scplotMulti(respi, type = "series", joint = TRUE)
