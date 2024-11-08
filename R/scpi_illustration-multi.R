################################################################################
## SCPI R Package
## R-file for Empirical Illustration - Multiple Treated Units
## Authors: Matias D. Cattaneo, Yingjie Feng, Filippo Palomba, and Rocio Titiunik
################################################################################
### Clear R environment
rm(list = ls(all = TRUE))

### Install R library
#install.packages("scpi")

### Load SCPI package
library(scpi)

set.seed(8894)

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
ggsave("germany_est_multi_1.png", height=6, width=9, dpi="retina")


respi <- scpi(df, w.constr = list("name" = "simplex"), cores = 1, sims = 50,
               e.method = "gaussian")

# plot series
scplotMulti(respi, type = "series")
ggsave("germany_est_multi_2.png", height=6, width=9, dpi="retina")

# plot treatment
scplotMulti(respi, type = "series", joint = TRUE)
ggsave("germany_est_multi_3.png", height=6, width=9, dpi="retina")


######################################################
# average unit treatment effect (\tau_{i.})
######################################################

df <- scdataMulti(data, id.var = "country", outcome.var = "gdp",
                  treatment.var = "treatment", time.var = "year", constant = TRUE,
                  cointegrated.data = TRUE, features = list(c("gdp", "trade")),
                  cov.adj = list(c("constant", "trend")), effect = "unit")

res <- scest(df, w.constr = list("name" = "simplex"))
scplotMulti(res)
ggsave("germany_est_multi_4.png", height=6, width=9, dpi="retina")

respi <- scpi(df, w.constr = list("name" = "simplex"), cores = 1, sims = 50,
               e.method = "gaussian")

# plot series
scplotMulti(respi, type = "series")
ggsave("germany_est_multi_5.png", height=6, width=9, dpi="retina")

# plot treatment
scplotMulti(respi, type = "series", joint = TRUE)
ggsave("germany_est_multi_6.png", height=6, width=9, dpi="retina")

######################################################
# average treatment effect on the treated (\tau_{.k})
######################################################

df <- scdataMulti(data, id.var = "country", outcome.var = "gdp",
                  treatment.var = "treatment", time.var = "year", constant = TRUE,
                  cointegrated.data = TRUE, features = list(c("gdp", "trade")),
                  cov.adj = list(c("constant", "trend")), effect = "time")

res <- scest(df, w.constr = list("name" = "simplex"))
scplotMulti(res)
ggsave("germany_est_multi_7.png", height=6, width=9, dpi="retina")

respi <- scpi(df, w.constr = list("name" = "simplex"), cores = 1, sims = 50,
               e.method = "gaussian")

# plot series
scplotMulti(respi, type = "series")
ggsave("germany_est_multi_8.png", height=6, width=9, dpi="retina")

# plot treatment
scplotMulti(respi, type = "series", joint = TRUE)
ggsave("germany_est_multi_9.png", height=6, width=9, dpi="retina")
