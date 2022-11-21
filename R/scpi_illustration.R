################################################################################
## SCPI R Package
## R-file for Empirical Illustration - Single Treated Unit
## Authors: Matias D. Cattaneo, Yingjie Feng, Filippo Palomba and Rocio Titiunik  
################################################################################
### Clear R environment
rm(list=ls(all=TRUE))

### Install R library
#install.packages('scpi')

### Load SCPI package
library(scpi)
library(ggplot2)
library(latex2exp)

###############################################################################
# SINGLE TREATED UNIT
###############################################################################

### Load data
data <- scpi_germany

####################################
### Set options for data preparation
id.var      <- "country"                              # ID variable
time.var    <- "year"                                 # Time variable
period.pre  <- seq(from = 1960, to = 1990, by = 1)    # Pre-treatment period
period.post <- (1991:2003)                            # Post-treatment period
unit.tr     <- "West Germany"                         # Treated unit (in terms of id.var)
unit.co     <- setdiff(unique(data$country), unit.tr) # Donors pool
outcome.var <- "gdp"                                  # Outcome variable
cov.adj     <- NULL                                   # Covariates for adjustment
features    <- NULL                                   # No features other than outcome
constant    <- FALSE                                  # No constant term
report.missing <- FALSE                               # To check where missing values are
cointegrated.data <- TRUE                             # Belief that the data are cointegrated


####################################
### Data preparation
df  <-   scdata(df = data, id.var = id.var, time.var = time.var, outcome.var = outcome.var,
                period.pre = period.pre, period.post = period.post,
                unit.tr = unit.tr, unit.co = unit.co, cov.adj = cov.adj, features = features,
                constant = constant, cointegrated.data = cointegrated.data)

####################################
### SC - point estimation with simplex
est.si  <- scest(data = df, w.constr = list(name="simplex"))
# Use print or summary methods to check results
print(est.si)
summary(est.si)
est.si2 <- scest(data = df, w.constr = list(p = "L1", dir = "==", Q = 1, lb = 0))
summary(est.si2)

####################################
### SC - plot results
scplot(result = est.si, fig.path = ".",
       fig.name = "germany_est", fig.format = "png", plot.range = (1960:2003),
       label.xy = list(x.lab = "Year", y.lab = "GDP per capita (thousand US dollars)"),
       event.label = list(lab = "Reunification", height = 10))


####################################
### SC - point estimation with lasso
est.lasso <- scest(data = df, w.constr = list(name="lasso"))
summary(est.lasso)
est.lasso2 <- scest(data = df, w.constr = list(p = "L1", dir = "<=", Q = 1, lb = -Inf))
summary(est.lasso2)


####################################
### SC - point estimation with ridge
est.ridge <- scest(data = df, w.constr = list(name="ridge"))
summary(est.ridge)
Qest <- est.ridge$est.results$w.constr$Q
est.ridge2 <- scest(data = df, w.constr = list(p = "L2", dir = "<=", Q = Qest, lb = -Inf))
summary(est.ridge2)


####################################
### SC - point estimation with L1-L2
est.l1l2 <- scest(data = df, w.constr = list(name="L1-L2"))
summary(est.l1l2)
est.l1l2.2 <- scest(data = df, w.constr = list(p = "L1-L2", dir = "==/<=", Q = 1, 
                                             Q2 = Qest, lb = 0))
summary(est.l1l2.2)


####################################
### SC - point estimation with least squares
est.ls <- scest(data = df, w.constr = list(name="ols"))
summary(est.ls)
est.ls2 <- scest(data = df, w.constr = list(p = "no norm", dir = NULL, Q = NULL, lb = -Inf))
summary(est.ls2)


####################################
## Set options for inference
u.alpha  <- 0.05                         # Confidence level (in-sample uncertainty)
e.alpha  <- 0.05                         # Confidence level (out-of-sample uncertainty)
rho      <- NULL                         # Regularization parameter (if NULL it is estimated)
rho.max  <- 1                            # Maximum value attainable by rho
sims     <- 200                          # Number of simulations
u.order  <- 1                            # Degree of polynomial in B and C when modelling u
u.lags   <- 0                            # Lags of B to be used when modelling u
u.sigma  <- "HC1"                        # Estimator for the variance-covariance of u
u.missp  <- T                            # If TRUE then the model is treated as misspecified
e.lags   <- 0                            # Degree of polynomial in B and C when modelling e
e.order  <- 1                            # Lags of B to be used when modelling e
e.method <- "gaussian"                   # Estimation method for out-of-sample uncertainty
cores    <- 1                            # Number of cores to be used by scpi
w.constr <- list(name = "simplex")       # Simplex-type constraint set

set.seed(8894)
pi.si   <- scpi(data = df,u.order = u.order, u.lags = u.lags, u.sigma = u.sigma, 
                u.missp = u.missp, sims = sims, e.order = e.order, e.lags = e.lags,
                e.method = e.method, cores = cores, w.constr = w.constr, u.alpha = u.alpha,
                e.alpha = e.alpha, rho = rho, rho.max = rho.max) 

# Use print or summary methods to check results
print(pi.si)
summary(pi.si)

####################################
### SC - plot results
scplot(result = pi.si, fig.path = ".",
       fig.name = "germany_unc", fig.format = "png", plot.range = (1960:2003),
       label.xy = list(x.lab = "Year", y.lab = "GDP per capita (thousand US dollars)"),
       x.ticks = NULL, e.out = T, event.label = list(lab = "Reunification", height = 10))



########################################################
# Sensitivity Analysis for 1997 using subgaussian bounds
########################################################
set.seed(8894)
res.si  <- scpi(data = df, sims = sims, e.method = "gaussian", e.order = e.order, e.lags = e.lags,
                u.order = u.order, u.lags = u.lags, u.sigma = u.sigma, u.missp = u.missp,
                cores = cores, w.constr = list(name = "simplex"), lgapp = "linear") 

e.alpha <- 0.05  # default level in scpi
sens <- c(0.25, 0.5, 1, 1.5, 2)
time <- c(1997)
emean <- res.si$inference.results$e.mean
esig <- sqrt(res.si$inference.results$e.var)
sc.l.0 <- res.si$inference.results$CI.in.sample[,1,drop = F]
sc.r.0 <- res.si$inference.results$CI.in.sample[,2,drop = F]
y <- res.si$data$Y.post

for (l in 1:length(time)) {
  ssc.l.1 <- ssc.r.1 <- c()
  e.mean <- emean[time[l]-1990]
  sig <- esig[time[l]-1990]
  sig.seq <- sens*sig
  
  for (s in 1:length(sig.seq)) {
    eps  <- sqrt(-log(e.alpha/2)*2*(sig.seq[s]^2))
    ssc.l.1[s] <- sc.l.0[time[l]-1990] + e.mean - eps
    ssc.r.1[s] <- sc.r.0[time[l]-1990] + e.mean + eps
  }
  
  sen.dat <- data.frame(t=c(1:5), lb1=ssc.l.1, ub1=ssc.r.1,
                        lb=rep(sc.l.0[time[l]-1990], 5),
                        ub=rep(sc.r.0[time[l]-1990], 5),
                        lab=as.factor(sens))
  plot <- ggplot() + theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x="sd. of e", y="GDP per capita (thousand US dollars)")
  plot <- plot + geom_errorbar(data=sen.dat, aes(x=lab, ymin=lb1, ymax=ub1),
                               col="maroon", width=0.2, linetype=5) +
    geom_errorbar(data=sen.dat, aes(x=lab, ymin=lb, ymax=ub),
                  col="blue", width=0.2, linetype=1) +
    geom_hline(yintercept = y[time[l]-1990], linetype=1, linewidth=0.3, alpha=0.8) +
    annotate("text", x=5.4, y=y[time[l]-1990]-.1,label="Y(1)", size=3.5) +
    scale_x_discrete(labels=c(parse(text=TeX("$0.25\\hat{\\sigma}$")),
                              parse(text=TeX("$0.5\\hat{\\sigma}$")),
                              parse(text=TeX("$\\hat{\\sigma}$")),
                              parse(text=TeX("$1.5\\hat{\\sigma}$")),
                              parse(text=TeX("$2\\hat{\\sigma}$"))))
}



################################################################################
### Other examples of data preparation

## multiple features
df  <-   scdata(df = data, id.var = id.var, time.var = time.var, outcome.var = outcome.var,
                period.pre = period.pre, period.post = period.post,
                unit.tr = unit.tr, unit.co = unit.co, cov.adj = cov.adj, features = c("gdp", "trade"),
                constant = constant, cointegrated.data = cointegrated.data)

## multiple features and featuer-specific covariate adjustment
df  <-   scdata(df = data, id.var = id.var, time.var = time.var, outcome.var = outcome.var,
                period.pre = period.pre, period.post = period.post,
                unit.tr = unit.tr, unit.co = unit.co, features = c("gdp", "trade"), 
                cov.adj = list('gdp' = c("constant","trend"), 'trade' = c("constant")),
                constant = constant, cointegrated.data = cointegrated.data)
