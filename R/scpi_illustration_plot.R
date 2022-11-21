################################################################################
## SCPI R Package
## R-file for Empirical Illustration - Single Treated Unit
## Authors: Matias D. Cattaneo, Yingjie Feng, Filippo Palomba and Rocio Titiunik  
################################################################################
### Clear R environment
rm(list=ls(all=TRUE))

### Install R library
#install.packages('scpi')

### Load packages
library(scpi)
library(ggplot2)

theme_set(theme_bw())

##############################################################################
# SINGLE TREATED UNIT
##############################################################################

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
## Set options for inference
u.alpha  <- 0.05                         # Confidence level (in-sample uncertainty)
e.alpha  <- 0.05                         # Confidence level (out-of-sample uncertainty)
rho      <- NULL                         # Regularization parameter (if NULL it is estimated)
rho.max  <- 1                            # Maximum value attainable by rho
sims     <- 200                          # Number of simulations
V        <- NULL                         # Weighting matrix (if NULL it is the identity matrix)
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
result  <- scpi(data = df,u.order = u.order, u.lags = u.lags, u.sigma = u.sigma, 
                u.missp = u.missp, sims = sims, e.order = e.order, e.lags = e.lags,
                e.method = e.method, cores = cores, w.constr = w.constr, u.alpha = u.alpha,
                e.alpha = e.alpha, rho = rho, rho.max = rho.max) 

####################################
### SC - plot results
scplot(result = result, fig.path = ".",
       fig.name = "germany_unc", fig.format = "png", plot.range = (1960:2003),
       label.xy = list(x.lab = "Year", y.lab = "GDP per capita (thousand US dollars)"),
       x.ticks = NULL, e.out = T, event.label = list(lab = "Reunification", height = 10))


####################################
### SC - manually reproduce plot
# Store data on treated unit, synthetic unit, and prediction bars
y.fit <- rbind(result$est.results$Y.pre.fit, result$est.results$Y.post.fit)
y.act <- rbind(result$data$Y.pre, result$data$Y.post)

sc.l  <- result$inference.results$CI.all.gaussian[, 1, drop = FALSE]
sc.r  <- result$inference.results$CI.all.gaussian[, 2, drop = FALSE]

# Store other specifics
period.pre  <- result$data$specs$period.pre
period.post <- result$data$specs$period.post
T0          <- period.pre[length(period.pre)] # intercept
plot.range  <- c(period.pre, period.post)

# Actual data
dat    <- data.frame(t     = c(period.pre, period.post),
                     Y.act = c(y.act),
                     sname = "Treated")

# Fill with NAs Y.fit and confidence bounds where missing
Y.fit.na  <- matrix(NA, nrow = length(c(period.pre, period.post)))
sc.l.na   <- matrix(NA, nrow = length(c(period.pre, period.post)))
sc.r.na   <- matrix(NA, nrow = length(c(period.pre, period.post)))

names <- strsplit(rownames(y.fit), "\\.")
not.missing.plot <- c(period.pre,period.post) %in% unlist(lapply(names, "[[", 2))
names <- strsplit(rownames(sc.l), "\\.")
not.missing.ci   <- c(period.pre,period.post) %in% unlist(lapply(names, "[[", 2))

Y.fit.na[not.missing.plot, 1] <- y.fit
sc.l.na[not.missing.ci, 1]    <- sc.l
sc.r.na[not.missing.ci, 1]    <- sc.r

# Synthetic unit data
dat.sc <- data.frame(t        = c(period.pre, period.post),
                     Y.sc     = Y.fit.na,
                     lb       = c(sc.l.na), ub = c(sc.r.na),
                     sname    = "SC Unit")

# Set ticks, event label and merge
x.ticks <- c(seq(plot.range[1], plot.range[length(plot.range)], length.out = 5), T0)
x.ticks <- round(unique(x.ticks))


event.lab <- paste("\n", "Reunification", sep = "")
event.lab.height <- 10

dat.plot    <- subset(dat,    t %in% plot.range)
dat.sc.plot <- subset(dat.sc, t %in% plot.range)

plotdf <- dplyr::left_join(dat.plot, dat.sc.plot, by = 't')

## Plot specs
plot <- ggplot() + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  labs(x = "Year", y = "GDP per capita (thousand US dollars)") +
  theme(legend.position = "bottom", legend.box = "horizontal", legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = "black"))

## Add Series to plot
plot <- plot + 
  geom_line( data = plotdf, aes(x = t, y = Y.act, colour = sname.x), linetype = 'solid') +
  geom_point(data = plotdf, aes(x = t, y = Y.act, colour = sname.x), shape = 1) +
  geom_line( data = plotdf, aes(x = t, y = Y.sc,  colour = sname.y), linetype = 'dashed') +
  geom_point(data = plotdf, aes(x = t, y = Y.sc,  colour = sname.y), shape = 19) +
  geom_vline(xintercept = T0, linetype = "dashed") +
  geom_text(aes(x = T0, label = event.lab, y = event.lab.height), angle = 90, size = 4) +
  scale_x_continuous(breaks = x.ticks) + 
  scale_color_manual(name = "", values = c("mediumblue", "grey46"),
                     labels = c("Synthetic Control", "Treated"),
                     guide = guide_legend(override.aes = list(
                       linetype = c('dashed','solid'), shape = c(19, 1))))

## Add confidence bars and plot
plot + geom_errorbar(data = plotdf,
                    aes(x = t, ymin = lb, ymax = ub, colour = sname.y),
                    width = 0.5, linetype = 1) + ggtitle("In and Out of Sample Uncertainty")
