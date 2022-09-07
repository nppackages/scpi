################################################################################
## SCPI R Package
## R-file for Empirical Illustration - Multiple Treated Units
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

###############################################################################
# MULTIPLE TREATED UNITS
###############################################################################

### Load data
data <- scpi_germany

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

res.pi <- scpi(df, w.constr = list("name" = "simplex"), cores = 1, sims = 50,
               e.method = "gaussian")

# plot series
scplotMulti(res.pi, type = "series", joint = TRUE, save.data = '__scpi_data')

load('__scpi_data.RData')

plot <- ggplot(toplot) + xlab("Date") + ylab("Outcome") +
  geom_line(aes(x=Time, y=Y, colour=Type)) + 
  geom_point(aes(x=Time, y=Y, colour=Type), size=1.5) + 
  geom_vline(aes(xintercept=Tdate)) +
  facet_wrap(~ID, ncol = 2) + theme(legend.position="bottom") +
  scale_color_manual(name = "", values = c("black", "blue"),
                     labels = c("Treated", "Synthetic Control"))

plot.w1 <- plot + geom_errorbar(data = toplot,
                                aes(x = Time, ymin = lb.gaussian, ymax = ub.gaussian), colour = "blue",
                                width = 0.5, linetype = 1) + ggtitle("In and Out of Sample Uncertainty - Subgaussian Bounds")

plotdf <- subset(toplot, Type == "Synthetic")
plot.w1 + geom_ribbon(data=plotdf, aes(x=Time, ymin=lb.joint, ymax=ub.joint), fill="blue", alpha=0.1)

