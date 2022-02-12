###############################################################################
# SCPI Python Package
# Script for Empirical Illustration
# Authors: Matias D. Cattaneo, Yingjie Feng, Filippo Palomba and Rocio Titiunik
###############################################################################

########################################
# Load SCPI_PKG package
import pandas
import numpy
import random
import os
from plotnine import ggsave
from warnings import filterwarnings
from scpi_pkg.scdata import scdata
from scpi_pkg.scest import scest
from scpi_pkg.scpi import scpi
from scpi_pkg.scplot import scplot

os.chdir('YOUR_PATH_HERE')
filterwarnings("ignore")
########################################
# Load database
data = pandas.read_csv("scpi_germany.csv")

########################################
# Set options for data preparation
id_var = 'country'
outcome_var = 'gdp'
time_var = 'year'
features = None
cov_adj = None
period_pre = numpy.arange(1960, 1991)
period_post = numpy.arange(1991, 2004)
unit_tr = 'West Germany'
unit_co = list(set(data[id_var].to_list()))
unit_co = [cou for cou in unit_co if cou != 'West Germany']
constant = True
report_missing = False
cointegrated_data = True

data_prep = scdata(df=data, id_var=id_var, time_var=time_var,
                   outcome_var=outcome_var, period_pre=period_pre,
                   period_post=period_post, unit_tr=unit_tr,
                   unit_co=unit_co, features=features, cov_adj=cov_adj,
                   cointegrated_data=cointegrated_data, constant=constant,
                   report_missing=report_missing)


####################################
# SC - point estimation with simplex
est_si = scest(data_prep, w_constr={'name': "simplex"})
print(est_si)
est_si2 = scest(data_prep, w_constr={'p': 'L1', 'dir': '==', 'Q': 1, 'lb': 0})
print(est_si2)


####################################
# SC - plot results
plot = scplot(est_si)
ggsave(filename='germany_est.png', plot=plot)


####################################
# SC - point estimation with lasso
est_lasso = scest(data_prep, w_constr={'name': "lasso"})
print(est_lasso)
est_lasso2 = scest(data_prep, w_constr={'p': 'L1', 'dir': '<=', 'Q': 1, 'lb': -numpy.inf})
print(est_lasso2)

####################################
# SC - point estimation with ridge
est_ridge = scest(data_prep, w_constr={'name': "ridge"})
print(est_ridge)
Q_est = est_ridge.w_constr['Q']
est_ridge2 = scest(data_prep, w_constr={'p': 'L2', 'dir': '<=', 'Q': Q_est, 'lb': -numpy.inf})
print(est_ridge2)

####################################
# SC - point estimation with ols
est_ls = scest(data_prep, w_constr={'name': "ols"})
print(est_ls)
est_ls2 = scest(data_prep, w_constr={'p': 'no norm', 'dir': None, 'Q': None, 'lb': -numpy.inf})
print(est_ls2)


####################################
# Set options for inference
w_constr = {'name': 'simplex', 'Q': 1}
u_missp = True
u_sigma = "HC1"
u_order = 1
u_lags = 0
e_method = "qreg"
e_order = 1
e_lags = 0
e_alpha = 0.05
u_alpha = 0.05
sims = 500
cores = 1

random.seed(8894)
pi_si = scpi(data_prep, sims=sims, w_constr=w_constr, u_order=u_order, u_lags=u_lags,
             e_order=e_order, e_lags=e_lags, e_method=e_method, u_missp=u_missp,
             u_sigma=u_sigma, cores=cores, e_alpha=e_alpha, u_alpha=u_alpha)
print(pi_si)

####################################
# SC - plot results
plot = scplot(pi_si)
ggsave(filename='germany_unc.png', plot=plot)

################################################################################
# Other examples of data preparation

# multiple features
data_prep = scdata(df=data, id_var=id_var, time_var=time_var,
                   outcome_var=outcome_var, period_pre=period_pre,
                   period_post=period_post, unit_tr=unit_tr,
                   unit_co=unit_co, features=['gdp', 'trade'], cov_adj=cov_adj,
                   cointegrated_data=cointegrated_data, constant=constant,
                   report_missing=report_missing)

# multiple features and feature-specific covariate adjustment
data_prep = scdata(df=data, id_var=id_var, time_var=time_var,
                   outcome_var=outcome_var, period_pre=period_pre,
                   period_post=period_post, unit_tr=unit_tr,
                   unit_co=unit_co, features=['gdp', 'trade'], cov_adj=[['constant', 'trend'], ['constant']],
                   cointegrated_data=cointegrated_data, constant=constant,
                   report_missing=report_missing)
