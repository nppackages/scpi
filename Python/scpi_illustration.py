###############################################################################
# SCPI Python Package
# Script for Empirical Illustration - Single Treated Unit
# Authors: Matias D. Cattaneo, Yingjie Feng, Filippo Palomba and Rocio Titiunik
###############################################################################

########################################
# Load SCPI_PKG package
import pandas
import numpy
import random
import os
from copy import deepcopy
from warnings import filterwarnings
from scpi_pkg.scdata import scdata
from scpi_pkg.scdataMulti import scdataMulti
from scpi_pkg.scest import scest
from scpi_pkg.scpi import scpi
from scpi_pkg.scplot import scplot
from scpi_pkg.scplotMulti import scplotMulti

########################################
# Load database
data = pandas.read_csv("scpi_germany.csv")

filterwarnings("ignore")

##############################################################################
# SINGLE TREATED UNIT
##############################################################################

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
cointegrated_data = True

data_prep = scdata(df=data, id_var=id_var, time_var=time_var,
                   outcome_var=outcome_var, period_pre=period_pre,
                   period_post=period_post, unit_tr=unit_tr,
                   unit_co=unit_co, features=features, cov_adj=cov_adj,
                   cointegrated_data=cointegrated_data, constant=constant)


####################################
# SC - point estimation with simplex
est_si = scest(data_prep, w_constr={'name': "simplex"})
print(est_si)
est_si2 = scest(data_prep, w_constr={'p': 'L1', 'dir': '==', 'Q': 1, 'lb': 0})
print(est_si2)


####################################
# SC - plot results
plot = scplot(est_si)

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
Q_est = est_ridge.w_constr['West Germany']['Q']
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
e_method = "gaussian"
e_order = 1
e_lags = 0
e_alpha = 0.05
u_alpha = 0.05
sims = 200
cores = 1

random.seed(8894)
pi_si = scpi(data_prep, sims=sims, w_constr=w_constr, u_order=u_order, u_lags=u_lags,
             e_order=e_order, e_lags=e_lags, e_method=e_method, u_missp=u_missp,
             u_sigma=u_sigma, cores=cores, e_alpha=e_alpha, u_alpha=u_alpha)
print(pi_si)

####################################
# SC - plot results
plot = scplot(pi_si)

################################################################################
# Other examples of data preparation

# multiple features
data_prep = scdata(df=data, id_var=id_var, time_var=time_var,
                   outcome_var=outcome_var, period_pre=period_pre,
                   period_post=period_post, unit_tr=unit_tr,
                   unit_co=unit_co, features=['gdp', 'trade'], cov_adj=cov_adj,
                   cointegrated_data=cointegrated_data, constant=constant)

# multiple features and feature-specific covariate adjustment
data_prep = scdata(df=data, id_var=id_var, time_var=time_var,
                   outcome_var=outcome_var, period_pre=period_pre,
                   period_post=period_post, unit_tr=unit_tr,
                   unit_co=unit_co, features=['gdp', 'trade'], cov_adj=[['constant', 'trend'], []],
                   cointegrated_data=cointegrated_data, constant=constant)

################################################################################
# Specifying features for different pre-treatment periods or just use pre-treatment averages

# I) we want to include "trade" just for some selected periods, i.e., 1960, 1970, 1980, 1990
tradeUnselPeriods = [t for t in range(1960, 1991) if t not in [1960, 1970, 1980, 1990]]
data['tradeAux'] = data['trade']
data.loc[data['year'].isin(tradeUnselPeriods), 'tradeAux'] = numpy.nan

data_prep = scdata(df=data, id_var=id_var, time_var=time_var,
                   outcome_var=outcome_var, period_pre=period_pre,
                   period_post=period_post, unit_tr=unit_tr,
                   unit_co=unit_co, features=['gdp', 'tradeAux'],
                   cointegrated_data=cointegrated_data, constant=constant)

data_prep.B

# II) we want to include just the pre-treatment average of "infrate"
dataAux = data[data["year"] <= 1990].groupby(['country'])[['infrate']].mean()
dataAux.rename(columns={"infrate": "infrateAvg"}, inplace=True)
data = data.merge(dataAux, on="country")
data.loc[data['year'] != 1990, "infrateAvg"] = numpy.nan

data_prep = scdata(df=data, id_var=id_var, time_var=time_var,
                   outcome_var=outcome_var, period_pre=period_pre,
                   period_post=period_post, unit_tr=unit_tr,
                   unit_co=unit_co, features=['gdp', 'infrateAvg'],
                   cointegrated_data=cointegrated_data, constant=constant)

data_prep.B
