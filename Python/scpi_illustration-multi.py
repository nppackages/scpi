###############################################################################
# SCPI Python Package
# Script for Empirical Illustration - Multiple Treated Units
# Authors: Matias D. Cattaneo, Yingjie Feng, Filippo Palomba and Rocio Titiunik
###############################################################################

########################################
# Load SCPI_PKG package
import pandas
import numpy
import os
from warnings import filterwarnings
from scpi_pkg.scdataMulti import scdataMulti
from scpi_pkg.scest import scest
from scpi_pkg.scpi import scpi
from scpi_pkg.scplotMulti import scplotMulti

########################################
# Load database
os.chdir("YOUR_PATH_HERE")
data = pandas.read_csv("scpi_germany.csv")

filterwarnings("ignore")
numpy.random.seed(8894)

##############################################################################
# MULTIPLE TREATED UNITS
##############################################################################

# create treatment dummy
data['status'] = 0
data.loc[(data['country'] == "West Germany") & (data['year'] >= 1991), 'status'] = 1

# Create a second placebo treated unit
data.loc[(data['country'] == "Italy") & (data['year'] >= 1992), 'status'] = 1

id_var = 'country'
time_var = 'year'
treatment_var = 'status'
outcome_var = 'gdp'
covs_adj = {'Italy': ['constant', 'trend'],
            'West Germany': [['constant', 'trend'], ['constant', 'trend']]}


###############################################
# unit-time treatment effect
###############################################

aux = scdataMulti(df=data,
                  id_var=id_var,
                  treatment_var=treatment_var,
                  outcome_var=outcome_var,
                  time_var=time_var,
                  features={"Italy": ["gdp", "trade"],
                            "West Germany": ["gdp", "infrate"]},
                  constant={'Italy': True, 'West Germany': False},
                  cointegrated_data=True, cov_adj=covs_adj)

res = scest(aux, w_constr={'name': 'simplex'})
scplotMulti(res)

res_pi = scpi(aux, w_constr={'name': 'simplex'}, sims=50, cores=1)

# plot series
scplotMulti(res_pi, ptype="series")

# plot treatment effects
scplotMulti(res_pi, ptype="treatment", joint=True)

###############################################
# average unit post-treatment effect
###############################################

aux = scdataMulti(df=data,
                  id_var=id_var,
                  treatment_var=treatment_var,
                  outcome_var=outcome_var,
                  time_var=time_var,
                  features={"Italy": ["gdp", "trade"],
                            "West Germany": ["gdp", "infrate"]},
                  constant={'Italy': True, 'West Germany': False},
                  cointegrated_data=True,
                  cov_adj=covs_adj, effect="unit")

res = scest(aux, w_constr={'name': 'simplex'})
scplotMulti(res)

res_pi = scpi(aux, w_constr={'name': 'simplex'}, sims=50, cores=1)

# plot series
scplotMulti(res_pi, ptype="series")

# plot treatment effects
scplotMulti(res_pi, ptype="treatment", joint=True)

###############################################
# average post-treatment effect on the treated
###############################################

aux = scdataMulti(df=data,
                  id_var=id_var,
                  treatment_var=treatment_var,
                  outcome_var=outcome_var,
                  time_var=time_var,
                  features={"Italy": ["gdp", "trade"],
                            "West Germany": ["gdp", "infrate"]},
                  constant={'Italy': True, 'West Germany': False},
                  cointegrated_data=True,
                  cov_adj=covs_adj, effect="time")

res = scest(aux, w_constr={'name': 'simplex'})
scplotMulti(res)

res_pi = scpi(aux, w_constr={'name': 'simplex'}, sims=50, cores=1)

# plot series
scplotMulti(res_pi, ptype="series")

# plot treatment effects
scplotMulti(res_pi, ptype="treatment", joint=True)
