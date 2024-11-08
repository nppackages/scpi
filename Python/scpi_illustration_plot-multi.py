################################################################################
# SCPI Python Package
# Script for Visualization - Multiple Treated Units
# Authors: Matias D. Cattaneo, Yingjie Feng, Filippo Palomba and Rocio Titiunik
################################################################################

########################################
# Load SCPI_PKG package
import pandas
import numpy
import os
from copy import deepcopy
from plotnine import ggplot, aes, geom_point, geom_errorbar, geom_vline, geom_line, theme, theme_bw
from plotnine import element_blank, labs, guide_legend, scale_color_manual, ggtitle, facet_wrap, geom_ribbon
from warnings import filterwarnings
from scpi_pkg.scdata import scdata
from scpi_pkg.scdataMulti import scdataMulti
from scpi_pkg.scest import scest
from scpi_pkg.scpi import scpi
from scpi_pkg.scplot import scplot
from scpi_pkg.scplotMulti import scplotMulti


########################################
# Load database
os.chdir("YOUR_PATH_HERE")
data = pandas.read_csv("scpi_germany.csv")

numpy.random.seed(8894)
##############################################################################
# MULTIPLE TREATED UNITS
##############################################################################

########################################
# Load database
data = pandas.read_csv("scpi_germany.csv")

# create treatment dummy
data['status'] = 0
data.loc[(data['country'] == "West Germany") & (data['year'] >= 1991), 'status'] = 1

# fictitious second treated unit for the sake of the example
data.loc[(data['country'] == "Italy") & (data['year'] >= 1992), 'status'] = 1

id_var = 'country'
time_var = 'year'
treatment_var = 'status'
outcome_var = 'gdp'
covs_adj = {'Italy': ['constant', 'trend'],
            'West Germany': [['constant', 'trend'], ['constant', 'trend']]}

aux = scdataMulti(df=data,
                  id_var=id_var,
                  treatment_var=treatment_var,
                  outcome_var=outcome_var,
                  time_var=time_var,
                  features={"Italy": ["gdp", "trade"],
                            "West Germany": ["gdp", "infrate"]},
                  constant={'Italy': True, 'West Germany': False},
                  cointegrated_data=True,
                  cov_adj=covs_adj)

result = scpi(aux, w_constr={'name': 'simplex'}, sims=10, cores=1)

plot = scplotMulti(result, ptype="series", joint=True, save_data='__scpi_plot')

toplot = pandas.read_csv('__scpi_plot.csv')

Y_actual = deepcopy(result.Y_df)
Y_actual.columns = ['ID', 'Time', 'Treatment', 'Actual']
Y_actual.set_index(['ID', 'Time'], drop=True, inplace=True)
treated_periods = Y_actual.loc[Y_actual['Treatment'] == 1]
treated_periods.reset_index(drop=False, inplace=True)
treated_reception = treated_periods.groupby('ID').min()
treated_reception.columns = ["Tdate", "Treatment", "Actual"]
treated_reception['Tdate'] = treated_reception['Tdate'] - 1 / 2
treated_reception.reset_index(drop=False, inplace=True)

plot_struc = (ggplot(data=toplot) + theme_bw() +
              theme(panel_grid=element_blank(),
                    legend_position=(.5, .05),
                    legend_direction='horizontal',
                    legend_box='horizontal',
                    legend_title=element_blank(),
                    subplots_adjust={'bottom': 0.2}) +
              labs(x="Date", y="Outcome"))

plot = (plot_struc +
        geom_line(mapping=aes(x='Time', y='Outcome', colour='Type')) +
        geom_point(mapping=aes(x='Time', y='Outcome', colour='Type'), size=1.5) +
        geom_vline(data=treated_reception, mapping=aes(xintercept='Tdate', group='ID')) +
        facet_wrap('ID', ncol=2) +
        scale_color_manual(name="", values=["black", "blue"],
                           labels=["Treated", "Synthetic Control"]))

# add prediction intervals
plot_w1 = (plot +
           geom_errorbar(data=toplot,
                         mapping=aes(x='Time', ymin='Lower_gaussian', ymax='Upper_gaussian'),
                         colour='blue', width=0.5, linetype="solid") +
           ggtitle("In and Out of Sample Uncertianty - Subgaussian Bounds"))

# add shaded area
plotdf = toplot[toplot['Type'] == 'Synthetic']
plot_w1 = (plot_w1 +
           geom_ribbon(data=plotdf,
                       mapping=aes(x='Time', ymin='Lower_joint', ymax='Upper_joint'),
                       fill='blue', alpha=0.1))
