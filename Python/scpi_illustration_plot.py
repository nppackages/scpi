################################################################################
# SCPI Python Package
# Script for Visualization - Single Treated Unit
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

numpy.random.seed(8894)
result = scpi(data_prep, sims=sims, w_constr=w_constr, u_order=u_order, u_lags=u_lags,
              e_order=e_order, e_lags=e_lags, e_method=e_method, u_missp=u_missp,
              u_sigma=u_sigma, cores=cores, e_alpha=e_alpha, u_alpha=u_alpha)


####################################
# SC - plot results
plot = scplot(result)

####################################
# SC - manually reproduce plot

sc_l = result.CI_all_gaussian.iloc[:, [0]].to_numpy()
sc_r = result.CI_all_gaussian.iloc[:, [1]].to_numpy()

# Store data for treated unit
time = numpy.concatenate([period_pre, period_post])
y_act = pandas.concat([result.Y_pre, result.Y_post]).to_numpy().flatten()
data_points_act = pandas.DataFrame({'time': time,
                                    'y_act': y_act})

# Store data for synthetic control unit
y_sc_df = pandas.concat([result.Y_pre_fit, result.Y_post_fit])

y_sc_na = pandas.DataFrame(numpy.array([numpy.nan] * len(time)))
sc_l_na = pandas.DataFrame(numpy.array([numpy.nan] * len(time)))
sc_r_na = pandas.DataFrame(numpy.array([numpy.nan] * len(time)))

# Check if some periods have missing point estimate/missing CI
not_miss_plot = [t in y_sc_df.index.get_level_values(1).tolist() for t in time]
not_miss_ci = [t in result.CI_all_gaussian.index.get_level_values(1).tolist() for t in time]

y_sc_na.loc[not_miss_plot, ] = y_sc_df.iloc[:, [0]].to_numpy()
sc_l_na.loc[not_miss_ci, ] = sc_l
sc_r_na.loc[not_miss_ci, ] = sc_r

data_points_act = pandas.DataFrame({'time': time,
                                    'y_act': y_act
                                    })

data_points = pandas.concat([data_points_act, y_sc_na,
                            sc_l_na, sc_r_na], axis=1)
data_points.columns = ['time', 'y_act', 'y_sc', 'lb', 'ub']

data_points['tr'] = ['Treated'] * len(y_sc_na)
data_points['sc'] = ['Synthetic Control'] * len(y_sc_na)

T0 = period_pre[len(period_pre) - 1]
col_dots_t = 'black'
col_line_t = 'black'
col_dots_s = 'blue'
col_line_s = 'blue'

plot_struc = (ggplot(data_points) +
              theme_bw() +
              theme(panel_grid=element_blank(),
                    legend_position=(.5, .05),
                    legend_direction='horizontal',
                    legend_box='horizontal',
                    legend_title=element_blank(),
                    subplots_adjust={'bottom': 0.2}) +
              labs(x='Year', y='GDP per capita (thousand US dollars)')
              )

plot_lines = (plot_struc +
              geom_point(mapping=aes(x='time', y='y_act', color='tr'), shape='o', fill='white', na_rm=False) +
              geom_point(mapping=aes(x='time', y='y_sc', color='sc'), shape='o', na_rm=False) +
              geom_line(mapping=aes(x='time', y='y_act', color='tr'), na_rm=False) +
              geom_line(mapping=aes(x='time', y='y_sc', color='sc'), linetype='dashed', na_rm=False) +
              geom_vline(xintercept=T0, linetype='dotted') +
              scale_color_manual(name="", values=[col_dots_t, col_dots_s],
                                 labels=["Treated", "Synthetic Control"],
                                 guide=guide_legend(override_aes={'linetype': ['solid', 'dashed'],
                                                                  'shape': ['o', 'o']})))

title_plot = 'In and Out of Sample Uncertainty'
plot = plot_lines + geom_errorbar(mapping=aes(x='time', ymin='lb', ymax='ub', color='sc'),
                                  size=0.5, linetype='solid') + ggtitle(title_plot)
