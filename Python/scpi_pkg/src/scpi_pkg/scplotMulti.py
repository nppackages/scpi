# -*- coding: utf-8 -*-

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas
pandas.options.mode.chained_assignment = None
import numpy
from plotnine import ggplot, aes, geom_point, geom_errorbar, geom_vline, geom_line, geom_hline, theme, theme_bw, scale_x_datetime
from plotnine import element_blank, labs, scale_color_manual, ggtitle, facet_wrap, coord_flip, geom_ribbon
from copy import deepcopy
from .funs import CIrename, ix2rn
from math import ceil, floor
from mizani.breaks import date_breaks
from mizani.formatters import date_format


def scplotMulti(result,
                ptype="series",
                e_out=True,
                joint=False,
                col_treated="black",
                col_synth="mediumblue",
                scales="fixed",
                point_size=1.5,
                ncols=3,
                dateBreaks='10 years',
                dateFormat='%Y',
                e_method_input=None,
                save_data=None,
                verbose=True):

    """

    Parameters
    ----------
    result : scest_output/scpi_output
        a class 'scest_multi_output' object, obtained by calling scest, or a class
        'scpi_multi_output' object, obtained by calling scpi. The data object given as input to
        this command has to be processed with scdataMulti.

    ptype : str, default "series"
        a string that specifies the type of plot to be produced. If set to 'treatment', then treatment effects are
        plotted. If set to 'series' (default), the actual and synthetic time series are reported.

    e_out : bool, default True
        a logical specifying whether out-of-sample uncertainty should be included in the plot(s).

    joint : bool, default False
        a logical specifying whether simultaneous prediction intervals should be included in the plot(s).
        It requires e_out=True.

    col_treated : str, default "black"
        a string specifying the color for the treated unit series. Find the full list at
        http://sape.inf.usi.ch/quick-reference/ggplot2/colour.

    col_synth : str, default "mediumblue"
        a string specifying the color for the synthetic unit series. Find the full list at
        http://sape.inf.usi.ch/quick-reference/ggplot2/colour.

    scales : str, defaul "fixed"
        should axes scales be fixed ("fixed", the default), free ("free"),
        or free in one dimension ("free_x", "free_y")?

    point_size : float, default 1.5
        a scalar controlling the size of points in the scatter plot. Default is 1.5.

    ncols : int, default 3
        an integer controlling the number of columns in the plot.

    e_method_input : str, default None
        a string specifying the type of uncertainty estimation used for out-of-sample uncertainty quantification.
        To be used only when scpi received the option "e_method='all'" and the user wants to choose among the three
        techniques to quantify out-of-sample uncertainty.

    dateBreaks: str, default "10 years"
        a string specifying the breaks in the x-axis label. It is
        an interval specification, one of "sec", "min", "hour", "day", "week", "month", "year".
        Can be specified as an integer and a space, or followed by "s". Fractional seconds are supported.
        Some examples are "10 years", "2 months", etc. To be used only if time_var is not an integer!

    dateFormat: str, default "%Y"
        a string specifying the date/time format using standard POSIX specification.
        To be used only if time_var is not an integer!

    save_data : str, default None
        a string specifying the name (and the folder) of the saved dataframe containing the processed data used to
        produce the plot. The data is saved in .csv format and the folder specified.

    verbose : bool, default True
        if False prevents printing additional information in the console.

    Returns
    ----------
    plot : plotnine
        plotnine object that can be further modified.

    plotdata : dataframe
        dataframe object containing the processed data used to produce the plot.

    References
    ----------
    Abadie, A. (2021), “Using Synthetic Controls: Feasibility, Data Requirements, and Methodological
    Aspects,” Journal of Economic Literature, 59, 391-425.

    Cattaneo, M. D., Feng, Y., and Titiunik, R. (2021), “Prediction Intervals for Synthetic Control
    Methods,” Journal of the American Statistical Association, 116, 1865-1880.

    Cattaneo, M. D., Palomba, F., Feng, Y., and Titiunik, R. (2025), “scpi: Uncertainty Quantification for
    Synthetic Control Estimators”, Journal of Statistical Software, 113(2), 1-38.

    Cattaneo, M. D., Palomba, F., Feng, Y., and Titiunik, R. (2025), “Uncertainty Quantification in Synthetic
    Controls with Staggered Treatment Adoption”, Review of Economic Studies, doi:10.1080/01621459.2021.1979561.

    See Also
    --------
    scdata, scdataMulti, scest, scpi, scplot

    """

    class_input = result.__class__.__name__

    if class_input not in ['scest_multi_output', 'scpi_multi_output']:
        raise Exception("The object 'result' should be the output of scest or scpi " +
                        "when the data have been processed through scdataMulti!")

    if ptype not in ["series", "treatment"]:
        raise Exception("'type' should be either 'series' or 'treatment'!")

    if e_method_input not in [None, "gaussian", "ls", "qreg"]:
        raise Exception("'e_method_input' should be either None, 'gaussian', 'ls', or 'qreg'!")

    plot_type = deepcopy(result.effect)
    iota = deepcopy(result.iota)

    if plot_type == "time":
        iota = 1

    if iota > 20 and (plot_type != "unit" or ptype != "treatment") and verbose:
        warnings.warn(str(iota) + " treated units detected, therefore some graphs might be too crowded!" +
                      "We suggest saving the data with the option save_data, consult our replication " +
                      "files at https://nppackages.github.io/scpi/, and reproduce the same graph for just" +
                      " a fraction of the sample at a time!")

    Y_pre_fit = deepcopy(result.Y_pre_fit)
    Y_post_fit = deepcopy(result.Y_post_fit)

    Y_df = deepcopy(result.Y_df)
    Y_df.columns = ['ID', 'Time', 'Treatment', 'Actual']

    sel_units = [i in result.units_est for i in Y_df['ID']]
    res_df = deepcopy(Y_df[sel_units])

    # create hash maps if required
    if result.timeConvert is True:
        time_unique_ts = sorted(set(Y_df['Time'].tolist()))
        if plot_type != "time":
            int2ts = {i: time_unique_ts[i] for i in range(len(time_unique_ts))}
            ts2int = {time_unique_ts[i]: i for i in range(len(time_unique_ts))}
        else:
            ts2int = {time_unique_ts[i]: i + 2000 for i in range(len(time_unique_ts))}

        Y_df['Time'] = Y_df['Time'].map(ts2int)
        res_df = deepcopy(Y_df[sel_units])

        Y_pre_fit.reset_index(drop=False, inplace=True)
        Y_pre_fit['Time'] = Y_pre_fit['Time'].map(ts2int)
        Y_pre_fit.set_index(['ID', 'Time'], drop=True, inplace=True)

        if plot_type != "time":
            Y_post_fit.reset_index(drop=False, inplace=True)
            Y_post_fit['Time'] = Y_post_fit['Time'].map(ts2int)
            Y_post_fit.set_index(['ID', 'Time'], drop=True, inplace=True)

    synth_mat = pandas.concat([Y_pre_fit, Y_post_fit], axis=0)

    if plot_type != "time":
        synth_mat.index = synth_mat.index.rename(['ID', 'Time'])
        synth_mat.columns = ['Synthetic']

    if plot_type == "unit":
        Y_actual_pre = res_df[res_df['Treatment'] == 0]
        Y_actual_post = res_df[res_df['Treatment'] == 1]
        Y_actual_post_agg = Y_actual_post[['ID', 'Actual']].groupby(by='ID').mean()
        Y_actual_post_agg['Treatment'] = 1
        Y_actual_post_agg.set_index(Y_post_fit.index, inplace=True)
        Y_actual_pre.set_index(['ID', 'Time'], append=False, inplace=True)
        Y_actual = pandas.concat([Y_actual_pre, Y_actual_post_agg], axis=0)
    else:
        Y_actual = deepcopy(res_df)
        Y_actual.set_index(['ID', 'Time'], drop=True, inplace=True)

    treated_periods = Y_actual.loc[Y_actual['Treatment'] == 1]
    treated_periods.reset_index(drop=False, inplace=True)
    treated_reception = treated_periods.groupby('ID').min()
    treated_reception.columns = ["Tdate", "Treatment", "Actual"]
    ant_df = pandas.DataFrame.from_dict(result.anticipation, orient='index')
    ant_df.index.rename('ID', inplace=True)
    ant_df.columns = ["anticipation"]
    treated_reception = treated_reception.merge(ant_df, on="ID")
    treated_reception['Tdate'] = treated_reception['Tdate'] - treated_reception['anticipation'] - 1 / 2

    treated_reception.reset_index(drop=False, inplace=True)
    sel_units = [i in result.units_est for i in treated_reception['ID']]
    treated_reception = treated_reception[sel_units]

    if plot_type != "time":
        toplot = pandas.concat([Y_actual, synth_mat], axis=1, join='inner')

    elif plot_type == "time":
        res_df = res_df.merge(treated_reception[['ID', 'Tdate']], on="ID")
        Y_actual_pre = res_df[res_df['Time'] < res_df['Tdate']]
        Y_actual_post = res_df[res_df['Time'] > res_df['Tdate']]
        Y_actual_pre['Tdate'] = Y_actual_pre['Tdate'] + 1 / 2
        Y_actual_post['Tdate'] = Y_actual_post['Tdate'] + 1 / 2
        Y_actual_pre['tstd'] = Y_actual_pre['Time'] - Y_actual_pre['Tdate']
        Y_actual_post['tstd'] = Y_actual_post['Time'] - Y_actual_post['Tdate']

        names = synth_mat.index.values.tolist()
        names = [ix2rn(n).split(',') for n in names]
        unit = []
        unitagg = []
        time = []
        for n in names:
            if len(n) == 2:
                unit.append(n[0])
                time.append(n[1])
            else:
                unitagg.append(n[0])

        synth_pre = pandas.DataFrame({'ID': unit,
                                      'Time': time,
                                      'Synthetic': synth_mat.iloc[0:len(unit), 0].values})
        synth_pre = synth_pre.astype({'Time': Y_actual_pre['Time'].dtypes})
        Y_pre = Y_actual_pre.merge(synth_pre, on=['ID', 'Time'], how='left')

        max_pre = max(Y_actual_pre.groupby(['ID'])['tstd'].min()) + 1
        min_post = min([v for v in result.T1_outcome.values()]) - 1

        Y_pre_agg = Y_pre.groupby(['tstd'])[['Actual', 'Synthetic']].mean()
        Y_pre_agg.reset_index(inplace=True)
        Y_pre_agg.columns = ['Time', 'Actual', 'Synthetic']
        Y_pre_agg = Y_pre_agg[Y_pre_agg['Time'] >= max_pre]

        Y_post_agg = Y_actual_post.groupby(['tstd'])[['Actual']].mean()
        Y_post_agg.reset_index(inplace=True)
        Y_post_agg = Y_post_agg[Y_post_agg['tstd'] <= min_post]

        Y_post_agg = pandas.DataFrame({'ID': unitagg,
                                       'Actual': Y_post_agg['Actual'],
                                       'Synthetic': synth_mat.iloc[len(unit):, 0].values,
                                       'Time': range(0, len(unitagg))})

        Y_pre_agg['Treatment'] = 0
        Y_post_agg['Treatment'] = 1
        Y_pre_agg['ID'] = "aggregate"
        Y_post_agg['ID'] = "aggregate"

        Y_actual = pandas.concat([Y_pre_agg, Y_post_agg], axis=0)
        Y_actual['Tdate'] = 0

        plot_type = "unit-time"
        iota = 1
        treated_reception = pandas.DataFrame({'ID': 'aggregate',
                                              'Tdate': [0.5]})
        toplot = Y_actual
        toplot['Time'] = toplot['Time'] + 1

    toplot['Effect'] = toplot['Actual'] - toplot['Synthetic']

    if result.timeConvert is True and result.effect != "time":
        toplot.reset_index(drop=False, inplace=True)
        toplot['Time'] = toplot['Time'].map(int2ts)
        toplot.set_index(['ID', 'Time'], drop=True, inplace=True)
        tr_recp_ts = []
        for dd in treated_reception['Tdate'].values:
            avgtime = (int2ts[floor(dd)].asm8.astype(numpy.int64) + int2ts[ceil(dd)].asm8.astype(numpy.int64)) / 2
            tr_recp_ts.append(pandas.Timestamp(avgtime.astype('<M8[ns]')))
        treated_reception['Tdate'] = tr_recp_ts

    if plot_type != "time":
        toplot.reset_index(drop=False, inplace=True)

    toplot = toplot.merge(treated_reception[['ID', 'Tdate']], on='ID')

    if plot_type == 'unit-time' and ptype == "series":  # plot series

        toplot = toplot[["ID", "Time", "Actual", "Synthetic"]]
        toplot = pandas.melt(toplot,
                             id_vars=['ID', 'Time'],
                             value_vars=['Actual', 'Synthetic'],
                             var_name='Type',
                             value_name='Outcome')
        toplot = toplot.merge(treated_reception[['ID', 'Tdate']], on='ID')

        plot_struc = (ggplot(data=toplot) + theme_bw() +
                      theme(panel_grid=element_blank(),
                            legend_position=(.5, .05),
                            legend_direction='horizontal',
                            legend_box='horizontal',
                            legend_title=element_blank(),
                            subplots_adjust={'bottom': 0.2}) +
                      labs(x="Date", y="Outcome"))

        if result.timeConvert is True and result.effect != "time":
            plot_struc = (plot_struc + scale_x_datetime(breaks=date_breaks(dateBreaks), labels=date_format(dateFormat)))

        plot = (plot_struc +
                geom_line(mapping=aes(x='Time', y='Outcome', colour='Type')) +
                geom_point(mapping=aes(x='Time', y='Outcome', colour='Type'), size=point_size) +
                geom_vline(data=treated_reception, mapping=aes(xintercept='Tdate', group='ID')) +
                facet_wrap('ID', ncol=ncols, scales=scales) +
                scale_color_manual(name="", values=[col_treated, col_synth],
                                   labels=["Treated", "Synthetic Control"]))

    elif plot_type == "unit-time" and ptype == "treatment":

        plot_struc = (ggplot(data=toplot) + theme_bw() +
                      theme(panel_grid=element_blank(),
                            legend_position=(.5, .05),
                            legend_direction='horizontal',
                            legend_box='horizontal',
                            legend_title=element_blank(),
                            subplots_adjust={'bottom': 0.2}) +
                      labs(x="Date", y="Effect"))

        if result.timeConvert is True and result.effect != "time":
            plot_struc = (plot_struc + scale_x_datetime(breaks=date_breaks(dateBreaks), labels=date_format(dateFormat)))

        plot = (plot_struc +
                geom_line(mapping=aes(x='Time', y='Effect'), colour=col_synth) +
                geom_point(mapping=aes(x='Time', y='Effect'), colour=col_synth, size=point_size) +
                geom_vline(data=treated_reception, mapping=aes(xintercept='Tdate', group='ID')) +
                geom_hline(mapping=aes(yintercept=0, group='ID'), linetype="dashed") +
                facet_wrap('ID', ncol=ncols, scales=scales))

    elif plot_type == "unit" and ptype == "series":

        # add number of post-treatment periods that have been averaged
        auxdf = pandas.DataFrame({'ID': result.T1_outcome.keys(),
                                  'T1': result.T1_outcome.values()})
        treated_reception = treated_reception.merge(auxdf, on="ID")

        toplot = toplot[["ID", "Time", "Actual", "Synthetic"]]
        toplot = pandas.melt(toplot,
                             id_vars=['ID', 'Time'],
                             value_vars=['Actual', 'Synthetic'],
                             var_name='Type',
                             value_name='Outcome')
        toplot = toplot.merge(treated_reception[['ID', 'Tdate', 'T1']], on='ID')

        plot_struc = (ggplot(data=toplot) + theme_bw() +
                      theme(panel_grid=element_blank(),
                            legend_position=(.5, .05),
                            legend_direction='horizontal',
                            legend_box='horizontal',
                            legend_title=element_blank(),
                            subplots_adjust={'bottom': 0.2}) +
                      labs(x="Date", y="Outcome"))

        if result.timeConvert is True:
            plot_struc = (plot_struc + scale_x_datetime(breaks=date_breaks(dateBreaks), labels=date_format(dateFormat)))

        plot = (plot_struc +
                geom_line(data=toplot[toplot['Time'] < toplot['Tdate']],
                          mapping=aes(x='Time', y='Outcome', colour='Type')) +
                geom_point(mapping=aes(x='Time', y='Outcome', colour='Type'), size=point_size) +
                facet_wrap('ID', ncol=ncols, scales=scales) +
                scale_color_manual(name="", values=[col_treated, col_synth],
                                   labels=["Treated", "Synthetic Control"]))

    elif plot_type == "unit" and ptype == "treatment":

        # add number of post-treatment periods that have been averaged
        auxdf = pandas.DataFrame({'ID': result.T1_outcome.keys(),
                                  'T1': result.T1_outcome.values()})
        toplot = toplot.merge(auxdf, on="ID")

        plot_struc = (ggplot(data=toplot[toplot["Treatment"] == 1]) + theme_bw() +
                      theme(panel_grid=element_blank(),
                            legend_position=(.5, .05),
                            legend_direction='horizontal',
                            legend_box='horizontal',
                            legend_title=element_blank(),
                            subplots_adjust={'bottom': 0.2}) +
                      labs(x="Date", y="Effect"))

        plot = (plot_struc +
                geom_point(mapping=aes(x='ID', y='Effect', size='T1'), colour="#D67236", shape='o') +
                geom_hline(mapping=aes(yintercept=0), linetype="dashed") +
                coord_flip())

    #############################################################################
    # ADD UNCERTAINTY
    #############################################################################

    if class_input == "scpi_multi_output":
        e_method = deepcopy(result.e_method)
        if e_method == "all":
            e_method = "gaussian"
        if e_method_input is not None:
            e_method = deepcopy(e_method_input)

        if ptype == "series":

            toplot = toplot.merge(CIrename(result.CI_in_sample, 'insample'), on=['ID', 'Time'], how='left')

            if e_method in ["gaussian"]:
                toplot = toplot.merge(CIrename(result.CI_all_gaussian, 'gaussian'), on=['ID', 'Time'], how='left')

            if e_method in ["ls"]:
                toplot = toplot.merge(CIrename(result.CI_all_ls, 'ls'), on=['ID', 'Time'], how='left')

            if e_method in ["qreg"]:
                toplot = toplot.merge(CIrename(result.CI_all_qreg, 'qreg'), on=['ID', 'Time'], how='left')

            if joint is True:
                toplot = toplot.merge(CIrename(result.bounds['joint'], 'bd'), on=['ID', 'Time'], how='left')
                toplot['Lower_joint'] = toplot['Outcome'] + toplot['Lower_bd']
                toplot['Upper_joint'] = toplot['Outcome'] + toplot['Upper_bd']
                toplot.drop(['Lower_bd', 'Upper_bd'], axis=1, inplace=True)

        elif ptype == "treatment":

            toplot = toplot.merge(CIrename(result.bounds['insample'], 'bd'), on=['ID', 'Time'], how='left')
            toplot['Lower_insample'] = toplot['Effect'] - toplot['Lower_bd']
            toplot['Upper_insample'] = toplot['Effect'] - toplot['Upper_bd']
            toplot.drop(['Lower_bd', 'Upper_bd'], axis=1, inplace=True)

            if e_method in ["gaussian"]:
                toplot = toplot.merge(CIrename(result.bounds['subgaussian'], 'bd'), on=['ID', 'Time'], how='left')
                toplot['Lower_gaussian'] = toplot['Effect'] - toplot['Lower_bd']
                toplot['Upper_gaussian'] = toplot['Effect'] - toplot['Upper_bd']
                toplot.drop(['Lower_bd', 'Upper_bd'], axis=1, inplace=True)

            if e_method in ["ls"]:
                toplot = toplot.merge(CIrename(result.bounds['ls'], 'bd'), on=['ID', 'Time'], how='left')
                toplot['Lower_ls'] = toplot['Effect'] - toplot['Lower_bd']
                toplot['Upper_ls'] = toplot['Effect'] - toplot['Upper_bd']
                toplot.drop(['Lower_bd', 'Upper_bd'], axis=1, inplace=True)

            if e_method in ["qreg"]:
                toplot = toplot.merge(CIrename(result.bounds['qreg'], 'bd'), on=['ID', 'Time'], how='left')
                toplot['Lower_qreg'] = toplot['Effect'] - toplot['Lower_bd']
                toplot['Upper_qreg'] = toplot['Effect'] - toplot['Upper_bd']
                toplot.drop(['Lower_bd', 'Upper_bd'], axis=1, inplace=True)

            if joint is True:
                toplot = toplot.merge(CIrename(result.bounds['joint'], 'bd'), on=['ID', 'Time'], how='left')
                toplot['Lower_joint'] = toplot['Effect'] - toplot['Lower_bd']
                toplot['Upper_joint'] = toplot['Effect'] - toplot['Upper_bd']
                toplot.drop(['Lower_bd', 'Upper_bd'], axis=1, inplace=True)

        # In sample uncertainty
        if e_out is False:

            if ptype == "treatment" and plot_type == "unit":
                plot_w = (plot +
                          geom_errorbar(data=toplot,
                                        mapping=aes(x='ID', ymin='Lower_insample', ymax='Upper_insample'),
                                        colour="#FD6467", width=0.25, linetype="solid") +
                          ggtitle("In-sample Uncertainty"))

            else:
                plot_w = (plot +
                          geom_errorbar(data=toplot,
                                        mapping=aes(x='Time', ymin='Lower_insample', ymax='Upper_insample'),
                                        colour=col_synth, width=0.5, linetype="solid") +
                          ggtitle("In-sample Uncertainty"))

        # In and out of sample uncertainty: subgaussian
        if e_out is True and e_method in ["gaussian"]:

            if ptype == "treatment" and plot_type == "unit":
                plot_w1 = (plot +
                           geom_errorbar(data=toplot,
                                         mapping=aes(x='ID', ymin='Lower_gaussian', ymax='Upper_gaussian'),
                                         colour="#FD6467", width=0.25, linetype="solid") +
                           ggtitle("In and Out of Sample Uncertainty - Subgaussian Bounds"))

            else:
                plot_w1 = (plot +
                           geom_errorbar(data=toplot,
                                         mapping=aes(x='Time', ymin='Lower_gaussian', ymax='Upper_gaussian'),
                                         colour=col_synth, width=0.5, linetype="solid") +
                           ggtitle("In and Out of Sample Uncertainty - Subgaussian Bounds"))

            if joint is True and plot_type == "unit-time":
                if ptype == "treatment":
                    plotdf = toplot[toplot['Treatment'] == 1]
                else:
                    plotdf = toplot[toplot['Type'] == "Synthetic"]

                plot_w1 = (plot_w1 +
                           geom_ribbon(data=plotdf,
                                       mapping=aes(x='Time', ymin='Lower_joint', ymax='Upper_joint'),
                                       fill=col_synth, alpha=0.1))

            elif joint is True and plot_type == "unit":
                if ptype == "treatment":
                    plotdf = toplot[toplot['Treatment'] == 1]
                    plot_w1 = (plot_w1 +
                               geom_errorbar(data=plotdf,
                                             mapping=aes(x='ID', ymin='Lower_joint', ymax='Upper_joint'),
                                             colour="#5B1A18", width=0.25))
                elif ptype == "series":
                    plotdf = toplot[toplot['Type'] == 'Synthetic']
                    plot_w1 = (plot_w1 +
                               geom_errorbar(data=plotdf,
                                             mapping=aes(x='Time', ymin='Lower_joint', ymax='Upper_joint'),
                                             color="darkred", width=0.6))

        # In and out of sample uncertainty: location-scale

        if e_out is True and e_method in ["ls"]:
            if ptype == "treatment" and plot_type == "unit":
                plot_w2 = (plot +
                           geom_errorbar(data=toplot,
                                         mapping=aes(x='ID', ymin='Lower_ls', ymax='Upper_ls'),
                                         colour="#FD6467", width=0.25, linetype="solid") +
                           ggtitle("In and Out of Sample Uncertainty - Location-scale Model"))

            else:
                plot_w2 = (plot +
                           geom_errorbar(data=toplot,
                                         mapping=aes(x='Time', ymin='Lower_ls', ymax='Upper_ls'),
                                         colour=col_synth, width=0.5, linetype="solid") +
                           ggtitle("In and Out of Sample Uncertainty - Location-scale Model"))

            if joint is True and plot_type == "unit-time":
                if ptype == "treatment":
                    plotdf = toplot[toplot['Treatment'] == 1]
                else:
                    plotdf = toplot[toplot['Type'] == "Synthetic"]

                plot_w2 = (plot_w2 +
                           geom_ribbon(data=plotdf,
                                       mapping=aes(x='Time', ymin='Lower_joint', ymax='Upper_joint'),
                                       fill=col_synth, alpha=0.1))

            elif joint is True and plot_type == "unit":
                if ptype == "treatment":
                    plotdf = toplot[toplot['Treatment'] == 1]
                    plot_w2 = (plot_w2 +
                               geom_errorbar(data=plotdf,
                                             mapping=aes(x='ID', ymin='Lower_joint', ymax='Upper_joint'),
                                             colour="#5B1A18", width=0.25))
                elif ptype == "series":
                    plotdf = toplot[toplot['Synthetic'] == 1]
                    plot_w2 = (plot_w2 +
                               geom_errorbar(data=plotdf,
                                             mapping=aes(x='Time', ymin='Lower_joint', ymax='Upper_joint'),
                                             color="darkred", width=0.6))

        # In and out of sample uncertainty: quantile regression

        if e_out is True and e_method in ["qreg"]:
            if ptype == "treatment" and plot_type == "unit":
                plot_w3 = (plot +
                           geom_errorbar(data=toplot,
                                         mapping=aes(x='ID', ymin='Lower_qreg', ymax='Upper_qreg'),
                                         colour="#FD6467", width=0.25, linetype="solid") +
                           ggtitle("In and Out of Sample Uncertainty - Qunatile Regression"))

            else:
                plot_w3 = (plot +
                           geom_errorbar(data=toplot,
                                         mapping=aes(x='Time', ymin='Lower_qreg', ymax='Upper_qreg'),
                                         colour=col_synth, width=0.5, linetype="solid") +
                           ggtitle("In and Out of Sample Uncertainty - Qunatile Regression"))

            if joint is True and plot_type == "unit-time":
                if ptype == "treatment":
                    plotdf = toplot[toplot['Treatment'] == 1]
                else:
                    plotdf = toplot[toplot['Type'] == "Synthetic"]

                plot_w3 = (plot_w3 +
                           geom_ribbon(data=plotdf,
                                       mapping=aes(x='Time', ymin='Lower_joint', ymax='Upper_joint'),
                                       fill=col_synth, alpha=0.1))

            elif joint is True and plot_type == "unit":
                if ptype == "treatment":
                    plotdf = toplot[toplot['Treatment'] == 1]
                    plot_w3 = (plot_w3 +
                               geom_errorbar(data=plotdf,
                                             mapping=aes(x='ID', ymin='Lower_joint', ymax='Upper_joint'),
                                             colour="#5B1A18", width=0.25))
                elif ptype == "series":
                    plotdf = toplot[toplot['Synthetic'] == 1]
                    plot_w3 = (plot_w3 +
                               geom_errorbar(data=plotdf,
                                             mapping=aes(x='Time', ymin='Lower_joint', ymax='Upper_joint'),
                                             color="darkred", width=0.6))

    # Save data to reproduce plot
    if save_data is not None:
        data_name = save_data + ".csv"
        toplot.to_csv(data_name)

    if class_input == "scpi_multi_output":
        if e_out is False:
            return plot_w
        elif e_method == "gaussian":
            return plot_w1
        elif e_method == "ls":
            return plot_w2
        elif e_method == "qreg":
            return plot_w3

    else:
        return plot
