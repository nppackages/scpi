# -*- coding: utf-8 -*-

# Temporary code to suppress pandas FutureWarning
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas
pandas.options.mode.chained_assignment = None
import numpy
from copy import deepcopy
from .funs import w_constr_prep, b_est, b_est_multi, V_prep, mat2dict, ix2rn
from .scplot import scplot
from .scplotMulti import scplotMulti


def scest(df, w_constr=None, V="separate", Vmat=None, plot=False):

    """

    Parameters
    ----------
    df : scdata_output
        a class scdata_output object, obtained by calling scdata

    w_constr : dictionary, default {"name": "simplex"}
        a dictionary specifying the constraint set the estimated weights of the donors must belong to.
        w_constr can contain up to four objects:
        - p, a string indicating the norm to be used (p should be one of "no norm", "L1", and "L2")
        - dir, a string indicating whether the constraint on the norm is an equality ("==") or inequality ("<=")
        - Q, a scalar defining the value of the constraint on the norm
        - lb, a scalar defining the lower bound on the weights. It can be either 0 or -numpy.inf.
        - name, a character selecting one of the default proposals.

    V : str, default "separate"
        a weighting matrix to be used when minimizing the sum of squared residuals.
        The default is the identity matrix ("separate"), so equal weight is given to all observations.
        The other possibility is to specify V = "pooled" for the pooled fit.

    Vmat : numpy.array, defaul None
        a conformable weighting matrix to be used in the minimization of the sum of squared residuals. To check the proper
        dimensions, we suggest to check the output of scdata or scdataMulti and inspect the dimensions of B and C.

    plot : bool, default False
        a logical specifying whether scplot should be called and a plot saved in the current working directory. For more
        options see scplot.

    Returns
    -------
    The function returns an object of class 'scest_output' containing the following objects

    w : pandas.DataFrame
        a dataframe containing the weights of the donors.

    r : pandas.DataFrame
        a dataframe containing the values of the covariates used for adjustment.

    b : pandas.DataFrame
        a dataframe containing w and r.

    Y_pre_fit : pandas.DataFrame
        a dataframe containing the estimated pre-treatment outcome for the SC unit.

    Y_post_fit : pandas.DataFrame
        a dataframe containing the estimated post-treatment outcome for the SC unit.

    A_hat : pandas.DataFrame
        a dataframe containing the predicted values of the features of the treated unit.

    res : pandas.DataFrame
        a dataframe containing the residuals A - A_hat.

    V : numpy.array
        an array containing the weighting matrix used in estimation.

    w_constr : dictionary
        a dictionary containing the specifics of the constraint set used on the weights.

    A : pandas.DataFrame
        a dataframe containing pre-treatment features of the treated unit.

    B : pandas.DataFrame
        a dataframe containing pre-treatment features of the control units.

    C : pandas.DataFrame
        a dataframe containing covariates for adjustment.

    P : pandas.DataFrame
        a dataframe whose rows are the vectors used to predict the out-of-sample series for the synthetic unit(s).

    Y_pre : pandas.DataFrame
        a dataframe containing the pre-treatment outcome of the treated unit(s). If multiple treated units are present and the desired
        predictand involves aggregation (e.g., effect = "time" or effect = "unit) then it contains only the raw data before aggregation.
        For the aggregated data see 'Y_actual'.

    Y_post : pandas.DataFrame
        a dataframe containing the post-treatment outcome of the treated unit(s). If multiple treated units are present and the desired
        predictand involves aggregation (e.g., effect = "time" or effect = "unit) then it contains only the raw data before aggregation.
        For the aggregated data see 'Y_actual'.

    Y_actual : pandas.DataFrame
        a dataframe containing the pre- and post-treatment outcome of the treated unit(s). If the desired predictand
        involves aggregation (e.g., effect = "time" or effect = "unit) then it contains the data after aggregation.
        For the disaggregated data see 'Y_pre' and 'Y_post'.

    Y_donors : pandas.DataFrame
        a dataframe containing the pre-treatment outcome of the control units.

    J : int
        the number of control units

    K : array
        a numeric array with the number of covariates used for adjustment for each feature

    KM : int
        the total number of covariates used for adjustment for each treated unit

    KMI : int
        the total number of covariates used for adjustment

    M : int
        number of features

    period_pre : array
        a numeric array with the pre-treatment period

    period_post : array
        a numeric array with the post-treatment period

    T0_features : array
        a numeric array with the number of periods used in estimation for each feature

    T1_outcome : int
        the number of post-treatment periods

    glob_cons : bool
        for internal use only

    out_in_features : bool
        for internal use only

    cointegrated_data: bool
        logical indicating whether the data has been treated as cointegrated.

    features : list
        a list with the name of the features

    out_in_features : bool
        for internal use only

    effect : str
        for internal use only

    References
    ----------
    Abadie, A. (2021), “Using Synthetic Controls: Feasibility, Data Requirements, and Methodological
    Aspects,” Journal of Economic Literature, 59, 391-425.

    Cattaneo, M. D., Feng, Y., and Titiunik, R. (2021), “Prediction Intervals for Synthetic Control
    Methods,” Journal of the American Statistical Association, 116, 1865-1880.

    Cattaneo, M. D., Palomba, F., Feng, Y., and Titiunik, R. (2022), “scpi: Uncertainty Quantification for
    Synthetic Control Estimators”.

    Cattaneo, M. D., Palomba, F., Feng, Y., and Titiunik, R. (2022), “Uncertainty Quantification in Synthetic
    Controls with Staggered Treatment Adoption”.

    See Also
    --------
    scdata, scdataMulti, scpi, scplot, scplotMulti

    """

    if df.__class__.__name__ not in ['scdata_output', 'scdata_multi_output']:
        raise Exception("df should be the object returned by running scdata or scdataMulti!")

    w_names = []
    w_values = []
    if w_constr is not None:
        if not isinstance(w_constr, dict):
            raise Exception("w_constr should be an object of type dict!")

        for name, value in w_constr.items():
            w_names.append(name)
            w_values.append(value)

        if 'name' in w_names:

            if not w_constr['name'] in ["simplex", "lasso", "ridge", "ols", "L1-L2"]:
                raise Exception("If 'name' is specified in w_constr, then it should be chosen" +
                                " among 'simplex', 'lasso', 'ridge', 'ols', 'L1-L2'.")
        else:
            w_constr['name'] = None

    # Data matrices and specs
    A = deepcopy(df.A)
    B = deepcopy(df.B)
    C = deepcopy(df.C)
    P = deepcopy(df.P)
    Z = pandas.concat([B, C], axis=1)
    Y_donors = df.Y_donors

    if df.__class__.__name__ == "scdata_output":
        class_type = "scpi_data"
    else:
        class_type = "scpi_data_multi"

    V_type = deepcopy(V)
    J = df.J
    K = df.K
    KM = df.KM
    M = df.M

    if class_type == "scpi_data":
        iota = 1
        KMI = deepcopy(KM)
        Jtot = deepcopy(J)

    elif class_type == "scpi_data_multi":
        KMI = df.KMI
        iota = df.iota
        Jtot = sum(J.values())

    T0_features = df.T0_features
    out_in_features = df.out_in_features
    outcome_var = df.outcome_var

    ##########################
    # Set up the estimation problem

    # Create weighting matrix
    if not isinstance(V, str):
        raise Exception("The object V should be a string! If you want to manually specify the weighting matrix " +
                        "consider using the option 'Vmat'!")

    if Vmat is not None:
        if not isinstance(Vmat, (pandas.DataFrame, numpy.ndarray, str)):
            raise Exception("The object Vmat should a pandas.dataframe or a numpy.array!")

        else:
            V_shape = numpy.shape(Vmat)
            if V_shape[0] != len(B) or V_shape[1] != len(B):
                raise Exception("The object Vmat should be of shape (" + str(len(B)) +
                                "," + str(len(B)) + ")!")

            Vmat = pandas.DataFrame(Vmat, index=B.index,
                                    columns=B.index.get_level_values('ID'))
    else:
        Vmat = V_prep(V_type, B, T0_features, iota)

    V = Vmat

    # Estimate SC
    if class_type == "scpi_data":
        w_constr = w_constr_prep(w_constr, w_names, A, Z, V, J, KM)
        result = b_est(A=A, Z=Z, J=J, KM=KM, w_constr=w_constr, V=V)
        cnms = [c.split("_", 1)[1] for c in Z.columns.tolist()]
        idx = pandas.MultiIndex.from_product([df.treated_units, cnms], names=['ID', 'donor'])
        b = pandas.DataFrame(result, index=idx)
        w_constr_dict = {df.treated_units[0]: w_constr}

    elif class_type == "scpi_data_multi":
        A_dict = mat2dict(A, cols=False)
        B_dict = mat2dict(B)
        C_dict = mat2dict(C)
        V_dict = mat2dict(V)

        w_constr_dict = {}

        w_store = pandas.DataFrame(None)
        r_store = pandas.DataFrame(None)

        for tr in df.treated_units:
            A_i = A_dict[tr]
            Z_i = pandas.concat([B_dict[tr], C_dict[tr]], axis=1)
            V_i = V_dict[tr]

            w_constr_dict[tr] = w_constr_prep(w_constr, w_names, A_i, Z_i, V_i, J[tr], KM[tr])

            if V_type == "separate":
                result = b_est(A=A_i, Z=Z_i, J=J[tr], KM=KM[tr], w_constr=w_constr_dict[tr],
                               V=V_i)
                idx = pandas.MultiIndex.from_product([[tr], df.donors_dict[tr]])
                auxdf = pandas.DataFrame(result[0:J[tr]], index=idx)
                w_store = pandas.concat([w_store, auxdf], axis=0)

                if KM[tr] > 0:
                    cnm = [c.split("_", 1)[1] for c in C_dict[tr].columns.tolist()]
                    idx = pandas.MultiIndex.from_product([[tr], cnm])
                    auxdf = pandas.DataFrame(result[J[tr]:len(result)], index=idx)
                    r_store = pandas.concat([r_store, auxdf], axis=0)

        if V_type != "separate":
            w_constr_list = [v for v in w_constr_dict.values()]
            b = b_est_multi(A=A, Z=Z, J=J, KM=KM, iota=iota, w_constr=w_constr_list,
                            V=V)

            j_lb = 0
            for tr in df.treated_units:
                j_ub = j_lb + J[df.treated_units[0]]
                idx = pandas.MultiIndex.from_product([[tr], df.donors_dict[tr]])
                auxdf = pandas.DataFrame(b[j_lb:j_ub], index=idx)
                w_store = pandas.concat([w_store, auxdf], axis=0)
                j_lb = j_ub

            k_lb = Jtot
            for tr in df.treated_units:
                k_ub = k_lb + KM[df.treated_units[0]]
                if KM[tr] > 0:
                    cnm = [c.split("_", 1)[1] for c in C_dict[tr].columns.tolist()]
                    idx = pandas.MultiIndex.from_product([[tr], cnm])
                    auxdf = pandas.DataFrame(b[k_lb:k_ub], index=idx)
                    r_store = pandas.concat([r_store, auxdf], axis=0)
                    k_lb = k_ub

        b = pandas.concat([w_store, r_store], axis=0)
        b.index.rename(['ID', 'donor'], inplace=True)
    w_constr = w_constr_dict

    ##########################
    # Create useful objects
    # Store point estimates
    if KMI == 0:
        w = b
        r = pandas.DataFrame(None)
    else:
        w = b[0:Jtot]
        r = b[Jtot:]

    # Fitted values and residuals
    A_hat = Z.dot(numpy.array(b))
    A_hat.columns = A.columns
    res = A - A_hat

    # Pre-treatment fit of outcome of interest
    if class_type == "scpi_data":
        if out_in_features is True:   # outcome is among features
            fit_pre = A_hat.iloc[:T0_features[outcome_var], ].reset_index(level='feature', drop=True)
        else:  # outcome is not among features
            fit_pre = Y_donors.dot(w)

    elif class_type == "scpi_data_multi":
        if out_in_features[tr] is True:
            fit_pre = A_hat.loc[pandas.IndexSlice[:, outcome_var, :]]
        else:
            fit_pre = df.Y_donors.dot(w)

    fit_post = P.dot(numpy.array(b))
    fit_pre.columns = A.columns
    fit_post.columns = A.columns

    Y_pre_fit = deepcopy(fit_pre)
    Y_post_fit = deepcopy(fit_post)

    ############################################################################
    # Store actual values for the outcome pre and post when scdataMulti is used
    if class_type == "scpi_data_multi":
        Y_df = deepcopy(df.Y_df)
        Y_df.columns = ['ID', 'Time', 'Treatment', 'Actual']

        sel_units = [i in df.units_est for i in Y_df['ID']]
        res_df = deepcopy(Y_df[sel_units])

        # create hash maps if required
        if df.timeConvert is True:
            time_unique_ts = sorted(set(Y_df['Time'].tolist()))
            if df.effect != "time":
                int2ts = {i: time_unique_ts[i] for i in range(len(time_unique_ts))}
                ts2int = {time_unique_ts[i]: i for i in range(len(time_unique_ts))}
            else:
                ts2int = {time_unique_ts[i]: i + 2000 for i in range(len(time_unique_ts))}

            Y_df['Time'] = Y_df['Time'].map(ts2int)
            res_df = deepcopy(Y_df[sel_units])

            fit_pre.reset_index(drop=False, inplace=True)
            fit_pre['Time'] = fit_pre['Time'].map(ts2int)
            fit_pre.set_index(['ID', 'Time'], drop=True, inplace=True)

            if df.effect != "time":
                fit_post.reset_index(drop=False, inplace=True)
                fit_post['Time'] = fit_post['Time'].map(ts2int)
                fit_post.set_index(['ID', 'Time'], drop=True, inplace=True)

        synth_mat = pandas.concat([fit_pre, fit_post], axis=0)

        if df.effect != "time":
            synth_mat.index = synth_mat.index.rename(['ID', 'Time'])
            synth_mat.columns = ['Synthetic']

        if df.effect == "unit":
            Y_actual_pre = res_df[res_df['Treatment'] == 0]
            Y_actual_post = res_df[res_df['Treatment'] == 1]
            Y_actual_post_agg = Y_actual_post[['ID', 'Actual']].groupby(by='ID').mean()
            Y_actual_post_agg['Treatment'] = 1
            Y_actual_post_agg.set_index(fit_post.index, inplace=True)
            Y_actual_pre.set_index(['ID', 'Time'], append=False, inplace=True)
            Y_actual = pandas.concat([Y_actual_pre, Y_actual_post_agg], axis=0)
        else:
            Y_actual = deepcopy(res_df)
            Y_actual.set_index(['ID', 'Time'], drop=True, inplace=True)

        treated_periods = Y_actual.loc[Y_actual['Treatment'] == 1]
        treated_periods.reset_index(drop=False, inplace=True)
        treated_reception = treated_periods.groupby('ID').min()
        treated_reception.columns = ["Tdate", "Treatment", "Actual"]
        ant_df = pandas.DataFrame.from_dict(df.anticipation, orient='index')
        ant_df.index.rename('ID', inplace=True)
        ant_df.columns = ["anticipation"]
        treated_reception = treated_reception.merge(ant_df, on="ID")
        treated_reception['Tdate'] = treated_reception['Tdate'] - treated_reception['anticipation'] - 1 / 2

        treated_reception.reset_index(drop=False, inplace=True)
        sel_units = [i in df.units_est for i in treated_reception['ID']]
        treated_reception = treated_reception[sel_units]

        if df.effect != "time":
            toplot = pandas.concat([Y_actual, synth_mat], axis=1, join='inner')

        elif df.effect == "time":
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
            min_post = min([v for v in df.T1_outcome.values()]) - 1

            Y_pre_agg = Y_pre.groupby(['tstd'])[['Actual', 'Synthetic']].mean()
            Y_pre_agg.reset_index(inplace=True, drop=False)
            Y_pre_agg.columns = ['Time', 'Actual', 'Synthetic']
            Y_pre_agg = Y_pre_agg[Y_pre_agg['Time'] >= max_pre]

            Y_post_agg = Y_actual_post.groupby(['tstd'])[['Actual']].mean()
            Y_post_agg.reset_index(inplace=True, drop=False)
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

        if df.effect != "time" and df.timeConvert is True:
            Y_actual.reset_index(drop=False, inplace=True)
            Y_actual['Time'] = Y_actual['Time'].map(int2ts)
            Y_actual.set_index(['ID', 'Time'], drop=True, inplace=True)

    ##################################################
    # Plot

    if plot is True:
        if class_type == "scpi_data":
            to_plot = scest_output(b=b, w=w, r=r, Y_pre_fit=Y_pre_fit,
                                   Y_post_fit=Y_post_fit, A_hat=A_hat, res=res,
                                   V=V, w_constr=w_constr, A=A, B=B, C=C,
                                   P=P, P_diff=None, Y_pre=df.Y_pre, Y_post=df.Y_post,
                                   Y_donors=Y_donors, J=J, K=K, KM=KM,
                                   M=M, iota=iota, KMI=df.KMI,
                                   cointegrated_data=df.cointegrated_data,
                                   period_pre=df.period_pre, period_post=df.period_post,
                                   T0_features=T0_features, T1_outcome=df.T1_outcome,
                                   outcome_var=df.outcome_var, features=df.features,
                                   glob_cons=df.glob_cons, out_in_features=out_in_features,
                                   plotres=None, treated_units=df.treated_units,
                                   donors_dict={df.treated_units[0]: df.donors_units},
                                   units_est=df.units_est, anticipation=df.anticipation,
                                   effect="unit-time", timeConvert=df.timeConvert)
            plotres = scplot(result=to_plot)
        else:
            to_plot = scest_multi_output(b=b, w=w, r=r, Y_pre_fit=Y_pre_fit, Y_post_fit=Y_post_fit,
                                         Y_pre=df.Y_pre, Y_post=df.Y_post, Y_actual=Y_actual,
                                         A_hat=A_hat, res=res, V=V, w_constr=w_constr, A=A,
                                         B=B, C=C, P=P, P_diff=df.P_diff, Y_df=df.Y_df, Y_donors=df.Y_donors,
                                         J=J, K=K, KM=KM, M=M, iota=df.iota, KMI=df.KMI,
                                         cointegrated_data=df.cointegrated_data,
                                         period_pre=df.period_pre,
                                         period_post=df.period_post,
                                         T0_features=df.T0_features,
                                         T1_outcome=df.T1_outcome,
                                         features=df.features,
                                         outcome_var=df.outcome_var,
                                         glob_cons=df.glob_cons,
                                         out_in_features=df.out_in_features,
                                         plotres=None,
                                         donors_dict=df.donors_dict,
                                         treated_units=df.treated_units,
                                         effect=df.effect,
                                         units_est=df.units_est,
                                         anticipation=df.anticipation, timeConvert=df.timeConvert)
            plotres = scplotMulti(result=to_plot)
    else:
        plotres = None

    if class_type == "scpi_data":
        return scest_output(b=b,
                            w=w,
                            r=r,
                            Y_pre_fit=Y_pre_fit,
                            Y_post_fit=Y_post_fit,
                            A_hat=A_hat,
                            res=res,
                            V=V,
                            w_constr=w_constr,
                            A=A,
                            B=B,
                            C=C,
                            P=P,
                            P_diff=None,
                            Y_pre=df.Y_pre,
                            Y_post=df.Y_post,
                            Y_donors=Y_donors,
                            J=J,
                            K=K,
                            KM=KM,
                            M=M,
                            KMI=df.KMI,
                            iota=iota,
                            cointegrated_data=df.cointegrated_data,
                            period_pre=df.period_pre,
                            period_post=df.period_post,
                            T0_features=T0_features,
                            T1_outcome=df.T1_outcome,
                            outcome_var=df.outcome_var,
                            features=df.features,
                            glob_cons=df.glob_cons,
                            out_in_features=out_in_features,
                            plotres=plotres,
                            treated_units=df.treated_units,
                            donors_dict={df.treated_units[0]: df.donors_units},
                            units_est=df.units_est,
                            anticipation=df.anticipation,
                            effect="unit-time",
                            timeConvert=df.timeConvert)
    else:
        return scest_multi_output(b=b, w=w, r=r, Y_pre_fit=Y_pre_fit, Y_post_fit=Y_post_fit, Y_pre=df.Y_pre,
                                  Y_post=df.Y_post, Y_actual=Y_actual, A_hat=A_hat,
                                  res=res, V=V, w_constr=w_constr,
                                  A=A, B=B, C=C, P=P,
                                  P_diff=df.P_diff, Y_df=df.Y_df, Y_donors=df.Y_donors,
                                  J=J, K=K, KM=KM, M=M, iota=df.iota, KMI=df.KMI,
                                  cointegrated_data=df.cointegrated_data,
                                  period_pre=df.period_pre,
                                  period_post=df.period_post,
                                  T0_features=df.T0_features,
                                  T1_outcome=df.T1_outcome,
                                  features=df.features,
                                  outcome_var=df.outcome_var,
                                  glob_cons=df.glob_cons,
                                  out_in_features=df.out_in_features,
                                  plotres=plotres,
                                  donors_dict=df.donors_dict,
                                  treated_units=df.treated_units,
                                  effect=df.effect,
                                  units_est=df.units_est,
                                  anticipation=df.anticipation, timeConvert=df.timeConvert)


class scest_output:
    def __init__(self, b, w, r, Y_pre_fit, Y_post_fit, A_hat, res, V, w_constr,
                 A, B, C, P, P_diff, Y_pre, Y_post, Y_donors, J, K, KM, M, iota, KMI,
                 cointegrated_data, period_pre, period_post, T0_features,
                 T1_outcome, features, outcome_var, glob_cons, out_in_features,
                 plotres, treated_units, donors_dict, units_est, anticipation, effect, timeConvert):

        self.b = b
        self.w = w
        self.r = r
        self.Y_pre_fit = Y_pre_fit
        self.Y_post_fit = Y_post_fit
        self.A_hat = A_hat
        self.res = res
        self.V = V
        self.w_constr = w_constr
        self.A = A
        self.B = B
        self.C = C
        self.P = P
        self.P_diff = P_diff
        self.Y_pre = Y_pre
        self.Y_post = Y_post
        self.Y_donors = Y_donors
        self.J = J
        self.K = K
        self.KM = KM
        self.M = M
        self.iota = iota
        self.KMI = KMI
        self.cointegrated_data = cointegrated_data
        self.period_pre = period_pre
        self.period_post = period_post
        self.T0_features = T0_features
        self.T1_outcome = T1_outcome
        self.features = features
        self.outcome_var = outcome_var
        self.glob_cons = glob_cons
        self.out_in_features = out_in_features
        self.plotres = plotres
        self.treated_units = treated_units
        self.donors_dict = donors_dict
        self.units_est = units_est
        self.anticipation = anticipation
        self.effect = effect
        self.timeConvert = timeConvert

    def __repr__(self):

        # Prepare objects to print
        fw = 30
        if self.M == 1:
            fw_r = 35
        if self.M > 1:
            fw_r = 40
        w_constr = self.w_constr[self.treated_units[0]]

        constr = w_constr['name']
        if w_constr['Q'] is not None:
            Qsize = round(w_constr['Q'], 3)
        else:
            Qsize = "-"
        tr_unit = self.A.columns.values.tolist()
        pt_in = self.period_pre[0]
        pt_fi = self.period_pre[len(self.period_pre) - 1]
        ppre = str(pt_in) + '-' + str(pt_fi)

        Weights = self.w.rename(columns={0: 'Weights'}, inplace=False)
        Weights = round(Weights, 3)
        idx1 = Weights.index.get_level_values(0).tolist()
        idx2 = Weights.index.get_level_values(1).tolist()
        idx = pandas.MultiIndex.from_product([[idx1[0]], idx2], names=['Treated Unit', 'Donor'])
        Weights.set_index(idx, inplace=True)
        activew = len(Weights.loc[abs(Weights['Weights']) > 0])

        if self.KM > 0:
            Covariates = self.r.rename(columns={0: 'Covariates'}, inplace=False)
            Covariates = round(Covariates, 3)
            idx1 = Covariates.index.get_level_values(0).tolist()
            toparse = Covariates.index.get_level_values(1).tolist()
            parsed = [don.split("_", 1)[1] for don in toparse]
            idx = pandas.MultiIndex.from_product([[idx1[0]], parsed], names=['', ''])
            Covariates.set_index(idx, inplace=True)

        if constr is None:
            constr = "User Provided"

        # Print stuff
        print('-----------------------------------------------------------------------')
        print('Call: scest')
        print('Synthetic Control Estimation - Setup')
        print('')

        print('Constraint Type:'.ljust(fw), str(constr).rjust(fw_r))
        print('Constraint Size (Q):'.ljust(fw), str(Qsize).rjust(fw_r))
        print('Treated Unit:'.ljust(fw), str(idx1[0]).rjust(fw_r))
        print('Size of the donor pool:'.ljust(fw), str(self.J).rjust(fw_r))
        print('Features'.ljust(fw), str(self.M).rjust(fw_r))
        print('Pre-treatment period'.ljust(fw), str(ppre).rjust(fw_r))

        if self.M == 1:
            T0 = [c for c in self.T0_features.values()]

            print('Pre-treatment periods used in estimation:'.ljust(fw),
                  str(T0[0]).rjust(fw_r - 11))
            print('Covariates used for adjustment:'.ljust(fw),
                  str(self.KM).rjust(fw_r - 1))
        else:
            toprint = pandas.DataFrame.from_dict(self.T0_features, orient='index').reset_index()
            toprint.columns = ['Feature', 'Observations']
            print('Pre-treatment periods used in estimation per feature:'.ljust(fw))
            print(toprint.to_string(index=False))

            toprint = pandas.DataFrame(self.K)
            print('Covariates used for adjustment per feature:'.ljust(fw))
            print(toprint.to_string(index=False))

        print('')
        print('Synthetic Control Estimation - Results')
        print('')
        print('Active donors:', activew)
        print('')
        print('Coefficients:')
        print(Weights)
        if self.KM > 0:
            print('')
            print(Covariates)

        return ''

class scest_multi_output:
    def __init__(self, b, w, r, Y_pre_fit, Y_post_fit, Y_pre, Y_post, Y_actual, A_hat, res, V, w_constr,
                 A, B, C, P, P_diff, Y_df, Y_donors, J, K, KM, M, iota, KMI,
                 cointegrated_data, period_pre, period_post, T0_features,
                 T1_outcome, features, outcome_var, glob_cons, out_in_features, plotres,
                 donors_dict, treated_units, effect, units_est, anticipation, timeConvert):

        self.b = b
        self.w = w
        self.r = r
        self.Y_pre_fit = Y_pre_fit
        self.Y_post_fit = Y_post_fit
        self.Y_pre = Y_pre
        self.Y_post = Y_post
        self.Y_actual = Y_actual
        self.A_hat = A_hat
        self.res = res
        self.V = V
        self.w_constr = w_constr
        self.A = A
        self.B = B
        self.C = C
        self.P = P
        self.P_diff = P_diff
        self.Y_df = Y_df
        self.Y_donors = Y_donors
        self.J = J
        self.K = K
        self.KM = KM
        self.M = M
        self.iota = iota
        self.KMI = KMI
        self.cointegrated_data = cointegrated_data
        self.period_pre = period_pre
        self.period_post = period_post
        self.T0_features = T0_features
        self.T1_outcome = T1_outcome
        self.features = features
        self.outcome_var = outcome_var
        self.glob_cons = glob_cons
        self.out_in_features = out_in_features
        self.plotres = plotres
        self.donors_dict = donors_dict
        self.treated_units = treated_units
        self.effect = effect
        self.units_est = units_est
        self.anticipation = anticipation
        self.timeConvert = timeConvert

    def __repr__(self):

        return ''
