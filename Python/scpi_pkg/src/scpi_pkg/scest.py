# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 11:08:57 2021

@author: Filippo Palomba
"""
# Temporary code to suppress pandas FutureWarning
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas
import numpy
from copy import deepcopy
from .funs import w_constr_prep, b_est, b_est_multi, V_prep, mat2dict
from .scplot import scplot


def scest(df, w_constr=None, V="separate", opt_dict=None, plot=False):

    '''

    Parameters
    ----------
    df : scdata_output
        a class scdata_output object, obtained by calling scdata

    w_constr : dictionary
        a dictionary specifying the constraint set the estimated weights of the donors must belong to.
        w_constr can contain up to four objects:
        - p, a string indicating the norm to be used (p should be one of "no norm", "L1", and "L2")
        - dir, a string indicating whether the constraint on the norm is an equality ("==") or inequality ("<=")
        - Q, a scalar defining the value of the constraint on the norm
        - lb, a scalar defining the lower bound on the weights. It can be either 0 or -numpy.inf.
        - name, a character selecting one of the default proposals.

    V : str/numpy.array, default "separate"
        a weighting matrix to be used when minimizing the sum of squared residuals.
        The default is the identity matrix ("separate"), so equal weight is given to all observations.
        The other possibility is to specify V = "pooled" for the pooled fit.

    opt_dict : dictionary
        a dictionary specifying the stopping criteria used by the underling optimizer (nlopt) for point estimation.
        The default is a sequential quadratic programming (SQP) algorithm for nonlinearly constrained gradient-based
        optimization ('SLSQP'). In case a lasso-type constraint is implemented, cvxpy is used for optimization.
        More information on the stopping criteria can be obtained reading the official documentation at
        https://www.cvxpy.org/. The default values are 'maxeval = 5000', 'xtol_rel = 1e-8', 'xtol_abs = 1e-8',
        'ftol_rel = 1e-8', 'ftol_abs = 1e-8', 'tol_eq = 1e-8', and 'tol_ineq = 1e-8'.

    plot : bool, default False
        a logical specifying whether scplot should be called and a plot saved in the current working directory. For more
        options see scplot.

    Returns
    -------
    The function returns an object of class `scest_output' containing the following objects

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
        a dataframe whose rows are the vectors used to predict the out-of-sample series for the synthetic unit.

    Y_pre : pandas.DataFrame
        a dataframe containing the pre-treatment outcome of the treated unit.

    Y_post : pandas.DataFrame
        a dataframe containing the post-treatment outcome of the treated unit.

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

    See Also
    --------
    scdata, scdataMulti, scpi, scplot, scplotMulti

    '''

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

            if not w_constr['name'] in ["simplex", "lasso", "ridge", "ols", "L1/L2"]:
                raise Exception("If 'name' is specified in w_constr, then it should be chosen" +
                                " among 'simplex', 'lasso', 'ridge', 'ols', 'L1/L2'.")
        else:
            w_constr['name'] = None

    if opt_dict is not None:
        if not isinstance(opt_dict, dict):
            raise Exception("opt_dict should be an object of type dict!")

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
    if not isinstance(V, (pandas.DataFrame, numpy.ndarray, str)):
        raise Exception("The object V should be a string, a dataframe or a matrix!")

    if isinstance(V, (pandas.DataFrame, numpy.ndarray)):
        V_shape = numpy.shape(V)
        if V_shape[0] != len(B) or V_shape[1] != len(B):
            raise Exception("The object V should be of shape (" + str(len(B)) +
                            "," + str(len(B)) + ")!")
        V_type = "separate"
        V = pandas.DataFrame(V, index=B.index,
                             columns=B.index.get_level_values('treated_unit'))
    else:
        V = V_prep(V_type, B, T0_features, iota)

    # Estimate SC
    if class_type == "scpi_data":
        w_constr = w_constr_prep(w_constr, w_names, A, Z, V, J, KM)
        result = b_est(A=A, Z=Z, J=J, KM=KM, w_constr=w_constr, V=V, opt_dict=opt_dict)
        cnms = [c.split("_", 1)[1] for c in Z.columns.tolist()]
        idx = pandas.MultiIndex.from_product([df.treated_units, cnms], names=['treated_unit', 'donor'])
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
                               V=V_i, opt_dict=opt_dict)
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
                            V=V, opt_dict=opt_dict)

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
        b.index.rename(['treated_unit', 'donor'], inplace=True)
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

    ##################################################
    # Plot

    if plot is True:
        if class_type == "scpi_data":
            to_plot = scest_output(b=b, w=w, r=r, Y_pre_fit=fit_pre,
                                   Y_post_fit=fit_post, A_hat=A_hat, res=res,
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
                                   effect="unit-time")
        else:
            top_plot = scest_multi_output(b=b, w=w, r=r, Y_pre_fit=fit_pre, Y_post_fit=fit_post,
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
                                          anticipation=df.anticipation)
        plotres = scplot(result=to_plot)

    else:
        plotres = None

    if class_type == "scpi_data":
        return scest_output(b=b,
                            w=w,
                            r=r,
                            Y_pre_fit=fit_pre,
                            Y_post_fit=fit_post,
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
                            effect="unit-time")
    else:
        return scest_multi_output(b=b,
                                  w=w,
                                  r=r,
                                  Y_pre_fit=fit_pre,
                                  Y_post_fit=fit_post,
                                  A_hat=A_hat,
                                  res=res,
                                  V=V,
                                  w_constr=w_constr,
                                  A=A,
                                  B=B,
                                  C=C,
                                  P=P,
                                  P_diff=df.P_diff,
                                  Y_df=df.Y_df,
                                  Y_donors=df.Y_donors,
                                  J=J,
                                  K=K,
                                  KM=KM,
                                  M=M,
                                  iota=df.iota,
                                  KMI=df.KMI,
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
                                  anticipation=df.anticipation)


class scest_output:
    def __init__(self, b, w, r, Y_pre_fit, Y_post_fit, A_hat, res, V, w_constr,
                 A, B, C, P, P_diff, Y_pre, Y_post, Y_donors, J, K, KM, M, iota, KMI,
                 cointegrated_data, period_pre, period_post, T0_features,
                 T1_outcome, features, outcome_var, glob_cons, out_in_features,
                 plotres, treated_units, donors_dict, units_est, anticipation, effect):

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
    def __init__(self, b, w, r, Y_pre_fit, Y_post_fit, A_hat, res, V, w_constr,
                 A, B, C, P, P_diff, Y_df, Y_donors, J, K, KM, M, iota, KMI,
                 cointegrated_data, period_pre, period_post, T0_features,
                 T1_outcome, features, outcome_var, glob_cons, out_in_features, plotres,
                 donors_dict, treated_units, effect, units_est, anticipation):

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

    def __repr__(self):

        return ''
