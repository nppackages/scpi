# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 14:46:14 2021

@author: Filippo Palomba
"""
# Temporary code to suppress pandas FutureWarning
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas
import numpy
import warnings
from copy import deepcopy
import multiprocessing as mp
from sklearn.linear_model import LinearRegression
from .funs import local_geom, u_des_prep, e_des_prep, complete_cases, df_EST, u_sigma_est, scpi_in, scpi_out
from .funs import executionTime, createPoly, DUflexGet
from .scest import scest
from .scplot import scplot
from scipy.linalg import sqrtm

def scpi(data,
         w_constr=None,
         V=None,
         P=None,
         u_missp=True,
         u_sigma="HC1",
         u_order=1,
         u_lags=0,
         u_design=None,
         u_alpha=0.05,
         e_method="all",
         e_order=1,
         e_lags=0,
         e_design=None,
         e_alpha=0.05,
         sims=200,
         rho=None,
         rho_max=None,
         cores=None,
         plot=False,
         w_bounds=None,
         e_bounds=None,
         opt_dict_est=None,
         opt_dict_inf=None,
         pass_stata=False):

    '''
    Parameters
    ----------
    data : scdata_output
        a class scdata_output object, obtained by calling scdata

    w_constr : dictionary
        a dictionary specifying the constraint set the estimated weights of the donors must belong to.
        w_constr can contain up to five objects:
        1. p, a scalar indicating the norm to be used (p should be one of "no norm", "L1", and "L2")
        2. dir, a string indicating whether the constraint on the norm is an equality ("==") or inequality ("<=")
        3. Q, a scalar defining the value of the constraint on the norm
        4. lb, a scalar defining the lower bound on the weights. It can be either 0 or -numpy.inf.
        5. name, a character selecting one of the default proposals.

    V : numpy.array, default numpy.identity
        an array specifying the weighting matrix to be used when minimizing the sum of squared residuals.
        The default is the identity matrix, so equal weight is given to all observations.

    P : numpy.array, default None
        a T_1 x (J+K_1) array containing the design matrix to be used to obtain the predicted.
        post-intervention outcome of the synthetic control unit. T_1 is the number of post-treatment periods,
        J is the size of the donor pool, and K_1 is the number of covariates used for adjustment in
        the outcome equation.

    u_missp : bool, default True
        a logical indicating if misspecification should be taken into account when dealing with u.

    u_sigma : str, default "HC1"
        a string specifying the type of variance-covariance estimator to be used when estimating
        the conditional variance of u. Available choices are HC0, HC1, HC2, and HC3.

    u_order : int, default 1
        an integer that sets the order of the polynomial in B when predicting moments of u.

    u_lags : int, default 0
        an integer that sets the number of lags of B when predicting moments of u.

    u_design : numpy.array, default None
        an array with the same number of rows of A and B and whose columns specify the design matrix
        to be used when modeling the estimated pseudo-true residuals u.

    u_alpha : float, default 0.05
        the confidence level for in-sample uncertainty.

    e_method : str, default "all"
        a string selecting the method to be used in quantifying out-of-sample uncertainty among:
        "gaussian" which uses conditional subgaussian bounds; "ls" which specifies a location-scale model for u; "qreg"
        which employs a quantile regressions to get the conditional bounds; "all" uses each one of the previous methods.

    e_order : int, default 1
        an integer that sets the order of the polynomial in B when predicting moments of e.

    e_lags: int, default 0
        a scalar that sets the number of lags of B when predicting moments of e.

    e_design : numpy.array, default None
        an array with the same number of rows of A and B and whose columns specify the design matrix
        to be used when modeling the estimated out-of-sample residuals e.

    e_alpha : float, default 0.05
        an integer specifying the confidence level for out-of-sample uncertainty.

    sims : int, default 200
        an integer providing the number of simulations to be used in quantifying in-sample uncertainty.

    rho : float/str, default 'type-1'
        a float specifying the regularizing parameter that imposes sparsity on the estimated vector of weights. If
        rho = 'type-1', then the tuning parameter is computed based on optimization inequalities. Other options are
        'type-2', and 'type-3'. See the software article for more information.

    rho_max : float, default 1
        a float indicating the maximum value attainable by the tuning parameter rho.

    cores : integer, default multiprocessing.cpu_count() - 1
        number of cores to be used by the command. The default is half the cores available.

    plot : bool, default False
        a logical specifying whether scplot should be called and a plot saved in the current working directory.
        For more options see scplot.

    w_bounds : numpy.array
        a T1 x 2 array with the user-provided bounds on beta. If w_bounds is provided, then
        the quantification of in-sample uncertainty is skipped. It is possible to provide only the lower bound or the
        upper bound by filling the other column with NAs.

    e_bounds : numpy.array
        a T1 x 2 array with the user-provided bounds on e. If e_bounds is provided, then
        the quantification of out-of-sample uncertainty is skipped. It is possible to provide only the lower bound or
        the upper bound by filling the other column with NAs.

    opt_dict_est : dictionary
        a dictionary specifying the stopping criteria used by the underling optimizer (nlopt) for point estimation.
        The default is a sequential quadratic programming (SQP) algorithm for nonlinearly constrained gradient-based
        optimization ('SLSQP'). In case a lasso-type constraint is implemented, cvxpy is used for optimization.
        More information on the stopping criteria can be obtained reading the official documentation at
        https://www.cvxpy.org/. The default values are 'maxeval = 5000', 'xtol_rel = 1e-8', 'xtol_abs = 1e-8',
        'ftol_rel = 1e-4', 'ftol_abs = 1e-4', 'tol_eq = 1e-8', and 'tol_ineq = 1e-8'.

    opt_dict_inf : dictionary
        same as above but for inference purposes. The default values are 'maxeval = 5000', 'xtol_rel = 1e-8',
        'xtol_abs = 1e-8', 'ftol_rel = 1e-4', 'ftol_abs = 1e-4', 'tol_eq = 1e-8', and 'tol_ineq = 1e-8'.

    pass_stat : bool
        for internal use only.

    Returns
    -------
    The function returns an object of class `scpi_output' containing the following objects

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

    CI_in_sample : pandas.DataFrame
        a dataframe containing the prediction intervals taking only in-sample uncertainty in to account.

    CI_all_gaussian : pandas.DataFrame
        a dataframe containing the prediction intervals taking out-of-sample uncertainty in to account.

    CI_all_ls : pandas.DataFrame
        a dataframe containing the prediction intervals taking out-of-sample uncertainty in to account.

    CI_all_qreg : pandas.DataFrame
        a dataframe containing the prediction intervals taking out-of-sample uncertainty in to account.

    Sigma : numpy.array
        an array containing the estimated variance-covariance Sigma.

    u_mean : numpy.array
        an array containing the estimated conditional mean of the the pseudo-residuals u.

    u_var : numpy.array
        an array containing the estimated conditional variance-covariance of the pseudo-residuals u.

    e_mean : numpy.array
        an array containing the estimated conditional mean of the out-of-sample error e.

    e_var : numpy.array
        an array containing the estimated conditional variance of the out-of-sample error e.

    u_missp : bool
        a logical indicating whether the model has been treated as misspecified or not.

    u_lags : int
        an integer containing the number of lags in B used in predicting moments of the pseudo-residuals u.

    u_order : int
        an integer containing the order of the polynomial in B used in predicting moments of the pseudo-residuals u.

    u_sigma : str
        a string indicating the estimator used for Sigma.

    u_user : bool
        a logical indicating whether the design matrix to predict moments of u was user-provided.

    u_alpha : float
        a float indicating the confidence level used for in-sample uncertainty.

    e_method : str
        a string indicating the specification used to predict moments of the out-of-sample error e.

    e_lags : int
        an integer containing the number of lags in B used in predicting moments of the pseudo-residuals u.

    e_order : int
        an integer containing the number of lags in B used in predicting moments of the pseudo-residuals u.

    e_user : bool
        a logical indicating whether the design matrix to predict moments of e was user-provided.

    e_alpha : float
        a float indicating the confidence level used for out-of-sample uncertainty.

    rho : str/float
        an integer specifying the estimated regularizing parameter that imposes sparsity on
        the estimated vector of weights.

    sims : int
        an integer indicating the number of simulations.

    failed_sims : numpy.array
        an array containing the number of failed simulations per post-treatment period to estimate
        lower and upper bounds.


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
    scdata, scest, splot

    '''

    if data.__class__.__name__ != 'scdata_output':
        raise Exception("data should be the object returned by running scdata!")

    ######################################
    # Estimation of synthetic weights
    if pass_stata is False:
        print("-----------------------------------------------")
        print("Estimating Weights...")
    sc_pred = scest(df=data, w_constr=w_constr, V=V, opt_dict=opt_dict_est)

    ######################################
    # Retrieve processed data from scest

    A = sc_pred.A                           # Features of treated unit
    B = sc_pred.B                           # Features of control units
    C = sc_pred.C                           # Covariates for adjustment
    Z = pandas.concat([B, C], axis=1)       # B and C column-bind
    Y_donors = sc_pred.Y_donors             # Outcome variable of control units
    K = sc_pred.K                           # Number of covs for adjustment per feature
    KM = sc_pred.KM                         # Dimension of r (total number of covs for adj)
    J = sc_pred.J                           # Number of donors
    M = sc_pred.M                           # Number of features
    T0 = sc_pred.T0_features                # Time periods used per feature
    T1 = sc_pred.T1_outcome                 # Number of out-of-sample periods
    outcome_var = sc_pred.outcome_var       # Name of outcome variable
    constant = sc_pred.glob_cons            # Logical indicating whether a constant is included
    out_feat = sc_pred.out_in_features      # Logical indicating whether the outcome variable is among features
    coig_data = sc_pred.cointegrated_data   # Logical indicating whether B is cointegrated
    w_constr = deepcopy(sc_pred.w_constr)   # Constraints on w
    V = sc_pred.V                           # Weighting matrix
    w = deepcopy(sc_pred.w)                 # Estimated vector of weights
    r = sc_pred.r                           # Estimated coefficients of covariates
    Y_post_fit = sc_pred.Y_post_fit         # Estimated post-treatment outcome for SC unit
    res = sc_pred.res                       # Residuals from estimation
    T0_tot = sum(T0.values())               # Total number of observations used in estimation

    if P is None:           # Matrix for out-of-sample prediction
        P = sc_pred.P
    else:                   # User-provided prediction matrix P (should be T1 by (J+KM))
        if not isinstance(P, (pandas.DataFrame, numpy.ndarray)):
            raise Exception("The object P should be a dataframe or an array!")

        if isinstance(P, numpy.ndarray):
            P = pandas.DataFrame(P,
                                 index=Y_post_fit.index,
                                 columns=Y_post_fit.columns)

        P_shape = numpy.shape(P)
        if P_shape[0] != T1:
            raise Exception("The object P has " + str(P_shape[0]) + " rows, when" +
                            " instead " + str(T1) + " were expected (i.e. the number of post-intervention periods)!")

        if P_shape[1] != (J + int(K.iloc[0, 0])):
            raise Exception("The object P has " + str(P_shape[1]) + " columns, when" +
                            " instead " + str(J + int(K.iloc[0, 0])) + " were expected (i.e. the size of the donor " +
                            "pool plus the number of covariates used in adjustment in the outcome equation)!")

        # Add zeros to avoid loading the coefficient of covs used for adj in other eqs
        if (KM - int(K.iloc[0, 0])) > 0:
            zeros = numpy.zeros((T1, sum(K.iloc[1:, ])))
            P = pandas.DataFrame(numpy.c_[P, zeros])

    if not isinstance(u_sigma, str):
        raise Exception("The object u_sigma should be of type character!")
    else:
        if u_sigma not in ['HC0', 'HC1', 'HC2', 'HC3']:
            raise Exception("Supported variance estimators are 'HC0','HC1','HC2','HC3'.")

    if w_bounds is not None:
        if not isinstance(w_bounds, (pandas.DataFrame, numpy.ndarray)):
            raise Exception("The object w_bounds should be a dataframe or an array!")

        if isinstance(w_bounds, numpy.ndarray):
            w_bounds = pandas.DataFrame(w_bounds,
                                        index=Y_post_fit.index,
                                        columns=Y_post_fit.columns)

        w_bounds_shape = numpy.shape(w_bounds)
        if w_bounds_shape[0] != len(Y_post_fit):
            raise Exception("w.bounds should be a matrix with two columns: the first column for the " +
                            "lower bound, the second for the upper. In case you don't want to specify " +
                            "the lower or the upper bound just fill the specific column with NAs.")

        if w_bounds_shape[1] != len(Y_post_fit):
            raise Exception("w.bounds should be a matrix with " + str(len(Y_post_fit))
                            + " rows (i.e. the number of post-intervention periods).")

    if sims < 10:
        raise Exception("The number of simulations needs to be larger or equal than 10!")

    if e_bounds is not None:
        if not isinstance(e_bounds, (pandas.DataFrame, numpy.ndarray)):
            raise Exception("The object e_bounds should be a dataframe or an array!")

        if isinstance(e_bounds, numpy.ndarray):
            e_bounds = pandas.DataFrame(e_bounds,
                                        index=Y_post_fit.index,
                                        columns=Y_post_fit.columns)

        e_bounds_shape = numpy.shape(e_bounds)
        if e_bounds_shape[0] != len(Y_post_fit):
            raise Exception("w.e_bounds should be a matrix with two columns: the first column for the " +
                            "lower bound, the second for the upper. In case you don't want to specify " +
                            "the lower or the upper bound just fill the specific column with NAs.")

        if e_bounds_shape[1] != len(Y_post_fit):
            raise Exception("e_bounds should be a matrix with " + str(len(Y_post_fit))
                            + " rows (i.e. the number of post-intervention periods).")

    # Check rho
    if rho is not None:
        if isinstance(rho, str):
            if rho not in ['type-1', 'type-2', 'type-3']:
                raise Exception("When not a scalar, 'rho' must be 'type-1', 'type-2', or 'type-3'.")
    else:
        rho = 'type-1'

    # Check on number of cores
    if cores is not None:
        n_cores = mp.cpu_count()
        if cores > n_cores:
            raise Exception("You selected " + str(cores) + " cores, but only " + str(n_cores) +
                            " are available on your machine!")
    else:
        cores = mp.cpu_count() - 1

    if pass_stata is False:
        print("Quantifying Uncertainty")
        executionTime(T0, T1, J, cores, sims, w_constr['name'])
        print(" ")

    ######################################
    ######################################
    # Estimate In-Sample Uncertainty
    ######################################
    ######################################

    # Regularize w
    w_constr_inf, w_star, index_w, rho, Q_star = local_geom(w_constr, rho,
                                                            rho_max, res, B, C, coig_data, T0_tot, J, w)

    # Create an index that selects all non-zero weights and additional covariates
    if isinstance(index_w, pandas.core.indexes.base.Index):
        index = deepcopy(index_w.tolist())
    else:
        index = deepcopy(index_w)
    index.extend(C.columns.tolist())
    beta = pandas.concat([w_star, r], axis=0)

    ##########################################################################
    # Prepare design matrix for in-sample uncertainty
    u_des_0 = u_des_prep(B, C, Z, u_order, u_lags, coig_data, T0_tot, M, constant, index,
                         index_w, u_design, res)

    ##########################################################################
    # Prepare design matrix for out_of-sample uncertainty
    e_res, e_des_0, e_des_1 = e_des_prep(B, C, Z, P, e_order, e_lags, res, sc_pred, Y_donors, out_feat, J, M, index,
                                         index_w, coig_data, T0, T0_tot, T1, outcome_var, constant, e_design)

    ###########################################################################
    # Remove NA - In sample uncertainty
    X = pandas.concat([A, res, u_des_0, Z], axis=1)
    tosel = complete_cases(X)
    X_na = X.loc[tosel, ]

    j2 = 1
    j3 = j2 + 1
    j4 = j3 + len(u_des_0.columns)

    res_na = X_na.iloc[:, j2:j3]
    u_des_0_na = X_na.iloc[:, j3:j4]
    Z_na = X_na.iloc[:, j4:]

    V_na = V[numpy.ix_(tosel, tosel)]

    # Effective number of observation used for inference (not yet adjusted for df used)
    TT = len(Z_na)

    # Remove NA - In sample uncertainty
    X = pandas.concat([e_res, e_des_0], axis=1)
    tosel = complete_cases(X)
    X_na = X.loc[tosel, ]

    e_res_na = X_na.iloc[:, 0]
    e_des_0_na = X_na.iloc[:, 1:]

    # Proceed cleaning missing data in the post-treatment period
    tosel = complete_cases(P)
    P_na = P.loc[tosel, ]

    ############################################################
    # Augment H with powers and interactions of B (not of C!!!)
    # (PolynomialFeatures is not designed to handle nans)
    u_des_0_na, e_des_0_na, e_des_1 = createPoly(u_order, e_order, index_w, u_des_0_na, e_des_0_na, e_des_1, out_feat)

    ###########################################################################
    # Estimate E[u|H], V[u|H], and Sigma

    # If the model is thought to be misspecified then E[u|H] is estimated
    if u_missp is True:
        u_des_0_na = DUflexGet(u_des_0_na, C)
        if len(u_des_0_na) <= len(u_des_0_na.columns):
            warnings.warn("Consider specifying a less complicated model for u. The number of observations used " +
                          "to parametrically predict moments is smaller than the number of covariates used." +
                          " Consider reducing either the number " +
                          "of lags (u_lags) or the order of the polynomial (u_order)!")
        u_mean = LinearRegression().fit(u_des_0_na, res_na).predict(u_des_0_na)

    elif u_missp is False:
        u_mean = 0

    # Use HC inference to estimate V[u|H]
    df = df_EST(w_constr=w_constr, w=w, A=A, B=B, J=J, KM=KM)

    Sigma, Omega = u_sigma_est(u_mean=u_mean, u_sigma=u_sigma,
                               res=res_na, Z=Z_na, V=V_na,
                               index=index, TT=TT, M=M, df=df)

    Sigma_root = sqrtm(Sigma).real

    # Auxiliary logical values to estimate bounds for w
    w_lb_est = True
    w_ub_est = True

    if w_bounds is not None:
        if sum(numpy.isnan(w_bounds[0])) == len(Y_post_fit):
            w_lb_est = False
        if sum(numpy.isnan(w_bounds[1])) == len(Y_post_fit):
            w_ub_est = False

    # Define constrained problem to be simulated
    if w_lb_est is True or w_ub_est is True:
        Q = numpy.array(Z_na).T.dot(V_na).dot(numpy.array(Z_na)) / TT
        b_arr = numpy.array(beta).flatten()

        lb = w_constr_inf['lb']
        dire = w_constr_inf['dir']
        p = w_constr_inf['p']
        QQ = w_constr_inf['Q']

        if p == "no norm":
            p_int = 0
        if p == "L1":
            p_int = 1
        if p == "L2":
            p_int = 2

        vsig = scpi_in(sims, b_arr, Sigma_root, Q, P_na, J, KM, w_lb_est,
                       w_ub_est, p, p_int, QQ, dire, lb, cores, opt_dict_inf, pass_stata)

    if w_lb_est is True:
        w_lb = numpy.nanquantile(vsig[:, :len(P_na)], q=u_alpha / 2, axis=0)
        w_lb = pandas.DataFrame(w_lb,
                                index=Y_post_fit.index,
                                columns=Y_post_fit.columns)
        fail_lb = 100 * numpy.sum(numpy.isnan(vsig[:, :len(P_na)]), axis=0) / sims

    else:
        w_lb = w_bounds[0]
        fail_lb = numpy.array([0] * len(w_bounds[0]))

    if w_ub_est is True:
        w_ub = numpy.nanquantile(vsig[:, len(P_na):], q=(1 - u_alpha / 2), axis=0)
        w_ub = pandas.DataFrame(w_ub,
                                index=Y_post_fit.index,
                                columns=Y_post_fit.columns)
        fail_ub = 100 * numpy.sum(numpy.isnan(vsig[:, len(P_na):]), axis=0) / sims
    else:
        w_ub = w_bounds[1]
        fail_ub = numpy.array([0] * len(w_bounds[1]))

    fail = numpy.c_[fail_lb, fail_ub]
    fail = numpy.array([fail_lb, fail_ub])
    failed_sims = pandas.DataFrame(fail,
                                   index=['lb', 'ub'])

    if fail.sum() > 0.1 * sims * vsig.shape[1]:
        warnings.warn("For some of the simulations used to quantify in-sample uncertainty the solution of " +
                      "the optimization problem was not found! We suggest inspecting the magnitude of this issue " +
                      "by consulting the percentage of simulations that failed contained in " +
                      "YOUR_SCPI_OBJECT_NAME.failed_sims." +
                      "In case the number of unsuccessful simulations is high, you might want to consider " +
                      "changing the stopping criteria of the algorithm through the option 'opt_list_inf'.")

    # PIs for w
    sc_l_0 = Y_post_fit + w_lb        # Left bound
    sc_r_0 = Y_post_fit + w_ub        # Left bound
    len_0 = sc_r_0 - sc_l_0           # Length

    ######################################
    ######################################
    # Estimate out-of-sample uncertainty
    ######################################
    ######################################

    sc_l_1 = sc_l_2 = sc_l_3 = sc_r_1 = sc_r_2 = sc_r_3 = None
    e_mean = e_var = numpy.nan
    # Auxiliary logical values to estimate bounds for e
    e_lb_est = True
    e_ub_est = True

    if e_bounds is not None:
        if sum(numpy.isnan(e_bounds[0])) == len(Y_post_fit):
            e_lb_est = False
        if sum(numpy.isnan(e_bounds[1])) == len(Y_post_fit):
            e_ub_est = False

    if e_lb_est is True or e_ub_est is True:
        if len(e_des_0_na) <= len(e_des_0_na.columns):
            warnings.warn("Consider specifying a less complicated model for e. The number of observations used " +
                          "to parametrically predict moments is smaller than the number of covariates used. Consider " +
                          "reducing either the number of lags (e_lags) or the order of the polynomial (e_order)!")

    if e_method == 'gaussian' or e_method == 'all':
        e_lb, e_ub, e_1, e_2 = scpi_out(y=e_res_na, x=e_des_0_na, preds=e_des_1,
                                        e_method="gaussian", alpha=e_alpha / 2,
                                        e_lb_est=e_lb_est, e_ub_est=e_ub_est)

        # Overwrite with user's input
        if e_lb_est is False:
            e_lb = e_bounds[0]
        if e_ub_est is False:
            e_ub = e_bounds[1]

        sc_l_1 = sc_l_0 + e_lb
        sc_r_1 = sc_r_0 + e_ub
        len_1 = sc_r_1 - sc_l_1

        e_mean = e_1
        e_var = e_2

    if e_method == 'ls' or e_method == 'all':
        e_lb, e_ub, e_1, e_2 = scpi_out(y=e_res_na, x=e_des_0_na, preds=e_des_1,
                                        e_method="ls", alpha=e_alpha / 2,
                                        e_lb_est=e_lb_est, e_ub_est=e_ub_est)

        # Overwrite with user's input
        if e_lb_est is False:
            e_lb = e_bounds[0]
        if e_ub_est is False:
            e_ub = e_bounds[1]

        sc_l_2 = sc_l_0 + e_lb
        sc_r_2 = sc_r_0 + e_ub
        len_2 = sc_r_2 - sc_l_2

    if e_method == 'qreg' or e_method == 'all':
        if e_order == 0:
            e_lb = numpy.quantile(e_res_na, q=e_alpha / 2)
            e_ub = numpy.quantile(e_res_na, q=1 - e_alpha / 2)
        else:
            e_lb, e_ub, e_1, e_2 = scpi_out(y=e_res_na, x=e_des_0_na, preds=e_des_1,
                                            e_method="qreg", alpha=e_alpha / 2,
                                            e_lb_est=e_lb_est, e_ub_est=e_ub_est)
        # Overwrite with user's input
        if e_lb_est is False:
            e_lb = e_bounds[0]
        if e_ub_est is False:
            e_ub = e_bounds[1]

        sc_l_3 = sc_l_0 + e_lb
        sc_r_3 = sc_r_0 + e_ub
        len_3 = sc_r_3 - sc_l_3

    ###############################################
    # Return objects
    CI_0 = pandas.concat([sc_l_0, sc_r_0, len_0], axis=1)
    CI_0.columns = ['Left Bound', 'Right Bound', 'Length']
    CI_0.index.rename('Time', inplace=True)

    if sc_l_1 is not None:
        CI_1 = pandas.concat([sc_l_1, sc_r_1, len_1], axis=1)
        CI_1.columns = ['Left Bound', 'Right Bound', 'Length']
        CI_1.index.rename('Time', inplace=True)
    else:
        CI_1 = None

    if sc_l_2 is not None:
        CI_2 = pandas.concat([sc_l_2, sc_r_2, len_2], axis=1)
        CI_2.columns = ['Left Bound', 'Right Bound', 'Length']
        CI_2.index.rename('Time', inplace=True)
    else:
        CI_2 = None

    if sc_l_3 is not None:
        CI_3 = pandas.concat([sc_l_3, sc_r_3, len_3], axis=1)
        CI_3.columns = ['Left Bound', 'Right Bound', 'Length']
        CI_3.index.rename('Time', inplace=True)
    else:
        CI_3 = None

    u_user = u_design is not None  # True if user provided the design matrix for in-sample inference
    e_user = e_design is not None  # True if user provided the design matrix for out-of-sample inference

    ##################################################
    # Plot

    if plot is True:
        to_plot = scpi_output(b=sc_pred.b,
                              w=sc_pred.w,
                              r=sc_pred.r,
                              Y_pre_fit=sc_pred.Y_pre_fit,
                              Y_post_fit=sc_pred.Y_post_fit,
                              A_hat=sc_pred.A_hat,
                              res=sc_pred.res,
                              V=sc_pred.V,
                              w_constr=sc_pred.w_constr,
                              w_constr_inf=w_constr_inf,
                              A=sc_pred.A,
                              B=sc_pred.B,
                              C=sc_pred.C,
                              P=P_na,
                              Y_pre=sc_pred.Y_pre,
                              Y_post=sc_pred.Y_post,
                              Y_donors=sc_pred.Y_donors,
                              J=sc_pred.J,
                              K=sc_pred.K,
                              KM=sc_pred.KM,
                              M=sc_pred.M,
                              cointegrated_data=sc_pred.cointegrated_data,
                              period_pre=sc_pred.period_pre,
                              period_post=sc_pred.period_post,
                              T0_features=sc_pred.T0_features,
                              T1_outcome=sc_pred.T1_outcome,
                              outcome_var=sc_pred.outcome_var,
                              features=sc_pred.features,
                              glob_cons=sc_pred.glob_cons,
                              out_in_features=sc_pred.out_in_features,
                              CI_in_sample=CI_0,
                              CI_all_gaussian=CI_1,
                              CI_all_ls=CI_2,
                              CI_all_qreg=CI_3,
                              Sigma=Sigma,
                              u_mean=u_mean,
                              u_var=Omega,
                              e_mean=e_mean,
                              e_var=e_var,
                              u_missp=u_missp,
                              u_lags=u_lags,
                              u_order=u_order,
                              u_sigma=u_sigma,
                              u_user=u_user,
                              e_method=e_method,
                              e_lags=e_lags,
                              e_order=e_order,
                              e_user=e_user,
                              rho=rho,
                              u_alpha=u_alpha,
                              e_alpha=e_alpha,
                              sims=sims,
                              failed_sims=failed_sims,
                              plotres=None)

        plotres = scplot(result=to_plot)
    else:
        plotres = None

    return scpi_output(b=sc_pred.b,
                       w=sc_pred.w,
                       r=sc_pred.r,
                       Y_pre_fit=sc_pred.Y_pre_fit,
                       Y_post_fit=sc_pred.Y_post_fit,
                       A_hat=sc_pred.A_hat,
                       res=sc_pred.res,
                       V=sc_pred.V,
                       w_constr=sc_pred.w_constr,
                       w_constr_inf=w_constr_inf,
                       A=sc_pred.A,
                       B=sc_pred.B,
                       C=sc_pred.C,
                       P=P_na,
                       Y_pre=sc_pred.Y_pre,
                       Y_post=sc_pred.Y_post,
                       Y_donors=sc_pred.Y_donors,
                       J=sc_pred.J,
                       K=sc_pred.K,
                       KM=sc_pred.KM,
                       M=sc_pred.M,
                       cointegrated_data=sc_pred.cointegrated_data,
                       period_pre=sc_pred.period_pre,
                       period_post=sc_pred.period_post,
                       T0_features=sc_pred.T0_features,
                       T1_outcome=sc_pred.T1_outcome,
                       outcome_var=sc_pred.outcome_var,
                       features=sc_pred.features,
                       glob_cons=sc_pred.glob_cons,
                       out_in_features=sc_pred.out_in_features,
                       CI_in_sample=CI_0,
                       CI_all_gaussian=CI_1,
                       CI_all_ls=CI_2,
                       CI_all_qreg=CI_3,
                       Sigma=Sigma,
                       u_mean=u_mean,
                       u_var=Omega,
                       e_mean=e_mean,
                       e_var=e_var,
                       u_missp=u_missp,
                       u_lags=u_lags,
                       u_order=u_order,
                       u_sigma=u_sigma,
                       u_user=u_user,
                       e_method=e_method,
                       e_lags=e_lags,
                       e_order=e_order,
                       e_user=e_user,
                       rho=rho,
                       u_alpha=u_alpha,
                       e_alpha=e_alpha,
                       sims=sims,
                       failed_sims=failed_sims,
                       plotres=plotres)


# Define class

class scpi_output:
    def __init__(self, b, w, r, Y_pre_fit, Y_post_fit, A_hat, res, V, w_constr, w_constr_inf,
                 A, B, C, P, Y_pre, Y_post, Y_donors, J, K, KM, M,
                 cointegrated_data, period_pre, period_post, T0_features,
                 T1_outcome, features, outcome_var, glob_cons, out_in_features,
                 CI_in_sample, CI_all_gaussian, CI_all_ls, CI_all_qreg, Sigma,
                 u_mean, u_var, e_mean, e_var, u_missp, u_lags, u_order,
                 u_sigma, u_user, e_method, e_lags, e_order, e_user,
                 rho, u_alpha, e_alpha, sims, failed_sims, plotres):

        self.b = b
        self.w = w
        self.r = r
        self.Y_pre_fit = Y_pre_fit
        self.Y_post_fit = Y_post_fit
        self.A_hat = A_hat
        self.res = res
        self.V = V
        self.w_constr = w_constr
        self.w_constr_inf = w_constr_inf
        self.A = A
        self.B = B
        self.C = C
        self.P = P
        self.Y_pre = Y_pre
        self.Y_post = Y_post
        self.Y_donors = Y_donors
        self.J = J
        self.K = K
        self.KM = KM
        self.M = M
        self.cointegrated_data = cointegrated_data
        self.period_pre = period_pre
        self.period_post = period_post
        self.T0_features = T0_features
        self.T1_outcome = T1_outcome
        self.features = features
        self.outcome_var = outcome_var
        self.glob_cons = glob_cons
        self.out_in_features = out_in_features
        self.CI_in_sample = CI_in_sample
        self.CI_all_gaussian = CI_all_gaussian
        self.CI_all_ls = CI_all_ls
        self.CI_all_qreg = CI_all_qreg
        self.Sigma = Sigma
        self.u_mean = u_mean
        self.u_var = u_var
        self.e_mean = e_mean
        self.e_var = e_var
        self.u_missp = u_missp
        self.u_lags = u_lags
        self.u_order = u_order
        self.u_sigma = u_sigma
        self.u_user = u_user
        self.e_method = e_method
        self.e_lags = e_lags
        self.e_order = e_order
        self.e_user = e_user
        self.rho = rho
        self.u_alpha = u_alpha
        self.e_alpha = e_alpha
        self.sims = sims
        self.failed_sims = failed_sims
        self.plotres = plotres

    def __repr__(self):

        # Prepare objects to print (estimation)
        fw = 30
        if self.M == 1:
            fw_r = 14
        if self.M > 1:
            fw_r = 40

        constr = self.w_constr['name']
        if self.w_constr['Q'] is not None:
            Qsize = round(self.w_constr['Q'], 3)
        else:
            Qsize = "-"
        tr_unit = self.A.columns.values.tolist()
        pt_in = self.period_pre[0]
        pt_fi = self.period_pre[len(self.period_pre) - 1]
        ppre = str(pt_in) + '-' + str(pt_fi)

        Weights = self.w.rename(columns={0: 'Weights'}, inplace=False)
        Weights = round(Weights, 3)
        activew = len(Weights.loc[abs(Weights['Weights']) > 0])

        if self.KM > 0:
            Covariates = self.r.rename(columns={0: 'Covariates'}, inplace=False)
            Covariates = round(Covariates, 3)

        if constr is None:
            constr = "User Provided"

        # Print stuff
        print('-----------------------------------------------------------------------')
        print('Call: scpi')
        print('Synthetic Control Estimation - Setup')
        print('')

        print('Constraint Type:'.ljust(fw), str(constr).rjust(fw_r))
        print('Constraint Size (Q):'.ljust(fw), str(Qsize).rjust(fw_r))
        print('Treated Unit:'.ljust(fw), str(tr_unit[0]).rjust(fw_r))
        print('Size of the donor pool:'.ljust(fw), str(self.J).rjust(fw_r))
        print('Features'.ljust(fw), str(self.M).rjust(fw_r))
        print('Pre-treatment period'.ljust(fw), str(ppre).rjust(fw_r))

        if self.M == 1:
            T0 = [c for c in self.T0_features.values()]

            print('Pre-treatment periods (used):'.ljust(fw), str(T0[0]).rjust(fw_r))

            print('Adjustment Covariates:'.ljust(fw), str(self.KM).rjust(fw_r))
        else:
            T0v = [c for c in self.T0_features.values()]
            T0n = [c for c in self.T0_features]
            Kv = [k for k in self.K.values]
            Kn = [k for k in self.K.index.get_level_values(0)]

            print('Pre-treatment periods used in estimation per feature:'.ljust(fw))
            for cov in range(len(T0n)):
                toprint = str(T0n[cov]) + ' ' + str(T0v[cov])
                print(toprint.rjust(fw_r))

            print('Covariates used for adjustment per feature:'.ljust(fw))
            for cov in range(len(Kn)):
                toprint = str(Kn[cov]) + ' ' + str(Kv[cov])
                print(toprint.rjust(fw_r))

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

        # Prepare objects to print (inference)
        e_method = self.e_method
        Y_tr_post = round(self.Y_post, 2)
        Y_tr_post = self.Y_post.astype(float).round(2)
        Y_sc_post = round(self.Y_post_fit, 2)
        Y_tr_post.columns = pandas.Index(['Treated'])
        Y_sc_post.columns = pandas.Index(['Synthetic'])

        if e_method == 'gaussian':
            CI = round(self.CI_all_gaussian.iloc[:, 0:2], 2)
        elif e_method == 'ls':
            CI = round(self.CI_all_ls.iloc[:, 0:2], 2)
        elif e_method == 'qreg':
            CI = self.CI_all_qreg.iloc[:, 0:2].astype(float).round(2)
        elif e_method == 'all':
            CI1 = round(self.CI_all_gaussian.iloc[:, 0:2], 2)
            CI2 = round(self.CI_all_ls.iloc[:, 0:2], 2)
            CI3 = self.CI_all_qreg.iloc[:, 0:2].astype(float).round(2)

        print('')
        print('-----------------------------------------------------------------------')
        print('Synthetic Control Inference - Setup')
        print('')

        if (self.u_user is False):
            print("In-sample Inference:                          ")
            print("     Misspecified model                       " + str(self.u_missp))
            print("     Order of polynomial (B)                  " + str(self.u_order))
            print("     Lags (B)                                 " + str(self.u_lags))
            print("     Variance-Covariance Estimator            " + self.u_sigma)
        else:
            print("In-sample Inference:")
            print("     User provided")

        if (self.e_user is False):
            print("Out-of-sample Inference:                          ")
            print("     Method                                   " + self.e_method)
            print("     Order of polynomial (B)                  " + str(self.e_order))
            print("     Lags (B)                                 " + str(self.e_lags))
        else:
            print("Out-of-sample Inference:")
            print("     User provided")

        if e_method == 'gaussian':
            print('   Inference with subgaussian bounds')
            dfprint = pandas.concat([Y_tr_post, Y_sc_post, CI], axis=1)
            dfprint.index.rename('Time', inplace=True)
            print(dfprint)
        elif e_method == 'ls':
            print('   Inference with location-scale model')
            dfprint = pandas.concat([Y_tr_post, Y_sc_post, CI], axis=1)
            dfprint.index.rename('Time', inplace=True)
            print(dfprint)
        elif e_method == 'qreg':
            print('   Inference with quantile regression')
            dfprint = pandas.concat([Y_tr_post, Y_sc_post, CI], axis=1)
            dfprint.index.rename('Time', inplace=True)
            print(dfprint)
        elif e_method == 'all':
            print('                              Subgaussian')
            dfprint = pandas.concat([Y_tr_post, Y_sc_post, CI1], axis=1)
            dfprint.index.rename('Time', inplace=True)
            print(dfprint)
            print('            Location Scale          Quantile Reg')
            print(pandas.concat([CI2, CI3], axis=1))

        else:
            print(pandas.concat([Y_tr_post, Y_sc_post, CI], axis=1).index.rename('Time'))

        if constr == "ridge":
            warnings.warn("The current version of the package does not take into account a small probability loss"
                          " in prediction intervals when using" +
                          "'ridge' as estimation method! An updated version will be rolled out soon!!")

        return ''
