# -*- coding: utf-8 -*-

# Temporary code to suppress pandas FutureWarning
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from statsmodels.tools.sm_exceptions import IterationLimitWarning
warnings.filterwarnings("ignore", category=IterationLimitWarning)

import pandas
pandas.options.mode.chained_assignment = None
import numpy
from copy import deepcopy
import multiprocessing as mp
from scipy.linalg import block_diag
from sklearn.linear_model import LinearRegression
from .funs import local_geom, u_des_prep, e_des_prep, complete_cases
from .funs import df_EST, u_sigma_est, scpi_in, scpi_out, epskappaGet
from .funs import executionTime, createPoly, DUflexGet, mat2dict
from .funs import localgeom2step, detectConstant, simultaneousPredGet
from .scest import scest
from .scplot import scplot
from .scplotMulti import scplotMulti
from scipy.linalg import sqrtm

def scpi(data,
         w_constr=None,
         V="separate",
         Vmat=None,
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
         lgapp="generalized",
         cores=None,
         plot=False,
         w_bounds=None,
         e_bounds=None,
         verbose=True,
         pass_stata=False):

    """
    Parameters
    ----------
    data : scdata_output
        a class scdata_output object, obtained by calling scdata

    w_constr : dict, default {"name": "simplex"}
        a dictionary specifying the constraint set the estimated weights of the donors must belong to.
        w_constr can contain up to five objects:
        1. p, a scalar indicating the norm to be used (p should be one of "no norm", "L1", and "L2")
        2. dir, a string indicating whether the constraint on the norm is an equality ("==") or inequality ("<=")
        3. Q, a scalar defining the value of the constraint on the norm
        4. lb, a scalar defining the lower bound on the weights. It can be either 0 or -numpy.inf.
        5. name, a character selecting one of the default proposals.

    V : str, default "separate"
        a weighting matrix to be used when minimizing the sum of squared residuals.
        The default is the identity matrix ("separate"), so equal weight is given to all observations.
        The other possibility is to specify V = "pooled" for the pooled fit.

    Vmat : numpy.array, defaul None
        a conformable weighting matrix to be used in the minimization of the sum of squared residuals. To check the proper
        dimensions, we suggest to check the output of scdata or scdataMulti and inspect the dimensions of B and C.

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
        If there is risk of over-fitting the command automatically sets u_order = 0.

    u_lags : int, default 0
        an integer that sets the number of lags of B when predicting moments of u. If there is risk of over-fitting the
        command automatically sets u_lags=0.

    u_design : numpy.array, default None
        an array with the same number of rows of A and B and whose columns specify the design matrix
        to be used when modeling the estimated pseudo-true residuals u.

    u_alpha : float, default 0.05
        a float specifying the confidence level for in-sample uncertainty, i.e. 1 - u_alpha is the confidence level.

    e_method : str, default "all"
        a string selecting the method to be used in quantifying out-of-sample uncertainty among:
        "gaussian" which uses conditional subgaussian bounds; "ls" which specifies a location-scale model for u; "qreg"
        which employs a quantile regressions to get the conditional bounds; "all" uses each one of the previous methods.

    e_order : int, default 1
        an integer that sets the order of the polynomial in B when predicting moments of e.
        If there is risk of over-fitting the command automatically sets e_order=0.

    e_lags: int, default 0
        a scalar that sets the number of lags of B when predicting moments of e. If there is risk of over-fitting the
        command automatically sets e_lags=0.

    e_design : numpy.array, default None
        an array with the same number of rows of A and B and whose columns specify the design matrix
        to be used when modeling the estimated out-of-sample residuals e.

    e_alpha : float, default 0.05
        an float specifying the confidence level for out-of-sample uncertainty, i.e. 1- e_alpha is the confidence level.

    sims : int, default 200
        an integer providing the number of simulations to be used in quantifying in-sample uncertainty.

    rho : float/str, default 'type-1'
        a float specifying the regularizing parameter that imposes sparsity on the estimated vector of weights. If
        rho = 'type-1', then the tuning parameter is computed based on optimization inequalities. Other options are
        'type-2', and 'type-3'. See the software article for more information.

    rho_max : float, default 1
        a float indicating the maximum value attainable by the tuning parameter rho.

    lgapp : str, default "generalized"
        selects the way local geometry is approximated in simulation. The options are "generalized"
        and "linear". The first one accommodates for possibly non-linear constraints, whilst the second one is valid
        with linear constraints only.

    cores : integer, default multiprocessing.cpu_count() - 1
        number of cores to be used by the command. The default is half the cores available.

    plot : bool, default False
        a logical specifying whether scplot should be called and a plot saved in the current working directory.
        For more options see scplot.

    w_bounds : numpy.array, default None
        a T1 x 2 array with the user-provided bounds on beta. If w_bounds is provided, then
        the quantification of in-sample uncertainty is skipped. It is possible to provide only the lower bound or the
        upper bound by filling the other column with NAs.

    e_bounds : numpy.array, default None
        a T1 x 2 array with the user-provided bounds on e. If e_bounds is provided, then
        the quantification of out-of-sample uncertainty is skipped. It is possible to provide only the lower bound or
        the upper bound by filling the other column with NAs.

    verbose : bool, default True
        if False prevents printing additional information in the console.

    pass_stat : bool, default False
        for internal use only.

    Returns
    -------
    The function returns an object of class 'scpi_output' containing the following objects

    w : pandas.DataFrame
        a dataframe containing the weights of the donors.

    r : pandas.DataFrame
        a dataframe containing the values of the covariates used for adjustment.

    b : pandas.DataFrame
        a dataframe containing w and r.

    Y_pre_fit : pandas.DataFrame
        a dataframe containing the estimated pre-treatment outcome for the SC unit(s).

    Y_post_fit : pandas.DataFrame
        a dataframe containing the estimated post-treatment outcome for the SC unit(s).

    Y_pre: pandas.DataFrame
        a dataframe containing the actual pre-treatment outcome for the treated unit(s).

    Y_post: pandas.DataFrame
        a dataframe containing the actual post-treatment outcome for the treated unit(s).

    A_hat : pandas.DataFrame
        a dataframe containing the predicted values of the features of the treated unit(s).

    res : pandas.DataFrame
        a dataframe containing the residuals A - A_hat.

    V : numpy.array
        an array containing the weighting matrix used in estimation.

    w_constr : dictionary
        a dictionary containing the specifics of the constraint set used on the weights.

    A : pandas.DataFrame
        a dataframe containing pre-treatment features of the treated unit(s).

    B : pandas.DataFrame
        a dataframe containing pre-treatment features of the control units.

    C : pandas.DataFrame
        a dataframe containing covariates for adjustment.

    P : pandas.DataFrame
        a dataframe whose rows are the vectors used to predict the out-of-sample series for the synthetic unit(s).

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

    bounds : dict
        a dictionary containing all the estimated bounds.

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
        a float specifying the confidence level used for in-sample uncertainty, i.e. 1-u_alpha is the confidence level.

    u_T : int
        an integer indicating the number of observations used to estimate (conditional) moments of the pseudo-residuals.

    u_params : int
        an integer indicating the number of parameters used to estimate (conditional) moments of the pseudo-residuals.

    u_D : array
        the design matrix used to predict moments of the pseudo-residuals.

    e_method : str
        a string indicating the specification used to predict moments of the out-of-sample error e.

    e_lags : int
        an integer containing the number of lags in B used in predicting moments of the pseudo-residuals u.

    e_order : int
        an integer containing the number of lags in B used in predicting moments of the pseudo-residuals u.

    e_user : bool
        a logical indicating whether the design matrix to predict moments of e was user-provided.

    e_T : int
        an integer indicating the number of observations used to estimate (conditional) moments of the out-of-sample
        error.

    e_params : int
        an integer indicating the number of parameters used to estimate (conditional) moments of the out-of-sample
        error.

    e_D : array
        the design matrix used to predict moments of the out-of-sample error.

    e_alpha : float
        a float indicating the confidence level used for out-of-sample uncertainty, i.e. 1-e_alpha is the confidence
        level.

    rho : str/float
        an integer specifying the estimated regularizing parameter that imposes sparsity on
        the estimated vector of weights.

    Q_star : dict
        a dictionary containing the regularized constraint on the norm of the weights.

    epskappa : pandas.DataFrame
        a dataframe containing the estimates for epsilon_kappa.

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

    Cattaneo, M. D., Palomba, F., Feng, Y., and Titiunik, R. (2022), “Uncertainty Quantification in Synthetic
    Controls with Staggered Treatment Adoption”.

    See Also
    --------
    scdata, scdataMulti, scest, splot, scplotMulti

    """

    if data.__class__.__name__ not in ['scdata_output', 'scdata_multi_output']:
        raise Exception("data should be the object returned by running scdata or scdataMulti!")

    if data.__class__.__name__ == "scdata_output":
        class_type = "scpi_data"
    else:
        class_type = "scpi_data_multi"

    ######################################
    # Estimation of synthetic weights
    if pass_stata is False and verbose:
        print("-----------------------------------------------")
        print("Estimating Weights...")
    sc_pred = scest(df=data, w_constr=w_constr, V=V, Vmat=Vmat)

    ######################################
    # Retrieve processed data from scest

    A = deepcopy(sc_pred.A)                 # Features of treated unit
    B = deepcopy(sc_pred.B)                 # Features of control units
    C = deepcopy(sc_pred.C)                 # Covariates for adjustment
    Z = pandas.concat([B, C], axis=1)       # B and C column-bind
    Y_donors = sc_pred.Y_donors             # Outcome variable of control units
    K = sc_pred.K                           # Number of covs for adjustment per feature
    KM = sc_pred.KM                         # Total number of covariates per treated unit
    J = sc_pred.J                           # Number of donors
    M = sc_pred.M                           # Number of features
    T0 = sc_pred.T0_features                # Time periods used per feature
    T1 = sc_pred.T1_outcome                 # Number of out-of-sample periods
    outcome_var = sc_pred.outcome_var       # Name of outcome variable
    w_constr = deepcopy(sc_pred.w_constr)   # Constraints on w
    V = sc_pred.V                           # Weighting matrix
    w = deepcopy(sc_pred.w)                 # Estimated vector of weights
    r = sc_pred.r                           # Estimated coefficients of covariates
    Y_post_fit = sc_pred.Y_post_fit         # Estimated post-treatment outcome for SC unit
    res = sc_pred.res                       # Residuals from estimation
    sc_effect = sc_pred.effect              # Causal quantity of interest
    tr_units = sc_pred.treated_units

    if class_type == "scpi_data":
        Jtot = J
        KMI = KM
        iota = 1
        T0_tot = sum(T0.values())
        T0_M = {tr_units[0]: T0_tot}
        T1_tot = T1
        features = {tr_units[0]: sc_pred.features}
        coig_data = {tr_units[0]: sc_pred.cointegrated_data}
        out_feat = {tr_units[0]: sc_pred.out_in_features}
        constant = {tr_units[0]: sc_pred.glob_cons}
        anticipation = {tr_units[0]: sc_pred.anticipation}
        T0 = {tr_units[0]: T0}
        T1 = {tr_units[0]: T1}
        J = {tr_units[0]: J}
        KM = {tr_units[0]: KM}
        M = {tr_units[0]: M}

    elif class_type == "scpi_data_multi":
        Jtot = sum(J.values())
        T0_M = {}
        for n, v in sc_pred.T0_features.items():
            T0_M[n] = sum(v.values())
        KMI = sc_pred.KMI
        iota = sc_pred.iota
        T0_tot = sum(T0_M.values())
        T1_tot = sum(T1.values())
        out_feat = sc_pred.out_in_features
        coig_data = sc_pred.cointegrated_data
        constant = sc_pred.glob_cons
        anticipation = sc_pred.anticipation

    if sc_effect == "unit":
        T1_tot = iota
        T1 = dict(zip(tr_units, [1] * iota))

    if P is None:           # Matrix for out-of-sample prediction
        P = deepcopy(sc_pred.P)
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

    # if sims < 10:
        # raise Exception("The number of simulations needs to be larger or equal than 10!")

    if w_bounds is not None:
        if not isinstance(w_bounds, (pandas.DataFrame, numpy.ndarray)):
            raise Exception("The object w_bounds should be a dataframe or an array!")

        if not isinstance(w_bounds, numpy.ndarray):
            w_bounds = numpy.array(w_bounds)

        w_bounds_shape = w_bounds.shape
        if w_bounds_shape[1] != 2:
            raise Exception("w.bounds should be a matrix with two columns: the first column for the " +
                            "lower bound, the second for the upper. In case you don't want to specify " +
                            "the lower or the upper bound just fill the specific column with NAs.")

        if w_bounds_shape[0] != len(Y_post_fit):
            raise Exception("w.bounds should be a matrix with " + str(len(Y_post_fit))
                            + " rows (i.e. the number of post-intervention periods).")

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

    if pass_stata is False and verbose:
        constr_type = w_constr[tr_units[0]]['name']
        print("Quantifying Uncertainty")
        executionTime(T0_tot, T1_tot, Jtot, iota, cores, sims, constr_type)
        print(" ")

    # create lists of matrices
    A_dict = mat2dict(A, cols=False)
    B_dict = mat2dict(B)
    Z_dict = mat2dict(Z)
    if len(C.columns) > 0:
        C_dict = mat2dict(C)
    else:
        C_dict = dict(zip(tr_units, [pandas.DataFrame(None)] * iota))

    if sc_effect == "time":
        P_dict = {}
        for tr in tr_units:
            csel = [c.split("_")[0] == tr for c in P.columns.tolist()]
            X_rc = P.loc[:, numpy.array(csel)]
            P_dict[tr] = X_rc
    else:
        P_dict = mat2dict(P)

    if sc_pred.P_diff is not None:
        Pd_dict = mat2dict(sc_pred.P_diff)
    else:
        Pd_dict = dict(zip(tr_units, [None] * iota))

    V_dict = mat2dict(V)
    w_dict = mat2dict(w, cols=False)
    res_dict = mat2dict(res, cols=False)
    Yd_dict = mat2dict(Y_donors)

    ############################################################################
    ############################################################################
    # Estimate In-Sample Uncertainty
    ############################################################################
    ############################################################################

    w_star = pandas.DataFrame(None)
    iw_dict = {}
    rho_dict = {}
    Q_star = {}
    Q2_star = {}
    w_constr_inf = {}
    u_des_0_na_dict = {}
    e_des_0_na_dict = {}
    e_des_1_dict = {}
    e_res_na_dict = {}
    res_na_dict = {}
    Z_na_dict = {}
    P_na_dict = {}
    w_constr_aux = deepcopy(w_constr)
    col_order = Z.columns.tolist()

    for tr in tr_units:
        # Regularize w
        w_constr_inf[tr], ws, iw_dict[tr], rho_dict[tr], Q_star[tr], Q2_star[tr] = local_geom(w_constr_aux[tr], rho,
                                                                                              rho_max, res_dict[tr],
                                                                                              B_dict[tr], C_dict[tr],
                                                                                              coig_data[tr], T0_M[tr],
                                                                                              J[tr], w_dict[tr],
                                                                                              verbose)

        w_star = pandas.concat([w_star, ws.reset_index(drop=True)], axis=0, ignore_index=True)
        index_i = iw_dict[tr].tolist() + [True] * KM[tr]
        lenfeat = [obs <= u_lags for obs in T0[tr].values()]
        if (sum(lenfeat) and u_lags > 0):
            if verbose:
                warnings.warn("At least one of your features is observed for less periods than the number of lags," +
                              " u.lags reverted to 0.")
            u_lags = 0

        ##########################################################################
        # Prepare design matrix for in-sample uncertainty
        ud0 = u_des_prep(B_dict[tr], C_dict[tr], u_order, 1, coig_data[tr],
                         T0_M[tr], M[tr], constant[tr], index_i, iw_dict[tr],
                         u_design, res_dict[tr])

        ##########################################################################
        # Prepare design matrix for out_of-sample uncertainty
        er, ed0, ed1 = e_des_prep(B_dict[tr], C_dict[tr], P_dict[tr], e_order,
                                  e_lags, res_dict[tr], sc_pred, Yd_dict[tr],
                                  out_feat[tr], J[tr], index_i, iw_dict[tr],
                                  coig_data[tr], T0[tr][outcome_var], T1[tr],
                                  constant[tr], e_design, outcome_var, Pd_dict[tr],
                                  sc_effect, iota)

        ###########################################################################
        # Remove NA - In sample uncertainty
        X = pandas.concat([A_dict[tr], res_dict[tr], ud0, Z_dict[tr]], axis=1)
        tosel = complete_cases(X)
        X_na = X.loc[tosel, ]

        j2 = 1
        j3 = j2 + 1
        j4 = j3 + len(ud0.columns)

        r_na = X_na.iloc[:, j2:j3]
        ud0_na = X_na.iloc[:, j3:j4]
        Z_na = X_na.iloc[:, j4:]

        # Remove NA - Out of sample uncertainty
        X = pandas.concat([er, ed0], axis=1)
        tosel = complete_cases(X)
        X_na = X.loc[tosel, ]

        er_na = X_na.iloc[:, [0]]
        ed0_na = X_na.iloc[:, 1:]

        # Proceed cleaning missing data in the post-treatment period
        tosel = complete_cases(P_dict[tr])
        P_na = P_dict[tr].loc[tosel, ]

        ############################################################
        # Augment H with powers and interactions of B (not of C!!!)
        # (PolynomialFeatures is not designed to handle nans)
        ud0_na, ed0_na, ed1 = createPoly(u_order, e_order, iw_dict[tr],
                                         ud0_na, ed0_na, ed1, out_feat[tr])

        ud0_na.insert(0, "ID", tr)
        r_na.insert(0, "ID", tr)
        er_na.insert(0, "ID", tr)
        ed0_na.insert(0, "ID", tr)
        if sc_effect == "time":
            idx = pandas.MultiIndex.from_product([[tr], ed1.index.get_level_values(0).tolist()],
                                                 names=['ID', 'Time'])
            ed1.set_index(idx, inplace=True, append=False)
        else:
            ed1.insert(0, "ID", tr)
            ed1.set_index('ID', append=True, drop=True, inplace=True)
            P_na.insert(0, "ID", tr)
            P_na.set_index('ID', append=True, drop=True, inplace=True)

        Z_na.insert(0, "ID", tr)
        ud0_na.set_index('ID', append=True, drop=True, inplace=True)
        r_na.set_index('ID', append=True, drop=True, inplace=True)
        er_na.set_index('ID', append=True, drop=True, inplace=True)
        ed0_na.set_index('ID', append=True, drop=True, inplace=True)
        Z_na.set_index('ID', append=True, drop=True, inplace=True)
        u_des_0_na_dict[tr] = ud0_na
        res_na_dict[tr] = r_na
        e_res_na_dict[tr] = er_na
        e_des_0_na_dict[tr] = ed0_na
        e_des_1_dict[tr] = ed1
        Z_na_dict[tr] = Z_na
        P_na_dict[tr] = P_na

    # Transform dictionaries to block diagonal matrices and fill NaNs
    tr_count = 1
    for tr in tr_units:
        if tr_count == 1:
            Z_na = deepcopy(Z_na_dict[tr])
            P_na = deepcopy(P_na_dict[tr])
            u_des_0_na = deepcopy(u_des_0_na_dict[tr])
            e_des_0_na = deepcopy(e_des_0_na_dict[tr])
            e_des_1 = deepcopy(e_des_1_dict[tr])
            res_na = deepcopy(res_na_dict[tr])
            e_res_na = deepcopy(e_res_na_dict[tr])
            index_w = deepcopy(iw_dict[tr])
        else:
            Z_na = pandas.concat([Z_na, Z_na_dict[tr]], axis=1)
            P_na = pandas.concat([P_na, P_na_dict[tr]], axis=1)
            u_des_0_na = pandas.concat([u_des_0_na, u_des_0_na_dict[tr]], axis=1)
            e_des_0_na = pandas.concat([e_des_0_na, e_des_0_na_dict[tr]], axis=1)
            e_des_1 = pandas.concat([e_des_1, e_des_1_dict[tr]], axis=1)
            e_res_na = pandas.concat([e_res_na, e_res_na_dict[tr]], axis=0)
            res_na = pandas.concat([res_na, res_na_dict[tr]], axis=0)
            index_w = numpy.append(index_w, iw_dict[tr])
        tr_count = tr_count + 1

    Z_na.fillna(0, inplace=True)
    P_na.fillna(0, inplace=True)
    u_des_0_na.fillna(0, inplace=True)
    e_des_0_na.fillna(0, inplace=True)
    e_des_1.fillna(0, inplace=True)
    e_res_na.fillna(0, inplace=True)
    res_na.fillna(0, inplace=True)
    u_des_0_na = u_des_0_na.reorder_levels(['ID', 'feature', 'Time'])
    e_des_0_na = e_des_0_na.reorder_levels(['ID', 'Time'])
    e_des_1 = e_des_1.reorder_levels(['ID', 'Time'])
    e_res_na = e_res_na.reorder_levels(['ID', 'Time'])
    res_na = res_na.reorder_levels(['ID', 'feature', 'Time'])
    Z_na = Z_na.reorder_levels(['ID', 'feature', 'Time'])
    if sc_effect != "time":
        P_na = P_na.reorder_levels(['ID', 'Time'])

    Z_na = Z_na[col_order]
    P_na = P_na[col_order]

    # Effective number of observation used for inference (not yet adjusted for df used)
    TT = len(Z_na)

    V_na = V.loc[V.index.isin(Z_na.index), V.index.isin(Z_na.index)]

    # Create an index that selects all non-zero weights and additional covariates
    if isinstance(index_w, numpy.ndarray):
        index = deepcopy(index_w.tolist())
    else:
        index = deepcopy(index_w)
    index.extend([True] * KMI)

    Q = {}
    for tr in tr_units:
        Q[tr] = w_constr_aux[tr]['Q']

    if lgapp == "generalized":  # we use rho only to impose sparsity on B when predicting moments
        beta = sc_pred.b
        Q_star, lb = localgeom2step(w, r, rho_dict, w_constr, Q, tr_units)

    elif lgapp == "linear":  # we use rho to regularize w too
        beta = pandas.concat([w_star, r], axis=0).set_index(sc_pred.b.index)
        Q_star = deepcopy(Q)
        lb = []
        for tr in tr_units:
            lb = lb + [w_constr_inf[tr]['lb']] * J[tr]

    ###########################################################################
    # Estimate E[u|H], V[u|H], and Sigma

    # If the model is thought to be misspecified then E[u|H] is estimated
    if u_missp is True:
        T_u = len(u_des_0_na)
        u_des_dict = mat2dict(u_des_0_na)
        u_des_0_flex = pandas.DataFrame(None)
        u_des_0_noflex = pandas.DataFrame(None)

        for tr in tr_units:
            udict = DUflexGet(u_des_dict[tr], C_dict[tr])
            udict.insert(0, "ID", tr)
            udict.set_index('ID', append=True, drop=True, inplace=True)
            u_des_0_flex = pandas.concat([u_des_0_flex, udict], axis=1)
            udict = u_des_dict[tr]
            udict.insert(0, "ID", tr)
            udict.set_index('ID', append=True, drop=True, inplace=True)
            u_des_0_noflex = pandas.concat([u_des_0_noflex, udict], axis=1)

        df_U = T_u - 10

        u_simple = df_U <= len(u_des_0_noflex.columns)
        u_noflex = (df_U > len(u_des_0_noflex.columns)) and (df_U <= len(u_des_0_flex.columns))
        u_flex = df_U > len(u_des_0_flex.columns)

        if u_simple is True:
            if verbose and (u_order > 0 or u_lags > 0):
                warnings.warn("One of u_order > 0 and  u_lags > 0 was specified, however the current number of " +
                              "observations (" + str(T_u) + ") used to estimate conditional moments of the " +
                              "pseudo-residuals is not larger than the number of parameters used in " +
                              "estimation (" + str(len(u_des_0_flex.columns)) + ") plus 10. " +
                              "To avoid over-fitting issues u_order and u_lags were set to 0.")
            u_des_0_na = pandas.DataFrame([1] * len(u_des_0_na), index=u_des_0_na.index)
            u_order = 0
            u_lags = 0

        elif u_noflex is True:
            if verbose and (u_order > 0 or u_lags > 0):
                warnings.warn("The current number of observations (" + str(T_u) + ") used to estimate conditional " +
                              "moments of the pseudo-residuals is not larger than the number of parameters used in " +
                              "estimation (" + str(len(u_des_0_flex.columns)) + ") plus 10 when allowing for a " +
                              "feature specific model. To avoid over-fitting issues, the conditional moments of " +
                              "the pseudo-residuals are estimated with the same model across features")
            u_des_0_na = u_des_0_noflex
            u_des_0_na = u_des_0_na.reorder_levels(['ID', 'feature', 'Time'])
            u_des_0_na.fillna(0, inplace=True)

        elif u_flex is True:
            u_des_0_na = u_des_0_flex
            u_des_0_na = u_des_0_na.reorder_levels(['ID', 'feature', 'Time'])
            u_des_0_na.fillna(0, inplace=True)

        u_mean = LinearRegression().fit(u_des_0_na, res_na).predict(u_des_0_na)
        params_u = len(u_des_0_na.columns)

    elif u_missp is False:
        u_mean = 0
        params_u = 0
        T_u = 0

    # Use HC inference to estimate V[u|H] (note that w is pre-regularization)
    df = df_EST(w_constr=w_constr[tr], w=w, B=B, J=Jtot, KM=KMI)
    if df >= TT:
        df = TT - 1
        warnings.warn("The current specification uses more degrees of freedom than observations. We suggest to increase the level of " +
                      "sparsity or consider using a smaller donor pool.")

    Sigma, Omega = u_sigma_est(u_mean=u_mean, u_sigma=u_sigma,
                               res=res_na, Z=Z_na, V=V_na,
                               index=index, TT=TT, df=df)

    Sigma_root = sqrtm(Sigma).real

    # Auxiliary logical values to estimate bounds for w
    w_lb_est = True
    w_ub_est = True

    if w_bounds is not None:
        if sum(numpy.isnan(w_bounds[:, 0])) == 0:
            w_lb_est = False
        if sum(numpy.isnan(w_bounds[:, 1])) == 0:
            w_ub_est = False

    vsig = numpy.empty((2 * T1_tot, sims))
    vsig[:] = numpy.nan

    # Define constrained problem to be simulated
    if w_lb_est is True or w_ub_est is True:
        Q = numpy.array(Z_na).T.dot(V_na).dot(numpy.array(Z_na)) / TT
        b_arr = numpy.array(beta).flatten()
        dire = w_constr_inf[tr]['dir']
        p = w_constr_inf[tr]['p']

        if p == "no norm":
            p_int = 0
        elif p == "L1":
            p_int = 1
        elif p == "L2":
            p_int = 2
        elif p == "L1-L2":
            p_int = None

        vsig = scpi_in(sims, b_arr, Sigma_root, Q, P_na, J, KM, iota, w_lb_est,
                       w_ub_est, p, p_int, Q_star, Q2_star, dire, lb, cores, pass_stata, verbose)

    if w_lb_est is True:
        w_lb = numpy.nanquantile(vsig[:, :len(P_na)], q=u_alpha / 2, axis=0)
        w_lb = pandas.DataFrame(w_lb,
                                index=Y_post_fit.index,
                                columns=Y_post_fit.columns)
        fail_lb = 100 * numpy.sum(numpy.isnan(vsig[:, :len(P_na)]), axis=0) / sims

    else:
        w_lb = w_bounds[:, 0]
        fail_lb = numpy.array([0] * len(w_bounds[0]))
        w_lb = pandas.DataFrame(w_lb,
                                index=Y_post_fit.index,
                                columns=Y_post_fit.columns)

    if w_ub_est is True:
        w_ub = numpy.nanquantile(vsig[:, len(P_na):], q=(1 - u_alpha / 2), axis=0)
        w_ub = pandas.DataFrame(w_ub,
                                index=Y_post_fit.index,
                                columns=Y_post_fit.columns)
        fail_ub = 100 * numpy.sum(numpy.isnan(vsig[:, len(P_na):]), axis=0) / sims
    else:
        w_ub = w_bounds[:, 1]
        fail_ub = numpy.array([0] * len(w_bounds[1]))
        w_ub = pandas.DataFrame(w_ub,
                                index=Y_post_fit.index,
                                columns=Y_post_fit.columns)

    fail = numpy.c_[fail_lb, fail_ub]
    fail = numpy.array([fail_lb, fail_ub])
    failed_sims = pandas.DataFrame(fail,
                                   index=['lb', 'ub'])

    if verbose and (fail.sum() > 0.1 * sims * 2 * T1_tot):
        warnings.warn("For some of the simulations used to quantify in-sample uncertainty the solution of " +
                      "the optimization problem was not found! We suggest inspecting the magnitude of this issue " +
                      "by consulting the percentage of simulations that failed contained in " +
                      "YOUR_SCPI_OBJECT_NAME.failed_sims.")

    # PIs for w
    sc_l_0 = Y_post_fit + w_lb        # Left bound
    sc_r_0 = Y_post_fit + w_ub        # Left bound
    len_0 = sc_r_0 - sc_l_0           # Length

    ######################################
    ######################################
    # Estimate out-of-sample uncertainty
    ######################################
    ######################################

    if w_constr_inf[tr_units[0]]['p'] in ["L1-L2", "L2"]:
        epsk = epskappaGet(P, rho_dict, beta, tr_units, effect=sc_effect)
        epsk = pandas.DataFrame(epsk, index=w_lb.index)
        epsk_j = epskappaGet(P, rho_dict, beta, tr_units, effect=sc_effect, joint=True)
    else:
        epsk = pandas.DataFrame([0] * len(w_lb), index=w_lb.index)
        epsk_j = 0

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

    er_dict = mat2dict(e_res_na, cols=False)
    ed0_dict = mat2dict(e_des_0_na)
    ed1_dict = mat2dict(e_des_1)

    T_e = len(e_des_0_na)
    params_e = len(e_des_0_na.columns)
    if (T_e - 10) <= params_e:
        ix0 = e_des_0_na.index
        ix1 = e_des_1.index
        e_des_0_na = pandas.DataFrame(None)
        e_des_1 = pandas.DataFrame(None)

        for tr in tr_units:
            edict0 = pandas.DataFrame(numpy.ones(len(ed0_dict[tr])), columns=[tr])
            edict1 = pandas.DataFrame(numpy.ones(len(ed1_dict[tr])), columns=[tr])
            edict0.insert(0, "ID", tr)
            edict0.set_index('ID', append=False, drop=True, inplace=True)
            edict1.insert(0, "ID", tr)
            edict1.set_index('ID', append=False, drop=True, inplace=True)
            e_des_0_na = pandas.concat([e_des_0_na, edict0], axis=0)
            e_des_1 = pandas.concat([e_des_1, edict1], axis=0)

        if verbose and (e_order > 0 or e_lags > 0):
            warnings.warn("One of e_order > 0 and e_lags > 0 was specified, however the current number of " +
                          "observations (" + str(T_e) + ") used to estimate conditional moments of the " +
                          "out-of-sample error is not larger than the number of parameters used in " +
                          "estimation (" + str(params_e) + ") plus 10. To avoid over-fitting issues " +
                          "e_order and e_lags were set to 0.")
        e_order = 0
        e_lags = 0
        params_e = len(e_des_0_na.columns)
        e_des_0_na.set_index(ix0, append=False, inplace=True)
        e_des_1.set_index(ix1, append=False, inplace=True)
        e_des_0_na.fillna(0, inplace=True)
        e_des_1.fillna(0, inplace=True)

    e_des_0_na = pandas.DataFrame(None)
    e_des_1 = pandas.DataFrame(None)

    if sc_effect == "time":
        scale_x = sc_pred.iota
    else:
        scale_x = 1

    for tr in tr_units:
        ed0_dict[tr] = detectConstant(ed0_dict[tr], tr)
        ed1_dict[tr] = detectConstant(ed1_dict[tr], tr, scale_x)
        ed0_dict[tr].insert(0, "ID", tr)
        ed0_dict[tr].set_index('ID', append=False, drop=True, inplace=True)
        ed1_dict[tr].insert(0, "ID", tr)
        ed1_dict[tr].set_index('ID', append=False, drop=True, inplace=True)
        e_des_0_na = pandas.concat([e_des_0_na, ed0_dict[tr]], axis=0)
        e_des_1 = pandas.concat([e_des_1, ed1_dict[tr]], axis=0)

    e_des_0_na.set_index(e_res_na.index, append=False, inplace=True)
    if sc_effect == "time":
        idx = pandas.MultiIndex.from_product([e_des_1.index.unique('ID').tolist(),
                                             [i for i in range(1, len(P_na) + 1)]],
                                             names=['ID', 'Time'])
        e_des_1.set_index(idx, append=False, inplace=True)
    else:
        e_des_1.set_index(P_na.index, append=False, inplace=True)
    e_des_0_na.fillna(0, inplace=True)
    e_des_1.fillna(0, inplace=True)

    if e_method == 'gaussian' or e_method == 'all':
        e_lb_gau, e_ub_gau, e_1, e_2 = scpi_out(y=e_res_na, x=e_des_0_na, preds=e_des_1,
                                                e_method="gaussian", alpha=e_alpha / 2,
                                                e_lb_est=e_lb_est, e_ub_est=e_ub_est,
                                                effect=sc_effect, out_feat=out_feat[tr])

        # Overwrite with user's input
        if e_lb_est is False:
            e_lb_gau = e_bounds[0]
        if e_ub_est is False:
            e_ub_gau = e_bounds[1]

        sc_l_1 = sc_l_0.iloc[:, 0] + e_lb_gau.iloc[:, 0] - epsk.iloc[:, 0]
        sc_r_1 = sc_r_0.iloc[:, 0] + e_ub_gau.iloc[:, 0] + epsk.iloc[:, 0]
        len_1 = sc_r_1 - sc_l_1

        e_mean = e_1
        e_var = e_2

    if e_method == 'ls' or e_method == 'all':
        e_lb_ls, e_ub_ls, e_1, e_2 = scpi_out(y=e_res_na, x=e_des_0_na, preds=e_des_1,
                                              e_method="ls", alpha=e_alpha / 2,
                                              e_lb_est=e_lb_est, e_ub_est=e_ub_est,
                                              effect=sc_effect, out_feat=out_feat[tr])

        # Overwrite with user's input
        if e_lb_est is False:
            e_lb_ls = e_bounds[0]
        if e_ub_est is False:
            e_ub_ls = e_bounds[1]

        sc_l_2 = sc_l_0.iloc[:, 0] + e_lb_ls.iloc[:, 0] - epsk.iloc[:, 0]
        sc_r_2 = sc_r_0.iloc[:, 0] + e_ub_ls.iloc[:, 0] + epsk.iloc[:, 0]
        len_2 = sc_r_2 - sc_l_2

    if e_method == 'qreg' or e_method == 'all':
        if e_order == 0:
            e_lb_qreg = numpy.quantile(e_res_na, q=e_alpha / 2)
            e_lb_qreg = pandas.DataFrame([e_lb_qreg] * len(w_lb), index=w_lb.index)
            e_ub_qreg = numpy.quantile(e_res_na, q=1 - e_alpha / 2)
            e_ub_qreg = pandas.DataFrame([e_ub_qreg] * len(w_lb), index=w_ub.index)

        else:
            e_lb_qreg, e_ub_qreg, e_1, e_2 = scpi_out(y=e_res_na, x=e_des_0_na, preds=e_des_1,
                                                      e_method="qreg", alpha=e_alpha / 2,
                                                      e_lb_est=e_lb_est, e_ub_est=e_ub_est,
                                                      effect=sc_effect, out_feat=out_feat[tr])
        # Overwrite with user's input
        if e_lb_est is False:
            e_lb_qreg = e_bounds[0]
        if e_ub_est is False:
            e_ub_qreg = e_bounds[1]

        sc_l_3 = sc_l_0.iloc[:, 0] + e_lb_qreg.iloc[:, 0] - epsk.iloc[:, 0]
        sc_r_3 = sc_r_0.iloc[:, 0] + e_ub_qreg.iloc[:, 0] + epsk.iloc[:, 0]
        len_3 = sc_r_3 - sc_l_3

    ####################################################
    # Simultaneous Prediction Intervals (for each unit)

    if sc_effect == "unit-time":  # joint within unit
        T1val = [v for v in T1.values()]
        ML, MU = simultaneousPredGet(vsig, T1val, len(P_na), iota, u_alpha,
                                     e_alpha, e_res_na, e_des_0_na, e_des_1,
                                     w_lb_est, w_ub_est, w_bounds,
                                     w_constr_aux[tr_units[0]]['name'], sc_effect, out_feat)

    elif sc_effect == "unit":  # joint across units
        ML, MU = simultaneousPredGet(vsig, [len(P_na)], len(P_na), 1, u_alpha,
                                     e_alpha, e_res_na, e_des_0_na, e_des_1,
                                     w_lb_est, w_ub_est, w_bounds,
                                     w_constr_aux[tr_units[0]]['name'], sc_effect, out_feat)

    elif sc_effect == "time":  # joint within aggregate unit
        ML, MU = simultaneousPredGet(vsig, [len(P_na)], len(P_na), 1, u_alpha,
                                     e_alpha, e_res_na, e_des_0_na, e_des_1,
                                     w_lb_est, w_ub_est, w_bounds,
                                     w_constr_aux[tr_units[0]]['name'], sc_effect, out_feat)

    ML.set_index(P_na.index, inplace=True, append=False)
    MU.set_index(P_na.index, inplace=True, append=False)

    ###############################################
    # Store all bounds
    if sc_effect == "time":
        idx = pandas.MultiIndex.from_product([["aggregate"], w_lb.index.get_level_values(0).tolist()],
                                             names=['ID', 'Time'])
    else:
        idx = w_lb.index
        idx.rename(["ID", "Time"], inplace=True)

    df_insample = pandas.DataFrame(numpy.c_[w_lb, w_ub], columns=["Lower", "Upper"], index=idx)
    df_insample = df_insample.apply(pandas.to_numeric, errors='coerce', axis=1)

    if e_method == "all" or e_method == "gaussian":
        sub_lb = w_lb.iloc[:, 0] + e_lb_gau.iloc[:, 0] - epsk.iloc[:, 0]
        sub_ub = w_ub.iloc[:, 0] + e_ub_gau.iloc[:, 0] + epsk.iloc[:, 0]
        df_subgaussian = pandas.DataFrame(numpy.c_[sub_lb, sub_ub], columns=["Lower", "Upper"], index=idx)
        df_subgaussian = df_subgaussian.apply(pandas.to_numeric, errors='coerce', axis=1)
    else:
        df_subgaussian = None

    if e_method == "all" or e_method == "ls":
        ls_lb = w_lb.iloc[:, 0] + e_lb_ls.iloc[:, 0] - epsk.iloc[:, 0]
        ls_ub = w_ub.iloc[:, 0] + e_lb_ls.iloc[:, 0] + epsk.iloc[:, 0]
        df_ls = pandas.DataFrame(numpy.c_[ls_lb, ls_ub], columns=["Lower", "Upper"], index=idx)
        df_ls = df_ls.apply(pandas.to_numeric, errors='coerce', axis=1)
    else:
        df_ls = None

    if e_method == "all" or e_method == "qreg":
        qreg_lb = w_lb.iloc[:, 0] + e_lb_qreg.iloc[:, 0] - epsk.iloc[:, 0]
        qreg_ub = w_ub.iloc[:, 0] + e_ub_qreg.iloc[:, 0] + epsk.iloc[:, 0]
        df_qreg = pandas.DataFrame(numpy.c_[qreg_lb, qreg_ub], columns=["Lower", "Upper"], index=idx)
        df_qreg = df_qreg.apply(pandas.to_numeric, errors='coerce', axis=1)
    else:
        df_qreg = None

    Mlb = ML.iloc[:, 0] - epsk_j
    Mub = MU.iloc[:, 0] + epsk_j
    df_joint = pandas.DataFrame(numpy.c_[Mlb, Mub], columns=["Lower", "Upper"], index=idx)
    df_joint = df_joint.apply(pandas.to_numeric, errors='coerce', axis=1)

    bounds = {"insample": df_insample,
              "subgaussian": df_subgaussian,
              "ls": df_ls,
              "qreg": df_qreg,
              "joint": df_joint}

    ###############################################
    # Return objects

    CI_0 = pandas.concat([sc_l_0, sc_r_0, len_0], axis=1)
    CI_0.columns = ['Lower', 'Upper', 'Length']
    if sc_effect == "time":
        CI_0.set_index(idx, inplace=True)
    CI_0.index.rename(['ID', 'Time'], inplace=True)
    CI_0 = CI_0.apply(pandas.to_numeric, errors='coerce', axis=1)

    if sc_l_1 is not None:
        CI_1 = pandas.concat([sc_l_1, sc_r_1, len_1], axis=1)
        CI_1.columns = ['Lower', 'Upper', 'Length']
        if sc_effect == "time":
            CI_1.set_index(idx, inplace=True)
        CI_1.index.rename(['ID', 'Time'], inplace=True)
        CI_1 = CI_1.apply(pandas.to_numeric, errors='coerce', axis=1)
    else:
        CI_1 = None

    if sc_l_2 is not None:
        CI_2 = pandas.concat([sc_l_2, sc_r_2, len_2], axis=1)
        CI_2.columns = ['Lower', 'Upper', 'Length']
        if sc_effect == "time":
            CI_2.set_index(idx, inplace=True)
        CI_2.index.rename(['ID', 'Time'], inplace=True)
        CI_2 = CI_2.apply(pandas.to_numeric, errors='coerce', axis=1)
    else:
        CI_2 = None

    if sc_l_3 is not None:
        CI_3 = pandas.concat([sc_l_3, sc_r_3, len_3], axis=1)
        CI_3.columns = ['Lower', 'Upper', 'Length']
        if sc_effect == "time":
            CI_3.set_index(idx, inplace=True)
        CI_3.index.rename(['ID', 'Time'], inplace=True)
        CI_3 = CI_3.apply(pandas.to_numeric, errors='coerce', axis=1)
    else:
        CI_3 = None

    u_user = u_design is not None  # True if user provided the design matrix for in-sample inference
    e_user = e_design is not None  # True if user provided the design matrix for out-of-sample inference

    ##################################################
    # Plot

    if plot is True:
        if class_type == "scpi_data":

            toplot = scpi_output(b=sc_pred.b,
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
                                 KMI=sc_pred.KMI,
                                 M=sc_pred.M,
                                 iota=1,
                                 cointegrated_data=sc_pred.cointegrated_data,
                                 anticipation=sc_pred.anticipation,
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
                                 bounds=bounds,
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
                                 u_T=T_u,
                                 u_params=params_u,
                                 u_D=u_des_0_na,
                                 e_method=e_method,
                                 e_lags=e_lags,
                                 e_order=e_order,
                                 e_user=e_user,
                                 e_T=T_e,
                                 e_params=params_e,
                                 e_D=e_des_0_na,
                                 rho=rho_dict,
                                 Q_star=Q_star,
                                 u_alpha=u_alpha,
                                 e_alpha=e_alpha,
                                 epskappa=epsk,
                                 sims=sims,
                                 failed_sims=failed_sims,
                                 plotres=None,
                                 donors_dict=sc_pred.donors_dict,
                                 treated_units=sc_pred.treated_units,
                                 units_est=sc_pred.units_est,
                                 timeConvert=sc_pred.timeConvert)
            plotres = scplot(result=toplot)
        else:
            to_plot = scpi_multi_output(b=sc_pred.b,
                                        w=sc_pred.w,
                                        r=sc_pred.r,
                                        Y_pre_fit=sc_pred.Y_pre_fit,
                                        Y_post_fit=sc_pred.Y_post_fit,
                                        Y_pre=sc_pred.Y_pre,
                                        Y_post=sc_pred.Y_post,
                                        Y_actual=sc_pred.Y_actual,
                                        A_hat=sc_pred.A_hat,
                                        res=sc_pred.res,
                                        V=sc_pred.V,
                                        w_constr=sc_pred.w_constr,
                                        w_constr_inf=w_constr_inf,
                                        A=sc_pred.A,
                                        B=sc_pred.B,
                                        C=sc_pred.C,
                                        P=P_na,
                                        Y_df=sc_pred.Y_df,
                                        Y_donors=sc_pred.Y_donors,
                                        J=sc_pred.J,
                                        K=sc_pred.K,
                                        KM=sc_pred.KM,
                                        M=sc_pred.M,
                                        iota=sc_pred.iota,
                                        KMI=sc_pred.KMI,
                                        cointegrated_data=sc_pred.cointegrated_data,
                                        anticipation=sc_pred.anticipation,
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
                                        bounds=bounds,
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
                                        u_T=T_u,
                                        u_params=params_u,
                                        u_D=u_des_0_na,
                                        e_method=e_method,
                                        e_lags=e_lags,
                                        e_order=e_order,
                                        e_user=e_user,
                                        e_T=T_e,
                                        e_params=params_e,
                                        e_D=e_des_0_na,
                                        rho=rho_dict,
                                        Q_star=Q_star,
                                        u_alpha=u_alpha,
                                        e_alpha=e_alpha,
                                        epskappa=epsk,
                                        sims=sims,
                                        failed_sims=failed_sims,
                                        plotres=plotres,
                                        effect=sc_pred.effect,
                                        donors_dict=sc_pred.donors_dict,
                                        treated_units=sc_pred.treated_units,
                                        units_est=sc_pred.units_est,
                                        timeConvert=sc_pred.timeConvert)

            plotres = scplotMulti(result=to_plot)
    else:
        plotres = None

    if class_type == "scpi_data":

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
                           KMI=sc_pred.KMI,
                           M=sc_pred.M,
                           iota=1,
                           cointegrated_data=sc_pred.cointegrated_data,
                           anticipation=sc_pred.anticipation,
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
                           bounds=bounds,
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
                           u_T=T_u,
                           u_params=params_u,
                           u_D=u_des_0_na,
                           e_method=e_method,
                           e_lags=e_lags,
                           e_order=e_order,
                           e_user=e_user,
                           e_T=T_e,
                           e_params=params_e,
                           e_D=e_des_0_na,
                           rho=rho_dict,
                           Q_star=Q_star,
                           u_alpha=u_alpha,
                           e_alpha=e_alpha,
                           epskappa=epsk,
                           sims=sims,
                           failed_sims=failed_sims,
                           plotres=plotres,
                           donors_dict=sc_pred.donors_dict,
                           treated_units=sc_pred.treated_units,
                           units_est=sc_pred.units_est,
                           timeConvert=sc_pred.timeConvert)
    else:
        return scpi_multi_output(b=sc_pred.b,
                                 w=sc_pred.w,
                                 r=sc_pred.r,
                                 Y_pre_fit=sc_pred.Y_pre_fit,
                                 Y_post_fit=sc_pred.Y_post_fit,
                                 Y_pre=sc_pred.Y_pre,
                                 Y_post=sc_pred.Y_post,
                                 Y_actual=sc_pred.Y_actual,
                                 A_hat=sc_pred.A_hat,
                                 res=sc_pred.res,
                                 V=sc_pred.V,
                                 w_constr=sc_pred.w_constr,
                                 w_constr_inf=w_constr_inf,
                                 A=sc_pred.A,
                                 B=sc_pred.B,
                                 C=sc_pred.C,
                                 P=P_na,
                                 Y_df=sc_pred.Y_df,
                                 Y_donors=sc_pred.Y_donors,
                                 J=sc_pred.J,
                                 K=sc_pred.K,
                                 KM=sc_pred.KM,
                                 M=sc_pred.M,
                                 iota=sc_pred.iota,
                                 KMI=sc_pred.KMI,
                                 cointegrated_data=sc_pred.cointegrated_data,
                                 anticipation=sc_pred.anticipation,
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
                                 bounds=bounds,
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
                                 u_T=T_u,
                                 u_params=params_u,
                                 u_D=u_des_0_na,
                                 e_method=e_method,
                                 e_lags=e_lags,
                                 e_order=e_order,
                                 e_user=e_user,
                                 e_T=T_e,
                                 e_params=params_e,
                                 e_D=e_des_0_na,
                                 rho=rho_dict,
                                 Q_star=Q_star,
                                 u_alpha=u_alpha,
                                 e_alpha=e_alpha,
                                 epskappa=epsk,
                                 sims=sims,
                                 failed_sims=failed_sims,
                                 plotres=plotres,
                                 effect=sc_pred.effect,
                                 donors_dict=sc_pred.donors_dict,
                                 treated_units=sc_pred.treated_units,
                                 units_est=sc_pred.units_est,
                                 timeConvert=sc_pred.timeConvert)


# Define class

class scpi_output:
    def __init__(self, b, w, r, Y_pre_fit, Y_post_fit, A_hat, res, V, w_constr, w_constr_inf,
                 A, B, C, P, Y_pre, Y_post, Y_donors, J, K, KM, KMI, M, iota, anticipation,
                 cointegrated_data, period_pre, period_post, T0_features,
                 T1_outcome, features, outcome_var, glob_cons, out_in_features,
                 CI_in_sample, CI_all_gaussian, CI_all_ls, CI_all_qreg, bounds, Sigma,
                 u_mean, u_var, e_mean, e_var, u_missp, u_lags, u_order,
                 u_sigma, u_user, u_T, u_params, u_D, e_method, e_lags, e_order, e_user, e_T,
                 e_params, e_D, rho, Q_star, u_alpha, e_alpha, epskappa, sims, failed_sims, plotres,
                 donors_dict, treated_units, units_est, timeConvert):

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
        self.KMI = KMI
        self.M = M
        self.iota = iota
        self.cointegrated_data = cointegrated_data
        self.anticipation = anticipation
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
        self.bounds = bounds
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
        self.u_T = u_T
        self.u_params = u_params
        self.u_D = u_D
        self.e_method = e_method
        self.e_lags = e_lags
        self.e_order = e_order
        self.e_user = e_user
        self.e_T = e_T
        self.e_params = e_params
        self.e_D = e_D
        self.rho = rho
        self.Q_star = Q_star
        self.u_alpha = u_alpha
        self.e_alpha = e_alpha
        self.epskappa = epskappa
        self.sims = sims
        self.failed_sims = failed_sims
        self.plotres = plotres
        self.donors_dict = donors_dict
        self.treated_units = treated_units
        self.units_est = units_est
        self.timeConvert = timeConvert

    def __repr__(self):

        # Prepare objects to print (estimation)
        fw = 30
        if self.M == 1:
            fw_r = 14
        if self.M > 1:
            fw_r = 40

        w_constr = deepcopy(self.w_constr[self.treated_units[0]])
        constr = w_constr['name']
        if w_constr['Q'] is not None:
            Qsize = round(w_constr['Q'], 3)
        else:
            Qsize = "-"
        tr_unit = self.treated_units[0]
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
            print("     Parameters used to estimate moments      " + str(self.u_params))
        else:
            print("In-sample Inference:")
            print("     User provided")

        if (self.e_user is False):
            print("Out-of-sample Inference:                          ")
            print("     Method                                   " + self.e_method)
            print("     Order of polynomial (B)                  " + str(self.e_order))
            print("     Lags (B)                                 " + str(self.e_lags))
            print("     Parameters used to estimate moments      " + str(self.e_params))
        else:
            print("Out-of-sample Inference:")
            print("     User provided")

        if e_method == 'gaussian':
            print('   Inference with subgaussian bounds')
            dfprint = pandas.concat([Y_tr_post, Y_sc_post, CI], axis=1)
            dfprint.index.rename(['Treated Unit', 'Time'], inplace=True)
            print(dfprint)
        elif e_method == 'ls':
            print('   Inference with location-scale model')
            dfprint = pandas.concat([Y_tr_post, Y_sc_post, CI], axis=1)
            dfprint.index.rename(['Treated Unit', 'Time'], inplace=True)
            print(dfprint)
        elif e_method == 'qreg':
            print('   Inference with quantile regression')
            dfprint = pandas.concat([Y_tr_post, Y_sc_post, CI], axis=1)
            dfprint.index.rename(['Treated Unit', 'Time'], inplace=True)
            print(dfprint)
        elif e_method == 'all':
            print('                              Subgaussian')
            dfprint = pandas.concat([Y_tr_post, Y_sc_post, CI1], axis=1)
            dfprint.index.rename(['Treated Unit', 'Time'], inplace=True)
            print(dfprint)
            print('            Location Scale          Quantile Reg')
            print(pandas.concat([CI2, CI3], axis=1))

        else:
            print(pandas.concat([Y_tr_post, Y_sc_post, CI], axis=1).index.rename(['Treated Unit', 'Time']))

        return ''


class scpi_multi_output:
    def __init__(self, b, w, r, Y_pre_fit, Y_post_fit, Y_pre, Y_post, Y_actual, A_hat, res, V, w_constr, w_constr_inf,
                 A, B, C, P, Y_df, Y_donors, J, K, KM, M, iota, KMI,
                 cointegrated_data, anticipation, period_pre, period_post, T0_features,
                 T1_outcome, features, outcome_var, glob_cons, out_in_features,
                 CI_in_sample, CI_all_gaussian, CI_all_ls, CI_all_qreg, bounds, Sigma,
                 u_mean, u_var, e_mean, e_var, u_missp, u_lags, u_order,
                 u_sigma, u_user, u_T, u_params, u_D, e_method, e_lags, e_order, e_user, e_T,
                 e_params, e_D, rho, Q_star, u_alpha, e_alpha, epskappa, sims, failed_sims, plotres,
                 effect, donors_dict, treated_units, units_est, timeConvert):

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
        self.w_constr_inf = w_constr_inf
        self.A = A
        self.B = B
        self.C = C
        self.P = P
        self.Y_df = Y_df
        self.Y_donors = Y_donors
        self.J = J
        self.K = K
        self.KM = KM
        self.M = M
        self.iota = iota
        self.KMI = KMI
        self.cointegrated_data = cointegrated_data
        self.anticipation = anticipation
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
        self.bounds = bounds
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
        self.u_T = u_T
        self.u_params = u_params
        self.u_D = u_D
        self.e_method = e_method
        self.e_lags = e_lags
        self.e_order = e_order
        self.e_user = e_user
        self.e_T = e_T
        self.e_params = e_params
        self.e_D = e_D
        self.rho = rho
        self.Q_star = Q_star
        self.u_alpha = u_alpha
        self.e_alpha = e_alpha
        self.epskappa = epskappa
        self.sims = sims
        self.failed_sims = failed_sims
        self.plotres = plotres
        self.effect = effect
        self.donors_dict = donors_dict
        self.treated_units = treated_units
        self.units_est = units_est
        self.timeConvert = timeConvert

    def __repr__(self):

        return ''
