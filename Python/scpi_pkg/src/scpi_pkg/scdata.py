# -*- coding: utf-8 -*-

# Temporary code to suppress pandas FutureWarning
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas
pandas.options.mode.chained_assignment = None
import pandas
import numpy
import scipy.linalg
from copy import deepcopy
from collections import Counter
from .funs import complete_cases


def scdata(df,
           id_var,
           time_var,
           outcome_var,
           period_pre,
           period_post,
           unit_tr,
           unit_co,
           features=None,
           cov_adj=None,
           cointegrated_data=False,
           anticipation=0,
           constant=False,
           verbose=True,
           report_missing=False):

    """
    Parameters
    ----------
    df : pandas.DataFrame
        a dataframe object containing the data to be processed

    id_var : str
        a character with the name of the variable containing units IDs

    time_var : str
        a character with the name of the time variable. The time variable has to be numpy.int64,
        numpy.datetime64, or pandas.Timestamp. Input a numeric time variable is suggested when working with
        yearly data, whereas for all other frequencies numpy.datetime64 type is preferred.

    outcome_var : str
        a character with the name of the outcome variable. The outcome variable has to be numeric.

    period_pre : array
        a numeric vector that identifies the pre-treatment period in time_var.

    period_post : array
        a numeric vector that identifies the post-treatment period in time_var.

    unit_tr : int
        a scalar that identifies the treated unit in id_var.

    unit_co : array
         a numeric vector that identifies the donor pool in id_var.

    features : list, default None
        a character list containing the name of the feature variables used for estimation.
        If this option is not specified the default is features = outcome_var.

    cov_adj : list, default None
        a list specifying the names of the covariates to be used for adjustment for each feature. If the user wants
        to specify the same set of covariates for all features, a single list should be provided. If instead a
        different set of covariates per feature has to be specified, then a list of lists should be provided. Note that
        in this latter case the number of sub-lists must be equal to the number of features. Moreover, the order of the
        sub-lists matters, in the sense that the first sub-list is interpreted as the set of covariates for the first
        feature, and so on. Finally, the user can specify 'constant' and 'trend' as covariates even if they are not
        present in the loaded dataframe.

    constant : bool, default False
         a logical which controls the inclusion of a constant term across features. The default value is False.

    cointegrated_data : bool, default False
        a logical that indicates if there is a belief that the data is cointegrated or not. The default value is False.
        See the Details section for more.

    anticipation : int, default 0
        a scalar that indicates the number of periods of potential anticipation effects. Default is 0.

    verbose : bool, default True
        a logical to print additional information in the console.

    report_missing : bool, default False
        a logical which prints the location of missing values if present. The default value is False.

    Returns
    -------
    The function returns an object of class 'scdata_output' containing the following objects

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

    features : list
        a list with the name of the features

    out_in_features : bool
        for internal use only

    effect : str
        for internal use only

    glob_cons : bool
        for internal use only

    out_in_features : bool
        for internal use only

    timeConvert : bool
        for internal use only

    References
    ----------
    Abadie, A. (2021), “Using Synthetic Controls: Feasibility, Data Requirements, and Methodological
    Aspects,” Journal of Economic Literature, 59, 391-425.

    Cattaneo, M. D., Feng, Y., and Titiunik, R. (2021), “Prediction Intervals for Synthetic Control
    Methods,” Journal of the American Statistical Association, 116, 1865-1880.

    Cattaneo, M. D., Palomba, F., Feng, Y., and Titiunik, R. (2022), “scpi: Uncertainty Quantification for
    Synthetic Control Estimators”.

    Cattaneo, M. D., Palomba, F., Feng, Y., and Titiunik, R. (2023), “Uncertainty Quantification in Synthetic
    Controls with Staggered Treatment Adoption”.

    See Also
    --------
    scdataMulti, scest, scpi, scplot, scplotMulti

    """

    # Check main input is a dataframe
    if not isinstance(df, pandas.DataFrame):
        raise Exception('Data input should be a dataframe object!')

    data = deepcopy(df)

    # Store variable names and indexes
    var_names = data.columns
    indexes = data.index.names

    # Check inputs are strings
    if not isinstance(id_var, str):
        raise Exception("You should specify the name of id_var as a string! (eg. id_var = 'ID')")

    if not isinstance(outcome_var, str):
        raise Exception("You should specify the name of outcome_var as a string! (eg. outcome_var = 'outcome')")

    if not isinstance(time_var, str):
        raise Exception("You should specify the name of time_var as a string! (eg. time_var = 'time')")

    if features is None:
        features = [outcome_var]
    elif not isinstance(features, list):
        raise Exception("The object 'features' should be a list!")
    else:
        features.sort()

    # Check covariates for adjustment are in dataframe
    if cov_adj is not None:
        if not isinstance(cov_adj, list):
            raise Exception("The argument cov_adj should be a list!")

        if isinstance(cov_adj[0], list):  # feature-specific cov adj
            if len(cov_adj) != len(features):
                raise Exception("When specifying covariate adjustment separately for each feature make sure " +
                                "to do it for all features! You specified covariate adjustment for " +
                                str(len(cov_adj)) + " features when you currently have " +
                                str(len(features)) + " features!")
            unique_covs = [cov for sublist in cov_adj for cov in sublist]

        else:
            unique_covs = cov_adj
        # Check all covariates other than constant and trend are in dataframe

        unique_cov = list(set(unique_covs))
        unique_covs = [x for x in unique_cov if x not in ['constant', 'trend']]

        if len(unique_covs) > 0:
            if not all(cov in var_names for cov in unique_covs):
                cov_adj_not_found = [cov for cov in unique_covs if cov in var_names]
                cov_print = ' '.join(str(cov) for cov in cov_adj_not_found)
                raise Exception("Some of the covariates in cov_adj are not in the input dataframe!"
                                + cov_print)

    # Check variables are in dataframe (id and time can be either variables or indices)
    if (id_var not in var_names) and (id_var not in indexes):
        raise Exception("ID variable (id_var) not found in the input dataframe neither as a variable nor as an index!")

    if (time_var not in var_names) and (time_var not in indexes):
        raise Exception("Time variable (time_var) not found in the input df neither as a variable nor as an index!")

    # Make time and id columns if these variables are indexes of the dataframe
    if id_var in indexes:
        data['__ID'] = data.index.get_level_values(id_var)
    else:
        data.rename(columns={id_var: '__ID'}, inplace=True)

    if time_var in indexes:
        data['__time'] = data.index.get_level_values(time_var)
    else:
        data.rename(columns={time_var: '__time'}, inplace=True)

    timeConvert = False
    dd = data.iloc[0, data.columns.get_loc('__time')]
    if not isinstance(dd, (numpy.int64, numpy.int32, numpy.int16, pandas.Timestamp, numpy.datetime64)):
        raise Exception("The object time_var should be of type int, numpy.datetime64, or pandas.Timestamp!")

    elif isinstance(dd, (pandas.Timestamp, numpy.datetime64)):
        time_unique_ts = sorted(set(data['__time'].tolist()))
        int2ts = {i: time_unique_ts[i] for i in range(len(time_unique_ts))}
        ts2int = {time_unique_ts[i]: i for i in range(len(time_unique_ts))}
        data['__time'] = data['__time'].map(ts2int)
        period_pre = pandas.Series(period_pre).map(ts2int).to_numpy()
        period_post = pandas.Series(period_post).map(ts2int).to_numpy()
        timeConvert = True

    if outcome_var not in var_names:
        raise Exception("Outcome variable (outcome_var) not found in the input dataframe!")

    # CVXR does not like _ symbols
    if pandas.api.types.is_string_dtype(data['__ID']) is False:  # convert to string if not string
        data['__ID'] = data['__ID'].astype(str)
        unit_co = [str(co) for co in unit_co]
        unit_tr = str(unit_tr)
    data['__ID'] = data['__ID'].str.replace('_', ' ')
    unit_co = [s.replace('_', ' ') for s in unit_co]
    unit_tr = unit_tr.replace('_', ' ')

    data['treated_unit'] = unit_tr

    # Create index
    data.set_index(['treated_unit', '__ID', '__time'], drop=False, inplace=True)

    if not isinstance(period_pre, numpy.ndarray):
        raise Exception("The object period_pre should be of type numpy.ndarray (eg. use numpy.arange or numpy.array)!")

    if not isinstance(period_post, numpy.ndarray):
        raise Exception("The object period_post should be of type numpy.ndarray!")

    if not isinstance(period_pre[0], (numpy.int64, numpy.int32, numpy.int16, numpy.datetime64, pandas.Timestamp)):
        raise Exception("Elements of period_pre should either be of type int, pandas.Timestamp, or numpy.datetime64!")

    if not isinstance(period_post[0], (numpy.int64, numpy.int32, numpy.int16, numpy.datetime64, pandas.Timestamp)):
        raise Exception("Elements of period_post should either be of type int pandas.Timestamp, or numpy.datetime64!")

    if not isinstance(unit_co, list):
        raise Exception("The object unit_co should be of type list!")

    if not isinstance(unit_tr, list):
        unit_tr_list = [unit_tr]
    else:
        unit_tr_list = unit_tr

    if not data[outcome_var].dtype in ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']:
        raise Exception("Outcome variable (outcome_var) must be numeric!")

    if not all(feature in var_names for feature in features):
        fe_not_found = [feature for feature in features if feature not in var_names]
        fe_print = ' '.join(str(feature) for feature in fe_not_found)
        raise Exception("The following features are not in the input dataframe:"
                        + fe_print)

    if outcome_var in features:
        out_in_features = True
    else:
        out_in_features = False

    if cov_adj is not None:
        if len(features) == 1 and constant is True and 'constant' in unique_cov:
            raise Exception("When specifying just one feature you either specify constant == True" +
                            " or include 'constant' in cov_adj!")

    # Sort period variables
    period_pre = numpy.sort(period_pre)
    period_post = numpy.sort(period_post)

    # Create ID and time variables
    idd = list(set(data['__ID'].to_list()))   # unique IDs

    time = data['__time'].unique()             # unique periods

    # Check that specified units are in dataframe
    if not all(tr in idd for tr in unit_tr_list):
        raise Exception("There is no treated unit with the specified ID (unit_tr) in the specified" +
                        "ID variable (id_var)!")

    if not all(donor in idd for donor in unit_co):
        donors_not_found = [dnr for dnr in unit_co if dnr not in idd]
        donors_print = ' '.join(str(donor) for donor in donors_not_found)

        raise Exception("The following control unit(s) are not in the input dataframe:" +
                        donors_print)

    if any(co in unit_tr_list for co in unit_co):
        raise Exception("The treated unit is also contained in the donor pool!")

    if len(unit_co) < 2:
        raise Exception("Please provide at least two control units!")

    # Check specified time periods are in dataframe
    if not all(period in time for period in period_pre):
        period_not_found = [period for period in period_pre if period not in time]
        period_print = ' '.join(str(period) for period in period_not_found)
        raise Exception("The following pre-treatment period(s) are not in the input dataframe:" +
                        period_print)

    if not all(period in time for period in period_post):
        period_not_found = [period for period in period_post if period not in time]
        period_print = ' '.join(str(period) for period in period_not_found)
        raise Exception("The following post-treatment period(s) are not in the input dataframe:" +
                        period_print)

    if any(period in period_post for period in period_pre):
        raise Exception("There is an overlap between the pre-treatment period and post-treatment period!")

    # Consider eventual anticipation effect
    if anticipation > 0:
        t0 = len(period_pre)
        d = anticipation
        period_post = numpy.insert(period_post, 0, period_pre[(t0 - d):t0])
        period_pre = period_pre[:(t0 - d)]

    ############################################################################
    # Data preparation

    # Make the panel balanced
    data_bal = deepcopy(data.unstack().stack(dropna=False))
    data_bal['__ID'] = data_bal.index.get_level_values('__ID')
    data_bal['__time'] = data_bal.index.get_level_values('__time')

    rows_tr_pre = data_bal.loc[(unit_tr, unit_tr_list, period_pre), ]
    rows_tr_post = data_bal.loc[(unit_tr, unit_tr_list, period_post), ]
    rows_co_pre = data_bal.loc[(unit_tr, unit_co, period_pre), ]
    rows_co_post = data_bal.loc[(unit_tr, unit_co, period_post), ]

    ############################################################################
    # Estimation Data

    # Actual Pre-treatment Series
    Y_pre = rows_tr_pre[[outcome_var]].reset_index(level='__ID', drop=True)
    Y_pre.columns = unit_tr_list

    if len(unit_tr_list) == 1:
        label = unit_tr_list[0]
    else:
        label = "Treated"

    # Create A
    A_df = pandas.melt(rows_tr_pre,
                       id_vars=['treated_unit', '__ID', '__time'],
                       value_vars=features,
                       var_name='feature',
                       value_name=label)

    A = A_df.loc[:, [unit_tr]]

    # Create B

    # Stack features one on top of the other
    B_df = pandas.melt(rows_co_pre,
                       id_vars=['treated_unit', '__ID', '__time'],
                       value_vars=features,
                       var_name='feature',
                       value_name='value')

    # make the df wide so that countries are one next to the other
    B = B_df.pivot(index=['treated_unit', 'feature', '__time'],   # index to keep fixed
                   columns='__ID',                                # dimension to wide
                   values='value')                                # values to be widen

    if numpy.nan in B.index.get_level_values(0).unique().tolist():
        unitnonan = [un for un in B.index.get_level_values(0).unique().tolist() if un is not numpy.nan]
        idx = pandas.IndexSlice
        B = B.loc[idx[unitnonan[0], :, :]]
        B.insert(0, "treated_unit", unitnonan[0])
        B.set_index('treated_unit', append=True, drop=True, inplace=True)
        B = B.reorder_levels(['treated_unit', 'feature', '__time'])

    # Create list that dictates order of donors
    donor_order = B.columns.values.tolist()

    # Create matrix with outcome for the donors
    sel = rows_co_pre[[outcome_var, 'treated_unit', '__time', '__ID']]
    Y_donors = sel.pivot(index=['treated_unit', '__time'],
                         columns='__ID', values=outcome_var)
    Y_donors = Y_donors[donor_order]

    # Create C
    C = pandas.DataFrame(None)
    C_names = []

    if cov_adj is not None:

        # Fixed covariate adjustment across features
        if not isinstance(cov_adj[0], list):
            covs_adj = deepcopy(cov_adj)

            # Check that constant/time trend are required by the user
            if 'constant' in cov_adj:
                covs_adj = [cov for cov in covs_adj if not cov == 'constant']
                C['constant'] = numpy.ones(len(rows_tr_pre))

            if 'trend' in cov_adj:
                covs_adj = [cov for cov in covs_adj if not cov == 'trend']
                C['trend'] = period_pre - period_pre[0] + 1

            # Add other covariates from dataframe
            rows_C = data_bal.loc[(unit_tr, unit_co[0], period_pre), covs_adj].reset_index()
            C[covs_adj] = rows_C[covs_adj]

            for num in range(1, len(features) + 1):
                for cov in C.columns:
                    C_names.append(str(num) + '_' + cov)

            # Create block diagonal matrix
            C_arr = numpy.kron(numpy.identity(len(features)), numpy.array(C))
            C = pandas.DataFrame(data=C_arr, columns=C_names)

        else:
            C_names = []
            for m in range(len(features)):
                C_m = pandas.DataFrame(None)
                covs_adj = cov_adj[m]

                # Check that constant/time trend are required by the user
                if 'constant' in covs_adj:
                    covs_adj = [cov for cov in covs_adj if not cov == 'constant']
                    C_m['constant'] = numpy.ones(len(rows_tr_pre))

                if 'trend' in covs_adj:
                    covs_adj = [cov for cov in covs_adj if not cov == 'trend']
                    C_m['trend'] = period_pre - period_pre[0] + 1

                # Add other covariates from dataframe
                rows_C = data_bal.loc[(unit_tr, unit_co[0], period_pre), covs_adj].reset_index()
                C_m[covs_adj] = rows_C[covs_adj]

                if m == 0:
                    C = C_m
                else:
                    C = scipy.linalg.block_diag(numpy.array(C), numpy.array(C_m))

                feat_num = m + 1

                for cov in cov_adj[m]:
                    C_names.append(str(feat_num) + '_' + cov)

            C = pandas.DataFrame(data=C, columns=C_names)

            # take into account empty covariate adjustemnt
            for m in range(len(features)):
                if len(cov_adj[m]) == 0:
                    i_start = len(rows_tr_pre) * m
                    zeros = pandas.DataFrame(numpy.zeros((len(rows_tr_pre), len(C.columns))),
                                             columns=C_names)
                    if i_start == 0:
                        C_up = pandas.DataFrame(None)
                    else:
                        C_up = C.iloc[:i_start]

                    C_bo = C.iloc[i_start:]
                    C = pandas.concat([C_up, zeros, C_bo]).reset_index(drop=True)

    if constant is True:
        glob_c = pandas.DataFrame(data=numpy.ones(len(B)),
                                  columns=['0_constant'])
        if len(C.columns) == 0:
            C = glob_c
        else:
            C.insert(0, '0_constant', glob_c)

        C_names.insert(0, '0_constant')

    ##############################################################################
    # Prediction Data

    # Actual post-treatment series
    Y_post = rows_tr_post[[outcome_var]].reset_index(level='__ID', drop=True)
    Y_post.columns = unit_tr_list

    # Prediction matrix
    P = rows_co_post[[outcome_var, 'treated_unit',
                      '__ID', '__time']].pivot(index=['treated_unit', '__time'],
                                               columns='__ID',
                                               values=outcome_var)

    P = P[donor_order]

    # If the outcome variable is within the specified features then we need to
    # augment P with the corresponding (eventual) covariates used for adjustment,
    # If instead the outcome variable is not within the specified features
    # P is composed by the outcome variable of the donors only

    if out_in_features is True:
        # Check that global constant is required by the user
        if constant is True:
            P['0_constant'] = numpy.ones(len(rows_tr_post))

        # Add covariates used for adjustment in outcome variable equation (if present)
        if cov_adj is not None:

            for m in range(1, len(features) + 1):
                if not isinstance(cov_adj[0], list):
                    covs_adj = deepcopy(cov_adj)
                else:
                    covs_adj = deepcopy(cov_adj[m - 1])

                if m == 1:

                    if 'constant' in covs_adj:
                        const = numpy.ones(len(rows_tr_post))
                        const_n = str(m) + '_constant'
                        P.insert(len(P.columns), const_n, const)
                        covs_adj.remove('constant')

                    if 'trend' in covs_adj:
                        trend = period_post - period_pre[0] + 1
                        trend_n = str(m) + '_trend'
                        P.insert(len(P.columns), trend_n, trend)
                        covs_adj.remove('trend')

                    # Add other covariates from dataframe
                    rows_P = data_bal.loc[(unit_tr, unit_co[0], period_post), covs_adj].reset_index()
                    rows_P.set_index(P.index, inplace=True)
                    coln = [str(m) + '_' + cov_n for cov_n in covs_adj]
                    P[coln] = rows_P[covs_adj]
                else:
                    coln = [str(m) + '_' + cov_n for cov_n in covs_adj]
                    zeros = numpy.zeros((len(P.index), len(covs_adj)))
                    P[coln] = zeros

                m = m + 1

    T1 = len(period_post)

    ############################################################################
    # Proceed cleaning missing data in the pre-treatment period
    # Check if there are annoying donors with ALL missing values in the pre-treatment period
    empty_cols = numpy.all(numpy.isnan(B), axis=0)

    if sum(empty_cols) > 0:
        names = B.columns[empty_cols].tolist()
        names_print = ' '.join(str(name) for name in names)
        if verbose is True:
            warnings.warn("The following donors have no observations in the pre-treatment" +
                          " period, hence they have been removed! " + names_print)
        unit_co_eff = [co for co in unit_co if co not in names]
        B = B[unit_co_eff]
    else:
        unit_co_eff = unit_co

    A.set_index(B.index, drop=True, inplace=True)
    if len(C.columns) > 0:
        C.set_index(B.index, drop=True, inplace=True)

    X = A.join([B, C], how='outer', sort=True)
    # X = pandas.concat([A, B, C], axis=1)
    X_na = X.loc[complete_cases(X), ]

    A_na = X_na.loc[:, [unit_tr]]

    # take care of donors that have been dropped because of absent pre-treatment observations
    if len(unit_co_eff) < len(donor_order):
        C_cols_bool = [col not in donor_order for col in P.columns.tolist()]  # P has C columns to be taken care of
        C_cols = P.columns[C_cols_bool].tolist()

        donor_order = unit_co_eff

        Y_donors = Y_donors[donor_order]

        P_col_sel = donor_order + C_cols
        P = P[P_col_sel]

    B_na = X_na[donor_order]
    C_na = pandas.DataFrame(None)
    if len(C.columns) == 1:
        C_na = X_na.loc[:, C_names]
    elif len(C.columns) > 1:
        C_na = X_na.loc[:, C_names]

    # Store effective number of observations per feature
    T0_features = Counter(A_na.index.get_level_values('feature'))

    ############################################################################
    # Proceed cleaning missing data in the post-treatment period
    # P_na = P.loc[complete_cases(P), ]

    ############################################################################
    # Throw warnings for missing values
    # Report missing values in pre-treatment period

    if report_missing is True and len(A) != len(A_na):
        warnings.warn('Missing values detected in the pre-treatment period!')

        # Report missing values in A
        A_missing = complete_cases(A) is False  # identify missing rows
        if numpy.sum(A_missing) > 0:
            A_rows = A.loc[A_missing, ].index.to_frame()
            A_rows.columns = ['Feature', 'Time']
            print('Missing values detected in the following feature(s) of the treated unit:')
            print(A_rows.to_string(index=False))
        else:
            print('The feature(s) of the treated unit do not contain missing values')

        # Report missing values in B
        B_missing = complete_cases(B) is False  # identify missing rows
        if numpy.sum(B_missing) > 0:
            B_rows = B.loc[B_missing, ].index.to_frame()
            B_rows.columns = ['Feature', 'Time']
            print('Missing values detected in the following feature(s) of the donor pool:')
            print(B_rows.to_string(index=False))
        else:
            print('The feature(s) of the control units do not contain missing values')

        # Report missing values in C
        if C is not None:
            C_missing = complete_cases(C) is False  # identify missing rows
            if numpy.sum(C_missing) > 0:
                C_rows = C.loc[C_missing, ].index.to_frame()
                C_rows.columns = ['Feature', 'Time']
                print('Missing values detected in the following feature(s) of the treated unit:')
                print(C_rows.to_string(index=False))
            else:
                print('The covariate(s) used for adjustment do not contain missing values.')

    # elif report_missing == False and len(A) != len(A_na):
    #     warnings.warn('Missing values detected in the pre-treatment period!' +
    #                   'If you want to see where the missing values are located re-run the command' +
    #                   'with the option report_missing = True.')
    #     suggest = 0

    # Report missing values in post-treatment period
    P_missing = complete_cases(P) is False
    if numpy.sum(P_missing) > 0:
        warnings.warn('Missing values detected in the post-treatment period! Point estimate and prediction interval' +
                      'will not be computed for some of the required periods!')

        if report_missing is True:
            P_rows = P.loc[P_missing, ].index.to_frame()
            P_rows.columns = ['Feature', 'Time']
            print('Missing values detected in the data for post-treatment prediction in the following periods:')
            print(P_rows.to_string(index=False))
        # elif suggest == 1:
        #     print('If you want to see where the missing values are located re-run the command' +
        #               'with the option report_missing = True.')

    ############################################################################
    # Store objects

    # Size of donor pool
    J = len(unit_co_eff)

    # Total number of covariates used for adjustment
    if C is not None:
        KM = len(C.columns)
    else:
        KM = 0

    # Numer of features
    M = len(features)

    # Array containing number of covariates used for adjustment in each equation

    if cov_adj is None:
        K = pandas.DataFrame(numpy.zeros(M), index=features)
    elif not isinstance(cov_adj[0], list):
        K = pandas.DataFrame([len(cov_adj)])
    elif isinstance(cov_adj[0], list):
        k = [len(cov) for cov in cov_adj]
        K = pandas.DataFrame(k)

    if constant is True:
        K = K + 1

    K.reset_index(inplace=True)
    K.columns = ['Feature', 'Num of Covariates']
    K = K.astype({'Feature': str, 'Num of Covariates': int})

    A_na.columns = ['A']

    B_na.columns = unit_tr + '_' + B_na.columns
    if len(C_na.columns) > 0:
        C_na.columns = unit_tr + '_' + C_na.columns

    P.columns = unit_tr + '_' + P.columns
    Y_donors.columns = unit_tr + '_' + Y_donors.columns

    # convert indices from integer to ts if necessary
    if timeConvert is True:
        A_na.reset_index(drop=False, inplace=True)
        A_na['__time'] = A_na['__time'].map(int2ts)
        A_na.set_index(['treated_unit', 'feature', '__time'], drop=True, inplace=True)

        B_na.set_index(A_na.index, inplace=True)
        if len(C_na.columns) > 0:
            C_na.set_index(A_na.index, inplace=True)

        P.reset_index(drop=False, inplace=True)
        P['__time'] = P['__time'].map(int2ts)
        P.set_index(['treated_unit', '__time'], drop=True, inplace=True)

        Y_pre.reset_index(drop=False, inplace=True)
        Y_pre['__time'] = Y_pre['__time'].map(int2ts)
        Y_pre.set_index(['treated_unit', '__time'], drop=True, inplace=True)

        Y_post.reset_index(drop=False, inplace=True)
        Y_post['__time'] = Y_post['__time'].map(int2ts)
        Y_post.set_index(['treated_unit', '__time'], drop=True, inplace=True)

        Y_donors.reset_index(drop=False, inplace=True)
        Y_donors['__time'] = Y_donors['__time'].map(int2ts)
        Y_donors.set_index(['treated_unit', '__time'], drop=True, inplace=True)

        period_pre = pandas.Series(period_pre).map(int2ts).to_numpy()
        period_post = pandas.Series(period_post).map(int2ts).to_numpy()

    # make all index names homogenous
    A_na.index.names = ['ID', 'feature', 'Time']
    B_na.index.names = ['ID', 'feature', 'Time']
    if len(C_na.columns) > 0:
        C_na.index.names = ['ID', 'feature', 'Time']
    P.index.names = ['ID', 'Time']
    Y_pre.index.names = ['ID', 'Time']
    Y_post.index.names = ['ID', 'Time']
    Y_donors.index.names = ['ID', 'Time']

    if cointegrated_data is True:
        if any([t == 1 for t in T0_features.values()]) is True:
            if verbose is True:
                warnings.warn("You have at least one feature with only one pre-treatment period, " +
                              "thus you cannot specify cointegrated_data = True! Remember that this " +
                              "option uses the difference of pre-treatment residuals rather than " +
                              "their levels. We set cointegrated_data = False.")
            cointegrated_data = False

    return scdata_output(A=A_na, B=B_na, C=C_na, P=P, Y_pre=Y_pre,
                         Y_post=Y_post, Y_donors=Y_donors, J=J, K=K,
                         KM=KM, KMI=KM, M=M, iota=1, cointegrated_data=cointegrated_data,
                         period_pre=period_pre, period_post=period_post,
                         T0_features=T0_features, T1_outcome=T1,
                         outcome_var=outcome_var, features=features,
                         glob_cons=constant, out_in_features=out_in_features,
                         treated_units=[unit_tr], donors_units=unit_co_eff, units_est=[unit_tr],
                         anticipation={unit_tr: anticipation}, effect="unit-time",
                         timeConvert=timeConvert)


class scdata_output:
    def __init__(self, A, B, C, P, Y_pre, Y_post, Y_donors, J, K, KM, KMI, M, iota,
                 cointegrated_data, period_pre, period_post, T0_features,
                 T1_outcome, outcome_var, features, glob_cons, out_in_features,
                 treated_units, donors_units, units_est, anticipation, effect, timeConvert):

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
        self.period_pre = period_pre
        self.period_post = period_post
        self.T0_features = T0_features
        self.T1_outcome = T1_outcome
        self.outcome_var = outcome_var
        self.features = features
        self.glob_cons = glob_cons
        self.out_in_features = out_in_features
        self.treated_units = treated_units
        self.donors_units = donors_units
        self.units_est = units_est
        self.anticipation = anticipation
        self.effect = effect
        self.timeConvert = timeConvert
