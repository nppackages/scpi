# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 09:59:24 2021

@author: Filippo Palomba
"""
# Temporary code to suppress pandas FutureWarning
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas
import numpy
import scipy.linalg
import warnings
from copy import deepcopy
from collections import Counter
from .funs import complete_cases


class scdata_output:
    def __init__(self, A, B, C, P, Y_pre, Y_post, Y_donors, J, K, KM, M,
                 cointegrated_data, period_pre, period_post, T0_features,
                 T1_outcome, outcome_var, features, glob_cons, out_in_features):

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
        self.outcome_var = outcome_var
        self.features = features
        self.glob_cons = glob_cons
        self.out_in_features = out_in_features


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
           report_missing=False):

    '''

    The command prepares the data to be used by scest or scpi for point estimation and inference procedures using
    Synthetic Control.
    It allows the user to specify the outcome variable, the features of the treated unit to be
    matched, and covariate-adjustment feature by feature. The names of the output matrices
    follow the notation proposed in Cattaneo, Feng, and Titiunik (2021).

    Parameters
    ----------
     df
     a pandas.dataframe object.

     id_var
     a character with the name of the variable containing units' IDs. The ID variable can be numeric or string.

     time_var
     a character with the name of the time variable. The time variable has to be numpy.int64, or one of
     pandas.Timestamp and numpy.datetime64. Input a numeric time variable is suggested when working with yearly
     data, whereas for all other frequencies either pandas.Timestamp or numpy.datetime64 types are preferred.

     outcome_var
     a character with the name of the outcome variable. The outcome variable has to be numeric.

     period_pre
     a numeric vector that identifies the pre-treatment period in time_var.

     period_post
     a numeric vector that identifies the post-treatment period in time_var.

     unit_tr
     a scalar that identifies the treated unit in id_var.

     unit_co
     a numeric vector that identifies the donor pool in id_var.

     features
     a character lsit containing the name of the feature variables used for estimation.
     If this option is not specified the default is features = outcome_var.

     cov_adj
     a list specifying the names of the covariates to be used for adjustment for each feature. If the user wants
     to specify the same set of covariates for all features, a single list should be provided. If instead a
     different set of covariates per feature has to be specified, then a list of lists should be provided. Note that
     in this latter case the number of sub-lists must be equal to the number of features. Moreover, the order of the
     sub-lists matters, in the sense that the first sub-list is interpreted as the set of covariates for the first
     feature, and so on. Finally, the user can specify 'constant' and 'trend' as covariates even if they are not
     present in the loaded dataframe.

     constant
     a logical which controls the inclusion of a constant term across features. The default value is False.

     cointegrated_data
     a logical that indicates if there is a belief that the data is cointegrated or not. The default value is False.
     See the Details section for more.

     anticipation
     a scalar that indicates the number of periods of potential anticipation effects. Default is 0.

     report_missing
     a logical which prints the location of missing values if present. The default value is False.

    Returns
    -------

    A
    a dataframe containing pre-treatment features of the treated unit.

    B
    a dataframe containing pre-treatment features of the control units.

    C
    a dataframe containing covariates for adjustment.

    P
    a dataframe whose rows are the vectors used to predict the out-of-sample series for the synthetic unit.

    Y_pre
    a dataframe containing the pre-treatment outcome of the treated unit.

    Y_post
    a dataframe containing the post-treatment outcome of the treated unit.

    Y_donors
    a dataframe containing the pre-treatment outcome of the control units.

    J
    the number of control units

    K
    a numeric array with the number of covariates used for adjustment for each feature

    KM
    the total number of covariates used for adjustment

    M
    number of features

    period_pre
    a numeric array with the pre-treatment period

    period_post
    a numeric array with the post-treatment period

    T0_features
    a numeric array with the number of periods used in estimation for each feature

    T1_outcome
    the number of post-treatment periods

    glob_cons
    for internal use only

    out_in_features
    for internal use only

    References
    ----------
    Cattaneo, M. D., Feng, Y., and Titiunik, R. (2021), “Prediction Intervals for Synthetic Control
    Methods,” Journal of the American Statistical Association, 116, 1865-1880.

    Cattaneo, M. D., Palomba, F., Feng, Y., and Titiunik, R. (2022), “Uncertainty Quantification in
    Synthetic Controls with Staggered Treatment Adoption,” working paper.

    See Also
    --------
    scest, scpi, scplot

    '''
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

    if outcome_var not in var_names:
        raise Exception("Outcome variable (outcome_var) not found in the input dataframe!")

    # Make time and id columns if these variables are indexes of the dataframe
    if id_var in indexes:
        data['__ID'] = data.index.get_level_values(id_var)
    else:
        data.rename(columns={id_var: '__ID'}, inplace=True)

    if time_var in indexes:
        data['__time'] = data.index.get_level_values(time_var)
    else:
        data.rename(columns={time_var: '__time'}, inplace=True)

    # Create index
    data.set_index(['__ID', '__time'], drop=False, inplace=True)

    dd = data.iloc[0, data.columns.get_loc('__time')]
    if not isinstance(dd, (numpy.int64, numpy.int32, numpy.int16, pandas.Timestamp, numpy.datetime64)):
        raise Exception("The object time_var should be of type int, numpy.datetime64, or pandas.Timestamp!")

    if not isinstance(period_pre, numpy.ndarray):
        raise Exception("The object period_pre should be of type numpy.ndarray (eg. use numpy.arange or numpy.array)!")

    if not isinstance(period_post, numpy.ndarray):
        raise Exception("The object period_post should be of type numpy.ndarray!")

    if not isinstance(period_pre[0], (numpy.int64, numpy.int32, numpy.int16, numpy.datetime64, pandas.Timestamp)):
        raise Exception("Elements of period_pre should either be of type int or numpy.datetime64!")

    if not isinstance(period_post[0], (numpy.int64, numpy.int32, numpy.int16, numpy.datetime64, pandas.Timestamp)):
        raise Exception("Elements of period_post should either be of type int or numpy.datetime64!")

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

    rows_tr_pre = data_bal.loc[(unit_tr_list, period_pre), ]
    rows_tr_post = data_bal.loc[(unit_tr_list, period_post), ]
    rows_co_pre = data_bal.loc[(unit_co, period_pre), ]
    rows_co_post = data_bal.loc[(unit_co, period_post), ]

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
                       id_vars=['__ID', '__time'],
                       value_vars=features,
                       var_name='feature',
                       value_name=label)

    A = A_df.loc[:, [unit_tr]]

    # Create B

    # Stack features one on top of the other
    B_df = pandas.melt(rows_co_pre,
                       id_vars=['__ID', '__time'],
                       value_vars=features,
                       var_name='feature',
                       value_name='value')

    # make the df wide so that countries are one next to the other
    B = B_df.pivot(index=['feature', '__time'],  # index to keep fixed
                   columns='__ID',                # dimension to wide
                   values='value')               # values to be widen

    # Create list that dictates order of donors
    donor_order = B.columns.values.tolist()

    # Create matrix with outcome for the donors
    sel = rows_co_pre[[outcome_var, '__time', '__ID']]
    Y_donors = sel.pivot(index='__time', columns='__ID', values=outcome_var)
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
            rows_C = data_bal.loc[(unit_co[0], period_pre), covs_adj].reset_index()
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
                rows_C = data_bal.loc[(unit_co[0], period_pre), covs_adj].reset_index()
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
                    i_start = len(rows_tr_pre) * m - 1
                    zeros = pandas.DataFrame(numpy.zeros((len(rows_tr_pre), len(C.columns))),
                                             columns=C_names)
                    if i_start < 0:
                        C_up = pandas.DataFrame(None)
                    else:
                        C_up = C.iloc[:i_start]

                    C_bo = C.iloc[(i_start + 1):]
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
    P = rows_co_post[[outcome_var, '__ID', '__time']].pivot(index='__time',
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

            if not isinstance(cov_adj[0], list):
                covs_adj = deepcopy(cov_adj)
                if 'constant' in covs_adj:
                    const = numpy.ones(len(rows_tr_post))
                    const_n = str(1) + '_constant'
                    P.insert(len(P.columns), const_n, const)
                    covs_adj.remove('constant')

                if 'trend' in covs_adj:
                    trend = period_post - period_pre[0] + 1
                    trend_n = str(1) + '_trend'
                    P.insert(len(P.columns), trend_n, trend)
                    covs_adj.remove('trend')

                # Add other covariates from dataframe
                rows_P = data_bal.loc[(unit_co[0], period_post), covs_adj].reset_index()
                rows_P.set_index(P.index, inplace=True)
                coln = [str(1) + '_' + cov_n for cov_n in covs_adj]
                P[coln] = rows_P[covs_adj]

            else:
                m = 1
                for covs_adj in cov_adj:
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
                        rows_P = data_bal.loc[(unit_co[0], period_post), covs_adj].reset_index()
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

    A.set_index(B.index, drop=True, inplace=True)
    if len(C.columns) > 0:
        C.set_index(B.index, drop=True, inplace=True)

    X = pandas.concat([A, B, C], axis=1)
    X_na = X.loc[complete_cases(X), ]

    A_na = X_na.loc[:, [unit_tr]]
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
    J = len(unit_co)

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

    return scdata_output(A=A_na, B=B_na, C=C_na, P=P, Y_pre=Y_pre,
                         Y_post=Y_post, Y_donors=Y_donors, J=J, K=K,
                         KM=KM, M=M, cointegrated_data=cointegrated_data,
                         period_pre=period_pre, period_post=period_post,
                         T0_features=T0_features, T1_outcome=T1,
                         outcome_var=outcome_var, features=features,
                         glob_cons=constant, out_in_features=out_in_features)
