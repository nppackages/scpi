# -*- coding: utf-8 -*-

# Temporary code to suppress pandas FutureWarning
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas
pandas.options.mode.chained_assignment = None
import pandas
import numpy
from math import ceil
from copy import deepcopy
from .scdata import scdata


def scdataMulti(df,
                id_var,
                time_var,
                outcome_var,
                treatment_var,
                features=None,
                cov_adj=None,
                cointegrated_data=False,
                post_est=None,
                units_est=None,
                donors_est=None,
                anticipation=0,
                effect="unit-time",
                constant=False,
                verbose=True):

    """
    Parameters
    ----------
    df : pandas.DataFrame
        a dataframe object containing the data to be processed.

    id_var : str
        a character with the name of the variable containing units' IDs.

    time_var : str
        a character with the name of the time variable. The time variable has to be numpy.int64,
        numpy.datetime64, or pandas.Timestamp. Input a numeric time variable is suggested when working with
        yearly data, whereas for all other frequencies numpy.datetime64 type is preferred.

    outcome_var : str
        a character with the name of the outcome variable. The outcome variable has to be numeric.

    treatment_var : str
        a character with the name of the treatment variable. The treatment variable has to be 1 only in periods
        where a unit is treated.

    features : dict, default None
        a dictionary whose elements are lists containing the name of the feature variables used for estimation.
        If a dictionary with a single key
        is provided as input, then the same features are used for all treated units. Alternatively, if the user wants to
        specify different features for different treated units, the dictionary must contain as many keys as the number
        of treated units in the data. Each key must correspond to the identifier (id_var) of one treated unit.
        If this option is not specified the default is features = outcome_var.

    cov_adj : dict, default None
        a dictionary whose elements are lists containing the name of the covariates used for adjustment.
        If a dictionary with a single key
        is provided as input, then the same covariates are used for adjustment for all treated units.
        Alternatively, if the user wants to specify different covariates for different treated units,
        the dictionary must contain as many keys as the number of treated units in the data. Each key must
        correspond to the identifier (id_var) of one treated unit.

        More in detail, if the user wants
        to specify the same set of covariates for all features, a single list should be provided. If instead a
        different set of covariates per feature has to be specified, then a list of lists should be provided. Note that
        in this latter case the number of sub-lists must be equal to the number of features. Moreover, the order of the
        sub-lists matters, in the sense that the first sub-list is interpreted as the set of covariates for the first
        feature, and so on. Finally, the user can specify 'constant' and 'trend' as covariates even if they are not
        present in the loaded dataframe.

    post_est : int/str, default None
        an integer or string specifying the number of post-treatment periods for which treatment effects have to be estimated for each
        treated unit. It must be an integer when time_var is integer, otherwise it must be a string of the form "10 years", "2 months",
        "1 day" and so on. Possible options are: 'year(s)', 'month(s)', 'week(s)', 'day(s), and 'hour(s)'.
        It is only effective when effect = "unit-time".

    units_est : list, default None
        a list specifying the treated units for which treatment effects have to be estimated.

    donors_est : dict, default None
        a dictionary specifying the donors units to be used. If the dictionary has length 1, then all treated units share the same
        potential donors. Otherwise, if the user requires different donor pools for different treated units, the dictionary must
        be of the same length of the number of treated units and each element has to be named with one treated unit's name as
        specified in id_var.

    constant : bool/dict, default False
        a logical which controls the inclusion of a constant term across features. If the user wants to specify this
        option indipendently for each treated unit, a dictionary must be provided instead of a boolean value.
        Specifically, the dictionary must contain as many keys as the number of treated units in the data.
        Each key must correspond to the identifier (id_var) of one treated unit.

    cointegrated_data : bool/dict, default False
        a logical that indicates if there is a belief that the data is cointegrated or not. If the user wants to specify
        this option indipendently for each treated unit, a dictionary must be provided instead of a boolean value.
        Specifically, the dictionary must contain as many keys as the number of treated units in the data.
        Each key must correspond to the identifier (id_var) of one treated unit.

    effect : str, default "unit-time"
        a string indicating the type of treatment effect to be estimated. Options are: 'unit-time', which estimates
        treatment effects for each treated unit- post treatment period combination; 'unit', which estimates the
        treatment effect for each unit by averaging post-treatment features over time; 'time', which estimates the
        average treatment effect on the treated at various horizons.

    anticipation : int/dict, default 0
        a scalar that indicates the number of periods of potential anticipation effects. If the user wants to specify
        this option indipendently for each treated unit, a dictionary must be provided instead of an integer value.
        Specifically, the dictionary must contain as many keys as the number of treated units in the data.
        Each key must correspond to the identifier (id_var) of one treated unit.

    Returns
    -------
    The function returns an object of class 'scdata_output' containing the following objects

    A : pandas.DataFrame
        a dataframe containing pre-treatment features of the treated units.

    B : pandas.DataFrame
        a dataframe containing pre-treatment features of the control units.

    C : pandas.DataFrame
        a dataframe containing covariates for adjustment.

    P : pandas.DataFrame
        a dataframe whose rows are the vectors used to predict the out-of-sample series for the synthetic units.

    Y_df : pandas.DataFrame
        a dataframe containing the outcome variable for all units.

    Y_pre: pandas.DataFrame
        a dataframe containing the actual pre-treatment outcome for the treated unit(s). Note that this is the raw data,
        therefore if effect is specified, it will not contain the aggregated data.

    Y_post: pandas.DataFrame
        a dataframe containing the actual post-treatment outcome for the treated unit(s). Note that this is the raw data,
        therefore if effect is specified, it will not contain the aggregated data.

    Y_donors : pandas.DataFrame
        a dataframe containing the pre-treatment outcome of the control units.

    J : dict
        a dictionary containing the number of donors for each treated unit

    K : dict
        a dictionary containing the number of covariates used for adjustment for each feature for each treated unit

    KM : dict
        a dictionary containing the total number of covariates used for adjustment for each treated unit

    M : dict
        a dictionary containing number of features used for each treated unit

    iota : int
        number of treated units

    KMI : int
        overall number of covariates used for adjustment

    period_pre : dict
        a dictionary containing a numeric vector with the pre-treatment period for each treated unit

    period_post : dict
        a dictionary containing a numeric vector with the post-treatment period for each treated unit

    T0_features : dict
        a dictionary containing a numeric vector with the number of periods used in estimation for each feature for each
        treated unit

    T1_outcome : dict
        a dictionary containing the number of post-treatment periods for each treated unit

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

    Cattaneo, M. D., Palomba, F., Feng, Y., and Titiunik, R. (2025), “scpi: Uncertainty Quantification for
    Synthetic Control Estimators”, Journal of Statistical Software, 113(2), 1-38.

    Cattaneo, M. D., Palomba, F., Feng, Y., and Titiunik, R. (2025), “Uncertainty Quantification in Synthetic
    Controls with Staggered Treatment Adoption”, Review of Economic Studies, doi:10.1080/01621459.2021.1979561.

    See Also
    --------
    scdata, scest, scpi, scplot, scplotMulti

    """

    # Error Checking

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

    if not isinstance(treatment_var, str):
        raise Exception("You should specify the name of treatment_var as a string! (eg. treatment_var = 'treatment')")

    if features is not None:
        if not isinstance(features, dict):
            raise Exception("The object 'features' should be a dictionary!")
        features = {k.replace('_', ' '): v for k, v in features.items()}

    else:
        features = {'features': [outcome_var]}

    if cov_adj is not None:
        if not isinstance(cov_adj, dict):
            raise Exception("The object 'cov_adj' should be a dictionary!")

    if effect not in ["unit", "unit-time", "time"]:
        raise Exception("The object 'effect' should be either 'unit', 'time', or 'unit-time'!")

    # Check variables are in dataframe (id and time can be either variables or indices)
    if (id_var not in var_names) and (id_var not in indexes):
        raise Exception("ID variable (id_var) not found in the input dataframe neither as a variable nor as an index!")

    if (time_var not in var_names) and (time_var not in indexes):
        raise Exception("Time variable (time_var) not found in the input df neither as a variable nor as an index!")

    if outcome_var not in var_names:
        raise Exception("Outcome variable (outcome_var) not found in the input dataframe!")

    if treatment_var not in var_names:
        raise Exception("Treatment variable (treatment_var) not found in the input dataframe!")

    # Make time and id columns if these variables are indexes of the dataframe and rename variables
    if id_var in indexes:
        data['__ID'] = data.index.get_level_values(id_var)
    else:
        data.rename(columns={id_var: '__ID'}, inplace=True)

    # CVXR does not like '_' symbols
    if pandas.api.types.is_string_dtype(data['__ID']) is False:  # convert to string if not string
        data['__ID'] = data['__ID'].astype(str)
    data['__ID'] = data['__ID'].str.replace('_', ' ')

    if time_var in indexes:
        data['__time'] = data.index.get_level_values(time_var)
    else:
        data.rename(columns={time_var: '__time'}, inplace=True)

    data.rename(columns={treatment_var: '__Treatment'}, inplace=True)

    Y_df = data[['__ID', '__time', '__Treatment', outcome_var]]

    # Check time_var type and eventually convert it
    timeConvert = False
    dd = data.iloc[0, data.columns.get_loc('__time')]
    if not isinstance(dd, (numpy.int64, numpy.int32, numpy.int16, numpy.datetime64, pandas.Timestamp)):
        raise Exception("The object time_var should be of type int, pandas.Timestamp, or numpy.datetime64!")

    elif isinstance(dd, (numpy.datetime64, pandas.Timestamp)):
        time_unique_ts = sorted(set(data['__time'].tolist()))
        int2ts = {i: time_unique_ts[i] for i in range(len(time_unique_ts))}
        ts2int = {time_unique_ts[i]: i for i in range(len(time_unique_ts))}
        data['__time'] = data['__time'].map(ts2int)
        timeConvert = True

    if post_est is not None:
        if not isinstance(post_est, (int, str)):
            raise Exception("You should specify post_est as an integer or a string!")

        if not isinstance(dd, (numpy.datetime64, pandas.Timestamp)):
            if not isinstance(post_est, int):
                raise Exception("You should specify post_est as an integer!")
        else:
            if not isinstance(post_est, str):
                raise Exception("You should specify post_est as a string!")
            aux = post_est.split(' ')
            if len(aux) != 2:
                raise Exception("You should specify post_est as a string of the form (e.g.) '10 years'!")
            post_est_delta = aux[0]
            post_est_freq = aux[1][0].upper()
            if post_est_freq == "H":
                post_est_freq = "h"

    # Identify treated units
    periods_treated = data[['__Treatment', '__ID']].groupby('__ID').sum()
    treated_units = periods_treated.loc[periods_treated.values > 0, ].index.values.tolist()
    treated_post = deepcopy(treated_units)

    if units_est is not None:
        if not isinstance(units_est, list):
            raise Exception("The object 'units_est' should be a list!")

        if not all(tr in treated_units for tr in units_est):
            raise Exception("The object 'units_est' must contain the identifiers (id_var) of the treated" +
                            " units for which treatment effects have to be estimated!")

        treated_units = [tr for tr in treated_units if tr in units_est]

    else:
        units_est = treated_units

    # Control covariates for adjustment and matching features
    if cov_adj is not None:
        if len(cov_adj) > 1:
            if len(cov_adj) != len(treated_units):
                raise Exception("If you want to specify covariate adjustment separately" +
                                " for each treated unit, make sure that 'cov.adj' has" +
                                " the same number of elements as there are treated" +
                                " units (" + str(len(treated_units)) + ")!")

            cov_adj = {k.replace('_', ' '): v for k, v in cov_adj.items()}

            names_dict = []
            for n, v in cov_adj.items():
                names_dict.append(n)
            if not all(tr in names_dict for tr in treated_units):
                tr_not_found = [tr for tr in treated_units if tr not in names_dict]
                tr_print = ' '.join(str(tr) for tr in tr_not_found)
                raise Exception("There is no match in the object 'cov_adj' for the " +
                                "following treated units: " + tr_print)

    if len(features) > 1:
        if len(features) != len(treated_units):
            raise Exception("If you want to specify features separately" +
                            " for each treated unit, make sure that 'features' has" +
                            " the same number of elements as there are treated" +
                            " units (" + str(len(treated_units)) + ")!")
        names_dict = []
        for n, v in features.items():
            names_dict.append(n)
        if not all(tr in names_dict for tr in treated_units):
            tr_not_found = [tr for tr in treated_units if tr not in names_dict]
            tr_print = ' '.join(str(tr) for tr in tr_not_found)
            raise Exception("There is no match in the object 'features' for the " +
                            "following treated units: " + tr_print)

    if not isinstance(constant, bool):
        if not isinstance(constant, dict):
            raise Exception("If you want to specify the presence of a constant separately" +
                            " for each treated unit then 'constant' has to be a dictionary!")

        if len(constant) != len(treated_units):
            raise Exception("If you want to specify the presence of a constant separately" +
                            " for each treated unit, make sure that 'constant' has" +
                            " the same number of elements as there are treated" +
                            " units (" + str(len(treated_units)) + ")!")
        names_dict = []
        for n, v in constant.items():
            names_dict.append(n)
        if not all(tr in names_dict for tr in treated_units):
            tr_not_found = [tr for tr in treated_units if tr not in names_dict]
            tr_print = ' '.join(str(tr) for tr in tr_not_found)
            raise Exception("There is no match in the object 'constant' for the " +
                            "following treated units: " + tr_print)

    if not isinstance(cointegrated_data, bool):
        if not isinstance(cointegrated_data, dict):
            raise Exception("If you want to specify the presence of cointegration separately" +
                            " for each treated unit then 'cointegrated_data' has to be a dictionary!")

        if len(cointegrated_data) != len(treated_units):
            raise Exception("If you want to specify the presence of cointegration separately" +
                            " for each treated unit, make sure that 'cointegrated_data' has" +
                            " the same number of elements as there are treated" +
                            " units (" + str(len(treated_units)) + ")!")
        names_dict = []
        for n, v in cointegrated_data.items():
            names_dict.append(n)
        if not all(tr in names_dict for tr in treated_units):
            tr_not_found = [tr for tr in treated_units if tr not in names_dict]
            tr_print = ' '.join(str(tr) for tr in tr_not_found)
            raise Exception("There is no match in the object 'cointegrated_data' for the " +
                            "following treated units: " + tr_print)

    if not isinstance(anticipation, int):
        if not isinstance(anticipation, dict):
            raise Exception("If you want to specify the presence of anticipation effects separately" +
                            " for each treated unit then 'anticipation' has to be a dictionary!")

        if len(anticipation) != len(treated_units):
            raise Exception("If you want to specify the presence of anticipation effects separately" +
                            " for each treated unit, make sure that 'anticipation' has" +
                            " the same number of elements as there are treated" +
                            " units (" + str(len(treated_units)) + ")!")
        names_dict = []
        for n, v in anticipation.items():
            names_dict.append(n)
        if not all(tr in names_dict for tr in treated_units):
            tr_not_found = [tr for tr in treated_units if tr not in names_dict]
            tr_print = ' '.join(str(tr) for tr in tr_not_found)
            raise Exception("There is no match in the object 'anticipation' for the " +
                            "following treated units: " + tr_print)

    if donors_est is not None:
        if not isinstance(donors_est, dict):
            raise Exception("The option 'donors_est' must be of type dictionary!")

        if len(donors_est) != 1 and len(donors_est) != len(treated_units):
            raise Exception("The option 'donors_est' must be a dictionary of either length 1" +
                            " or " + str(len(treated_units)) + " (the number of treated units" +
                            " for which treatment effects have to be computed)!")

        if len(donors_est) > 1:
            names_dict = []
            for n, v in donors_est.items():
                names_dict.append(n)
            if not all(tr in names_dict for tr in treated_units):
                raise Exception("If len(donors.est) > 1, all the names of the elements have to be values of 'id_var'!")

    # Data preparation
    # Get first treated period of each treated unit
    aux = data.loc[data['__Treatment'] == 1, ['__ID', '__time']]
    treated_periods = aux.groupby('__ID').min()

    tr_count = 1
    for treated_unit in treated_units:

        # parse options
        if cov_adj is not None:
            if len(cov_adj) == 1:
                cov_adj_tr = list(cov_adj.values())[0]
            else:
                cov_adj_tr = cov_adj[treated_unit]
        else:
            cov_adj_tr = None

        if len(features) == 1:
            features_tr = list(features.values())[0]
        else:
            features_tr = features[treated_unit]

        if isinstance(constant, bool):
            constant_tr = constant
        else:
            constant_tr = constant[treated_unit]

        if isinstance(cointegrated_data, bool):
            cointegrated_data_tr = cointegrated_data
        else:
            cointegrated_data_tr = cointegrated_data[treated_unit]

        if isinstance(anticipation, int):
            anticipation_tr = anticipation
        else:
            anticipation_tr = anticipation[treated_unit]

        # parse data
        treated_unit_T0 = treated_periods.loc[treated_unit, ][0]  # treatment date

        # get all other units before treatment of treated unit
        donors = data[numpy.invert(data['__ID'].isin(treated_post)) &
                      (data['__time'] < treated_unit_T0)]

        if len(donors) == 0:
            raise Exception("The current specification for " + treated_unit + " does not have observations!")

        # number of periods units have been treated before treatment of treated unit
        donors_count = donors[['__ID', '__Treatment']].groupby('__ID').sum()
        donors_units = donors_count[donors_count['__Treatment'] == 0].index.values.tolist()

        if post_est is not None:
            if timeConvert is False:
                T1_last = treated_unit_T0 + post_est - anticipation_tr
            else:
                T1_last = treated_unit_T0 + numpy.timedelta64(post_est_delta, post_est_freq) - anticipation_tr

            treated_donors = data[(data['__ID'].isin(treated_post)) &
                                  (data['__time'] < T1_last)]
            tr_donors_count = treated_donors[['__ID', '__Treatment']].groupby('__ID').sum()
            tr_donors_units = tr_donors_count[tr_donors_count['__Treatment'] == 0].index.values.tolist()

            # not-yet-treated are used as donors only with unit-time predictands
            if effect in ["unit-time", "unit"]:
                donors_units.extend(tr_donors_units)
                donors_units = (list(set(donors_units)))
                donors_units.sort()

        if donors_est is not None:
            if len(donors_est) == 1:
                donors_filter = donors_est[list(donors_est.keys())[0]]
            else:
                donors_filter = donors_est[treated_unit]

            donors_units = list(set(donors_units) & set(donors_filter))

        # subset dataset selecting treated units and proper donors
        df_aux = data[data['__ID'].isin(donors_units + [treated_unit])]

        # create time arrays
        time_array = df_aux['__time'].unique()
        period_pre = time_array[time_array < treated_unit_T0]
        period_post = time_array[time_array >= treated_unit_T0]

        if post_est is not None:
            sel_post = period_post < T1_last
            period_post = period_post[sel_post]

        if timeConvert is True:  # convert back to time series format to trigger hash maps in scdata
            df_aux['__time'] = df_aux['__time'].map(int2ts)
            period_pre = pandas.Series(period_pre).map(int2ts).to_numpy()
            period_post = pandas.Series(period_post).map(int2ts).to_numpy()

        try:
            scdata_out = scdata(df=df_aux,
                                id_var='__ID',
                                time_var='__time',
                                outcome_var=outcome_var,
                                period_pre=period_pre,
                                period_post=period_post,
                                unit_tr=treated_unit,
                                unit_co=donors_units,
                                features=features_tr,
                                cov_adj=cov_adj_tr,
                                constant=constant_tr,
                                cointegrated_data=cointegrated_data_tr,
                                anticipation=anticipation_tr)
        except Exception as e:
            str1 = "There is a problem with your specification for the treated unit: " + treated_unit
            str2 = ". Here is the original error message: "
            raise Exception(str1 + str2 + str(e))

        # Store data matrices
        A_tr = scdata_out.A
        B_tr = scdata_out.B
        C_tr = scdata_out.C

        if len(C_tr) == 0:  # to avoid bad behavior of block_diag when C is null
            C_tr = A_tr[[]]

        Y_donors_tr = scdata_out.Y_donors
        P_tr = scdata_out.P

        if effect == "time":
            time = P_tr.index.get_level_values('Time').tolist()
            if timeConvert is True:
                time = [ts2int[t] for t in time]
            time = [t - min(time) + 1 for t in time]
            P_tr = pandas.DataFrame(P_tr.values,
                                    index=time,
                                    columns=P_tr.columns)

        if not (effect == "unit" and cointegrated_data_tr is True):
            P_diff = None

        if effect == "unit":  # average within unit
            if cointegrated_data_tr is True:  # differentiate the data if cointegration
                JJ = scdata_out.J
                if scdata_out.out_in_features is False:
                    P_first = P_tr.iloc[[0], :JJ] - Y_donors_tr.iloc[[len(Y_donors_tr) - 1], :].values
                    P_diff = P_tr.iloc[:, :JJ].diff()
                    P_diff.iloc[0, :] = P_first

                elif scdata_out.out_in_features is True:
                    # Remove last observation of first feature from first period of P
                    P_first = P_tr.iloc[[0], :JJ] - B_tr.iloc[[scdata_out.T0_features[outcome_var] - 1], :].values
                    P_diff = P_tr.iloc[:, :JJ].diff()
                    P_diff.iloc[0, :] = P_first
                    P_diff = pandas.concat([P_diff.iloc[:, :JJ],
                                            P_tr.iloc[1:, JJ:]], axis=1)
                    aux = numpy.array([P_diff.mean(axis=0)])
                    time = scdata_out.period_post[ceil(scdata_out.T1_outcome / 2) - 1]
                    idx = pandas.MultiIndex.from_product([[treated_unit], [time]],
                                                         names=['ID', 'Time'])
                    P_diff = pandas.DataFrame(aux,
                                              index=idx,
                                              columns=P_diff.columns)

            aux = numpy.array([P_tr.mean(axis=0)])
            time = scdata_out.period_post[ceil(scdata_out.T1_outcome / 2) - 1]
            idx = pandas.MultiIndex.from_product([[treated_unit], [time]],
                                                 names=['ID', 'Time'])
            P_tr = pandas.DataFrame(aux,
                                    index=idx,
                                    columns=P_tr.columns)

        if tr_count == 1:
            A_stacked = deepcopy(A_tr)
            B_stacked = deepcopy(B_tr)
            C_stacked = deepcopy(C_tr)
            P_stacked = deepcopy(P_tr)
            Pd_stacked = deepcopy(P_diff)
            Y_donors_stacked = deepcopy(Y_donors_tr)

            J_dict = {treated_unit: scdata_out.J}
            K_dict = {treated_unit: scdata_out.K}
            KM_dict = {treated_unit: scdata_out.KM}
            M_dict = {treated_unit: scdata_out.M}
            Y_pre_dict = {treated_unit: scdata_out.Y_pre}
            Y_post_dict = {treated_unit: scdata_out.Y_post}
            period_pre_dict = {treated_unit: scdata_out.period_pre}
            period_post_dict = {treated_unit: scdata_out.period_post}
            T0_features_dict = {treated_unit: scdata_out.T0_features}
            T1_dict = {treated_unit: scdata_out.T1_outcome}
            out_in_features_dict = {treated_unit: scdata_out.out_in_features}
            constant_dict = {treated_unit: scdata_out.glob_cons}
            cointegrated_data_dict = {treated_unit: scdata_out.cointegrated_data}
            donors_dict = {treated_unit: scdata_out.donors_units}
            anticipation_dict = {treated_unit: anticipation_tr}

        else:
            A_stacked = pandas.concat([A_stacked, A_tr], axis=0)
            B_stacked = pandas.concat([B_stacked, B_tr], axis=0)
            C_stacked = pandas.concat([C_stacked, C_tr], axis=0)
            if effect == "time":
                P_stacked = pandas.concat([P_stacked, P_tr], axis=1, join='inner')  # stack horizontally
            else:
                P_stacked = pandas.concat([P_stacked, P_tr], axis=0)  # stack diagonally
            if Pd_stacked is not None:
                Pd_stacked = pandas.concat([Pd_stacked, P_diff], axis=0)
            Y_donors_stacked = pandas.concat([Y_donors_stacked, Y_donors_tr], axis=0)

            J_dict[treated_unit] = scdata_out.J
            K_dict[treated_unit] = scdata_out.K
            KM_dict[treated_unit] = scdata_out.KM
            M_dict[treated_unit] = scdata_out.M
            Y_pre_dict[treated_unit] = scdata_out.Y_pre
            Y_post_dict[treated_unit] = scdata_out.Y_post
            period_pre_dict[treated_unit] = scdata_out.period_pre
            period_post_dict[treated_unit] = scdata_out.period_post
            T0_features_dict[treated_unit] = scdata_out.T0_features
            T1_dict[treated_unit] = scdata_out.T1_outcome
            out_in_features_dict[treated_unit] = scdata_out.out_in_features
            constant_dict[treated_unit] = scdata_out.glob_cons
            cointegrated_data_dict[treated_unit] = scdata_out.cointegrated_data
            donors_dict[treated_unit] = scdata_out.donors_units
            anticipation_dict[treated_unit] = anticipation_tr

        tr_count = tr_count + 1

    B_stacked.set_index(A_stacked.index, inplace=True)
    C_stacked.set_index(A_stacked.index, inplace=True)

    B_stacked.fillna(0, inplace=True)
    C_stacked.fillna(0, inplace=True)
    P_stacked.fillna(0, inplace=True)
    if Pd_stacked is not None:
        Pd_stacked.fillna(0, inplace=True)
    Y_donors_stacked.fillna(0, inplace=True)

    # Rearrange P so that order of columns coincides with (B,C)
    bcols = B_stacked.columns.tolist()
    ccols = C_stacked.columns.tolist()
    P_stacked = P_stacked[bcols + ccols]

    if effect == "time":
        P_stacked = P_stacked / len(treated_units)
        # overwrite T1_dict because now T1 is min(T1_i) where i are treated units
        T1min = min(T1_dict.values())
        for tr in treated_units:
            T1_dict[tr] = T1min

    # number of treated units and number of covariates used for adjustment
    iota = len(treated_units)
    KMI = len(C_stacked.columns)

    # transform outcome dictionaries in dataframes
    Y_pre_df = pandas.DataFrame(columns=['Actual'])
    Y_post_df = pandas.DataFrame(columns=['Actual'])
    ix_pre = []
    ix_post = []
    for tr in treated_units:
        temp_df_pre = Y_pre_dict[tr]
        temp_df_pre.rename(columns={tr: "Actual"}, inplace=True)
        temp_df_post = Y_post_dict[tr]
        temp_df_post.rename(columns={tr: "Actual"}, inplace=True)
        ix_pre = temp_df_pre.index.union(ix_pre)
        ix_post = temp_df_post.index.union(ix_post)
        Y_pre_df = pandas.concat([Y_pre_df, temp_df_pre], axis='index')
        Y_post_df = pandas.concat([Y_post_df, temp_df_post], axis='index')

    Y_pre_df.set_index(ix_pre, inplace=True)
    Y_post_df.set_index(ix_post, inplace=True)

    return scdata_multi_output(A=A_stacked, B=B_stacked, C=C_stacked, P=P_stacked, P_diff=Pd_stacked,
                               Y_df=Y_df, Y_donors=Y_donors_stacked, Y_pre=Y_pre_df, Y_post=Y_post_df,
                               J=J_dict, K=K_dict, KM=KM_dict, M=M_dict, iota=iota, KMI=KMI,
                               cointegrated_data=cointegrated_data_dict,
                               period_pre=period_pre_dict, period_post=period_post_dict,
                               T0_features=T0_features_dict, T1_outcome=T1_dict,
                               outcome_var=outcome_var, features=features,
                               glob_cons=constant_dict, out_in_features=out_in_features_dict,
                               donors_dict=donors_dict, treated_units=treated_units,
                               effect=effect, units_est=units_est, anticipation=anticipation_dict,
                               timeConvert=timeConvert)


class scdata_multi_output:
    def __init__(self, A, B, C, P, P_diff, Y_df, Y_donors, Y_pre, Y_post, J, K, KM, M, iota, KMI,
                 cointegrated_data, period_pre, period_post, T0_features,
                 T1_outcome, outcome_var, features, glob_cons, out_in_features,
                 donors_dict, treated_units, effect, units_est, anticipation, timeConvert):

        self.A = A
        self.B = B
        self.C = C
        self.P = P
        self.P_diff = P_diff
        self.Y_df = Y_df
        self.Y_donors = Y_donors
        self.Y_pre = Y_pre
        self.Y_post = Y_post
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
        self.outcome_var = outcome_var
        self.features = features
        self.glob_cons = glob_cons
        self.out_in_features = out_in_features
        self.donors_dict = donors_dict
        self.treated_units = treated_units
        self.effect = effect
        self.units_est = units_est
        self.anticipation = anticipation
        self.timeConvert = timeConvert
