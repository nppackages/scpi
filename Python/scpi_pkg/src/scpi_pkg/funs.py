# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 09:58:31 2021

@author: Filippo Palomba
"""
# Temporary code to suppress pandas FutureWarning
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings(action='ignore', message='All-NaN slice encountered')
warnings.filterwarnings(action='ignore', message='Mean of empty slice')

import numpy
import pandas
import nlopt
import os
import cvxpy
import shutil
import statsmodels.api as sm
from math import sqrt, log, ceil, floor
from copy import deepcopy
from sklearn.preprocessing import PolynomialFeatures

def complete_cases(x):
    return numpy.all(numpy.invert(numpy.isnan(x)), axis=1)


def crossprod(x, y=None):
    if y is None:
        return x.T @ x
    else:
        return x.T @ y


# Auxiliary functions for estimation
def obj_fun_est(x, grad, Z, V, A):
    Gram = Z.T.dot(V).dot(Z)
    a = 2 * numpy.array(A).T @ numpy.array(V) @ numpy.array(Z)

    if grad.size > 0:
        grad[:] = 2 * Gram.dot(x) - a

    f = x.T.dot(Gram).dot(x) - a.dot(x)

    return f[0]


def norm_co_est(x, grad, p, J, KM, QQ, dire):
    if p == 1:
        if grad.size > 0 and dire == "==":
            grad[:] = [1] * J + [0] * KM
        elif grad.size > 0:
            aux = [1 if xi > 0 else -1 for xi in x[0:J]]
            grad[:] = numpy.append(aux, [0] * KM)
        co = numpy.sum(abs(x[0:J])) - QQ
    else:
        if grad.size > 0:
            grad[:] = numpy.append(2 * x[0:J], [0] * KM)
        co = numpy.sum(x[0:J]**2) - QQ**2

    return co


def norm_co_est_multi(result, x, grad, p, J, KM, Dtot, iota, QQ, dire):
    j_lb = 0

    for i in range(iota):
        j_ub = j_lb + J[i]
        if p == 1:
            if grad.size > 0 and dire == "==":
                grad[i, :] = [0] * j_lb + [1] * J[i] + [0] * (Dtot - J[i] - j_lb)
            elif grad.size > 0:
                grad[i, :] = [0] * j_lb + [1 if xi > 0 else -1 for xi in x[j_lb:j_ub]] + [0] * (Dtot - J[i] - j_lb)
            result[i] = numpy.sum(abs(x[j_lb:j_ub])) - QQ[i]
        else:
            if grad.size > 0:
                grad[i, :] = [0] * j_lb + (2 * x[j_lb:j_ub]).tolist() + [0] * (Dtot - J[i] - j_lb)
            result[i] = numpy.sum(x[j_lb:j_ub]**2) - QQ[i]**2
        j_lb = j_ub


# Auxiliary functions for inference

def obj_fun_min(x, grad, xt, beta):
    f = -numpy.sum(xt * (x - beta))
    if grad.size > 0:
        grad[:] = -xt
    return f


def obj_fun_max(x, grad, xt, beta):
    f = numpy.sum(xt * (x - beta))
    if grad.size > 0:
        grad[:] = xt
    return f


# Loss function constraint
def single_ineq(x, grad, beta, Q, G):
    a = -2 * G - 2 * beta.T.dot(Q)
    d = 2 * G.dot(beta) + beta.dot(Q.dot(beta))
    f = x.T.dot(Q).dot(x) + a.dot(x) + d
    if grad.size > 0:
        grad[:] = 2 * Q.dot(x) + a

    return f


def w_constr_prep(cons, names, A, Z, V, J, KM):
    # Default method to estimate weights as in Abadie et al. (2010)
    if cons is None:
        cons = {'lb': 0,
                'p': 'L1',
                'dir': '==',
                'Q': 1,
                'name': 'simplex'}

    elif cons['name'] == 'simplex':
        if 'Q' not in names:
            cons['Q'] = 1
        cons = {'lb': 0,
                'p': 'L1',
                'dir': '==',
                'Q': cons['Q'],
                'name': 'simplex'}

    elif cons['name'] == 'ols':
        cons = {'lb': -numpy.inf,
                'p': 'no norm',
                'dir': None,
                'Q': None,
                'name': 'ols'}

    elif cons['name'] == 'lasso':
        if 'Q' not in names:
            cons['Q'], lambd = shrinkage_EST("lasso", A, Z, V, J, KM)
        cons = {'lb': -numpy.inf,
                'p': 'L1',
                'dir': '<=',
                'Q': cons['Q'],
                'name': 'lasso'}

    elif cons['name'] == 'ridge':

        features = A.index.get_level_values(1).unique().tolist()
        Qfeat = []
        for feat in features:
            Af = A.iloc[A.index.get_level_values(1) == feat]
            Zf = Z.iloc[Z.index.get_level_values(1) == feat]
            Vf = V.iloc[V.index.get_level_values(1) == feat, V.index.get_level_values(1) == feat]

            if len(Af) >= 5:
                try:
                    Qest, lambd = shrinkage_EST("ridge", Af, Zf, Vf, J, KM)
                    Qfeat.append(Qest)
                except:
                    pass

        if len(Qfeat) == 0:
            Qest, lambd = shrinkage_EST("ridge", A, Z, V, J, KM)

        Qest = min(Qfeat)

        if 'Q' not in names:
            cons['Q'] = Qest

        cons = {'lb': -numpy.inf,
                'p': 'L2',
                'dir': '<=',
                'Q': cons['Q'],
                'name': 'ridge',
                'lambda': lambd}

    elif cons['name'] == 'L1/L2':
        Qest, lambd = shrinkage_EST("ridge", A, Z, V, J, KM)
        if 'Q2' not in names:
            cons['Q'] = Qest
        cons = {'lb': -numpy.inf,
                'p': 'L1/L2',
                'dir': '==/<=',
                'Q': 1,
                'Q2': cons['Q'],
                'name': 'L1/L2',
                'lambda': lambd}

    else:
        if not all(req in names for req in ['p', 'dir', 'Q', 'lb']):
            raise Exception("If 'name' is not specified, w_constr should be" +
                            " a list whose elements  must be named 'p','dir','Q','lb'.")
        if cons['p'] not in ['no norm', 'L1', 'L2']:
            raise Exception("In w_constr specify either p = 'no norm' (no constraint on the norm of weights)," +
                            " p = 'L1' (L1-norm), p = 'L2' (L2-norm).")
        if cons['lb'] != 0 and cons['lb'] != -numpy.inf:
            raise Exception("In w_constr specify either lb = 0 or lb = -numpy.inf.")
        cons['name'] = "user provided"
    return cons


def shrinkage_EST(method, A, Z, V, J, KM):
    lamb = None
    if method == "lasso":
        Q = 1

    if method == "ridge":

        wls = sm.WLS(A, Z, weights=numpy.diag(V)).fit()
        sig = wls.scale
        L2 = sum(wls.params**2)
        lamb = sig * (J + KM) / L2            # rule of thumb
        Q = sqrt(L2) / (1 + lamb)                   # convert lambda into Q

        if len(Z) <= len(Z.columns) + 10:
            lasso_cols = b_est(A=A, Z=Z, J=J, KM=KM, V=V, opt_dict=None,
                               w_constr={'dir': "<=", 'lb': -numpy.inf, 'Q': 1, 'p': "L1"})
            active_cols = abs(lasso_cols) > 1e-8
            active_cols = active_cols[:, 0].tolist()

            if active_cols.sum() >= max(len(A) - 10, 2):
                lasso_colsdf = pandas.DataFrame(abs(lasso_cols))
                active_cols = lasso_colsdf.rank(ascending=False) <= max(len(A) - 10, 2)
                active_cols = active_cols[0].tolist()

            Z_sel = Z.loc[:, active_cols]
            wls = sm.WLS(A, Z, weights=numpy.diag(V)).fit()
            sig = wls.scale
            L2 = sum(wls.params**2)
            lamb = sig * (J + KM) / L2
            Q = sqrt(L2) / (1 + lamb)

    return Q, lamb


def b_est(A, Z, J, KM, w_constr, V, opt_dict):

    lb = w_constr['lb']
    dire = w_constr['dir']
    p = w_constr['p']

    if p == "no norm":
        pp = 0
    elif p == "L1":
        pp = 1
    elif p == "L2":
        pp = 2
    elif p == "L1/L2":
        pp = None

    Zarr = numpy.array(Z)

    opt_dict = prepareOptions(opt_dict, "scest")

    use_cvxpy = useCVXPY(w_constr)

    if use_cvxpy is True:
        x = cvxpy.Variable((J + KM, 1))

        Aarr = numpy.array(A)

        objective = cvxpy.Minimize(cvxpy.quad_form(Aarr - Zarr @ x, V))
        constraints = [cvxpy.norm1(x[0:J]) <= w_constr['Q'], x[0:J] >= lb]
        prob = cvxpy.Problem(objective, constraints)
        prob.solve()

        b = x.value
        alert = prob.status != 'optimal'

    else:

        if pp == 1 and lb == -numpy.inf:
            opt = nlopt.opt(nlopt.LD_MMA, J + KM)
        else:
            opt = nlopt.opt(nlopt.LD_SLSQP, J + KM)

        opt.set_min_objective(lambda x, grad: obj_fun_est(x, grad, Zarr, V, A))

        if dire == "==":
            opt.add_equality_constraint(lambda x, grad: norm_co_est(x, grad,
                                        pp, J, KM, w_constr['Q'], dire),
                                        opt_dict['tol_eq'])

        elif dire == "<=":
            opt.add_inequality_constraint(lambda x, grad: norm_co_est(x,
                                          grad, pp, J, KM, w_constr['Q'],
                                          dire), opt_dict['tol_ineq'])

        elif dire == "==/<=":
            opt_dict['ftol_rel'] = 1e-32
            opt_dict['ftol_abs'] = 1e-32
            opt.add_equality_constraint(lambda x, grad: norm_co_est(x, grad,
                                        1, J, KM, w_constr['Q'], "=="),
                                        opt_dict['tol_eq'])
            opt.add_inequality_constraint(lambda x, grad: norm_co_est(x,
                                          grad, 2, J, KM, w_constr['Q2'],
                                          "<="), opt_dict['tol_ineq'])

        opt.set_lower_bounds([lb] * J + [-numpy.inf] * KM)
        opt.set_ftol_rel(opt_dict['ftol_rel'])
        opt.set_ftol_abs(opt_dict['ftol_abs'])
        opt.set_maxeval(opt_dict['maxeval'])
        b = opt.optimize([0] * (J + KM))

        last_res = opt.last_optimize_result()
        alert = last_res < 0 or last_res >= 5
        if alert is True:
            raise Exception("Estimation algorithm not converged! The algorithm returned the value: " +
                            str(opt.last_optimize_result()) + ". To check to what errors it corresponds" +
                            "go to 'https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#return-values'.")

    return b

def b_est_multi(A, Z, J, KM, iota, w_constr, V, opt_dict):

    lb = w_constr[0]['lb']
    dire = w_constr[0]['dir']
    p = w_constr[0]['p']
    QQ = [co['Q'] for co in w_constr]
    if w_constr[0]['name'] == "L1/L2":
        QQ2 = [co['Q2'] for co in w_constr]

    if p == "no norm":
        pp = 0
    elif p == "L1":
        pp = 1
    elif p == "L2":
        pp = 2
    elif p == "L1/L2":
        pp = None

    Zarr = numpy.array(Z)
    opt_dict = prepareOptions(opt_dict, "scest")
    use_cvxpy = useCVXPY(w_constr[0])
    Jtot = sum(J.values())
    KMI = sum(KM.values())
    Jval = [v for v in J.values()]
    KMval = [v for v in KM.values()]

    if use_cvxpy is True:
        x = cvxpy.Variable((Jtot + KMI, 1))

        Aarr = numpy.array(A)
        Zarr = numpy.array(Z)
        Varr = numpy.array(V)

        objective = cvxpy.Minimize(cvxpy.quad_form(Aarr - Zarr @ x, Varr))
        constraints = [x[0:Jtot] >= lb]

        j_lb = 0
        for i in range(iota):
            j_ub = j_lb + Jval[i] + 1
            constraints.append(cvxpy.norm1(x[j_lb:j_ub]) <= QQ[i])
            j_lb = j_ub

        prob = cvxpy.Problem(objective, constraints)
        prob.solve()

        b = x.value
        alert = prob.status != 'optimal'

    else:
        Dtot = Jtot + KMI
        if pp == 1 and lb == -numpy.inf:
            opt = nlopt.opt(nlopt.LD_MMA, Dtot)
        else:
            opt = nlopt.opt(nlopt.LD_SLSQP, Dtot)

        opt.set_min_objective(lambda x, grad: obj_fun_est(x, grad, Zarr, V, A))

        if dire == "==":
            opt.add_equality_mconstraint(lambda result, x, grad: norm_co_est_multi(result, x,
                                         grad, pp, Jval, KMval, Dtot, iota, QQ, dire),
                                         [opt_dict['tol_eq']] * iota)

        elif dire == "<=":
            opt.add_inequality_mconstraint(lambda result, x, grad: norm_co_est_multi(result, x,
                                           grad, pp, Jval, KMval, Dtot, iota, QQ, dire),
                                           [opt_dict['tol_ineq']] * iota)

        elif dire == "==/<=":
            opt.add_equality_mconstraint(lambda result, x, grad: norm_co_est_multi(result, x,
                                         grad, 1, Jval, KMval, Dtot, iota, QQ, "=="),
                                         [opt_dict['tol_eq']] * iota)

            opt.add_inequality_mconstraint(lambda result, x, grad: norm_co_est_multi(result, x,
                                           grad, 2, Jval, KMval, Dtot, iota, QQ2, "<="),
                                           [opt_dict['tol_ineq']] * iota)

        opt.set_lower_bounds([lb] * Jtot + [-numpy.inf] * KMI)
        opt.set_xtol_rel(opt_dict['xtol_rel'])
        opt.set_xtol_abs(opt_dict['xtol_abs'])
        opt.set_ftol_rel(opt_dict['ftol_rel'])
        opt.set_ftol_abs(opt_dict['ftol_abs'])
        opt.set_maxeval(opt_dict['maxeval'])

        b = opt.optimize([0] * Dtot)

        last_res = opt.last_optimize_result()
        alert = last_res < 0 or last_res >= 5
        if alert is True:
            raise Exception("Estimation algorithm not converged! The algorithm returned the value: " +
                            str(opt.last_optimize_result()) + ". To check to what errors it corresponds" +
                            "go to 'https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#return-values'.")

    return b

def V_prep(type, B, T0_features, iota):
    if type == "separate":
        V = numpy.identity(len(B))

    elif type == "pooled":
        dim_V = [sum(v.values()) for v in T0_features.values()]
        max_dim = max(dim_V)
        ones = numpy.ones([iota, 1])
        eye = numpy.identity(max_dim)
        V = numpy.kron(ones @ ones.T, eye)

        sel = numpy.full(V.shape, True)

        for i in range(iota):
            if dim_V[i] < max_dim:
                shift = i * max_dim
                sel[(shift + dim_V[i]): (shift + max_dim), :] = False
                sel[:, (shift + dim_V[i]): (shift + max_dim)] = False

        row_trim = sel.sum(axis=1) != 0
        col_trim = sel.sum(axis=0) != 0

        V = V[:, col_trim][row_trim, :] / (iota ** 2)

    V = pandas.DataFrame(V, index=B.index,
                         columns=B.index.get_level_values('treated_unit'))

    return V

def u_des_prep(B, C, u_order, u_lags, coig_data, T0_tot, M, constant, index,
               index_w, u_design, res):

    Z = pandas.concat([B, C], axis=1)

    # Construct the polynomial terms in B
    if u_order == 0:       # Simple mean
        u_des_0 = pandas.DataFrame(numpy.ones(T0_tot), index=B.index)

    elif u_order > 0:  # Include covariates when predicting u_mean

        # Create first differences feature-by-feature of the matrix B (not of C!!)
        if coig_data is True:
            B_diff = B - B.groupby('feature').shift(1)
            u_des_0 = pandas.concat([B_diff, C], axis=1).loc[:, index]

        elif coig_data is False:
            u_des_0 = Z.loc[:, index]

        if constant is False:
            colname = B.columns[0].split("_", 1)[0] + '_0_constant'
            u_des_0.insert(loc=len(u_des_0.columns), column=colname,
                           value=numpy.ones(len(u_des_0)))

    # polynomials of second order are handled later in the code

    # Construct lags of B
    if u_lags > 0:
        B_lags = pandas.DataFrame(None, index=B.index)

        if coig_data is True:
            B_diff = B - B.groupby('feature').shift(1)

        for ll in range(1, u_lags + 1):
            if coig_data is False:
                B_l = B.loc[:, index_w].groupby('feature').shift(ll)
                B_l.columns = [b + '_l' + str(ll) for b in B_l.columns.tolist()]
                B_lags = pandas.concat([B_lags, B_l], axis=1)
            else:
                B_l = B_diff.loc[:, index_w].groupby('feature').shift(ll)
                B_l.columns = [b + '_l' + str(ll) for b in B_l.columns.tolist()]
                B_lags = pandas.concat([B_lags, B_l], axis=1)

        u_des_0 = pandas.concat([u_des_0, B_lags], axis=1)

    # User chose U design matrix
    if u_design is not None:
        if not isinstance(u_design, (pandas.DataFrame, numpy.ndarray)):
            raise Exception("The object u_design should be a dataframe or a matrix!")

        if isinstance(u_design, numpy.ndarray):
            u_design = pandas.DataFrame(u_design,
                                        index=res.index,
                                        columns=res.columns)

        if len(u_design) != len(res):
            raise Exception("The object u_design has " + len(u_design) + " rows when " +
                            len(res) + " where expected!")

        u_des_0 = u_design

    return u_des_0


def e_des_prep(B, C, P, e_order, e_lags, res, sc_pred, Y_donors, out_feat, J, index,
               index_w, coig_data, T0, T1, constant, e_design, outcome_var, P_diff_pre):

    Z = pandas.concat([B, C], axis=1)

    if out_feat is False:
        e_res = sc_pred.Y_pre - sc_pred.Y_pre_fit

        if coig_data is True:
            e_des_0 = Y_donors.loc[:, index_w] - Y_donors.loc[:, index_w].shift(1)
            P_first = P.iloc[[0], :J] - Y_donors.iloc[[len(Y_donors) - 1], :].values
            P_diff = P.iloc[:, :J].diff()
            P_diff.iloc[0, :] = P_first
            e_des_1 = P_diff.loc[:, index_w]
        else:
            e_des_0 = Y_donors.loc[:, index_w]
            e_des_1 = P.loc[:, index_w]

    elif out_feat is True:
        e_res = res.loc[(outcome_var,), ]

        # Construct the polynomial terms in B (e.des.0) and P (e.des.1)
        if e_order == 0:       # Simple mean
            ix = pandas.MultiIndex.from_product([[outcome_var], B.loc[(outcome_var,), :].index.tolist()])
            ix.rename(["feature", "__time"], inplace=True)

            e_des_0 = pandas.DataFrame(numpy.ones(T0), index=ix)
            e_des_1 = pandas.DataFrame(numpy.ones(T1), index=P.index)

        elif e_order > 0:  # Include covariates when predicting u_mean
            # Create first differences feature-by-feature of the matrix B (not of C!!)
            if coig_data is True:
                B_diff = B - B.groupby('feature').shift(1)
                e_des_0 = pandas.concat([B_diff, C], axis=1).loc[:, index]

                # Remove last observation of first feature from first period of P
                P_first = P.iloc[[0], :J] - B.iloc[T0 - 1, :].values
                P_diff = P.iloc[:, :J].diff()
                P_diff.iloc[0, :] = P_first
                e_des_1 = pandas.concat([P_diff.loc[:, index_w],
                                         P.iloc[:, J:]], axis=1)

            elif coig_data is False:
                e_des_0 = Z.loc[:, index]
                e_des_1 = P.loc[:, index]

            if constant is False:
                e_des_0.insert(loc=len(e_des_0.columns), column='0_constant',
                               value=numpy.ones(len(e_des_0)))
                e_des_1.insert(loc=len(e_des_1.columns), column='0_constant',
                               value=numpy.ones(len(e_des_1)))

        nolag = False

        if P_diff_pre is not None:
            e_des_1 = P_diff_pre.loc[:, index]
            nolag = True

        # Construct lags of B
        if e_lags > 0 and nolag is False:
            B_lags = pandas.DataFrame(None, index=B.index)
            P_lags = pandas.DataFrame(None, index=P.index)

            B_diff = B - B.groupby('feature').shift(1)
            aux = B.iloc[B.index.get_level_values('feature') == outcome_var]
            if coig_data is False:
                df1 = aux.iloc[(len(aux) - e_lags):len(aux)]
                df2 = P.iloc[:, 0:J]
                P_aug = pandas.concat([df1, df2], axis=0)
            else:
                df1 = aux.iloc[(len(aux) - e_lags - 1):len(aux)]
                df2 = P.iloc[:, 0:J]
                P_aug = pandas.concat([df1, df2], axis=0)

            P_aug_diff = P_aug - P_aug.shift(1)

            for ll in range(1, e_lags + 1):
                if coig_data is False:
                    B_lags = pandas.concat([B_lags,
                                           B.loc[:, index_w].groupby('feature').shift(ll)], axis=1)
                    P_lags = pandas.concat([P_lags,
                                           P_aug.iloc[e_lags:, ].loc[:, index_w].shift(ll)], axis=1)
                else:
                    B_lags = pandas.concat([B_lags,
                                           B_diff.loc[:, index_w].groupby('feature').shift(ll)], axis=1)
                    P_lags = pandas.concat([P_lags,
                                           P_aug_diff.shift(ll).iloc[(e_lags + 1):, ].loc[:, index_w]], axis=1)

            e_des_0 = pandas.concat([e_des_0, B_lags], axis=1)
            e_des_1 = pandas.concat([e_des_1, P_lags], axis=1)

        e_des_0 = e_des_0.loc[(outcome_var,), ]

        if e_design is not None:
            if not isinstance(e_design, (pandas.DataFrame, numpy.ndarray)):
                raise Exception("The object e_design should be a dataframe or a matrix!")

            if isinstance(e_design, numpy.ndarray):
                e_design = pandas.DataFrame(e_design,
                                            index=res.index,
                                            columns=res.columns)

            if len(e_design) != len(e_res):
                raise Exception("The object e_design has " + len(e_design) + " rows when " +
                                len(e_res) + " where expected!")

            e_des_0 = e_design

    return e_res, e_des_0, e_des_1

def createPoly(u_order, e_order, index_w, u_des_0_na, e_des_0_na, e_des_1, out_feat):
    if u_order > 1:
        poly = PolynomialFeatures(u_order, include_bias=False)     # create pipeline
        act_B = len(index_w)                   # count active B columns in u_des_0
        B_poly = poly.fit_transform(u_des_0_na.iloc[:, :act_B])  # evaluate poly with interactions
        B_poly = pandas.DataFrame(B_poly, index=u_des_0_na.index)
        u_des_0_na = pandas.concat([B_poly, u_des_0_na.iloc[:, act_B:]], axis=1)

    if e_order > 1 and out_feat is True:
        poly = PolynomialFeatures(e_order, include_bias=False)     # create pipeline

        # Model e
        act_B = len(index_w)                   # count active B columns in e_des_0
        B_poly = poly.fit_transform(e_des_0_na.iloc[:, 0:act_B])  # evaluate poly with interactions
        B_poly = pandas.DataFrame(B_poly, index=e_des_0_na.index)
        e_des_0_na = pandas.concat([B_poly, e_des_0_na.iloc[:, act_B:]], axis=1)

        # Predictors for e
        act_P = len(index_w)                   # count active B columns in e_des_1
        P_poly = poly.fit_transform(e_des_1.iloc[:, 0:act_P])  # evaluate poly with interactions
        P_poly = pandas.DataFrame(P_poly, index=e_des_1.index)
        e_des_1 = pandas.concat([P_poly, e_des_1.iloc[:, act_P:]], axis=1)

    return u_des_0_na, e_des_0_na, e_des_1

def DUflexGet(u_des_0, C):

    sel_u_des_B = [u not in C.columns.tolist() for u in u_des_0.columns]
    sel_u_des_C = [b is False for b in sel_u_des_B]
    D_b = u_des_0.loc[:, sel_u_des_B]
    D_c = u_des_0.loc[:, sel_u_des_C]
    f_id = pandas.get_dummies(u_des_0.index.get_level_values('feature'))
    f_id.set_index(D_b.index, inplace=True)
    D_bb = pandas.concat([D_b, f_id], axis=1)
    features = f_id.columns.tolist()
    tomult = D_b.columns.tolist()

    D_b_int = pandas.DataFrame(index=D_b.index)
    for f in features:
        aux = D_bb.loc[:, tomult].multiply(D_bb.loc[:, f], axis="index")
        D_b_int = pandas.concat([D_b_int, aux], axis=1)

    D_b_int.set_index(D_c.index, inplace=True)
    Du = pandas.concat([D_b_int, D_c], axis=1)

    return Du

def df_EST(w_constr, w, B, J, KM):
    if (w_constr['name'] == "ols") or (w_constr['p'] == "no norm"):
        df = J

    elif (w_constr['name'] == "lasso") or ((w_constr['p'] == "L1") and (w_constr['dir'] == "<=")):
        df = sum(abs(w.iloc[:, 0].values) >= 1e-6)

    elif (w_constr['name'] == "simplex") or ((w_constr['p'] == "L1") and (w_constr['dir'] == "==")):
        df = sum(abs(w.iloc[:, 0].values) >= 1e-6) - 1

    elif (w_constr['name'] == "ridge") or (w_constr['name'] == "L1/L2") or (w_constr["p"] == "L2"):
        d = numpy.linalg.svd(B)[1]
        d = d[d > 0]
        df = sum(d**2 / (d**2 + w_constr['lambda']))

    # add degrees of freedom coming from C block
    df = df + KM

    return df

def u_sigma_est(u_mean, u_sigma, res, Z, V, index, TT, df):
    ZZ = numpy.array(Z)

    if u_sigma == 'HC0':
        vc = 1

    elif u_sigma == 'HC1':
        vc = TT / (TT - df)

    elif u_sigma == 'HC2':
        ZVZinv = numpy.linalg.pinv(ZZ.T.dot(V).dot(ZZ))
        PP = ZZ.dot(ZVZinv).dot(ZZ.T).dot(V)
        vc = 1 / (1 - numpy.diag(PP))

    elif u_sigma == 'HC3':
        ZVZinv = numpy.linalg.pinv(ZZ.T.dot(V).dot(ZZ))
        PP = ZZ.dot(ZVZinv).dot(ZZ.T).dot(V)
        vc = 1 / (1 - numpy.diag(PP))**2

    Omega = numpy.diag(numpy.array((res - u_mean)**2).flatten() * vc)
    Sigma = ZZ.T.dot(V).dot(Omega).dot(V).dot(ZZ) / (TT**2)

    return Sigma, Omega


def cond_pred(y, x, xpreds, method, tau=None):

    if method == 'lm':
        pred = xpreds.dot(sm.OLS(y, x, missing='drop').fit().params)

    elif method == 'qreg':
        qr = sm.QuantReg(endog=y, exog=x)
        pred = numpy.empty(shape=(len(xpreds), 2), dtype='object')
        try:
            pred[:, 0] = qr.fit(q=tau[0]).predict(exog=xpreds)
            pred[:, 1] = qr.fit(q=tau[1]).predict(exog=xpreds)
        except ValueError:
            raise Exception("The design matrix used to model e is not of full-rank! Consider" +
                            "reducing either e_lags or e_order!")
    return pred


def scpi_in_simul(i, x0, beta, Sigma_root, Q, P, J, Jtot, KM, KMI, iota, w_lb_est,
                  w_ub_est, p, p_int, QQ, QQ2, dire, lb, cores, opt_dict, use_cvxpy):

    zeta = numpy.random.normal(loc=0, scale=1, size=len(beta))
    G = Sigma_root.dot(zeta)
    Dtot = Jtot + KMI

    res_ub = []
    res_lb = []

    for hor in range(0, len(P)):
        pt = numpy.array(P.iloc[hor, :])

        if w_lb_est is True:

            if use_cvxpy is True:
                x = cvxpy.Variable(Jtot + KMI)
                a = -2 * G - 2 * beta.T.dot(Q)
                d = 2 * G.dot(beta) + beta.dot(Q.dot(beta))

                objective = cvxpy.Minimize(-cvxpy.sum(cvxpy.multiply(pt, x - beta)))
                constraints = [cvxpy.quad_form(x, Q) + cvxpy.sum(cvxpy.multiply(a, x)) + d <= 0]

                if lb[0] > -numpy.inf:
                    constraints = constraints.append(x[0:Jtot] >= lb)

                j_lb = 0
                for i in range(iota):
                    j_ub = j_lb + J[i]
                    constraints.append(cvxpy.norm1(x[j_lb:j_ub]) <= QQ[i])
                    j_lb = j_ub

                prob = cvxpy.Problem(objective, constraints)
                sol = prob.solve()

                alert = prob.status not in ['optimal', 'optimal_inaccurate']

                if alert is True:
                    res_lb.append(numpy.nan)
                else:
                    res_lb.append(sol)

            else:
                opt = nlopt.opt(nlopt.LD_SLSQP, Jtot + KMI)
                opt.set_min_objective(lambda x, grad: obj_fun_min(x, grad, pt, beta))

                # add loss function constraint
                opt.add_inequality_constraint(lambda x, grad: single_ineq(x, grad, beta, Q, G), opt_dict['tol_ineq'])

                if dire == "==":  # equality constraint on norm
                    opt.add_equality_mconstraint(lambda result, x, grad: norm_co_est_multi(result, x,
                                                 grad, p_int, J, KM, Dtot, iota, QQ, dire),
                                                 [opt_dict['tol_eq']] * iota)

                elif dire == "<=":  # inequality constraint on norm
                    opt.add_inequality_mconstraint(lambda result, x, grad: norm_co_est_multi(result, x,
                                                   grad, p_int, J, KM, Dtot, iota, QQ, dire),
                                                   [opt_dict['tol_ineq']] * iota)

                elif dire == "==/<=":
                    opt.add_equality_mconstraint(lambda result, x, grad: norm_co_est_multi(result, x,
                                                 grad, 1, J, KM, Dtot, iota, QQ, "=="),
                                                 [opt_dict['tol_eq']] * iota)
                    opt.add_inequality_mconstraint(lambda result, x, grad: norm_co_est_multi(result, x,
                                                   grad, 2, J, KM, Dtot, iota, QQ2, "<="),
                                                   [opt_dict['tol_ineq']] * iota)

                opt.set_lower_bounds(lb + [-numpy.inf] * KMI)
                opt.set_xtol_rel(opt_dict['xtol_rel'])
                opt.set_xtol_abs(opt_dict['xtol_abs'])
                opt.set_ftol_abs(opt_dict['ftol_abs'])
                opt.set_maxeval(opt_dict['maxeval'])

                try:
                    x = opt.optimize(x0)
                    alert = opt.last_optimize_result() < 0 or opt.last_optimize_result() >= 5
                except nlopt.RoundoffLimited:
                    try:
                        x = opt.optimize(x0)
                        alert = opt.last_optimize_result() < 0 or opt.last_optimize_result() >= 5
                    except nlopt.RoundoffLimited:
                        alert = True

                if alert is True:
                    res_lb.append(numpy.nan)
                else:
                    res_lb.append(opt.last_optimum_value())

        else:
            res_lb.append(numpy.nan)

        if w_ub_est is True:

            if use_cvxpy is True:
                x = cvxpy.Variable(Jtot + KMI)
                a = -2 * G - 2 * beta.T.dot(Q)
                d = 2 * G.dot(beta) + beta.dot(Q.dot(beta))

                objective = cvxpy.Minimize(-cvxpy.sum(cvxpy.multiply(pt, x - beta)))
                constraints = [cvxpy.quad_form(x, Q) + cvxpy.sum(cvxpy.multiply(a, x)) + d <= 0]

                if lb[0] > -numpy.inf:
                    constraints = constraints.append(x[0:Jtot] >= lb)

                j_lb = 0
                for i in range(iota):
                    j_ub = j_lb + J[i]
                    constraints.append(cvxpy.norm1(x[j_lb:j_ub]) <= QQ[i])
                    j_lb = j_ub

                prob = cvxpy.Problem(objective, constraints)
                sol = prob.solve()

                alert = prob.status not in ['optimal', 'optimal_inaccurate']

                if alert is True:
                    res_ub.append(numpy.nan)
                else:
                    res_ub.append(-sol)

            else:

                opt = nlopt.opt(nlopt.LD_SLSQP, Jtot + KMI)

                opt.set_min_objective(lambda x, grad: obj_fun_max(x, grad, pt, beta))

                # add loss function constraint
                opt.add_inequality_constraint(lambda x, grad: single_ineq(x, grad, beta, Q, G), opt_dict['tol_ineq'])

                if dire == "==":  # equality constraint on norm
                    opt.add_equality_mconstraint(lambda result, x, grad: norm_co_est_multi(result, x,
                                                 grad, p_int, J, KM, Dtot, iota, QQ, dire),
                                                 [opt_dict['tol_eq']] * iota)

                elif dire == "<=":  # inequality constraint on norm
                    opt.add_inequality_mconstraint(lambda result, x, grad: norm_co_est_multi(result, x,
                                                   grad, p_int, J, KM, Dtot, iota, QQ, dire),
                                                   [opt_dict['tol_ineq']] * iota)

                elif dire == "==/<=":
                    opt.add_equality_mconstraint(lambda result, x, grad: norm_co_est_multi(result, x,
                                                 grad, 1, J, KM, Dtot, iota, QQ, "=="),
                                                 [opt_dict['tol_eq']] * iota)
                    opt.add_inequality_mconstraint(lambda result, x, grad: norm_co_est_multi(result, x,
                                                   grad, 2, J, KM, Dtot, iota, QQ2, "<="),
                                                   [opt_dict['tol_ineq']] * iota)

                opt.set_lower_bounds(lb + [-numpy.inf] * KMI)
                opt.set_xtol_rel(opt_dict['xtol_rel'])
                opt.set_xtol_abs(opt_dict['xtol_abs'])
                opt.set_ftol_abs(opt_dict['ftol_abs'])
                opt.set_maxeval(opt_dict['maxeval'])

                try:
                    x = opt.optimize(x0)
                    alert = opt.last_optimize_result() < 0 or opt.last_optimize_result() >= 5
                except nlopt.RoundoffLimited:
                    try:
                        x = opt.optimize(x0)
                        alert = opt.last_optimize_result() < 0 or opt.last_optimize_result() >= 5
                    except nlopt.RoundoffLimited:
                        alert = True

                if alert is True:
                    res_ub.append(numpy.nan)
                else:
                    res_ub.append(-opt.last_optimum_value())

        else:
            res_ub.append(numpy.nan)

    res_lb.extend(res_ub)

    return res_lb


def scpi_in(sims, beta, Sigma_root, Q, P, J, KM, iota, w_lb_est,
            w_ub_est, p, p_int, QQ, QQ2, dire, lb, cores, opt_dict, pass_stata, verbose):

    # Progress bar
    iters = round(sims / 10)
    perc = 0

    opt_dict = prepareOptions(opt_dict, "scpi")
    use_cvxpy = useCVXPY({'Q': QQ, 'dir': dire, 'p': p})
    Jtot = sum(J.values())
    KMI = sum(KM.values())
    Jval = [v for v in J.values()]
    KMval = [v for v in KM.values()]
    Qval = [v for v in QQ.values()]
    Q2val = [v for v in QQ2.values()]

    # Algorithm initial value is lower bound unless -Inf
    x0 = numpy.array(lb)
    infx0 = [numpy.isinf(x) for x in x0]
    x0[infx0] = 0
    x0 = numpy.append(x0, [0] * KMI)

    sims_res = []

    if cores == 1:
        for i in range(sims):
            rem = (i + 1) % iters
            if rem == 0:
                perc = perc + 10
                if pass_stata is False and verbose:
                    if any('SPYDER' in name for name in os.environ) and pass_stata is False:
                        print(str(i + 1) + "/" + str(sims) +
                              " iterations completed (" + str(perc) + "%)", end="\n")
                    else:
                        print('\x1b[1A\x1b[2K')
                        print(str(i + 1) + "/" + str(sims) +
                              " iterations completed (" + str(perc) + "%)", end="\r")

            res = scpi_in_simul(i, x0, beta, Sigma_root, Q, P, Jval, Jtot, KMval, KMI, iota, w_lb_est,
                                w_ub_est, p, p_int, Qval, Q2val, dire, lb, cores, opt_dict, use_cvxpy)
            sims_res.append(res)

        vsig = numpy.array(sims_res)

    if cores > 1:
        if verbose:
            print("Starting parallelization with " + str(cores) + " cores...")
        from dask import compute, delayed, config
        from dask.distributed import Client

        config.set(scheduler='multiprocessing')
        client = Client(n_workers=cores)
        for i in range(sims):
            res = delayed(scpi_in_simul)(i, x0, beta, Sigma_root, Q, P, Jval, Jtot, KMval, KMI, iota, w_lb_est,
                                         w_ub_est, p, p_int, Qval, Q2val, dire, lb, cores, opt_dict, use_cvxpy)
            sims_res.append(res)

        vs = compute(sims_res)
        vsig = numpy.array(vs[0])

        client.close()

        shutil.rmtree("dask-worker-space")

        if verbose:
            print("")
            print("Closing working clusters...")

    return vsig


def prepareOptions(opt_dict, prompt):

    opt_names = []
    opt_values = []

    if opt_dict is not None:
        for name, value in opt_dict.items():
            opt_names.append(name)
            opt_values.append(value)
    else:
        opt_dict = {}

    if 'xtol_rel' not in opt_names:
        opt_dict['xtol_rel'] = 1e-8
    if 'xtol_abs' not in opt_names:
        opt_dict['xtol_abs'] = 1e-8
    if 'maxeval' not in opt_names:
        opt_dict['maxeval'] = 5000
    if 'tol_eq' not in opt_names:
        opt_dict['tol_eq'] = 1e-8
    if 'tol_ineq' not in opt_names:
        opt_dict['tol_ineq'] = 1e-8

    if prompt == "scest":
        if 'ftol_rel' not in opt_names:
            opt_dict['ftol_rel'] = 1e-8
        if 'ftol_abs' not in opt_names:
            opt_dict['ftol_abs'] = 1e-8
    else:
        if 'ftol_rel' not in opt_names:
            opt_dict['ftol_rel'] = 1e-4
        if 'ftol_abs' not in opt_names:
            opt_dict['ftol_abs'] = 1e-4

    return opt_dict


def scpi_out(y, x, preds, e_method, alpha, e_lb_est, e_ub_est):
    e_1 = e_2 = None
    idx = preds.index
    y = deepcopy(numpy.array(y))[:, 0]
    x = deepcopy(numpy.array(x))
    preds = deepcopy(numpy.array(preds))

    if e_lb_est is True or e_ub_est is True:

        if e_method == 'gaussian':
            x_more = numpy.vstack((preds, x))
            fit = cond_pred(y=y, x=x, xpreds=x_more, method='lm')
            e_mean = fit[:len(preds)]
            y_fit = fit[len(preds):]
            y_var = numpy.log((y - y_fit)**2)
            var_pred = cond_pred(y=y_var, x=x, xpreds=x_more, method='lm')
            e_sig2 = numpy.exp(var_pred[:len(preds)])

            q_pred = cond_pred(y=y - y_fit, x=x, xpreds=x_more, method='qreg', tau=[0.25, 0.75])
            IQ_pred = q_pred[:len(preds), 1] - q_pred[:len(preds), 0]
            IQ_pred = numpy.absolute(IQ_pred)
            e_sig = numpy.c_[numpy.sqrt(e_sig2), IQ_pred / 1.34].min(axis=1)

            eps = numpy.sqrt(-numpy.log(alpha) * 2) * e_sig

            lb = e_mean - eps
            ub = e_mean + eps

            e_1 = e_mean
            e_2 = e_sig2

            lb = lb.reshape(len(lb), 1)
            ub = ub.reshape(len(ub), 1)

        elif e_method == 'ls':
            x_more = numpy.vstack((preds, x))
            fit = cond_pred(y=y, x=x, xpreds=x_more, method='lm')
            e_mean = fit[:len(preds)]
            y_fit = fit[len(preds):]
            y_var = numpy.log((y - y_fit)**2)
            var_pred = cond_pred(y=y_var, x=x, xpreds=x_more, method='lm')
            q_pred = cond_pred(y=y - y_fit, x=x, xpreds=preds, method='qreg', tau=[0.25, 0.75])
            IQ_pred = q_pred[:len(preds), 1] - q_pred[:len(preds), 0]
            IQ_pred = numpy.absolute(IQ_pred)
            e_sig = numpy.sqrt(numpy.exp(var_pred[:len(preds)]))
            e_sig = numpy.c_[e_sig, IQ_pred / 1.34].min(axis=1)
            y_st = (y - y_fit) / numpy.sqrt(numpy.exp(var_pred[len(preds):]))

            lb = e_mean + e_sig * numpy.quantile(y_st, q=alpha)
            ub = e_mean + e_sig * numpy.quantile(y_st, q=(1 - alpha))

            lb = lb.reshape(len(lb), 1)
            ub = ub.reshape(len(ub), 1)

            e_1 = e_mean
            e_2 = e_sig**2

        elif e_method == 'qreg':
            e_pred = cond_pred(y=y, x=x, xpreds=preds,
                               method='qreg', tau=[alpha, 1 - alpha])
            lb = e_pred[:, [0]]
            ub = e_pred[:, [1]]

        lb = pandas.DataFrame(lb, index=idx)
        ub = pandas.DataFrame(ub, index=idx)

    else:
        lb = None
        ub = None

    return lb, ub, e_1, e_2


def simultaneousPredGet(vsig, T1, T1_tot, iota, u_alpha, e_alpha, e_res_na, e_des_0_na,
                        e_des_1, w_lb_est, w_ub_est, w_bounds, w_name):

    vsigLB = vsig[:, :T1_tot]
    vsigUB = vsig[:, T1_tot:]

    e_lb, e_ub, e_1, e_2 = scpi_out(y=e_res_na, x=e_des_0_na, preds=e_des_1,
                                    e_method="gaussian", alpha=e_alpha / 2,
                                    e_lb_est=True, e_ub_est=True)

    jmin = 0
    w_lb_joint = []
    w_ub_joint = []
    for i in range(iota):
        jmax = T1[i] + jmin

        if w_name in ['simplex', 'ridge', 'ols', 'L1-L2']:
            lb_joint = numpy.nanquantile(numpy.nanquantile(vsigLB[:, jmin:jmax], q=0.05, axis=0),
                                         q=(u_alpha / 2), axis=0)
            ub_joint = numpy.nanquantile(numpy.nanquantile(vsigUB[:, jmin:jmax], q=0.95, axis=0),
                                         q=(1 - u_alpha / 2), axis=0)

        else:
            lb_joint = numpy.nanquantile(vsigLB[:, jmin:jmax].min(axis=0), q=(u_alpha / 2), axis=0)
            ub_joint = numpy.nanquantile(vsigUB[:, jmin:jmax].max(axis=0), q=(1 - u_alpha / 2), axis=0)

        w_lb_joint = w_lb_joint + [lb_joint] * T1[i]
        w_ub_joint = w_ub_joint + [ub_joint] * T1[i]
        jmin = jmax

    eps = 1
    if len(e_1) > 1:
        eps = []
        for i in range(iota):
            eps = eps + [sqrt(log(T1[i] + 1))] * T1[i]

    e_lb_joint = numpy.array(e_lb) * numpy.array(eps)
    e_ub_joint = numpy.array(e_ub) * numpy.array(eps)

    if w_lb_est is False:
        w_lb_joint = w_bounds[:, 0]

    if w_ub_est is False:
        w_ub_joint = w_bounds[:, 1]

    MU = e_ub_joint + w_ub_joint
    ML = e_lb_joint + w_lb_joint

    return pandas.DataFrame(ML), pandas.DataFrame(MU)

def epskappaGet(P, rho_dict, beta, tr_units, joint=False):
    P_dict = mat2dict(P)
    beta_dict = mat2dict(beta, cols=False)

    epskappa = []
    for tr in tr_units:
        pnorm = abs(P_dict[tr]).sum(axis=1)
        epskappai = pnorm * rho_dict[tr]**2 / (2 * sqrt(numpy.sum(beta_dict[tr]**2)[0]))
        epskappa = epskappa + epskappai.tolist()

    if joint is True:
        epskappa = max(epskappa)

    return epskappa

def regularize_w(rho, rho_max, res, B, C, coig_data, T0_tot):

    if rho == 'type-1':
        ssr = (res - res.mean())**2
        sigma_u = sqrt(ssr.mean())
        sigma_bj = min(B.std(axis=0))
        CC = sigma_u / sigma_bj

    elif rho == 'type-2':
        ssr = (res - res.mean())**2
        sigma_u = sqrt(ssr.mean())
        sigma_bj2 = min(B.var(axis=0))
        sigma_bj = max(B.var(axis=0))
        CC = sigma_u * sigma_bj / sigma_bj2

    elif rho == 'type-3':
        sigma_bj2 = min(B.var(axis=0))
        tempdf = pandas.concat([res, B], axis=1)
        sigma_bju = max(tempdf.cov().iloc[:, 0])
        CC = sigma_bju / sigma_bj2

    if coig_data is True:
        c = 1
    else:
        c = 0.5

    rho = (CC * (log(T0_tot))**c) / (sqrt(T0_tot))

    if rho_max is not None:
        rho = min(rho, rho_max)

    return rho


def regularize_check(w, index_w, rho, verbose):
    # If regularization is ill-behaved just select the first weight
    if sum(index_w) == 0:
        sel = w.rank(ascending=False) <= 1
        index_w = index_w | sel[0]
        if verbose:
            warnings.warn("Regularization paramater was too high (" + str(round(rho, 3)) + "). " +
                          "We set it so that at least one component in w is non-zero.")
    return(index_w)


def local_geom(w_constr, rho, rho_max, res, B, C, coig_data, T0_tot, J, w, verbose):
    Q = w_constr['Q']
    Q2_star = None

    if isinstance(rho, str):
        rho = regularize_w(rho, rho_max, res, B, C, coig_data, T0_tot)

    if (w_constr['name'] == "simplex") | ((w_constr['p'] == "L1") & (w_constr['dir'] == "==")):
        index_w = w[0] > rho
        index_w = regularize_check(w, index_w, rho, verbose)
        w_star = deepcopy(w)
        to_regularize = [col.split("_")[1] not in index_w for col in B.columns]
        w_star.loc[to_regularize, ] = 0
        Q_star = numpy.sum(w_star)[0]

    elif (w_constr['name'] == "lasso") | ((w_constr['p'] == "L1") & (w_constr['dir'] == "<=")):
        l1 = float(numpy.sum(abs(w)))
        if (l1 >= Q - rho * sqrt(J)) & (l1 <= Q):
            Q_star = l1
        else:
            Q_star = Q

        index_w = abs(w[0]) > rho
        index_w = regularize_check(w, index_w, rho, verbose)
        w_star = deepcopy(w)
        to_regularize = [col.split("_")[1] not in index_w for col in B.columns]
        w_star.loc[to_regularize, ] = 0

    elif (w_constr['name'] == "ridge") or (w_constr['p'] == "L2"):
        l2 = float(sqrt(numpy.sum(w**2)))
        if (l2 >= Q - rho) & (l2 <= Q):
            Q_star = l2
        else:
            Q_star = Q

        w_star = deepcopy(w)
        index_w = pandas.DataFrame([True] * len(B.columns), index=w.index)[0]

    elif (w_constr['name'] == "L1/L2"):
        index_w = w[0] > rho
        index_w = regularize_check(w, index_w, rho, verbose)
        w_star = deepcopy(w)
        to_regularize = [col.split("_")[1] not in index_w for col in B.columns]
        w_star.loc[to_regularize, ] = 0
        Q_star = numpy.sum(w_star)[0]

        l2 = float(sqrt(numpy.sum(w**2)))
        if (l2 >= Q - rho) & (l2 <= Q):
            Q2_star = l2
        else:
            Q2_star = Q

        w_constr['Q2'] = Q2_star

    else:
        Q_star = Q
        w_star = w
        index_w = pandas.DataFrame([True] * len(B.columns), index=w.index)[0]

    w_constr['Q'] = Q_star

    return w_constr, w_star, index_w.values, rho, Q_star, Q2_star

def localgeom2step(w, r, rho_dict, w_constr, Q, treated_units):
    beta = pandas.concat([w, r], axis=0)
    w_dict = mat2dict(w, cols=False)
    rhoj_dict = {}
    lb = []
    Q_star = deepcopy(Q)

    for tr in treated_units:
        # Constraint on the norm of the weights
        if w_constr[tr]['p'] == "no norm":  # Unconstrained problem
            rhoj_dict[tr] = rho_dict[tr]  # auxiliary dictionary never used in this case

        elif w_constr[tr]['p'] == "L1":
            rhoj_dict[tr] = rho_dict[tr]
            w_norm = numpy.sum(abs(w_dict[tr]))[0]

        elif w_constr[tr]['p'] in ["L1/L2", "L2"]:
            L1 = numpy.sum(abs(w_dict[tr]))[0]
            rhoj_dict[tr] = 2 * L1 * rho_dict[tr]
            w_norm = numpy.sum(w_dict[tr]**2)[0]

        # Check if constraint is equality or inequality
        if w_constr[tr]['dir'] in ["<=", "==/<="]:
            active = 1 * ((w_norm - Q[tr]) > - rhoj_dict[tr])
            Q_star[tr] = active * (w_norm - Q[tr]) + Q[tr]

        # Constraint on lower bound of the weights
        if w_constr[tr]['lb'] == 0:
            active = 1 * (w_dict[tr][0] < rhoj_dict[tr])
            lb = lb + (active * w_dict[tr][0]).tolist()
        else:
            lb = lb + [-numpy.inf] * len(w_dict[tr])

    return Q_star, lb

def executionTime(T0, T1, J, iota, cores, sims, name):

    tincr = T0 / 1000

    coefsJ = numpy.array([-0.54755616, 0.09985644])

    time = numpy.array([1, J]) @ coefsJ
    time = ceil(time) * sims / 10
    time = time / cores
    time = time * tincr
    time = time * T1
    time = time * iota
    time = time / 60
    time = ceil(time) * 2

    if name == "lasso":
        time = 8 * time

    if time < 60:
        if time < 1:
            toprint = "Maximum expected execution time: less than a minute."
        elif time == 1:
            toprint = "Maximum expected execution time: " + str(time) + " minute."
        else:
            toprint = "Maximum expected execution time: " + str(time) + " minutes."
    else:
        hours = floor(time / 60)
        if hours == 1:
            toprint = "Maximum expected execution time: " + str(hours) + " hour."
        else:
            toprint = "Maximum expected execution time: " + str(hours) + "hours."

    print(toprint)

def useCVXPY(w_constr):
    flag = False
    if w_constr['p'] == "L1" and w_constr['dir'] == "<=":
        flag = True

    return flag

def mat2dict(mat, cols=True):
    X = deepcopy(mat)
    tr_units = X.index.get_level_values('treated_unit').unique().tolist()
    matdict = {}

    if len(mat.columns) == 1:
        cols = False

    if cols is True:
        for tr in tr_units:
            X_r = X.loc[pandas.IndexSlice[tr, :, :]]
            csel = [c.split("_")[0] == tr for c in X_r.columns.tolist()]
            X_rc = X_r.loc[:, numpy.array(csel)]
            matdict[tr] = X_rc
    else:
        for tr in tr_units:
            X_r = X.loc[pandas.IndexSlice[tr, :, :]]
            matdict[tr] = X_r

    return matdict

def CIrename(ci, citype):
    CI = deepcopy(ci)
    CI.columns = [c + "_" + citype for c in ci.columns]

    return CI

def detectConstant(x, tr):
    n = len(x)
    col_keep = x.sum(axis=0) != n
    col_keep = numpy.logical_and(col_keep, (x.sum(axis=0) != 0))
    x = deepcopy(x.loc[:, col_keep])
    x.insert(0, tr + "_constant", 1, allow_duplicates=True)

    return x
