# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 09:58:31 2021

@author: Filippo Palomba
"""
# Temporary code to suppress pandas FutureWarning
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy
import pandas
import nlopt
import os
import cvxpy
import shutil
import statsmodels.api as sm
from math import sqrt, log, ceil
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


# Norm constraint
def norm_equal(x, grad, J, KM, QQ, p_int, dire):
    if p_int == 1:
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
        Qest, lambd = shrinkage_EST("ridge", A, Z, V, J, KM)
        if 'Q' not in names:
            cons['Q'] = Qest
        cons = {'lb': -numpy.inf,
                'p': 'L2',
                'dir': '<=',
                'Q': cons['Q'],
                'name': 'ridge',
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
        L2 = sqrt(sum(wls.params**2))
        Ahat = wls.predict(Z)
        sig = numpy.sum((A.iloc[:, 0] - Ahat)**2) / (len(A) - (J + KM))
        lamb = sig * (J + KM) / L2            # rule of thumb
        Q = L2 / (1 + lamb)                 # convert lambda into Q

    return Q, lamb


def b_est(A, Z, J, KM, w_constr, V, opt_dict):

    lb = w_constr['lb']
    dire = w_constr['dir']
    p = w_constr['p']

    if p == "no norm":
        pp = 0
    if p == "L1":
        pp = 1
    if p == "L2":
        pp = 2

    Zarr = numpy.array(Z)

    opt_dict = prepareOptions(opt_dict, "scest")

    use_cvxpy = useCVXPY(w_constr)

    if use_cvxpy is True:
        x = cvxpy.Variable((J + KM, 1))

        omega_inv = numpy.linalg.pinv(V)
        C = numpy.linalg.cholesky(omega_inv)
        C_inv = numpy.linalg.pinv(C)
        A_star = C_inv.dot(A)
        Z_star = C_inv.dot(Zarr)
        objective = cvxpy.Minimize(cvxpy.sum_squares(Z_star @ x - A_star))
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

        opt.set_lower_bounds([lb] * J + [0] * KM)
        opt.set_xtol_rel(1e-16)
        opt.set_xtol_abs(1e-16)
        opt.set_ftol_rel(1e-16)
        opt.set_ftol_abs(1e-16)
        opt.set_maxeval(opt_dict['maxeval'])
        b = opt.optimize([0] * (J + KM))

        last_res = opt.last_optimize_result()
        alert = last_res < 0 or last_res >= 5
        if alert is True:
            raise Exception("Estimation algorithm not converged! The algorithm returned the value: " +
                            str(opt.last_optimize_result()) + ". To check to what errors it corresponds" +
                            "go to 'https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#return-values'.")

    return b


def u_des_prep(B, C, Z, u_order, u_lags, coig_data, T0_tot, M, constant, index,
               index_w, u_design, res):

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
            u_des_0.insert(loc=len(u_des_0.columns), column='constant_term_u',
                           value=numpy.ones(len(u_des_0)))

    # polynomials of second order are handled later in the code

    # Construct lags of B
    if u_lags > 0:
        B_lags = pandas.DataFrame(None, index=B.index)

        if coig_data is True:
            B_diff = B - B.groupby('feature').shift(1)

        for ll in range(1, u_lags + 1):
            if coig_data is False:
                B_lags = pandas.concat([B_lags,
                                        B.loc[:, index_w].groupby('feature').shift(ll)], axis=1)
            else:
                B_lags = pandas.concat([B_lags,
                                        B_diff.loc[:, index_w].groupby('feature').shift(ll)], axis=1)

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


def e_des_prep(B, C, Z, P, e_order, e_lags, res, sc_pred, Y_donors, out_feat, J, M, index,
               index_w, coig_data, T0, T0_tot, T1, outcome_var, constant, e_design):
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
            e_des_0 = pandas.DataFrame(numpy.ones(T0_tot), index=B.index)
            e_des_1 = pandas.DataFrame(numpy.ones(T1), index=P.index)

        elif e_order > 0:  # Include covariates when predicting u_mean

            # Create first differences feature-by-feature of the matrix B (not of C!!)
            if coig_data is True:
                B_diff = B - B.groupby('feature').shift(1)
                e_des_0 = pandas.concat([B_diff, C], axis=1).loc[:, index]

                # Remove last observation of first feature from first period of P
                P_first = P.iloc[[0], :J] - B.iloc[[T0[outcome_var] - 1], :].values
                P_diff = P.iloc[:, :J].diff()
                P_diff.iloc[0, :] = P_first
                e_des_1 = pandas.concat([P_diff.loc[:, index_w], P.iloc[:, J:]], axis=1)

            elif coig_data is False:
                e_des_0 = Z.loc[:, index]
                e_des_1 = P.loc[:, index]

            if constant is False:
                e_des_0.insert(loc=len(e_des_0.columns), column='0_constant',
                               value=numpy.ones(len(e_des_0)))
                e_des_1.insert(loc=len(e_des_1.columns), column='0_constant',
                               value=numpy.ones(len(e_des_1)))

    # Construct lags of B
    if e_lags > 0:
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

def df_EST(w_constr, w, A, B, J, KM):
    if (w_constr['name'] == "ols") or (w_constr['p'] == "no norm"):
        df = J

    elif (w_constr['name'] == "lasso") or ((w_constr['p'] == "L1") and (w_constr['dir'] == "<=")):
        df = sum(w != 0)

    elif (w_constr['name'] == "simplex") or ((w_constr['p'] == "L1") and (w_constr['dir'] == "==")):
        df = sum(w != 0) - 1

    elif (w_constr['name'] == "ridge") or (w_constr["p"] == "L2"):
        d = numpy.linalg.svd(B)[1]
        d = d[d > 0]
        df = sum(d**2 / (d**2 + w_constr['lambda']))

    # add degrees of freedom coming from C block
    df = df + KM

    return df


def u_sigma_est(u_mean, u_sigma, res, Z, V, index, TT, M, df):
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

    if len(x) <= len(x.columns):
        warnings.warn("Consider specifying a less complicated model for e. The number of observations used " +
                      "to parametrically predict moments is smaller than the number of covariates used. Consider " +
                      "reducing either the number of lags (e_lags) or the order of the polynomial (e_order)!")

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


def scpi_in_simul(i, beta, Sigma_root, Q, P, J, KM, w_lb_est,
                  w_ub_est, p, p_int, QQ, dire, lb, cores, opt_dict):

    opt_dict = prepareOptions(opt_dict, "scpi")

    zeta = numpy.random.normal(loc=0, scale=1, size=len(beta))
    G = Sigma_root.dot(zeta)

    res_ub = []
    res_lb = []

    x0 = numpy.array([0] * len(beta))

    use_cvxpy = useCVXPY({'Q': QQ, 'dir': dire, 'p': p})

    for hor in range(0, len(P)):
        pt = numpy.array(P.iloc[hor, :])

        if w_lb_est is True:

            if use_cvxpy is True:
                x = cvxpy.Variable(J + KM)
                a = -2 * G - 2 * beta.T.dot(Q)
                d = 2 * G.dot(beta) + beta.dot(Q.dot(beta))

                objective = cvxpy.Minimize(-cvxpy.sum(cvxpy.multiply(pt, x - beta)))

                if lb == 0:
                    constraints = [cvxpy.norm1(x[0:J]) <= QQ,
                                   x[0:J] >= 0,
                                   cvxpy.quad_form(x, Q) + cvxpy.sum(cvxpy.multiply(a, x)) + d <= 0]
                else:
                    constraints = [cvxpy.norm1(x[0:J]) <= QQ,
                                   cvxpy.quad_form(x, Q) + cvxpy.sum(cvxpy.multiply(a, x)) + d <= 0]
                prob = cvxpy.Problem(objective, constraints)
                sol = prob.solve()

                alert = prob.status != 'optimal'

                if alert is True:
                    res_lb.append(numpy.nan)
                else:
                    res_lb.append(sol)

            else:
                opt = nlopt.opt(nlopt.LD_SLSQP, J + KM)
                opt.set_min_objective(lambda x, grad: obj_fun_min(x, grad, pt, beta))

                # add loss function constraint
                opt.add_inequality_constraint(lambda x, grad: single_ineq(x, grad, beta, Q, G), opt_dict['tol_ineq'])

                if dire == "==":  # equality constraint on norm
                    opt.add_equality_constraint(lambda x, grad: norm_equal(x, grad, J, KM, QQ, p_int, dire),
                                                opt_dict['tol_eq'])

                elif dire == "<=":  # inequality constraint on norm
                    opt.add_inequality_constraint(lambda x, grad: norm_equal(x, grad, J, KM, QQ, p_int, dire),
                                                  opt_dict['tol_ineq'])

                opt.set_lower_bounds([lb] * J + [0] * KM)
                opt.set_xtol_rel(opt_dict['xtol_rel'])
                opt.set_xtol_abs(opt_dict['xtol_abs'])
                opt.set_ftol_abs(opt_dict['ftol_abs'])
                opt.set_maxeval(opt_dict['maxeval'])
                x = opt.optimize(x0)

                alert = opt.last_optimize_result() < 0 or opt.last_optimize_result() >= 5
                flag = checkConstraints(x, dire, beta, Q, G, J, KM, QQ, p_int, 1e-2, 1e-2)

                if alert is True or flag is True:
                    res_lb.append(numpy.nan)
                else:
                    res_lb.append(opt.last_optimum_value())

        else:
            res_lb.append(numpy.nan)

        if w_ub_est is True:

            if use_cvxpy is True:
                x = cvxpy.Variable(J + KM)
                a = -2 * G - 2 * beta.T.dot(Q)
                d = 2 * G.dot(beta) + beta.dot(Q.dot(beta))

                objective = cvxpy.Minimize(sum(cvxpy.multiply(pt, x - beta)))

                if lb == 0:
                    constraints = [cvxpy.norm1(x[0:J]) <= QQ,
                                   x[0:J] >= 0,
                                   cvxpy.quad_form(x, Q) + sum(cvxpy.multiply(a, x)) + d <= 0]
                else:
                    constraints = [cvxpy.norm1(x[0:J]) <= QQ,
                                   cvxpy.quad_form(x, Q) + sum(cvxpy.multiply(a, x)) + d <= 0]
                prob = cvxpy.Problem(objective, constraints)
                sol = prob.solve()

                alert = prob.status != 'optimal'

                if alert is True:
                    res_ub.append(numpy.nan)
                else:
                    res_ub.append(-sol)

            else:

                opt = nlopt.opt(nlopt.LD_SLSQP, J + KM)

                opt.set_min_objective(lambda x, grad: obj_fun_max(x, grad, pt, beta))

                # add loss function constraint
                opt.add_inequality_constraint(lambda x, grad: single_ineq(x, grad, beta, Q, G), opt_dict['tol_ineq'])

                if dire == "==":  # equality constraint on norm
                    opt.add_equality_constraint(lambda x, grad: norm_equal(x, grad, J, KM, QQ, p_int, dire),
                                                opt_dict['tol_eq'])

                elif dire == "<=":  # inequality constraint on norm
                    opt.add_inequality_constraint(lambda x, grad: norm_equal(x, grad, J, KM, QQ, p_int, dire),
                                                  opt_dict['tol_ineq'])

                opt.set_lower_bounds([lb] * J + [0] * KM)
                opt.set_xtol_rel(opt_dict['xtol_rel'])
                opt.set_xtol_abs(opt_dict['xtol_abs'])
                opt.set_ftol_abs(opt_dict['ftol_abs'])
                opt.set_maxeval(opt_dict['maxeval'])
                x = opt.optimize(x0)

                alert = opt.last_optimize_result() < 0 or opt.last_optimize_result() >= 5
                flag = checkConstraints(x, dire, beta, Q, G, J, KM, QQ, p_int, 1e-2, 1e-2)

                if alert is True or flag is True:
                    res_ub.append(numpy.nan)
                else:
                    res_ub.append(-opt.last_optimum_value())

        else:
            res_ub.append(numpy.nan)

    res_lb.extend(res_ub)

    return res_lb


def scpi_in(sims, beta, Sigma_root, Q, P, J, KM, w_lb_est,
            w_ub_est, p, p_int, QQ, dire, lb, cores, opt_dict, pass_stata):

    # Progress bar
    iters = round(sims / 10)
    perc = 0

    sims_res = []

    if cores == 1:
        for i in range(sims):
            rem = (i + 1) % iters
            if rem == 0:
                perc = perc + 10
                if pass_stata is False:
                    if any('SPYDER' in name for name in os.environ) and pass_stata is False:
                        print(str(i + 1) + "/" + str(sims) +
                              " iterations completed (" + str(perc) + "%)", end="\n")
                    else:
                        print('\x1b[1A\x1b[2K')
                        print(str(i + 1) + "/" + str(sims) +
                              " iterations completed (" + str(perc) + "%)", end="\r")

            res = scpi_in_simul(i, beta, Sigma_root, Q, P, J, KM, w_lb_est,
                                w_ub_est, p, p_int, QQ, dire, lb, cores, opt_dict)
            sims_res.append(res)

        vsig = numpy.array(sims_res)

    if cores > 1:
        print("Starting parallelization with " + str(cores) + " cores...")
        from dask import compute, delayed, config
        from dask.distributed import Client

        config.set(scheduler='multiprocessing')
        client = Client(n_workers=cores)
        for i in range(sims):
            res = delayed(scpi_in_simul)(i, beta, Sigma_root, Q, P, J, KM, w_lb_est,
                                         w_ub_est, p, p_int, QQ, dire, lb, cores, opt_dict)
            sims_res.append(res)

        vs = compute(sims_res)
        vsig = numpy.array(vs[0])

        client.close()

        shutil.rmtree("dask-worker-space")

        print("")
        print("Closing working clusters...")

    return vsig


def checkConstraints(x, dire, beta, Q, G, J, KM, QQ, p_int, tol_eq, tol_ineq):

    grad = numpy.array([1] * len(x))

    if dire is None:
        flag = single_ineq(x, grad, beta, Q, G) > tol_ineq

    elif dire == "==":
        flag1 = single_ineq(x, grad, beta, Q, G) > tol_ineq
        flag2 = norm_equal(x, grad, J, KM, QQ, p_int, dire) > tol_eq
        flag = flag1 or flag2

    elif dire == "<=":
        flag1 = single_ineq(x, grad, beta, Q, G) > tol_ineq
        flag2 = norm_equal(x, grad, J, KM, QQ, p_int, dire) > tol_ineq
        flag = flag1 or flag2

    return flag


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

    if e_lb_est is True or e_ub_est is True:

        if e_method == 'gaussian':
            x_more = pandas.concat([preds, x], axis=0)
            fit = cond_pred(y=y, x=x, xpreds=x_more, method='lm')
            e_mean = fit[:len(preds)]
            y_fit = fit[len(preds):]
            y_var = numpy.log((y - y_fit)**2)
            var_pred = cond_pred(y=y_var, x=x, xpreds=x_more, method='lm')
            e_sig2 = numpy.exp(var_pred[:len(preds)])

            eps = numpy.sqrt(-numpy.log(alpha) * 2 * e_sig2)
            lb = e_mean - eps
            ub = e_mean + eps

            e_1 = e_mean
            e_2 = e_sig2

            lb = lb.to_numpy().reshape(len(lb), 1)
            ub = ub.to_numpy().reshape(len(ub), 1)

        elif e_method == 'ls':
            x_more = pandas.concat([preds, x], axis=0)
            fit = cond_pred(y=y, x=x, xpreds=x_more, method='lm')
            e_mean = fit[:len(preds)]
            y_fit = fit[len(preds):]
            y_var = numpy.log((y - y_fit)**2)
            var_pred = cond_pred(y=y_var, x=x, xpreds=x_more, method='lm')
            e_sig = numpy.sqrt(numpy.exp(var_pred[:len(preds)]))
            y_st = (y - y_fit) / numpy.sqrt(numpy.exp(var_pred[len(preds):]))
            lb = e_mean + e_sig * y_st.quantile(q=alpha)
            ub = e_mean + e_sig * y_st.quantile(q=(1 - alpha))

            lb = lb.to_numpy().reshape(len(lb), 1)
            ub = ub.to_numpy().reshape(len(ub), 1)

            e_1 = e_mean
            e_2 = e_sig**2

        elif e_method == 'qreg':
            e_pred = cond_pred(y=y, x=x, xpreds=preds,
                               method='qreg', tau=[alpha, 1 - alpha])
            lb = e_pred[:, [0]]
            ub = e_pred[:, [1]]

    else:
        lb = None
        ub = None

    return lb, ub, e_1, e_2


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


def regularize_check(w, index_w, rho):
    # If regularization is ill-behaved just select the first weight
    if len(index_w) == 0:
        sel = w.rank(ascending=False) <= 1
        index_w = w.loc[sel[0], ].index.get_level_values(0).tolist()
        print("Warning: regularization paramater was too high (" + str(round(rho, 3)) + "). "
              + "We set it so that at least one component in w is non-zero.")
    return(index_w)


def local_geom(w_constr, rho, rho_max, res, B, C, coig_data, T0_tot, J, w):
    Q = w_constr['Q']
    if isinstance(rho, str):
        rho = regularize_w(rho, rho_max, res, B, C, coig_data, T0_tot)

    if (w_constr['name'] == "simplex") | ((w_constr['p'] == "L1") & (w_constr['dir'] == "==")):
        index_w = w.loc[w[0] > rho, ].index.get_level_values(0).tolist()
        index_w = regularize_check(w, index_w, rho)
        w_star = deepcopy(w)
        to_regularize = [col for col in B.columns if col not in index_w]
        w_star.loc[to_regularize, ] = 0

        Q_star = float(numpy.sum(w_star))

    elif (w_constr['name'] == "lasso") | ((w_constr['p'] == "L1") & (w_constr['dir'] == "<=")):
        l1 = float(numpy.sum(abs(w)))
        if (l1 >= Q - rho * sqrt(J)) & (l1 <= Q):
            Q_star = l1
        else:
            Q_star = Q

        index_w = w.loc[abs(w[0]) > rho, ].index.get_level_values(0).tolist()
        index_w = regularize_check(w, index_w, rho)
        w_star = deepcopy(w)
        to_regularize = [col for col in B.columns if col not in index_w]
        w_star.loc[to_regularize, ] = 0

    elif (w_constr['name'] == "ridge") or (w_constr['p'] == "L2"):
        l2 = float(sqrt(numpy.sum(w**2)))
        if (l2 >= Q - rho) & (l2 <= Q):
            Q_star = l2
        else:
            Q_star = Q

        w_star = deepcopy(w)
        index_w = B.columns

    else:
        Q_star = Q
        w_star = w
        index_w = B.columns

    w_constr['Q'] = Q_star

    return w_constr, w_star, index_w, rho, Q_star

def executionTime(T0, T1, J, cores, sims, name):
    Ttot = sum(T0.values())
    tincr = Ttot / 1000

    coefsJ = numpy.array([-0.54755616, 0.09985644])

    time = numpy.array([1, J]) @ coefsJ
    time = time * sims / 10
    time = time / cores
    time = time * tincr
    time = time * T1
    time = time / 60
    time = 2 * ceil(time)

    if name == "lasso":
        time = 3 * time

    if time < 1:
        toprint = "Maximum expected execution time: less than a minute."
    elif time == 1:
        toprint = "Maximum expected execution time: " + str(time) + " minute."
    else:
        toprint = "Maximum expected execution time: " + str(time) + " minutes."

    print(toprint)

def useCVXPY(w_constr):
    flag = False
    if w_constr['p'] == "L1" and w_constr['dir'] == "<=":
        flag = True

    return flag
