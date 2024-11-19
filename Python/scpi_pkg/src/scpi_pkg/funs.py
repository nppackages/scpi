# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 09:58:31 2021

@author: Filippo Palomba
"""
# Temporary code to suppress pandas FutureWarning
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
from statsmodels.tools.sm_exceptions import IterationLimitWarning
warnings.filterwarnings("ignore", category=IterationLimitWarning)
warnings.filterwarnings(action='ignore', message='All-NaN slice encountered')
warnings.filterwarnings(action='ignore', message='Mean of empty slice')

import pandas
pandas.options.mode.chained_assignment = None
import numpy
import pandas
import ecos
import os
import cvxpy
import shutil
import statsmodels.api as sm
from math import sqrt, log, ceil, floor
from copy import deepcopy
from sklearn.preprocessing import PolynomialFeatures
from scipy import sparse
from scipy.linalg import eigh


def complete_cases(x):
    return numpy.all(numpy.invert(numpy.isnan(x)), axis=1)


def crossprod(x, y=None):
    if y is None:
        return x.T @ x
    else:
        return x.T @ y


def isDiagonal(mat):
    sumTot = abs(mat).sum()
    sumDiag = abs(numpy.diag(mat)).sum()
    boolean = sumTot == sumDiag  # if True, then mat is diagonal
    return boolean.item()  # item() converts boolean from numpy.bool_ to bool


def ECOS_get_n_slacks(p, dire, Jtot, iota):

    ns = 1  # simplex, ols

    if p == "L1" and dire == "<=":  # lasso
        ns = Jtot + ns

    elif p == "L2" and dire == "<=":  # ridge
        ns = iota + ns

    elif p == "L1-L2":  # L1-L2
        ns = iota + ns

    return ns


def ECOS_get_dims(Jtot, J, KMI, p, dire, iota, red):

    if p == "L1" and dire == "==":  # simplex
        dims = {'l': Jtot + 1, 'q': [Jtot + KMI + 2 - red]}

    elif p == "L1" and dire == "<=":  # lasso
        dims = {'l': iota + 2 * Jtot + 1, 'q': [Jtot + KMI + 2 - red]}

    elif p == "L2" and dire == "<=":  # ridge
        dims = {'l': iota + 1, 'q': [j + 2 for j in J] + [Jtot + KMI + 2 - red]}

    elif p == "L1-L2":  # L1-L2
        dims = {'l': Jtot + iota + 1, 'q': [j + 2 for j in J] + [Jtot + KMI + 2 - red]}

    elif p == "no norm":  # ols
        dims = {'l': 1, 'q': [Jtot + KMI + 2 - red]}

    return dims


def blockdiag(iota, Jtot, J, KMI, ns, slack=False):

    mat = numpy.zeros([iota, Jtot + KMI + ns])

    j_lb = 0
    j_ub = J[0]

    if slack is True:
        j_lb = j_lb + Jtot + KMI
        j_ub = j_ub + Jtot + KMI

    for i in range(iota):
        if i > 0:
            j_lb = j_ub + 1
            j_ub = j_lb + J[i]

        mat[i, j_lb:j_ub] = 1

    return mat


def blockdiagRidge(Jtot, J, KMI, iota):

    mat = numpy.zeros((Jtot + 2 * iota, Jtot + KMI + iota + 1))

    i_lb = 0 + 2
    i_ub = J[0] + 2
    j_lb = 0
    j_ub = J[0]

    for i in range(iota):
        if i > 0:
            j_lb = j_ub
            j_ub = j_lb + J[i]
            i_lb = i_ub + 2
            i_ub = i_lb + J[i]

        mat[(i_lb - 2):i_lb, Jtot + KMI + i] = [-1, 1]
        mat[i_lb:i_ub, j_lb:j_ub] = numpy.diag([-2] * J[i])

    return mat


def ECOS_get_A(J, Jtot, KMI, iota, p, dire, ns):

    if (p == "L1" and dire == "==") or (p == "L1-L2"):  # simplex, L1-L2
        A = blockdiag(iota, Jtot, J, KMI, ns)

    elif p == "L1" and dire == "<=":  # lasso
        A = numpy.zeros((1, Jtot + KMI + ns))

    elif p == "L2" and dire == "<=":  # ridge
        A = numpy.zeros((1, Jtot + KMI + ns))

    elif p == "no norm":  # ols
        A = numpy.zeros((1, Jtot + KMI + ns))

    return sparse.csc_matrix(A)


def ECOS_get_c(xt, ns):

    C = numpy.zeros(len(xt) + ns, dtype=numpy.float32)
    C[0:len(xt)] = xt

    return C


def ECOS_get_b(Q1, p, dire):

    if (p == "L1" and dire == "==") or (p == "L1-L2"):  # simplex, L1-L2
        b = numpy.array(Q1, dtype=numpy.float32)

    else:
        b = numpy.zeros(1, dtype=numpy.float32)

    return b


def ECOS_get_G(Jtot, KMI, J, iota, a, Q, p, dire, ns, red, scale):

    if p == "L1" and dire == "==":

        G = numpy.concatenate((numpy.reshape(numpy.append(a, scale), (1, len(a) + ns)),                                # linear part of QF
                               numpy.concatenate((-numpy.diag([1] * Jtot), numpy.zeros((Jtot, KMI + ns))), axis=1),    # lower bounds on w
                               numpy.array([[0] * (Jtot + KMI) + [-1]]),                                               # SOC definition (||sqrt(Q)beta|| <= t)
                               numpy.array([[0] * (Jtot + KMI) + [1]]),
                               -2 * numpy.concatenate((Q, numpy.zeros((Jtot + KMI - red, 1))), axis=1)),
                              axis=0)

    elif p == "L1" and dire == "<=":  # x= (beta, z, t)

        G = numpy.concatenate((numpy.reshape(numpy.concatenate((a, numpy.zeros(ns - 1), numpy.array([scale]))),      # linear part of QF
                                             (1, len(a) + ns)),
                               numpy.concatenate((numpy.diag([1] * Jtot), numpy.zeros((Jtot, KMI)),                  # lower bounds on w
                                                  numpy.diag([1] * Jtot), numpy.zeros((Jtot, 1))), axis=1),
                               numpy.concatenate((-numpy.diag([1] * Jtot), numpy.zeros((Jtot, KMI)),                 # lower bounds on w
                                                  numpy.diag([1] * Jtot), numpy.zeros((Jtot, 1))), axis=1),
                               -blockdiag(iota, Jtot, J, KMI, ns, True),                                             # norm-inequality constraint
                               numpy.array([[0] * (Jtot + KMI + Jtot) + [-1]]),                                      # SOC definition (||sqrt(Q)beta|| <= t)
                               numpy.array([[0] * (Jtot + KMI + Jtot) + [1]]),
                               -2 * numpy.concatenate((Q, numpy.zeros((Jtot + KMI - red, Jtot + 1))), axis=1)),
                              axis=0)

    elif p == "L2" and dire == "<=":  # x= (beta, s, t)

        G = numpy.concatenate((numpy.reshape(numpy.concatenate((a, numpy.zeros(iota), numpy.array([scale]))),      # linear part of QF
                                             (1, len(a) + ns)),
                               numpy.concatenate((numpy.zeros((iota, Jtot + KMI)), numpy.diag([1] * iota),         # s <= Q1^2
                                                  numpy.zeros((iota, 1))), axis=1),
                               blockdiagRidge(Jtot, J, KMI, iota),                                                 # SOC definition (||w|| <= s)
                               numpy.array([[0] * (Jtot + KMI + iota) + [-1]]),                                    # SOC definition (||sqrt(Q)beta|| <= t)
                               numpy.array([[0] * (Jtot + KMI + iota) + [1]]),
                               -2 * numpy.concatenate((Q, numpy.zeros((Jtot + KMI - red, ns))), axis=1)),
                              axis=0)

    elif p == "L1-L2":  # x= (beta, s, t)

        G = numpy.concatenate((numpy.reshape(numpy.concatenate((a, numpy.zeros(iota), numpy.array([scale]))),          # linear part of QF
                                             (1, len(a) + ns)),
                               numpy.concatenate((-numpy.diag([1] * Jtot), numpy.zeros((Jtot, KMI + ns))), axis=1),    # lower bounds on w
                               numpy.concatenate((numpy.zeros((iota, Jtot + KMI)), numpy.diag([1] * iota),             # s <= Q1^2
                                                  numpy.zeros((iota, 1))), axis=1),
                               blockdiagRidge(Jtot, J, KMI, iota),                                                     # SOC definition (||w|| <= s)
                               numpy.array([[0] * (Jtot + KMI + iota) + [-1]]),                                        # SOC definition (||sqrt(Q)beta|| <= t)
                               numpy.array([[0] * (Jtot + KMI + iota) + [1]]),
                               -2 * numpy.concatenate((Q, numpy.zeros((Jtot + KMI - red, ns))), axis=1)),
                              axis=0)

    elif p == "no norm":

        G = numpy.concatenate((numpy.reshape(numpy.append(a, scale), (1, len(a) + ns)),                                # linear part of QF
                               numpy.array([[0] * (Jtot + KMI) + [-1]]),                                               # SOC definition (||sqrt(Q)beta|| <= t)
                               numpy.array([[0] * (Jtot + KMI) + [1]]),
                               -2 * numpy.concatenate((Q, numpy.zeros((Jtot + KMI - red, 1))), axis=1)),
                              axis=0)

    return sparse.csc_matrix(numpy.real(G))


def ECOS_get_h(d, lb, J, Jtot, KMI, iota, p, dire, Q1, Q2, red):

    if p == "L1" and dire == "==":
        h = [-d] + [-ll for ll in lb] + [1] + [1] + [0] * (Jtot + KMI - red)

    elif p == "L1" and dire == "<=":
        h = [-d] + [0] * (2 * Jtot) + Q1 + [1] + [1] + [0] * (Jtot + KMI - red)

    elif p == "L2" and dire == "<=":

        aux = []
        for i in range(iota):
            aux = aux + [1, 1] + [0] * J[i]

        h = [-d] + [q**2 for q in Q1] + aux + [1] + [1] + [0] * (Jtot + KMI - red)

    elif p == "L1-L2":

        aux = []
        for i in range(iota):
            aux = aux + [1, 1] + [0] * J[i]

        h = [-d] + [-lL for lL in lb] + [q**2 for q in Q2] + aux + [1] + [1] + [0] * (Jtot + KMI - red)

    elif p == "no norm":
        h = [-d] + [1] + [1] + [0] * (Jtot + KMI - red)

    return numpy.array(h, dtype=numpy.float32)


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

        features = A.index.get_level_values(0).unique().tolist()
        Qfeat = []

        for feat in features:
            Af = A.loc[(feat,), ]
            Zf = Z.loc[(feat,), ]
            Vf = V.iloc[V.index.get_level_values(0) == feat, V.index.get_level_values(0) == feat]

            if len(Af) >= 5:
                try:
                    Qest, lambd = shrinkage_EST("ridge", Af, Zf, Vf, J, KM)
                    Qfeat.append(Qest)
                except:
                    pass

        if len(Qfeat) == 0:
            Qest, lambd = shrinkage_EST("ridge", A, Z, V, J, KM)

        Qest = max(numpy.nanmin(Qfeat), .5)

        if 'Q' not in names:
            cons['Q'] = Qest

        cons = {'lb': -numpy.inf,
                'p': 'L2',
                'dir': '<=',
                'Q': cons['Q'],
                'name': 'ridge',
                'lambda': lambd}

    elif cons['name'] == 'L1-L2':

        features = A.index.get_level_values(0).unique().tolist()
        Qfeat = []

        for feat in features:
            Af = deepcopy(A.loc[(feat,), ])
            Zf = deepcopy(Z.loc[(feat,), ])
            Vf = deepcopy(V.iloc[V.index.get_level_values(0) == feat, V.index.get_level_values(0) == feat])

            if len(Af) >= 5:
                try:
                    Qest, lambd = shrinkage_EST("ridge", Af, Zf, Vf, J, KM)
                    Qfeat.append(Qest)
                except:
                    pass

        if len(Qfeat) == 0:
            Qest, lambd = shrinkage_EST("ridge", A, Z, V, J, KM)

        Qest = max(numpy.nanmin(Qfeat), .5)

        if 'Q2' not in names:
            cons['Q'] = Qest
        cons = {'lb': 0,
                'p': 'L1-L2',
                'dir': '==/<=',
                'Q': 1,
                'Q2': cons['Q'],
                'name': 'L1-L2',
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


def shrinkage_EST(method, A, ZZ, V, J, KM):

    Z = deepcopy(ZZ)
    lamb = None
    if method == "lasso":
        Q = 1

    if method == "ridge":

        Z_columns = [c.split('_') for c in Z.columns.tolist()]
        selZ = []
        for z in Z_columns:
            if len(z) <= 2:
                selZ.append(True)
            else:
                if z[2] == "constant":
                    selZ.append(False)
                else:
                    selZ.append(True)
        Z = Z.loc[:, selZ]
        Z['constant'] = 1
        deltak = len(Z_columns) - sum(selZ) - 1
        wls = sm.WLS(A, Z, weights=numpy.diag(V)).fit()
        sig = wls.scale
        L2 = sum(wls.params**2)
        lamb = sig * (J + KM) / L2            # rule of thumb
        Q = sqrt(L2) / (1 + lamb)                   # convert lambda into Q

        if len(Z) <= len(Z.columns) + 10:
            lasso_cols = b_est(A=A, Z=Z, J=J, KM=(KM - deltak), V=V,
                               w_constr={'name': 'lasso', 'dir': "<=", 'lb': 0, 'Q': 1, 'p': "L1"})
            active_cols = abs(lasso_cols) > 1e-8
            sumac = active_cols.sum()
            active_cols = active_cols[:, 0].tolist()

            if sumac >= max(len(A) - 10, 2):
                lasso_colsdf = pandas.DataFrame(abs(lasso_cols))
                active_cols = lasso_colsdf.rank(ascending=False) <= max(len(A) - 10, 2)
                active_cols = active_cols[0].values.tolist()

            Z_sel = Z.loc[:, active_cols]
            wls = sm.WLS(A, Z_sel, weights=numpy.diag(V)).fit()
            sig = wls.scale
            L2 = sum(wls.params**2)
            lamb = sig * (J + KM) / L2
            Q = sqrt(L2) / (1 + lamb)

    return Q, lamb


def b_est(A, Z, J, KM, w_constr, V):

    lb = w_constr['lb']
    dire = w_constr['dir']
    p = w_constr['p']

    Zarr = numpy.array(Z)

    x = cvxpy.Variable((J + KM, 1))

    Aarr = numpy.array(A)

    objective = cvxpy.Minimize(cvxpy.quad_form(Aarr - Zarr @ x, V))

    if p == "no norm":
        constraints = []

    elif p == "L1":
        if dire == "==":
            constraints = [cvxpy.sum(x[0:J]) == w_constr['Q'], x[0:J] >= lb]
        elif dire == "<=":
            constraints = [cvxpy.norm1(x[0:J]) <= w_constr['Q'], x[0:J] >= lb]

    elif p == "L2":
        if dire == "==":
            constraints = [cvxpy.sum_squares(x[0:J]) == cvxpy.power(w_constr['Q'], 2)]
        elif dire == "<=":
            constraints = [cvxpy.sum_squares(x[0:J]) <= cvxpy.power(w_constr['Q'], 2)]

    elif p == "L1-L2":
        constraints = [cvxpy.sum(x[0:J]) == w_constr['Q'], x[0:J] >= lb,
                       cvxpy.sum_squares(x[0:J]) <= cvxpy.power(w_constr['Q2'], 2)]

    prob = cvxpy.Problem(objective, constraints)

    if w_constr['name'] == 'lasso' or (p == "L1" and dire == "<="):
        prob.solve(solver=cvxpy.OSQP)
    else:
        prob.solve(solver=cvxpy.ECOS)
    b = x.value
    alert = prob.status != 'optimal'

    if alert is True:
        raise Exception("Estimation algorithm not converged! The algorithm returned the value: " +
                        str(prob.status) + ". To check to what errors it corresponds" +
                        "go to 'https://www.cvxpy.org/tutorial/intro/index.html'. " +
                        "Typically, this occurs because the problem is badly-scaled." +
                        "If so, scaling the data fixes the issue.")

    return b


def b_est_multi(A, Z, J, KM, iota, w_constr, V):

    lb = w_constr[0]['lb']
    dire = w_constr[0]['dir']
    p = w_constr[0]['p']
    QQ = [co['Q'] for co in w_constr]

    if w_constr[0]['name'] == "L1-L2":
        QQ2 = [co['Q2'] for co in w_constr]

    Zarr = numpy.array(Z)
    Jtot = sum(J.values())
    KMI = sum(KM.values())
    Jval = [v for v in J.values()]

    x = cvxpy.Variable((Jtot + KMI, 1))

    Aarr = numpy.array(A)
    Zarr = numpy.array(Z)
    Varr = numpy.array(V)

    objective = cvxpy.Minimize(cvxpy.quad_form(Aarr - Zarr @ x, Varr))

    if lb != -numpy.inf:
        constraints = [x[0:Jtot] >= lb]
    else:
        constraints = []

    j_lb = 0
    for i in range(iota):
        j_ub = j_lb + Jval[i]

        if p == "L1":
            if dire == "==":
                constraints.append(cvxpy.sum(x[j_lb:j_ub]) == QQ[i])
            elif dire == "<=":
                constraints.append(cvxpy.norm1(x[j_lb:j_ub]) <= QQ[i])

        elif p == "L2":
            if dire == "==":
                constraints.append(cvxpy.sum_squares(x[j_lb:j_ub]) == cvxpy.power(QQ[i], 2))
            elif dire == "<=":
                constraints.append(cvxpy.sum_squares(x[j_lb:j_ub]) <= cvxpy.power(QQ[i], 2))

        elif p == "L1-L2":
            constraints.append(cvxpy.sum(x[j_lb:j_ub]) == QQ[i])
            constraints.append(cvxpy.sum_squares(x[j_lb:j_ub]) <= QQ2[i] ** 2)

        j_lb = j_ub

    prob = cvxpy.Problem(objective, constraints)
    if w_constr[0]['name'] == 'lasso' or (p == "L1" and dire == "<="):
        prob.solve(solver=cvxpy.OSQP)
    else:
        prob.solve(solver=cvxpy.ECOS)

    b = x.value
    alert = prob.status != 'optimal'

    if alert is True:
        raise Exception("Estimation algorithm not converged! The algorithm returned the value: " +
                        str(prob.status) + ". To check to what errors it corresponds" +
                        "go to 'https://www.cvxpy.org/tutorial/intro/index.html'.")

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
                         columns=B.index.get_level_values('ID'))

    return V


def u_des_prep(B, C, u_order, u_lags, coig_data, T0_tot, M, constant, index,
               index_w, u_design, res):

    Z = pandas.concat([B, C], axis=1)

    # Construct the polynomial terms in B
    if u_order == 0:       # Simple mean
        u_des_0 = pandas.DataFrame(numpy.ones(T0_tot), index=B.index)
        u_des_0.columns = [B.columns[0]]

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
               index_w, coig_data, T0, T1, constant, e_design, outcome_var, P_diff_pre,
               effect, iota, tr):

    P, selp = trendRemove(P)
    C, selc = trendRemove(C)
    index = numpy.array(index)[numpy.array(selp)].tolist()
    Z = pandas.concat([B, C], axis=1)

    if P_diff_pre is not None:
        P_diff_pre, auxx = trendRemove(P_diff_pre)

    if out_feat is False:
        e_res = sc_pred.Y_pre - sc_pred.Y_pre_fit

        if coig_data is True:
            e_des_0 = Y_donors.loc[:, index_w] - Y_donors.loc[:, index_w].shift(1)
            if effect == "time":
                P_first = (P.iloc[[0], :J] * iota - Y_donors.iloc[[len(Y_donors) - 1], :].values) / iota
            else:
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
            ix.rename(["feature", "Time"], inplace=True)

            aux_name = tr + "_constant"
            e_des_0 = pandas.DataFrame({aux_name: numpy.ones(T0)}, index=ix)
            if effect == "time":
                e_des_1 = pandas.DataFrame({aux_name: numpy.ones(T1) / iota}, index=P.index)
            else:
                e_des_1 = pandas.DataFrame({aux_name: numpy.ones(T1)}, index=P.index)

        elif e_order > 0:  # Include covariates when predicting u_mean
            # Create first differences feature-by-feature of the matrix B (not of C!!)
            if coig_data is True:
                B_diff = B - B.groupby('feature').shift(1)
                e_des_0 = pandas.concat([B_diff, C], axis=1).loc[:, index]

                # Remove last observation of first feature from first period of P
                if effect == "time":
                    P_first = (P.iloc[[0], :J] * iota - B.iloc[T0 - 1, :].values) / iota
                else:
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
                if effect == "time":
                    e_des_1.insert(loc=len(e_des_1.columns), column='0_constant',
                                   value=numpy.ones(len(e_des_1)) / iota)
                else:
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

    elif (w_constr['name'] == "ridge") or (w_constr['name'] == "L1-L2") or (w_constr["p"] == "L2"):
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


def scpi_in(sims, beta, Sigma_root, Q, P, J, KM, iota, w_lb_est,
            w_ub_est, p, p_int, QQ, QQ2, dire, lb, cores, pass_stata, verbose):

    # Progress bar
    iters = ceil(sims / 10)
    perc = 0

    Jtot = sum(J.values())
    KMI = sum(KM.values())
    Jval = [v for v in J.values()]
    KMval = [v for v in KM.values()]
    Qval = [v for v in QQ.values()]
    Q2val = [v for v in QQ2.values()]

    dataEcos = {}
    ns = ECOS_get_n_slacks(p, dire, Jtot, iota)

    scale, Qreg = matRegularize(Q)
    dimred = numpy.shape(Q)[0] - numpy.shape(Qreg)[0]

    dataEcos['dims'] = ECOS_get_dims(Jtot, Jval, KMI, p, dire, iota, dimred)
    dataEcos['A'] = ECOS_get_A(Jval, Jtot, KMI, iota, p, dire, ns)
    dataEcos['b'] = ECOS_get_b(Qval, p, dire)
    sims_res = []

    if cores == 1:
        for i in range(sims):
            rem = (i + 1) % iters
            perc = printIter(i, rem, perc, pass_stata, verbose, sims)

            res = scpi_in_simul(i, dataEcos, ns, beta, Sigma_root, Q, Qreg, scale, dimred, P, Jval, Jtot, KMval, KMI, iota,
                                w_lb_est, w_ub_est, p, p_int, Qval, Q2val, dire, lb)
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
            res = delayed(scpi_in_simul)(i, dataEcos, ns, beta, Sigma_root, Q, Qreg, scale, dimred, P, Jval, Jtot, KMval, KMI, iota, w_lb_est,
                                         w_ub_est, p, p_int, Qval, Q2val, dire, lb)
            sims_res.append(res)

        vs = compute(sims_res)
        vsig = numpy.array(vs[0])

        client.close()

        try:
            shutil.rmtree("dask-worker-space")
        except:
            pass

        if verbose:
            print("")
            print("Closing working clusters...")

    return vsig


def scpi_in_diag(sims, beta, Q, P, J, KM, iota, w_lb_est,
                 w_ub_est, p, p_int, QQ, QQ2, dire, lb, pass_stata,
                 verbose, tr_units, zeta, sc_effect):

    # prepare dictionaries
    Jdict = deepcopy(J)
    Qdict = mat2dict(Q)

    if sc_effect != "time":
        Pdict = mat2dict(P)
    else:
        Pdict = {}
        for tr in tr_units:
            csel = [c.split("_")[0] == tr for c in P.columns.tolist()]
            Pdict[tr] = P.loc[:, csel]

    betadict = mat2dict(beta)
    zetadict = mat2dict(zeta, cols=False)
    QQdict = deepcopy(QQ)
    QQ2dict = deepcopy(QQ2)
    lbdict = {}
    jmin = 0
    for tr in tr_units:
        jmax = jmin + Jdict[tr]
        lbdict[tr] = lb[jmin:jmax]
        jmin = jmax

    Ieff = deepcopy(iota)
    iota = 1

    j = 1
    for tr in tr_units:
        # Progress bar
        iters = ceil(sims / 10)
        perc = 0

        P = Pdict[tr]
        Q = Qdict[tr].to_numpy()
        KMI = KM[tr]
        KMval = KM[tr]
        beta = betadict[tr].to_numpy().flatten()
        Qval = [QQdict[tr]]
        Q2val = [QQ2dict[tr]]
        Jval = [Jdict[tr]]
        Jtot = deepcopy(Jdict[tr])
        zeta = zetadict[tr].to_numpy()
        lb = lbdict[tr]

        dataEcos = {}
        ns = ECOS_get_n_slacks(p, dire, Jtot, iota)

        scale, Qreg = matRegularize(Q)
        dimred = numpy.shape(Q)[0] - numpy.shape(Qreg)[0]

        dataEcos['dims'] = ECOS_get_dims(Jtot, Jval, KMI, p, dire, iota, dimred)
        dataEcos['A'] = ECOS_get_A(Jval, Jtot, KMI, iota, p, dire, ns)
        dataEcos['b'] = ECOS_get_b(Qval, p, dire)

        sims_res_lb = []
        sims_res_ub = []

        for i in range(sims):
            rem = (i + 1) % iters
            perc = printIter(i, rem, perc, pass_stata, verbose, sims, tr)

            res_lb, res_ub = scpi_in_simul_diag(i, dataEcos, ns, beta, Q, Qreg, scale, dimred,
                                                P, Jval, Jtot, KMval, KMI, iota,
                                                w_lb_est, w_ub_est, p, p_int, Qval, Q2val, dire, lb, zeta[:, i])
            sims_res_lb.append(res_lb)
            sims_res_ub.append(res_ub)

        vsig_lb_tr = numpy.array(sims_res_lb)  # the list of lists becomes a sims X horizon array
        vsig_ub_tr = numpy.array(sims_res_ub)

        if j == 1:
            if sc_effect != "time":
                vsig_lb = vsig_lb_tr
                vsig_ub = vsig_ub_tr
            else:
                vsig_lb = vsig_lb_tr / Ieff
                vsig_ub = vsig_ub_tr / Ieff
            j = 2
        else:
            if sc_effect != "time":
                vsig_lb = numpy.concatenate((vsig_lb, vsig_lb_tr), axis=1)
                vsig_ub = numpy.concatenate((vsig_ub, vsig_ub_tr), axis=1)
            else:
                vsig_lb = numpy.add(numpy.nan_to_num(vsig_lb, nan=0), numpy.nan_to_num(vsig_lb_tr, nan=0) / Ieff)
                vsig_ub = numpy.add(numpy.nan_to_num(vsig_ub, nan=0), numpy.nan_to_num(vsig_ub_tr, nan=0) / Ieff)

    vsig = numpy.concatenate((vsig_lb, vsig_ub), axis=1)

    return vsig


def scpi_in_simul(i, dataEcos, ns, beta, Sigma_root, Q, Qreg, scale, dimred, P, J, Jtot, KM, KMI, iota, w_lb_est,
                  w_ub_est, p, p_int, QQ, QQ2, dire, lb):

    zeta = numpy.random.normal(loc=0, scale=1, size=len(beta))
    G = Sigma_root.dot(zeta)

    a = -2 * G - 2 * beta.T.dot(Q)
    d = 2 * G.dot(beta) + beta.dot(Q.dot(beta))

    dataEcos['G'] = ECOS_get_G(Jtot, KMI, J, iota, a, Qreg, p, dire, ns, dimred, scale)
    dataEcos['h'] = ECOS_get_h(d, lb, J, Jtot, KMI, iota, p, dire, QQ, QQ2, dimred)

    res_ub = []
    res_lb = []

    # import pickle

    # Saving the objects:
    # with open('objs_l1l2.pkl', 'wb') as f:
    #     pickle.dump([dataEcos, ns, beta, Sigma_root, Q, Qreg, scale, dimred, P, J, Jtot, KM, KMI, iota,
    #                  p, p_int, QQ, QQ2, dire, lb], f)

    for hor in range(0, len(P)):
        pt = numpy.array(P.iloc[hor, :])

        # minimization
        if w_lb_est is True:
            dataEcos['c'] = ECOS_get_c(-pt, ns)

            solution = ecos.solve(c=dataEcos['c'],
                                  G=dataEcos['G'],
                                  h=dataEcos['h'],
                                  dims=dataEcos['dims'],
                                  A=dataEcos['A'],
                                  b=dataEcos['b'],
                                  verbose=False)

            if solution['info']['infostring'] in ['Optimal solution found', 'Close to optimal solution found']:
                xx = solution['x'][0:(Jtot + KMI)]
                sol = sum(pt * (beta - xx))
                res_lb.append(sol)
            else:
                res_lb.append(numpy.nan)
        else:
            res_lb.append(numpy.nan)

        # maximization
        if w_ub_est is True:
            dataEcos['c'] = ECOS_get_c(pt, ns)

            solution = ecos.solve(c=dataEcos['c'],
                                  G=dataEcos['G'],
                                  h=dataEcos['h'],
                                  dims=dataEcos['dims'],
                                  A=dataEcos['A'],
                                  b=dataEcos['b'],
                                  verbose=False)

            if solution['info']['infostring'] in ['Optimal solution found', 'Close to optimal solution found']:
                xx = solution['x'][0:(Jtot + KMI)]
                sol = sum(pt * (beta - xx))
                res_ub.append(sol)
            else:
                res_ub.append(numpy.nan)

        else:
            res_ub.append(numpy.nan)

    res_lb.extend(res_ub)

    return res_lb


def scpi_in_simul_diag(i, dataEcos, ns, beta, Q, Qreg, scale, dimred, P, J, Jtot,
                       KM, KMI, iota, w_lb_est, w_ub_est, p, p_int, QQ, QQ2, dire, lb, zeta):

    G = zeta
    a = -2 * G - 2 * beta.T.dot(Q)
    d = 2 * G.dot(beta) + beta.dot(Q.dot(beta))

    dataEcos['G'] = ECOS_get_G(Jtot, KMI, J, iota, a, Qreg, p, dire, ns, dimred, scale)
    dataEcos['h'] = ECOS_get_h(d, lb, J, Jtot, KMI, iota, p, dire, QQ, QQ2, dimred)

    res_ub = []
    res_lb = []

    for hor in range(0, len(P)):
        pt = numpy.array(P.iloc[hor, :])

        # minimization
        if w_lb_est is True:
            dataEcos['c'] = ECOS_get_c(-pt, ns)

            solution = ecos.solve(c=dataEcos['c'],
                                  G=dataEcos['G'],
                                  h=dataEcos['h'],
                                  dims=dataEcos['dims'],
                                  A=dataEcos['A'],
                                  b=dataEcos['b'],
                                  verbose=False)

            if solution['info']['infostring'] in ['Optimal solution found', 'Close to optimal solution found']:
                xx = solution['x'][0:(Jtot + KMI)]
                sol = sum(pt * (beta - xx))
                res_lb.append(sol)
            else:
                res_lb.append(numpy.nan)
        else:
            res_lb.append(numpy.nan)

        # maximization
        if w_ub_est is True:
            dataEcos['c'] = ECOS_get_c(pt, ns)

            solution = ecos.solve(c=dataEcos['c'],
                                  G=dataEcos['G'],
                                  h=dataEcos['h'],
                                  dims=dataEcos['dims'],
                                  A=dataEcos['A'],
                                  b=dataEcos['b'],
                                  verbose=False)

            if solution['info']['infostring'] in ['Optimal solution found', 'Close to optimal solution found']:
                xx = solution['x'][0:(Jtot + KMI)]
                sol = sum(pt * (beta - xx))
                res_ub.append(sol)
            else:
                res_ub.append(numpy.nan)
        else:
            res_ub.append(numpy.nan)

    return res_lb, res_ub


def scpi_out(y, x, preds, e_method, alpha, e_lb_est, e_ub_est, effect, out_feat):
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
            if effect == "time":
                e_mean = pandas.DataFrame(data=e_mean, index=idx)
                e_mean = e_mean.groupby(level='Time').mean().values[:, 0]

            y_fit = fit[len(preds):]
            y_var = numpy.log((y - y_fit)**2)
            var_pred = cond_pred(y=y_var, x=x, xpreds=x_more, method='lm')
            if effect == "time":
                var_pred = pandas.DataFrame(data=var_pred[:len(preds)], index=idx)
                var_pred = var_pred.groupby(level='Time').mean().values[:, 0]
            else:
                var_pred = var_pred[:len(preds)]
            e_sig2 = numpy.exp(var_pred)

            q_pred = cond_pred(y=y - y_fit, x=x, xpreds=x_more, method='qreg', tau=[0.25, 0.75])
            if effect == "time":
                q3_pred = pandas.DataFrame(data=q_pred[:len(preds), 1], index=idx)
                q1_pred = pandas.DataFrame(data=q_pred[:len(preds), 0], index=idx)
                q3_pred = q3_pred.groupby(level='Time').mean().values[:, 0]
                q1_pred = q1_pred.groupby(level='Time').mean().values[:, 0]
            else:
                q3_pred = q_pred[:len(preds), 1]
                q1_pred = q_pred[:len(preds), 0]
            IQ_pred = q3_pred - q1_pred
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
            if effect == "time":
                e_mean = pandas.DataFrame(data=e_mean, index=idx)
                e_mean = e_mean.groupby(level='Time').mean().values[:, 0]

            y_fit = fit[len(preds):]
            y_var = numpy.log((y - y_fit)**2)
            var_pred = cond_pred(y=y_var, x=x, xpreds=x_more, method='lm')
            res_var = var_pred[len(preds):]
            if effect == "time":
                var_pred = pandas.DataFrame(data=var_pred[:len(preds)], index=idx)
                var_pred = var_pred.groupby(level='Time').mean().values[:, 0]
            else:
                var_pred = var_pred[:len(preds)]

            q_pred = cond_pred(y=y - y_fit, x=x, xpreds=x_more, method='qreg', tau=[0.25, 0.75])
            if effect == "time":
                q3_pred = pandas.DataFrame(data=q_pred[:len(preds), 1], index=idx)
                q1_pred = pandas.DataFrame(data=q_pred[:len(preds), 0], index=idx)
                q3_pred = q3_pred.groupby(level='Time').mean().values[:, 0]
                q1_pred = q1_pred.groupby(level='Time').mean().values[:, 0]
            else:
                q3_pred = q_pred[:len(preds), 1]
                q1_pred = q_pred[:len(preds), 0]
            IQ_pred = q3_pred - q1_pred
            IQ_pred = numpy.absolute(IQ_pred)
            e_sig = numpy.sqrt(numpy.exp(var_pred[:len(preds)]))
            e_sig = numpy.c_[e_sig, IQ_pred / 1.34].min(axis=1)
            y_st = (y - y_fit) / numpy.sqrt(numpy.exp(res_var))

            lb = e_mean + e_sig * numpy.quantile(y_st, q=alpha)
            ub = e_mean + e_sig * numpy.quantile(y_st, q=(1 - alpha))

            lb = lb.reshape(len(lb), 1)
            ub = ub.reshape(len(ub), 1)

            e_1 = e_mean
            e_2 = e_sig**2

        elif e_method == 'qreg':
            e_pred = cond_pred(y=y, x=x, xpreds=preds,
                               method='qreg', tau=[alpha, 1 - alpha])
            if effect == "time":
                e_pred_lb = pandas.DataFrame(data=e_pred[:, 0], index=idx)
                e_pred_ub = pandas.DataFrame(data=e_pred[:, 1], index=idx)
                e_pred_lb = e_pred_lb.groupby(level='Time').mean().values[:, 0]
                e_pred_ub = e_pred_ub.groupby(level='Time').mean().values[:, 0]
            else:
                e_pred_lb = e_pred[:, [0]]
                e_pred_ub = e_pred[:, [1]]

            lb = e_pred_lb
            ub = e_pred_ub

        if effect == "time":
            idx = idx.unique('Time')

        lb = pandas.DataFrame(lb, index=idx)
        ub = pandas.DataFrame(ub, index=idx)

        # if model is heavily misspecified just give up on out-of-sample uncertainty
        if out_feat is False:
            if any(lb[0] > 0):
                lb = pandas.DataFrame([lb.min()] * len(lb), index=idx)
            if any(ub[0] < 0):
                ub = pandas.DataFrame([ub.max()] * len(ub), index=idx)

    else:
        lb = None
        ub = None

    return lb, ub, e_1, e_2


def simultaneousPredGet(vsig, T1, T1_tot, iota, u_alpha, e_alpha, e_res_na, e_des_0_na,
                        e_des_1, w_lb_est, w_ub_est, w_bounds, w_name, effect, out_feat):

    vsigLB = vsig[:, :T1_tot]
    vsigUB = vsig[:, T1_tot:]

    e_lb, e_ub, e_1, e_2 = scpi_out(y=e_res_na, x=e_des_0_na, preds=e_des_1,
                                    e_method="gaussian", alpha=e_alpha / 2,
                                    e_lb_est=True, e_ub_est=True, effect=effect, out_feat=out_feat)

    jmin = 0
    w_lb_joint = []
    w_ub_joint = []
    for i in range(iota):
        jmax = T1[i] + jmin

        lb_joint = numpy.nanquantile(numpy.nanmin(vsigLB[:, jmin:jmax], axis=0), q=(u_alpha / 2), axis=0)
        ub_joint = numpy.nanquantile(numpy.nanmax(vsigUB[:, jmin:jmax], axis=0), q=(1 - u_alpha / 2), axis=0)

        w_lb_joint = w_lb_joint + [lb_joint] * T1[i]
        w_ub_joint = w_ub_joint + [ub_joint] * T1[i]
        jmin = jmax

    eps = 1
    if len(e_1) > 1:
        eps = []
        for i in range(iota):
            eps = eps + [sqrt(log(T1[i] + 1))] * T1[i]

    e_lb_joint = numpy.multiply(numpy.array(e_lb[0].values), numpy.array(eps))
    e_ub_joint = numpy.multiply(numpy.array(e_ub[0].values), numpy.array(eps))

    if w_lb_est is False:
        w_lb_joint = w_bounds[:, 0]

    if w_ub_est is False:
        w_ub_joint = w_bounds[:, 1]

    MU = e_ub_joint + w_ub_joint
    ML = e_lb_joint + w_lb_joint

    return pandas.DataFrame(ML), pandas.DataFrame(MU)


def epskappaGet(P, rho_dict, beta, tr_units, effect, joint=False):

    if effect == "time":
        rho_avg = sum([r for r in rho_dict.values()]) / len(rho_dict)
        pnorm = abs(P).sum(axis=1)
        epskappai = pnorm * rho_avg**2 / (2 * sqrt(numpy.sum(beta**2)[0]))
        epskappa = epskappai.tolist()

        if joint is True:
            epskappa = max(epskappa)
    else:
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
        denomCheck(sigma_bj)
        CC = sigma_u / sigma_bj

    elif rho == 'type-2':
        ssr = (res - res.mean())**2
        sigma_u = sqrt(ssr.mean())
        sigma_bj2 = min(B.var(axis=0))
        sigma_bj = max(B.std(axis=0))
        denomCheck(sigma_bj2)
        CC = sigma_u * sigma_bj / sigma_bj2

    elif rho == 'type-3':
        sigma_bj2 = min(B.var(axis=0))
        denomCheck(sigma_bj2)
        tempdf = pandas.concat([res, B], axis=1)
        sigma_bju = max(tempdf.cov().iloc[1:, 0])
        CC = abs(sigma_bju) / sigma_bj2

    if coig_data is True:
        c = 1
    else:
        c = 0.5

    rho = (CC * (log(T0_tot))**c) / (sqrt(T0_tot))

    if rho_max is not None:
        rho = min(rho, rho_max)

    return rho


def denomCheck(den):
    if abs(den) == 0:
        raise Exception("One of your donors has no variation in the pre-treatment period!")


def regularize_check(w, index_w, rho, verbose, B):
    # If regularization is ill-behaved just select the first weight
    if sum(index_w) == 0:
        sel = w.rank(ascending=False) <= 1
        index_w = index_w | sel[0]
        if verbose:
            tr_id = B.columns.tolist()[0].split('_')[0]
            warnings.warn("Regularization paramater was too high (" + str(round(rho, 3)) + ") " +
                          "for the treated unit with id " + str(tr_id) + ". " +
                          "We set it so that at least one component in w is non-zero.")
    return index_w


def regularize_check_lb(w, rho, rho_max, res, B, C, coig_data, T0_tot, verbose):
    rho_old = rho

    if rho < 0.001:
        rho = max(regularize_w("type-1", 0.2, res, B, C, coig_data, T0_tot),
                  regularize_w("type-2", 0.2, res, B, C, coig_data, T0_tot),
                  regularize_w("type-3", 0.2, res, B, C, coig_data, T0_tot))

        # strong evidence in favor of collinearity, thus heavy shrinkage
        if rho < 0.05:
            rho = rho_max

            if verbose is True:
                tr_id = B.columns.tolist()[0].split('_')[0]
                warnings.warn("Regularization paramater was too low (" + str(round(rho_old, 5)) + ") " +
                              "for the treated unit with id " + str(tr_id) + ". "
                              "We increased it to " + str(round(rho, 3)) + " to favor shrinkage " +
                              "and avoid overfitting issues when computing out-of-sample uncertainty!" +
                              "Please check that there is no collinearity among the donors you used for " +
                              "this treated unit. To do so run a linear regression of the features (matrix A) onto B and C.")

    return rho


def local_geom(w_constr, rho, rho_max, res, B, C, coig_data, T0_tot, J, w, verbose):
    Q = w_constr['Q']
    Q2_star = None

    # if rho is not given by the user we use our own routine to compute it
    if isinstance(rho, str):
        rho = regularize_w(rho, rho_max, res, B, C, coig_data, T0_tot)
        # here we check that our rho is not too low to prevent overfitting later on
        rho = regularize_check_lb(w, rho, rho_max, res, B, C, coig_data, T0_tot, verbose)

    if (w_constr['name'] == "simplex") | ((w_constr['p'] == "L1") & (w_constr['dir'] == "==")):
        index_w = w[0] > rho
        index_w = regularize_check(w, index_w, rho, verbose, B)
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
        index_w = regularize_check(w, index_w, rho, verbose, B)
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

    elif (w_constr['name'] == "L1-L2"):
        index_w = w[0] > rho
        index_w = regularize_check(w, index_w, rho, verbose, B)
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


def localgeom2step(w, r, rho_dict, rho_max, w_constr, Q, treated_units):

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

        elif w_constr[tr]['p'] in ["L1-L2", "L2"]:
            L1 = numpy.sum(abs(w_dict[tr]))[0]
            rhoj_dict[tr] = min(2 * L1 * rho_dict[tr], rho_max)
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


def mat2dict(mat, cols=True):
    X = deepcopy(mat)
    tr_units = X.index.get_level_values('ID').unique().tolist()
    matdict = {}

    if len(mat.columns) == 1:
        cols = False

    if cols is True:
        for tr in tr_units:
            X_r = X.loc[pandas.IndexSlice[tr, :, :]]
            csel = [str(c).split("_")[0] == tr for c in X_r.columns.tolist()]
            X_rc = X_r.loc[:, numpy.array(csel)]
            # to ensure backward compatibility with Python 3.7 (IndexSlice did not remove the column automatically)
            if 'ID' in X_rc.index.names:
                X_rc = X_rc.reset_index(level=0, drop=True)
            matdict[tr] = X_rc
    else:
        for tr in tr_units:
            X_r = X.loc[pandas.IndexSlice[tr, :, :]]
            if 'ID' in X_r.index.names:  # to ensure backward compatibility with Python 3.7
                X_r = X_r.reset_index(level=0, drop=True)
            matdict[tr] = X_r

    return matdict


def CIrename(ci, citype):
    CI = deepcopy(ci)
    CI.columns = [c + "_" + citype for c in ci.columns]

    return CI


def detectConstant(x, tr, scale_x=1):
    x = x * scale_x
    n = len(x.loc[complete_cases(x), ])
    col_keep = x.sum(axis=0) != n
    col_keep = numpy.logical_and(col_keep, (x.sum(axis=0) != 0))
    x = deepcopy(x.loc[:, col_keep])
    x.insert(0, tr + "_constant", 1, allow_duplicates=True)
    x = x / scale_x
    return x


def trendRemove(xxx):
    x = deepcopy(xxx)
    sel = []
    for c in x.columns.tolist():
        cp = c.split('_')
        if len(cp) < 3:
            sel.append(True)
        else:
            if cp[2] == "trend":
                sel.append(False)
            else:
                sel.append(True)
    xx = x.loc[:, numpy.array(sel)]

    return xx, sel


def printIter(i, rem, perc, pass_stata, verbose, sims, tr=0):
    if rem == 0:
        perc = perc + 10
        if pass_stata is False and verbose:
            if any('SPYDER' in name for name in os.environ) and pass_stata is False:
                if tr == 0:
                    print(str(i + 1) + "/" + str(sims) +
                          " iterations completed (" + str(perc) + "%)", end="\n")
                else:
                    print("Treated unit " + str(tr) + ": " + str(i + 1) + "/" + str(sims) +
                          " iterations completed (" + str(perc) + "%)", end="\n")
            else:
                print('\x1b[1A\x1b[2K')
                print(str(i + 1) + "/" + str(sims) +
                      " iterations completed (" + str(perc) + "%)", end="\r")

    return perc


def matRegularize(P, cond=None):

    w, V = eigh(P)

    if cond is None:
        cond = 1e6 * numpy.finfo(float).eps

    scale = max(abs(w))

    if scale < cond:
        scale = 0
        Qreg = None

        return scale, Qreg

    w_scaled = w / scale
    maskp = w_scaled > cond
    maskn = w_scaled < -cond

    if maskp.any() and maskn.any():
        raise Exception("Forming a non-convex expression QuadForm(x, indefinite)")

    if sum(maskp) <= 1:
        M1 = V[:, maskp] * numpy.sqrt(w_scaled[maskp])
    else:
        M1 = numpy.dot(V[:, maskp], numpy.diag(numpy.sqrt(w_scaled[maskp])))

    if sum(maskn) <= 1:
        M2 = V[:, maskn] * numpy.sqrt(-w_scaled[maskn])
    else:
        M2 = numpy.dot(V[:, maskn], numpy.diag(numpy.sqrt(-w_scaled[maskn])))

    if numpy.prod(M1.shape) > 0:
        Qreg = M1.transpose()
    elif numpy.prod(M2.shape) > 0:
        scale = -scale
        Qreg = M2.transpose()

    return scale, Qreg


def ix2rn(s):
    return str(s).replace('(', '').replace(')', '').replace("'", '')
