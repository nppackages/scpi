Point estimation for Synthetic Control (:py:mod:`scest.scest`)
========================================================================

This page describes the function ``scest`` to implement point estimation for synthetic control methods
using least squares, lasso, ridge, or simplex-type constraints as developed in `Cattaneo, Feng, and Titiunik (2021) <https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf>`_.

Syntax
---------

.. currentmodule:: scpi_pkg.scest
.. autofunction:: scest


Dependencies
------------

**Python**: `cvxpy <https://www.cvxpy.org/>`_, 

            `dask <https://docs.dask.org/en/stable/>`_, 

            `nlopt <https://nlopt.readthedocs.io/en/latest/>`_,

            `numpy <https://numpy.org/>`_, 

            `pandas <https://pandas.pydata.org/>`_, 

            `plotnine <https://plotnine.readthedocs.io/en/stable/>`_,

            `scikit-learn <https://scikit-learn.org/stable/>`_,

            `scipy <https://scipy.org/>`_,

            `statsmodels <https://www.statsmodels.org/stable/index.html>`_

**scpi_pkg**: scplot

Usage
-----

This example shows how to prepare the data and conduct estimation for synthetic control methods. The raw data can be downloaded
`here <https://raw.githubusercontent.com/nppackages/scpi/main/Python/scpi_germany.csv>`_::

    import pandas
    from scpi_pkg.scdata import scdata
    from scpi_pkg.scest import scest

    data = pandas.read_csv("scpi_germany.csv")

    id_var = 'country'
    outcome_var = 'gdp'
    time_var = 'year'
    features = None
    cov_adj = None
    period_pre = numpy.arange(1960, 1991)
    period_post = numpy.arange(1991, 2004)
    unit_tr = 'West Germany'
    unit_co = list(set(data[id_var].to_list()))
    unit_co = [cou for cou in unit_co if cou != 'West Germany']
    constant = True
    report_missing = False
    cointegrated_data = True

    data_prep = scdata(df=data, id_var=id_var, time_var=time_var,
                    outcome_var=outcome_var, period_pre=period_pre,
                    period_post=period_post, unit_tr=unit_tr,
                    unit_co=unit_co, features=features, cov_adj=cov_adj,
                    cointegrated_data=cointegrated_data, constant=constant,
                    report_missing=report_missing)

    # SC - point estimation with simplex
    est_si = scest(data_prep, w_constr={'name': "simplex"})
    print(est_si)
    est_si2 = scest(data_prep, w_constr={'p': 'L1', 'dir': '==', 'Q': 1, 'lb': 0})
    print(est_si2)

    # SC - point estimation with lasso
    est_lasso = scest(data_prep, w_constr={'name': "lasso"})
    print(est_lasso)
    est_lasso2 = scest(data_prep, w_constr={'p': 'L1', 'dir': '<=', 'Q': 1, 'lb': -numpy.inf})
    print(est_lasso2)

    est_ridge = scest(data_prep, w_constr={'name': "ridge"})
    print(est_ridge)
    Q_est = est_ridge.w_constr['Q']
    est_ridge2 = scest(data_prep, w_constr={'p': 'L2', 'dir': '<=', 'Q': Q_est, 'lb': -numpy.inf})
    print(est_ridge2)

    # SC - point estimation with ols
    est_ls = scest(data_prep, w_constr={'name': "ols"})
    print(est_ls)
    est_ls2 = scest(data_prep, w_constr={'p': 'no norm', 'dir': None, 'Q': None, 'lb': -numpy.inf})
    print(est_ls2)
