Prediction Intervals for Synthetic Control (:py:mod:`scpi.scpi`)
========================================================================

This page describes the function ``scpi`` to implement point estimation and inference procedures for synthetic control methods
using least squares, lasso, ridge, or simplex-type constraints. Uncertainty is quantified using prediction intervals according to
`Cattaneo, Feng, and Titiunik (2021) <https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA.pdf>`_.
The names of the output matrices follow the notation proposed in `Cattaneo, Feng, Palomba, and Titiunik (2022) <https://arxiv.org/abs/2210.05026>`_

Functions
---------

.. currentmodule:: scpi_pkg.scpi
.. autofunction:: scpi


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

**scpi_pkg**: scest, scplot

Usage
-----

This example shows how to prepare the data and estimate prediction intervals for synthetic control methods. The raw data can be downloaded
`here <https://raw.githubusercontent.com/nppackages/scpi/main/Python/scpi_germany.csv>`_::

    import pandas
    from scpi_pkg.scdata import scdata
    from scpi_pkg.scpi import scpi

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

    # Set options for inference
    w_constr = {'name': 'simplex', 'Q': 1}
    u_missp = True
    u_sigma = "HC1"
    u_order = 1
    u_lags = 0
    e_method = "qreg"
    e_order = 1
    e_lags = 0
    e_alpha = 0.05
    u_alpha = 0.05
    sims = 500
    cores = 1

    pi_si = scpi(data_prep, sims=sims, w_constr=w_constr, u_order=u_order, u_lags=u_lags,
                e_order=e_order, e_lags=e_lags, e_method=e_method, u_missp=u_missp,
                u_sigma=u_sigma, cores=cores, e_alpha=e_alpha, u_alpha=u_alpha)
    print(pi_si)