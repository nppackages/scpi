Plots for Synthetic Control (:py:mod:`scplot.scplot`)
========================================================================

This page describes the function ``scplot`` to implement several Synthetic Control plots. 
The function is designed to be called after ``scest`` or ``scpi``
which implement estimation and inference procedures for Synthetic Control methods using least squares, lasso,
ridge, or simplex-type constraints according to
`Cattaneo, Feng, and Titiunik (2021) <https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf>`_ 
and `Cattaneo, Feng, Palomba, and Titiunik (2022) <https://arxiv.org/abs/2210.05026>`_

Syntax
---------

.. currentmodule:: scpi_pkg.scplot
.. autofunction:: scplot


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

**scpi_pkg**: none

Usage
-----

This example shows how to prepare the data and conduct estimation for synthetic control methods. The raw data can be downloaded
`here <https://raw.githubusercontent.com/nppackages/scpi/main/Python/scpi_germany.csv>`_::

    import pandas
    from scpi_pkg.scdata import scdata
    from scpi_pkg.scest import scest
    from scpi_pkg.scplot import scplot

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
    
    plot = scplot(est_si)
