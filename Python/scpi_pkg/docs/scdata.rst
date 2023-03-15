Data Preparation for Synthetic Control (:py:mod:`scdata.scdata`)
========================================================================

This page describes the function ``scdata`` to prepare data for synthetic control designs. The function
produces an object of class ``scdata_output`` to be passed to ``scest`` or ``scpi``.

The command prepares the data to be used by scest or scpi for point estimation and inference procedures using
Synthetic Control. It allows the user to specify the outcome variable, the features of the treated unit
to be matched, and covariate-adjustment feature by feature. The names of the output matrices
follow the notation proposed in `Cattaneo, Feng, and Titiunik (2021) <https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA.pdf>`_.

Syntax
---------

.. currentmodule:: scpi_pkg.scdata
.. autofunction:: scdata


Dependencies
------------

**Python**: `cvxpy <https://www.cvxpy.org/>`_, 

            `dask <https://docs.dask.org/en/stable/>`_, 

            `numpy <https://numpy.org/>`_, 

            `pandas <https://pandas.pydata.org/>`_, 

            `plotnine <https://plotnine.readthedocs.io/en/stable/>`_,

            `scikit-learn <https://scikit-learn.org/stable/>`_,

            `scipy <https://scipy.org/>`_,

            `statsmodels <https://www.statsmodels.org/stable/index.html>`_

**scpi_pkg**: none

Usage
-----

This example shows how to prepare the data for synthetic control methods. The raw data can be downloaded
`here <https://raw.githubusercontent.com/nppackages/scpi/main/Python/scpi_germany.csv>`_::

    import pandas
    from scpi_pkg.scdata import scdata

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
