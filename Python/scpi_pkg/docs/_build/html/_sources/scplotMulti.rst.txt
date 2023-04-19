Plots for Synthetic Control with Multiple Treated Units (:py:mod:`scplotMulti.scplotMulti`)
===============================================================================================

This page describes the function ``scplotMulti`` to implement several Synthetic Control plots when Multiple 
treated units and staggered adoption are features of the design. 
The function is designed to be called after ``scest`` or ``scpi``
which implement estimation and inference procedures for Synthetic Control methods using least squares, lasso,
ridge, or simplex-type constraints according to
`Cattaneo, Feng, and Titiunik (2021) <https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA.pdf>`_
and `Cattaneo, Feng, Palomba, and Titiunik (2022) <https://arxiv.org/abs/2210.05026>`_

Syntax
---------

.. currentmodule:: scpi_pkg.scplotMulti
.. autofunction:: scplotMulti


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
    from scdataMulti import scdataMulti
    from scest import scest
    from scplotMulti import scplotMulti

    df = pandas.read_csv("scpi_germany.csv")    

    # Create a second placebo treated unit
    df.loc[(df['country'] == "Italy") & (df['year'] >= 1992), 'status'] = 1

    # same features for all treated units
    features = {'anyname': [list of variables to match on]}

    # same covariate adjustment for all treated units for all features
    cov_adj = {'anyname': ['constant', 'trend']}

    # unit specific presence of global constant
    constant = {'treated_1': True, 'treated_2': False}

    # unit specific presence of cointegration
    cointegrated_data = {'treated_1': True, 'treated_2': False}

    # unit specific presence of anticipation effects
    anticipation = {'treated_1': 0, 'treated_2': 1}

    aux = scdataMulti(df=df, 
                    id_var='country',
                    treatment_var='treatment',
                    outcome_var='gdp',
                    time_var='year',
                    features=features,
                    constant=constant, 
                    anticipation=anticipation, 
                    cointegrated_data=cointegrated_data)

    res = scest(aux)

    scplotMulti(res)
