Data Preparation for Synthetic Control with Staggered Adoption (:py:mod:`scdataMulti.scdataMulti`)
====================================================================================================

This page describes the function ``scdataMulti`` to prepare data for synthetic control designs in the general
case of multiple treated units, possibly with staggered adoption. The function
produces an object of class ``scdataMulti_output`` to be passed to ``scest`` or ``scpi``.

The command prepares the data to be used by scest or scpi for point estimation and inference procedures using
Synthetic Control. It allows the user to specify the outcome variable and, for each treated unit, the features 
to be matched, and covariate-adjustment feature by feature. The names of the output matrices
follow the notation proposed in `Cattaneo, Feng, Palomba, and Titiunik (2022) <https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf>`_
[UPDATE LINK].

Syntax
---------

.. currentmodule:: scpi_pkg.scdataMulti
.. autofunction:: scdataMulti


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
Some examples regarding features and cov_adj::

    # same features for all treated units
    features = {'anyname': [list of variables to match on]}

    # specify features for each treated unit
    features = {'treated_1': [list of variables to match on],
                'treated_2': [(possibly different) list of variables to match on]}

    # same covariate adjustment for all treated units for all features
    cov_adj = {'anyname': ['constant', 'trend']}

    # same covariate adjustment for all treated units but feature-specific (say, three features)
    cov_adj = {'anyname': [['constant', 'trend'], [], ['constant']]}

    # different covariate adjustment for each treated unit and feature-specific (say, two features for the first unit
    # and three for the second treated unit)
    cov_adj = {'treated_1': [['constant'], [] ],
               'treated_2': [['constant', 'trend'], ['trend'], []]}

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

This example shows how to prepare the data for synthetic control methods. The raw data can be downloaded
`here <https://raw.githubusercontent.com/nppackages/scpi/main/Python/scpi_germany.csv>`_::

    import pandas
    from scdataMulti import scdataMulti

    df = pandas.read_csv("scpi_germany.csv") 

    # Create a second placebo treated unit
    df.loc[(df['country'] == "Italy") & (df['year'] >= 1992), 'status'] = 1

    aux = scdataMulti(df=df, 
                      id_var='country',
                      treatment_var='treatment',
                      outcome_var='gdp',
                      time_var='year')  