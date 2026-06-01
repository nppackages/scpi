.. scpi-pkg documentation master file, created by
   sphinx-quickstart on Tue Feb 22 10:21:54 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to scpi-pkg's documentation!
====================================
The scpi package provides Python, R, and Stata implementations of estimation and inference procedures
for synthetic control methods with multiple treated units and staggered adoption using least squares, 
lasso, ridge, or simplex-type constraints. Uncertainty is quantifed using
prediction intervals as developed in `Cattaneo, Feng, and Titiunik (2021) <https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA.pdf>`_
for a single treated unit and in `Cattaneo, Feng, Palomba, and Titiunik (2027) <https://nppackages.github.io/references/Cattaneo-Feng-Palomba-Titiunik_2027_RESTAT.pdf>`_ for
multiple treated units and staggered adoption.

Companion `R <www.r-project.org>`_ and `Stata <https://www.python.org/>`_ packages are described in 
`Cattaneo, Feng, Palomba, and Titiunik (2025) <https://nppackages.github.io/references/Cattaneo-Feng-Palomba-Titiunik_2025_JSS.pdf>`_.

Related Stata, R, and Python packages useful for inference in synthetic control methods are described 
in the following website:

`https://nppackages.github.io/scpi/ <https://nppackages.github.io/scpi/>`_

For an introduction to synthetic control methods, see `Abadie (2021) <https://economics.mit.edu/files/17847>`_ and 
references therein.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   scdata
   scdataMulti
   scest
   scpi
   scplot
   scplotMulti
   Legal


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


References
==========
Abadie, A. (2021), "Using Synthetic Controls: Feasibility, Data Requirements, and Methodological
Aspects," Journal of Economic Literature, 59(2), 391-425.

Cattaneo, M. D., Feng, Y., and Titiunik, R. (2021), "Prediction Intervals for Synthetic Control
Methods," Journal of the American Statistical Association, 116(536), 1865-1880.

Cattaneo, M. D., Feng, Y., Palomba, F., and Titiunik, R. (2025), "scpi: Uncertainty Quantification for
Synthetic Control Methods," Journal of Statistical Software, 113(1), 1-38.

Cattaneo, M. D., Feng, Y., Palomba, F., and Titiunik, R. (2027), "Uncertainty Quantification in Synthetic
Controls with Staggered Treatment Adoption," Review of Economics and Statistics, forthcoming.