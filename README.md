# SCPI

The `scpi` package provides Python, R and Stata implementations of estimation and inference procedures for synthetic control methods.

This work was supported by the National Science Foundation through grants [SES-1947805](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947805) and [SES-2019432](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2019432), and by the National Institutes of Health through grant [R01 GM072611-16](https://reporter.nih.gov/project-details/10093056).

## Queries and Requests

Please email: [scpi_pkg@googlegroups.com](mailto:scpi_pkg@googlegroups.com)

## Python Implementation

To install/update in Python type:
```
pip install scpi_pkg
```

- Help: [PyPI repository](https://pypi.org/project/scpi_pkg/).

- Replication: [py-script](Python/scpi_illustration.py), [plot illustration](Python/scpi_illustration_plot.py), [data](Python/scpi_germany.csv).

## R Implementation

To install/update in R from CRAN type:
```
install.packages('scpi')
````


- Help: [R Manual](R/scpi.pdf), [CRAN repository](https://cran.r-project.org/package=scpi).

- Replication: [R-script](R/scpi_illustration.R), [plot illustration](R/scpi_illustration_plot.R), [data](R/scpi_germany.csv).

## Stata Implementation

The Stata implementation relies on Python, which needs to be available in the system.

### How to install Python
There are at least two ways to install Python:
1. Download and install Python directly from [https://realpython.com/installing-python/](https://realpython.com/installing-python/).
2. Download and install Anaconda for [Windows](https://docs.anaconda.com/anaconda/install/windows/), [macOS](https://docs.anaconda.com/anaconda/install/mac-os/), or [Linux](https://docs.anaconda.com/anaconda/install/linux/).

After Python is installed, please install the `scpi` package in Python using the instructions above.

### How to link Stata and Python
Stata (16.0 or newer) and Python (>=3.8) can be linked following the [official tutorial](https://blog.stata.com/2020/08/18/stata-python-integration-part-1-setting-up-stata-to-use-python/) on the Stata blog.

### To install/update in Stata type:
```
net install scpi, from(https://raw.githubusercontent.com/nppackages/scpi/master/stata) replace
```

- Help: [scdata](stata/scdata.pdf), [scest](/stata/scest.pdf), [scpi](stata/scpi.pdf), [scplot](stata/scplot.pdf).

- Replication files: [do-file](stata/scpi_illustration.do), [plot illustration](stata/scpi_illustration_plot.do), [data](stata/scpi_germany.dta).


## References

### Software and Implementation

- Cattaneo, Feng, Palomba and Titiunik (2022): [scpi: Uncertainty Quantification for Synthetic Control Estimators](https://nppackages.github.io/references/Cattaneo-Feng-Palomba-Titiunik_2022_scpi.pdf).<br>
Working paper.

### Technical and Methodological

- Cattaneo, Feng, Palomba and Titiunik (2022): Uncertainty Quantification in Synthetic Controls with Staggered Treatment Adoption.<br>
Working paper (coming soon).

- Cattaneo, Feng and Titiunik (2021): [Prediction Intervals for Synthetic Control Methods](https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA.pdf).<br>
_Journal of the American Statistical Association_ 116(536): 1865-1880.<br>
[Supplemental](https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA--Supplement.pdf)<br>

<br><br>
