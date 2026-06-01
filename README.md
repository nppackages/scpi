# Synthetic Control Methods

The `scpi` package implements estimation, inference and graphical procedures
for synthetic control methods, including designs with a single treated unit,
multiple treated units and staggered treatment adoption.

- `scdata`: data preparation for synthetic control estimation with one treated unit.
- `scdataMulti`: data preparation for multiple treated units and staggered adoption.
- `scest`: synthetic control point estimation.
- `scpi`: synthetic control prediction intervals and related uncertainty summaries.
- `scplot`: plots for single-treated-unit synthetic control results.
- `scplotMulti`: plots for multiple-treated-unit and staggered-adoption results.

## Authors

- Matias D. Cattaneo, Princeton University, <matias.d.cattaneo@gmail.com>.
- Yingjie Feng, Tsinghua University, <fengyingjiepku@gmail.com>.
- Filippo Palomba, Princeton University, <filippo.palomba19@gmail.com>.
- Rocio Titiunik, Princeton University, <rocio.titiunik@gmail.com>.

## Python Implementation

To install/update in Python type:
```
pip install scpi_pkg
```

- Help: [PyPI repository](https://pypi.org/project/scpi_pkg/).

- Replication: [py-script](Python/scpi_illustration.py), [plot illustration](Python/scpi_illustration_plot.py), [data](Python/scpi_germany.csv).

- Illustration Staggered Adoption: [py-script](Python/scpi_illustration-multi.py), [plot illustration](Python/scpi_illustration_plot-multi.py).

## R Implementation

To install/update in R type:
```
install.packages('scpi')
```

- Help: [R Manual](https://cran.r-project.org/web/packages/scpi/scpi.pdf), [CRAN repository](https://cran.r-project.org/package=scpi).

- Replication: [R-script](R/scpi_illustration.R), [plot illustration](R/scpi_illustration_plot.R), [data](R/scpi_germany.csv).

- Illustration Staggered Adoption: [R-script](R/scpi_illustration-multi.R), [plot illustration](R/scpi_illustration_plot-multi.R).

## Stata Implementation

The Stata implementation relies on Python, which needs to be available in the
system. After Python is installed and linked to Stata, install the required
Python packages:

```
pip install luddite scpi_pkg
```

To install/update in Stata type:
```
net install grc1leg, from("http://www.stata.com/users/vwiggins/") replace force
net install scpi, from(https://raw.githubusercontent.com/nppackages/scpi/main/stata) replace force
```

- Help: [scdata](stata/scdata.pdf), [scdatamulti](stata/scdatamulti.pdf), [scest](stata/scest.pdf), [scpi](stata/scpi.pdf), [scplot](stata/scplot.pdf), [scplotmulti](stata/scplotmulti.pdf).

- Replication files: [do-file](stata/scpi_illustration.do), [plot illustration](stata/scpi_illustration_plot.do), [data](stata/scpi_germany.dta).

- Illustration Staggered Adoption: [do-file](stata/scpi_illustration-multi.do), [plot illustration](stata/scpi_illustration_plot-multi.do).


## References

### Software and Implementation

- Cattaneo, Feng, Palomba and Titiunik (2025): [scpi: Uncertainty Quantification for Synthetic Control Methods](https://nppackages.github.io/references/Cattaneo-Feng-Palomba-Titiunik_2025_JSS.pdf).<br>
_Journal of Statistical Software_ 113(1): 1-38.

### Technical and Methodological

- Cattaneo, Feng, Palomba and Titiunik (2027): [Uncertainty Quantification in Synthetic Controls with Staggered Treatment Adoption](https://nppackages.github.io/references/Cattaneo-Feng-Palomba-Titiunik_2027_RESTAT.pdf).<br>
_Review of Economics and Statistics_, forthcoming.<br>
[Supplemental](https://nppackages.github.io/references/Cattaneo-Feng-Palomba-Titiunik_2027_RESTAT--Supplement.pdf)

- Cattaneo, Feng and Titiunik (2021): [Prediction Intervals for Synthetic Control Methods](https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA.pdf).<br>
_Journal of the American Statistical Association_ 116(536): 1865-1880.<br>
[Supplemental](https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA--Supplement.pdf)


## Funding

This work was supported by the National Science Foundation through grants [SES-1947805](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947805), [SES-2019432](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2019432), and [SES-2241575](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2241575), and by the National Institutes of Health through grant [R01 GM072611-16](https://reporter.nih.gov/project-details/10093056).


<br><br>
