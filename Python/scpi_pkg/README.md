# SCPI_PKG

The `scpi_pkg` package provides Python implementations of estimation and inference procedures for Synthetic Control methods.

## Authors

Matias D. Cattaneo (<cattaneo@princeton.edu>)

Yingjie Feng (<fengyj@sem.tsinghua.edu.cn>)

Filippo Palomba (<fpalomba@princeton.edu>)

Rocio Titiunik (<titiunik@princeton.edu>)

## Website

https://nppackages.github.io/scpi/

## Installation

To install/update use pip

```
pip install scpi_pkg
```

# Usage

```
from from scpi_pkg.scdata import scdata
from from scpi_pkg.scdataMulti import scdataMulti
from scpi_pkg.scest import scest
from scpi_pkg.scpi import scpi
from scpi_pkg.scplot import scplot
from scpi_pkg.scplotMulti import scplotMulti
```

- Replication: [Germany reunification example](https://github.com/nppackages/scpi/blob/main/Python/scpi_illustration.py).

## Dependencies

- cvxpy          (>= 1.1.18)
- dask            (>= 2021.04.0)
- ecos            (>= 2.0.7)
- luddite         (>= 1.0.2)
- numpy         (>= 1.20.1)
- pandas        (>= 1.5.0)
- plotnine       (>= 0.8.0)
- scikit-learn  (>= 0.24.1)
- scipy            (>= 1.7.1)
- statsmodels (>= 0.12.2)

## References

For overviews and introductions, see [nppackages website](https://nppackages.github.io/).

### Software and Implementation

- Cattaneo, Feng, Palomba, and Titiunik (2024+) [scpi: Uncertainty Quantification for Synthetic Control Estimators](https://arxiv.org/abs/2202.05984). forthcoming at *Journal of Statistical Software.*

### Technical and Methodological

- Cattaneo, Feng, and Titiunik (2021): [Prediction Intervals for Synthetic Control Methods](https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf).
  _Journal of the American Statistical Association_.
- Cattaneo, Feng, Palomba, and Titiunik (2023): [Uncertainty Quantification in Synthetic Controls with Staggered Treatment Adoption](https://arxiv.org/abs/2210.05026), working paper.
