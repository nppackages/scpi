# Setup SCPI_PKG on Stata/Python
Python users would require just the first two steps described below, whereas Stata users would need to complete all the steps. Pypi page [here](https://test.pypi.org/project/scpi-pkg/).

## How to Install Python
There are at least two ways of installing Python:
1. Download and install Python directly ([see here](https://realpython.com/installing-python/))
2. Download and install Anaconda ([Windows](https://docs.anaconda.com/anaconda/install/windows/), [macOS](https://docs.anaconda.com/anaconda/install/mac-os/), [Linux](https://docs.anaconda.com/anaconda/install/linux/))


## How to Install SCPI_PKG in Python
Once you have installed the latest version of Python run

`pip install -i https://test.pypi.org/simple/ scpi_pkg==0.2.3`

or

`python -m pip install --index-url https://test.pypi.org/simple/ --no-deps scpi_pkg==0.2.3`


## Common errors in Python

### Import Errors
If python throws the following errors:
- Error (typically Anaconda):
```
ImportError:

IMPORTANT: PLEASE READ THIS FOR ADVICE ON HOW TO SOLVE THIS ISSUE!

Importing the numpy C-extensions failed. This error can happen for
many reasons, often due to issues with your setup or how NumPy was
installed.
```
Solution (in anaconda command prompt):
```
pip uninstall pandas && pip uninstall numpy && pip install pandas
```
- Error:
```
from . import _imaging as core
ImportError: DLL load failed: The specified procedure could not be found.
```
Solution:
```
pip install --upgrade pillow
```
- Error:
```
from matplotlib import ft2font
ImportError: DLL load failed: The specified procedure could not be found.
```
Solution:
```
pip install --upgrade matplotlib
```

## How to link Stata and Python
First of all, Stata and Python can be linked if the running version of Stata is >= 16.0 and the Python's one is >= 2.7. You can follow the [official tutorial](https://blog.stata.com/2020/08/18/stata-python-integration-part-1-setting-up-stata-to-use-python/) on the Stata blog on how to setting up Stata to use Python.

## How to Install SCPI in Stata

```net install scpi, from(ADD_LINK) replace```

## Warnings/Errors in Stata

### `pandas` warning
If you receive the following message when running `scdata`
```
 FutureWarning: pandas.Int64Index is deprecated and will be removed from pandas in a 
 future version. Use pandas.Index with the appropriate dtype instead.
  from pandas import Int64Index as NumericIndex
```
it is just a warning thrown by the Python library `pandas`. More information available [here](https://pandas.pydata.org/docs/reference/api/pandas.Int64Index.html).

### `pickle` warning
If Python does not find the module `pickle` just install it throught the command line (eg. `pip install pickle`) or follow the [official guidelines](https://docs.python.org/3/library/pickle.html). However,
for almost all Python versions, 'pickle' comes as a defaul module with the Python interpreter.

### `cannot erase, read permission only` error
In most of the cases this is because you are running you .do file in a Dropbox folder that Dropbox is trying to keep updated. The same happens if you `cd` is set to a Drobpox folder. An easy fix is to 
pause file syncronization in Drobpox (see [here](https://help.dropbox.com/installs-integrations/sync-uploads/pause-resume) for help).
