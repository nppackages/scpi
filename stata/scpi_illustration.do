******************************************************************************** 
** SCPI Stata Package
** Do-file for Empirical Illustration 
** Authors: Matias D. Cattaneo, Yingjie Feng, Filippo Palomba and Rocio Titiunik 
********************************************************************************
** hlp2winpdf, cdn(scdata) replace
** hlp2winpdf, cdn(scest) replace
** hlp2winpdf, cdn(scpi) replace
** hlp2winpdf, cdn(scplot) replace
********************************************************************************
** net install scpi, from([ADD LINK]) replace
********************************************************************************
clear all
set more off
set linesize 80

********************************************************************************
** load data
********************************************************************************
use "scpi_germany.dta", clear

********************************************************************************
** prepare data - one feature, cointegrated data
********************************************************************************
scdata gdp, dfname("python_scdata") id(country) outcome(gdp) time(year) ///
			treatment(status) cointegrated constant

********************************************************************************
** SC - point estimation with simplex
********************************************************************************
scest, dfname("python_scdata") name(simplex) 
scest, dfname("python_scdata") p(1) q(1) direc("==") lb(0)
********************************************************************************
** SC - plot results
********************************************************************************
scplot, scest gphoptions("ytitle(GDP per capita) xtitle(Year)")

********************************************************************************
** SC - point estimation with lasso
********************************************************************************
scest, dfname("python_scdata") name(lasso)
scest, dfname("python_scdata") p(1) direc("<=") q(1) lb("-inf")        // equivalent

********************************************************************************
** SC - point estimation with ridge
********************************************************************************
scest, dfname("python_scdata") name(ridge)
local Qridge = e(q)
scest, dfname("python_scdata") p(2) direc("<=") q(`Qridge') lb("-inf")  // equivalent

********************************************************************************
** SC - point estimation with least squares
********************************************************************************
scest, dfname("python_scdata") name(ols)
scest, dfname("python_scdata") p(0) lb("-inf") 					     // equivalent


********************************************************************************
** SC - Inference with simplex, misspecified model
********************************************************************************
set seed 8894
scpi, dfname("python_scdata") name(simplex) u_missp u_order(1) u_lags(0) ///
	  u_sigma("HC1") e_order(1) e_lags(0) e_method(qreg)


********************************************************************************
** SC - plot results
********************************************************************************
scplot, gphoptions("ytitle(GDP per capita) xtitle(Year)")


********************************************************************************
** other examples of data preparation
********************************************************************************

* multiple features
scdata gdp trade, dfname("python_scdata") id(country) outcome(gdp) time(year) ///
				  treatment(status) cointegrated				  
				  
* multiple features and feature-specific covariate adjustement 
scdata gdp infrate, dfname("python_scdata") id(country) outcome(gdp) time(year) ///
				  treatment(status) cointegrated covadj("constant, trend; constant")


erase "python_scdata.obj"
erase "__scest__output.obj"
erase "__scpi__output.obj"






