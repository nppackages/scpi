*****************************************************************************************
** SCPI Stata Package
** Do-file for Empirical Illustration - Single Treated Unit
** Authors: Matias D. Cattaneo, Yingjie Feng, Filippo Palomba and Rocio Titiunik 
*****************************************************************************************
** hlp2pdf scdata, replace
** hlp2pdf scdatamulti, replace
** hlp2pdf scest, replace
** hlp2pdf scpi, replace
** hlp2pdf scplot, replace
** hlp2pdf scplotmulti, replace
  
*****************************************************************************************
** net install scpi, from(https://raw.githubusercontent.com/nppackages/scpi/master/stata) replace
*****************************************************************************************
clear all
set more off
set linesize 80

*****************************************************************************************
*****************************************************************************************
** SINGLE TREATED UNIT
*****************************************************************************************
*****************************************************************************************

*****************************************************************************************
** load data
*****************************************************************************************
use "scpi_germany.dta", clear

*****************************************************************************************
** prepare data - one feature, cointegrated data
*****************************************************************************************
scdata gdp, dfname("python_scdata") id(country) outcome(gdp) time(year) ///
			treatment(status) cointegrated constant 

*****************************************************************************************
** SC - point estimation with simplex
*****************************************************************************************
scest, dfname("python_scdata") name(simplex) 
scest, dfname("python_scdata") p(1) q(1) direc("==") lb(0)

*****************************************************************************************
** SC - plot results
*****************************************************************************************
scplot, scest gphoptions("ytitle(GDP per capita) xtitle(Year)")

*****************************************************************************************
** SC - point estimation with lasso
*****************************************************************************************
scest, dfname("python_scdata") name(lasso)
scest, dfname("python_scdata") p(1) direc("<=") q(1) lb("-inf")  // equivalent

*****************************************************************************************
** SC - point estimation with ridge
*****************************************************************************************
scest, dfname("python_scdata") name(ridge)
local Qridge = e(Qstar)[1,1]
di `Qridge'
scest, dfname("python_scdata") p(2) direc("<=") q(`Qridge') lb("-inf") // equivalent

*****************************************************************************************
** SC - point estimation with L1-L2
*****************************************************************************************
scest, dfname("python_scdata") name(L1-L2)

*****************************************************************************************
** SC - point estimation with least squares
*****************************************************************************************
scest, dfname("python_scdata") name(ols)
scest, dfname("python_scdata") p(0) lb("-inf") 	// equivalent

*****************************************************************************************
** SC - Inference with simplex, misspecified model
*****************************************************************************************
set seed 8894
scpi, dfname("python_scdata") name(simplex) u_missp u_order(1) u_lags(0) ///
	  u_sigma("HC1") e_order(1) e_lags(0) e_method(gaussian)


*****************************************************************************************
** SC - plot results
*****************************************************************************************
scplot, gphoptions("ytitle(GDP per capita) xtitle(Year)")


*****************************************************************************************
** other examples of data preparation
*****************************************************************************************

* multiple features
scdata gdp trade, dfname("python_scdata") id(country) outcome(gdp) time(year) ///
				  treatment(status) cointegrated				  
				  
* multiple features and feature-specific covariate adjustement 
scdata gdp infrate, dfname("python_scdata") id(country) outcome(gdp) time(year) ///
				  treatment(status) cointegrated covadj("constant, trend; constant")

erase "python_scdata.obj"
erase "__scest__output.obj"
erase "__scpi__output.obj"
