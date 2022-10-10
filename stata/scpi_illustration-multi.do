******************************************************************************************
** SCPI Stata Package
** Do-file for Empirical Illustration - Multiple Treated Unit
** Authors: Matias D. Cattaneo, Yingjie Feng, Filippo Palomba and Rocio Titiunik 
******************************************************************************************
  
******************************************************************************************
** net install scpi, from(https://raw.githubusercontent.com/nppackages/scpi/master/stata) replace
******************************************************************************************
clear all
set more off
set linesize 80

******************************************************************************************
******************************************************************************************
** MULTIPLE TREATED UNITS
******************************************************************************************
******************************************************************************************

use "scpi_germany.dta", clear

* Create a second placebo treated unit
replace status = 1 if country == "Italy" & year >= 1992

******************************************************************************************
** prepare data - one feature, cointegrated data, unit-time treatment effect
******************************************************************************************

global coint "True"			
global consta "True"			
global covs "constant, trend"			
					
scdatamulti gdp trade, dfname("python_scdatamulti") ///
			id(country) outcome(gdp) time(year) treatment(status) cointegrated($coint)   ///
			constant($consta) covadj($covs)
			
scest, dfname("python_scdatamulti") name(simplex) 
scplotmulti, scest 

scpi, dfname("python_scdatamulti") name(simplex) sims(200) 

* plot series
scplotmulti, uncertainty("gaussian") ptype("series") 

* plot treatment effects
scplotmulti, uncertainty("gaussian") ptype("treatment") joint 


******************************************************************************************
** prepare data - one feature, cointegrated data, average unit post-treatment effect
******************************************************************************************

* note that you can flexibly specify options for different treated units!					

global coint "(Italy: True) (West Germany: False)"			
global consta "(Italy: True) (West Germany: False)"			
global covs "(Italy: constant, trend) (West Germany: constant, trend; constant, trend)"			
scdatamulti (Italy: gdp trade) (West Germany: gdp infrate), dfname("python_scdatamulti") ///
			id(country) outcome(gdp) time(year) treatment(status) cointegrated($coint)   ///
			constant($consta) covadj($covs) effect("unit") 
			
scest, dfname("python_scdatamulti") name(simplex) 
scplotmulti, scest 

scpi, dfname("python_scdatamulti") name(simplex) sims(200) 

* plot series
scplotmulti, uncertainty("gaussian") ptype("series") 

* plot treatment effects
scplotmulti, uncertainty("gaussian") ptype("treatment") joint 


******************************************************************************************
** prepare data - one feature, cointegrated data, average treatment effect on the treated
******************************************************************************************

global coint "(Italy: True) (West Germany: False)"			
global consta "(Italy: True) (West Germany: False)"			
global covs "(Italy: constant, trend) (West Germany: constant, trend; constant, trend)"			
					
scdatamulti (Italy: gdp trade) (West Germany: gdp infrate), dfname("python_scdatamulti") ///
			id(country) outcome(gdp) time(year) treatment(status) cointegrated($coint)   ///
			constant($consta) covadj($covs) effect("time") 
			
scest, dfname("python_scdatamulti") name(simplex) 
scplotmulti, scest 

scpi, dfname("python_scdatamulti") name(simplex) sims(200) 

* plot series
scplotmulti, uncertainty("gaussian") ptype("series") 

* plot treatment effects
scplotmulti, uncertainty("gaussian") ptype("treatment") joint 


erase "python_scdatamulti.obj"
erase "__scest__output.obj"
erase "__scpi__output.obj"
