******************************************************************************** 
** SCPI Stata Package
** Do-file for Visualization - Single Treated Unit
** Authors: Matias D. Cattaneo, Yingjie Feng, Filippo Palomba and Rocio Titiunik 
********************************************************************************
** net install scpi, from(https://raw.githubusercontent.com/nppackages/scpi/master/stata) replace
********************************************************************************

********************************************************************************
*** Single Treated Unit
********************************************************************************

** Load data
use "scpi_germany.dta", clear

** Prepare data
scdata gdp, dfname("python_scdata") id(country) outcome(gdp) time(year) ///
			treatment(status) cointegrated

** SC - Inference 
set seed 8894
scpi, dfname("python_scdata") name(simplex) u_missp sims(100) ///
	e_method(gaussian) set_seed(8894)

scplot , savedata("scplot_data") joint

** Load processed data
use "scplot_data.dta", clear

local tline = Tdate[1]

twoway (rarea lbj ubj time, color(blue%10) lcolor(blue%0))                  ///
	   (rcap lb0 ub0 time, lpattern(solid) lwidth(medium) lcolor(blue))     ///
	   (scatter y_act time, msymbol(Dh) msize(small) mcolor(black))  		///
	   (scatter y_sc time,  msymbol(Dh) msize(small) mcolor(blue))  		///
	   (line y_act time, lpattern(solid) lwidth(medium) lcolor(black))  	///
	   (line y_sc time,  lpattern(dash) lwidth(medium) lcolor(blue)), 	    ///
	   ytitle("Outcome Variable") xtitle("Time") ylabel(,nogrid)    		///
	   xline(`tline', lcolor(black) lpattern(dash))                 		///
	   legend(order(5 6) lab(5 "Treated") lab(6 "Synthetic Control")        ///
	   region(style(none)) nobox) graphregion(color(white)) 				///
	   plotregion(color(white)) scheme(s2manual) 							///
	   note("In-and-out of sample Uncertainty - sub-Gaussian bounds")  
	    
erase "scplot_data.dta"
erase "python_scdata.obj"
