******************************************************************************** 
** SCPI Stata Package
** Do-file for Visualization 
** Authors: Matias D. Cattaneo, Yingjie Feng, Filippo Palomba and Rocio Titiunik 
********************************************************************************
** net install scpi, from([ADD LINK]) replace
********************************************************************************
** Load data
use "scpi_germany.dta", clear

** Prepare data
scdata gdp, dfname("python_scdata") id(country) outcome(gdp) time(year) treatment(status) cointegrated

** SC - Inference 
set seed 8894
scpi, dfname("python_scdata") name(simplex) u_missp sims(200) e_method(qreg)

scplot , savedata("scplot_data")

** Load processed data
use "scplot_data.dta", clear

local tline = Tdate[1]

twoway (rcap lb0 ub0 time, lpattern(solid) lwidth(medium) lcolor(blue))     ///
	   (scatter y_act time, msymbol(Dh) msize(small) mcolor(gray))  		///
	   (scatter y_sc time,  msymbol(Dh) msize(small) mcolor(blue))  		///
	   (line y_act time, lpattern(solid) lwidth(medium) lcolor(gray))  	    ///
	   (line y_sc time,  lpattern(dash) lwidth(medium) lcolor(blue)), 	    ///
	   ytitle("Outcome Variable") xtitle("Time") ylabel(,nogrid)    		///
	   xline(`tline', lcolor(black) lpattern(dash))                 		///
	   legend(order(4 5) lab(4 "Treated") lab(5 "Synthetic Control")        ///
	   region(style(none)) nobox) graphregion(color(white)) 				///
	   plotregion(color(white)) scheme(s2manual) 							///
	   note("In-and-out of sample Uncertainty- Quantile Regression")  
	    
erase "scplot_data.dta"
erase "python_scdata.obj"