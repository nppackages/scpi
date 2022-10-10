******************************************************************************** 
** SCPI Stata Package
** Do-file for Visualization - Multiple Treated Units 
** Authors: Matias D. Cattaneo, Yingjie Feng, Filippo Palomba and Rocio Titiunik 
********************************************************************************
** net install scpi, from(https://raw.githubusercontent.com/nppackages/scpi/master/stata) replace
********************************************************************************

* ssc install grc1leg

********************************************************************************
*** Multiple Treated Unit
********************************************************************************

** Load data
use "scpi_germany.dta", clear

* suppose a second unit was treated
replace status = 1 if country == "Italy" & year >= 1994

** Prepare data
scdatamulti gdp, dfname("python_scdata") id(country) outcome(gdp) time(year) ///
				 treatment(status) cointegrated(True)

** SC - Inference 
set seed 8894
scpi, dfname("python_scdata") name(simplex) u_missp sims(100) e_method(gaussian)

scplotmulti , savedata("scplot_data") joint

use "scplot_data.dta", clear

drop __0000*
global graph_names
qui levelsof ID, local(tr_units)
local plotname = 0
foreach tr of local tr_units {
	
	local ++plotname 
	qui su Tdate if ID == "`tr'"
	local tline = r(mean)
	
	twoway (rarea Lower_joint Upper_joint Time if Type == "Synthetic" & ID == "`tr'", color(blue%10) lcolor(blue%0))	                                         ///
		   (rcap Lower_gaussian Upper_gaussian Time if Type == "Synthetic" & ID == "`tr'", lpattern(solid) lwidth(medium) lcolor(blue))  ///
		   (scatter Outcome Time if Type == "Actual" & ID == "`tr'", msymbol(Dh) msize(vsmall) mcolor(black))     		  			     ///
		   (scatter Outcome Time if Type == "Synthetic" & ID == "`tr'", msymbol(Dh) msize(vsmall) mcolor(blue))  		  			     ///
		   (line Outcome Time if Type == "Actual" & ID == "`tr'", lpattern(solid) lwidth(medium) lcolor(black))     		  			 ///
		   (line Outcome Time if Type == "Synthetic" & ID == "`tr'", lpattern(solid) lwidth(medium) lcolor(blue)), 		  			     ///
		   ylabel(,nogrid) xline(`tline', lcolor(black) lpattern(dash)) ytitle("") xtitle("") title("") 								 ///
		   legend(order(5 6) lab(5 "Treated") lab(6 "Synthetic Control") region(style(none)) nobox size(vsmall))   				  		 ///
		   graphregion(color(white)) plotregion(color(white)) scheme(s2manual) name("__gr__`plotname'", replace) nodraw 

    global graph_names "$graph_names __gr__`plotname'"
			
}

grc1leg $graph_names, l1title("Outcome Variable", size(small)) b1title("Time", size(small)) graphregion(color(white)) plotregion(color(white)) scheme(s2manual) ///
			  title("") ring(2) ycommon xcommon  

erase "scplot_data.dta"
erase "python_scdata.obj"
