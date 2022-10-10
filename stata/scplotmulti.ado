*! Date        : 07 Oct 2022
*! Version     : 2.0
*! Authors     : Filippo Palomba
*! Email       : fpalomba@princeton.edu
*! Description : Plot Synthetic Control


capture program drop scplotmulti
program define scplotmulti, eclass         
version 17.0           

	syntax , [scest uncertainty(string) ptype(string) joint yscalefree xscalefree    					  ///
										dots_tr_col(string) dots_tr_symb(string) dots_tr_size(string)     ///
										dots_sc_col(string) dots_sc_symb(string) dots_sc_size(string)     ///
										line_tr_col(string) line_tr_patt(string) line_tr_width(string)    ///
										line_sc_col(string) line_sc_patt(string) line_sc_width(string)    ///
										spike_sc_col(string) spike_sc_patt(string) spike_sc_width(string) ///
										gphoptions(string) gphcombineoptions(string) gphsave(string) savedata(string) keepsingleplots pypinocheck]

	tempvar outvar CIlb CIub idcode
										
	if mi("`pypinocheck'") & mi("$scpi_version_checked") {
		python: version_checker()
		if "`alert_version'" == "y" {
			di as error "The current version of scpi_pkg in Python is `python_local_version', but version `python_pypi_version' needed! Please update the package in Python and restart Stata!"
			exit 198
		}
		global scpi_version_checked "yes"
	}										
	
	if !mi("`uncertainty'") {	
		if !inlist("`uncertainty'", "insample", "gaussian", "ls", "qreg") {
			di as error "{err}The option 'uncertainty' should be one of 'insample', 'gaussian', 'ls', or 'qreg'!"
			exit 198
		}
	}
	else {
		local uncertainty "gaussian"
	}	
	
	if "`uncertainty'" == "insample" {
		local e_out = "False"
		local e_method = "gaussian"
	}
	else {
		local e_out = "True"
		local e_method = "`uncertainty'"
	}

	if !mi("`ptype'") {	
		if !inlist("`ptype'", "series", "treatment") {
			di as error "{err}The option 'ptype' should be one of 'series', or 'treatment'!"
			exit 198
		}
	}
	else {
		local ptype "series"
	}	

	if !mi("`joint'") {	
		local joint = "True"
	}
	else {
		local joint = "False"
	}
	
	if !mi("`yscalefree'") {
		local ycom
	}
	else {
		local ycom = "ycommon"
	}

	if !mi("`xscalefree'") {
		local xcom
	}
	else {
		local xcom = "xcommon"
	}	
	
	if mi("`dots_tr_col'") {
		local dots_tr_col "black"
	}
	if mi("`dots_tr_size'") {
		local dots_tr_size "vsmall"
	}
	if mi("`dots_tr_symb'") {
		local dots_tr_symb "Dh"
	}
	
	if mi("`dots_sc_col'") {
		local dots_sc_col "blue"
	}
	if mi("`dots_sc_size'") {
		local dots_sc_size "vsmall"
	}
	if mi("`dots_sc_symb'") {
		local dots_sc_symb "Dh"
	}
	
	if mi("`line_tr_col'") {
		local line_tr_col "black"
	}
	if mi("`line_tr_patt'") {
		local line_tr_patt "solid"
	}
	if mi("`line_tr_width'") {
		local line_tr_width "medium"
	}

	if mi("`line_sc_col'") {
		local line_sc_col "blue"
	}
	if mi("`line_sc_patt'") {
		local line_sc_patt "solid"
	}
	if mi("`line_sc_width'") {
		local line_sc_width "medium"
	}

	if mi("`spike_sc_col'") {
		local spike_sc_col "blue"
	}
	if mi("`spike_sc_patt'") {
		local spike_sc_patt "solid"
	}
	if mi("`spike_sc_width'") {
		local spike_sc_width "medium"
	}	
	
	
	local last_object "__scpi__output.obj"
	
	if !mi("`scest'") {
		local last_object "__scest__output.obj"
	}
	
	* Prepare data for plots using Python
	python: scplot_loader("`last_object'", "`ptype'", "`joint'", "`e_out'", "`e_method'")
	
	* Prepare plots
	preserve
	use "__scpi__stata_plot.dta", clear
	
	global graph_names
	qui levelsof ID, local(tr_units)
	
	if "`ptype'" == "series" {
		qui g `outvar' = Outcome
		local ytitle = "Outcome Variable"
	}
	else {
		qui g `outvar' = Effect
		local ytitle = "Effect"
	}

	if mi("`scest'") {
		if "`uncertainty'" == "insample" {
			qui g `CIlb' = Lower_insample
			qui g `CIub' = Upper_insample
			local graph_note "In-sample uncertainty"
		}

		if "`uncertainty'" == "gaussian" {
			qui g `CIlb' = Lower_gaussian
			qui g `CIub' = Upper_gaussian
			local graph_note "Out-of-sample uncertainty: subgaussian bounds"
		}

		if "`uncertainty'" == "ls" {
			qui g `CIlb' = Lower_ls
			qui g `CIub' = Upper_ls
			local graph_note "Out-of-sample uncertainty: location-scale model"
		}

		if "`uncertainty'" == "qreg" {
			qui g `CIlb' = Lower_qreg
			qui g `CIub' = Upper_qreg
			local graph_note "Out-of-sample uncertainty: quantile regression"
		}	
	}
	
	if "`effect'" == "time" {
		local effect = "unit-time"
	}
	
	local plotname = 0
	foreach tr of local tr_units {
		
		local ++plotname 
		
		qui su Tdate if ID == "`tr'" & Time < Tdate
		local tline = r(mean) + 1/2

		if !mi("`keepsingleplots'") {
			local savingplot "saving('__scplotmulti__`plotname'', replace)"
		}			

		
		*****************************************************************************
		* Scatter with Lines
		*****************************************************************************
		
		if "`last_object'" == "__scest__output.obj" {
		
			if "`ptype'" == "series" {
				
				twoway (scatter `outvar' Time if Type == "Actual" & ID == "`tr'", msymbol(`dots_tr_symb') msize(`dots_tr_size') mcolor(`dots_tr_col'))     ///
					   (scatter `outvar' Time if Type == "Synthetic" & ID == "`tr'", msymbol(`dots_sc_symb') msize(`dots_sc_size') mcolor(`dots_sc_col'))  ///
					   (line `outvar' Time if Type == "Actual" & ID == "`tr'", lpattern(`line_tr_patt') lwidth(`line_tr_width') lcolor(`line_tr_col'))     ///
					   (line `outvar' Time if Type == "Synthetic" & ID == "`tr'", lpattern(`line_sc_patt') lwidth(`line_tr_width') lcolor(`line_sc_col')), ///
					   ylabel(,nogrid) xline(`tline', lcolor(black) lpattern(dash)) ytitle("") xtitle("") title("`tr'", size(vsmall))    				   ///
					   legend(order(3 4) lab(3 "Treated") lab(4 "Synthetic Control") region(style(none)) nobox size(vsmall))   							   ///
					   graphregion(color(white)) plotregion(color(white)) scheme(s2manual) 						  										   ///
					   `gphoptions' name("__gr__`plotname'", replace) nodraw `savingplot'
				
			}

			if "`ptype'" == "treatment" {
				
				twoway (scatter `outvar' Time if ID == "`tr'", msymbol(`dots_tr_symb') msize(`dots_tr_size') mcolor(`dots_tr_col'))    ///
					   (line `outvar' Time if ID == "`tr'", lpattern(`line_tr_patt') lwidth(`line_tr_width') lcolor(`line_tr_col'))    ///
					   ylabel(,nogrid) xline(`tline', lcolor(black) lpattern(dash)) ytitle("") xtitle("") title("`tr'", size(vsmall))  ///
					   graphregion(color(white)) plotregion(color(white)) scheme(s2manual) legend(order(2) lab(2 "Effect"))            ///
					   yline(0, lp("dot") lw("thin")) `gphoptions' name("__gr__`plotname'", replace) nodraw `savingplot'
				
			}
			
			local graphnote ""
		}

		
		*****************************************************************************
		* Scatter with Lines and Uncertainty Intervals/Bands - Unit-time
		*****************************************************************************

		if "`last_object'" == "__scpi__output.obj" & "`effect'" == "unit-time" {

			if "`ptype'" == "series" {

				if "`joint'" == "False" {

					twoway (rcap `CIlb' `CIub' Time if Type == "Synthetic" & ID == "`tr'", lpattern(`spike_sc_patt') lwidth(`spike_sc_width') lcolor(`spike_sc_col')) ///
						   (scatter `outvar' Time if Type == "Actual" & ID == "`tr'", msymbol(`dots_tr_symb') msize(`dots_tr_size') mcolor(`dots_tr_col'))     		  ///
						   (scatter `outvar' Time if Type == "Synthetic" & ID == "`tr'", msymbol(`dots_sc_symb') msize(`dots_sc_size') mcolor(`dots_sc_col'))  		  ///
						   (line `outvar' Time if Type == "Actual" & ID == "`tr'", lpattern(`line_tr_patt') lwidth(`line_tr_width') lcolor(`line_tr_col'))     		  ///
						   (line `outvar' Time if Type == "Synthetic" & ID == "`tr'", lpattern(`line_sc_patt') lwidth(`line_tr_width') lcolor(`line_sc_col')), 		  ///
						   ylabel(,nogrid) xline(`tline', lcolor(black) lpattern(dash)) ytitle("") xtitle("") title("`tr'", size(vsmall)) 					  		  ///
						   legend(order(4 5) lab(4 "Treated") lab(5 "Synthetic Control") region(style(none)) nobox size(vsmall))   							  		  ///
						   graphregion(color(white)) plotregion(color(white)) scheme(s2manual) `gphoptions' name("__gr__`plotname'", replace) nodraw `savingplot'
				
				}

				
				if "`joint'" == "True" {

					twoway (rarea Lower_joint Upper_joint Time if Type == "Synthetic" & ID == "`tr'", color(`spike_sc_col'%10) lcolor(`spike_sc_col'%0))                         ///
						   (rcap `CIlb' `CIub' Time if Type == "Synthetic" & ID == "`tr'", lpattern(`spike_sc_patt') lwidth(`spike_sc_width') lcolor(`spike_sc_col'))            ///
						   (scatter `outvar' Time if Type == "Actual" & ID == "`tr'", msymbol(`dots_tr_symb') msize(`dots_tr_size') mcolor(`dots_tr_col'))     		  			 ///
						   (scatter `outvar' Time if Type == "Synthetic" & ID == "`tr'", msymbol(`dots_sc_symb') msize(`dots_sc_size') mcolor(`dots_sc_col'))  		  			 ///
						   (line `outvar' Time if Type == "Actual" & ID == "`tr'", lpattern(`line_tr_patt') lwidth(`line_tr_width') lcolor(`line_tr_col'))     		  			 ///
						   (line `outvar' Time if Type == "Synthetic" & ID == "`tr'", lpattern(`line_sc_patt') lwidth(`line_tr_width') lcolor(`line_sc_col')), 		  			 ///
						   ylabel(,nogrid) xline(`tline', lcolor(black) lpattern(dash)) ytitle("") xtitle("") title("`tr'", size(vsmall))       			  		  			 ///
						   legend(order(5 6) lab(5 "Treated") lab(6 "Synthetic Control") region(style(none)) nobox size(vsmall))   							  		  			 ///
						   graphregion(color(white)) plotregion(color(white)) scheme(s2manual) `gphoptions' name("__gr__`plotname'", replace) nodraw `savingplot'
				
					local graph_note "`graph_note' and joint prediction intervals"
				}				
			}
			

			if "`ptype'" == "treatment" {

				if "`joint'" == "False" {

					twoway (rcap `CIlb' `CIub' Time if ID == "`tr'", lpattern(`spike_sc_patt') lwidth(`spike_sc_width') lcolor(`spike_sc_col')) 				      ///
						   (scatter `outvar' Time if ID == "`tr'", msymbol(`dots_sc_symb') msize(`dots_sc_size') mcolor(`dots_sc_col'))  		                      ///
						   (line `outvar' Time if ID == "`tr'", lpattern(`line_sc_patt') lwidth(`line_sc_width') lcolor(`line_sc_col')), 		  					  ///
						   ylabel(,nogrid) xline(`tline', lcolor(black) lpattern(dash)) ytitle("") xtitle("") title("`tr'", size(vsmall))                             ///
						   yline(0, lp("dot") lw("thin")) legend(order(2 1) lab(1 "Prediction Intervals") lab(2 "Effect"))											  ///
						   graphregion(color(white)) plotregion(color(white)) scheme(s2manual) `gphoptions' name("__gr__`plotname'", replace) nodraw `savingplot'
				
				}

				
				if "`joint'" == "True" {
					
					twoway (rarea Lower_joint Upper_joint Time if ID == "`tr'", color(`spike_sc_col'%10) lcolor(`spike_sc_col'%0))	                             ///
						   (rcap `CIlb' `CIub' Time if ID == "`tr'", lpattern(`spike_sc_patt') lwidth(`spike_sc_width') lcolor(`spike_sc_col')) 				 ///
						   (scatter `outvar' Time if ID == "`tr'", msymbol(`dots_sc_symb') msize(`dots_sc_size') mcolor(`dots_sc_col'))     		  			 ///
						   (line `outvar' Time if ID == "`tr'", lpattern(`line_sc_patt') lwidth(`line_sc_width') lcolor(`line_sc_col')),     		  			 ///
						   ylabel(,nogrid) xline(`tline', lcolor(black) lpattern(dash)) ytitle("") xtitle("") title("`tr'", size(vsmall))                        ///
						   yline(0, lp("dot") lw("thin")) legend(order(3 2) lab(2 "Prediction Intervals") lab(3 "Effect"))										 ///
						   graphregion(color(white)) plotregion(color(white)) scheme(s2manual) `gphoptions' name("__gr__`plotname'", replace) nodraw `savingplot'
						   
					local graph_note "`graph_note' and joint prediction intervals"
				}				
			}			
			
		}				   		   


		*****************************************************************************
		* Scatter with Lines and Uncertainty Intervals/Bands - Unit
		*****************************************************************************

		if "`last_object'" == "__scpi__output.obj" & "`effect'" == "unit" {

		
			if "`ptype'" == "series" {

				if "`joint'" == "False" {

					twoway (rcap `CIlb' `CIub' Time if Type == "Synthetic" & ID == "`tr'", lpattern(`spike_sc_patt') lwidth(`spike_sc_width') lcolor(`spike_sc_col'))                     ///
						   (scatter `outvar' Time if Type == "Actual" & ID == "`tr'", msymbol(`dots_tr_symb') msize(`dots_tr_size') mcolor(`dots_tr_col'))     		  					  ///
						   (scatter `outvar' Time if Type == "Synthetic" & ID == "`tr'", msymbol(`dots_sc_symb') msize(`dots_sc_size') mcolor(`dots_sc_col'))  		  		              ///
						   (line `outvar' Time if Type == "Actual" & ID == "`tr'" & Time < Tdate, lpattern(`line_tr_patt') lwidth(`line_tr_width') lcolor(`line_tr_col'))     		      ///
						   (line `outvar' Time if Type == "Synthetic" & ID == "`tr'" & Time < Tdate, lpattern(`line_sc_patt') lwidth(`line_tr_width') lcolor(`line_sc_col')), 		      ///
						   ylabel(,nogrid) ytitle("") xtitle("") title("`tr'", size(vsmall))              				 									  		  					  ///
						   legend(order(4 5) lab(4 "Treated") lab(5 "Synthetic Control") region(style(none)) nobox size(vsmall))   							  		  					  ///
						   graphregion(color(white)) plotregion(color(white)) scheme(s2manual) `gphoptions' name("__gr__`plotname'", replace) nodraw `savingplot'
				
				}

				
				if "`joint'" == "True" {

					twoway (rcap Lower_joint Upper_joint Time if Type == "Synthetic" & ID == "`tr'", lcolor(`spike_sc_col'%20)) 												 ///
						   (rcap `CIlb' `CIub' Time if Type == "Synthetic" & ID == "`tr'", lpattern(`spike_sc_patt') lwidth(`spike_sc_width') lcolor(`spike_sc_col'))            ///
						   (scatter `outvar' Time if Type == "Actual" & ID == "`tr'", msymbol(`dots_tr_symb') msize(`dots_tr_size') mcolor(`dots_tr_col'))     		  			 ///
						   (scatter `outvar' Time if Type == "Synthetic" & ID == "`tr'", msymbol(`dots_sc_symb') msize(`dots_sc_size') mcolor(`dots_sc_col'))  		  			 ///
						   (line `outvar' Time if Type == "Actual" & ID == "`tr'" & Time < Tdate, lpattern(`line_tr_patt') lwidth(`line_tr_width') lcolor(`line_tr_col'))        ///
						   (line `outvar' Time if Type == "Synthetic" & ID == "`tr'" & Time < Tdate, lpattern(`line_sc_patt') lwidth(`line_tr_width') lcolor(`line_sc_col')), 	 ///
						   ylabel(,nogrid) ytitle("") xtitle("") title("`tr'", size(vsmall))																  		  			 ///
						   legend(order(5 6) lab(5 "Treated") lab(6 "Synthetic Control") region(style(none)) nobox size(vsmall))   							  		  			 ///
						   graphregion(color(white)) plotregion(color(white)) scheme(s2manual) `gphoptions' name("__gr__`plotname'", replace) nodraw `savingplot'
				
					local graph_note "`graph_note' and joint prediction intervals"
				}				
			}
			
		}				   		   
		
		global graph_names "$graph_names __gr__`plotname'"
	
	}
		
	if "`ptype'" == "treatment" & "`effect'" == "unit" {
		
		qui encode(ID), g(`idcode')
		qui ta ID
		local trvars = r(r)
		
		if "`joint'" == "False" {
			
			twoway (rcap `CIlb' `CIub' `idcode' if Treatment == 1, lpattern(`spike_sc_patt') lwidth(`spike_sc_width') lcolor(`spike_sc_col')) 					  ///
				   (scatter `outvar' `idcode' [aw=T1] if Treatment == 1, msymbol(`dots_sc_symb') msize(`dots_sc_size') mcolor(`dots_sc_col')), 				      ///
				   ylabel(,nogrid) ytitle("") xtitle("") title("") yline(0, lp("dot") lw("thin")) xlabel(1(1)`trvars', angle(45) valuelabel labsize(small))    	  ///
				   graphregion(color(white)) plotregion(color(white)) scheme(s2manual) legend(order(2) lab(2 "Effect")) `gphoptions'
		
		}

		
		if "`joint'" == "True" {
			
			twoway (rcap Lower_joint Upper_joint `idcode' if Treatment == 1, lcolor(`spike_sc_col'%20))							              							///
			       (rcap `CIlb' `CIub' `idcode' if Treatment == 1, lpattern(`spike_sc_patt') lwidth(`spike_sc_width') lcolor(`spike_sc_col')) 							///
				   (scatter `outvar' `idcode' [aw=T1] if Treatment == 1, msymbol(`dots_sc_symb') msize(`dots_sc_size') mcolor(`dots_sc_col')),  						///
				   ylabel(,nogrid) ytitle("") xtitle("") title("") yline(0, lp("dot") lw("thin")) xlabel(1(1)`trvars', angle(45) valuelabel labsize(small)) 			///
				   graphregion(color(white)) plotregion(color(white)) scheme(s2manual) legend(order(3 2) lab(2 "Prediction Intervals") lab(3 "Effect"))	`gphoptions' 
				   
			local graph_note "`graph_note' and joint prediction intervals"
			
		}				
	}			
	else {
		grc1leg $graph_names, l1title("`ytitle'", size(small)) b1title("Time", size(small)) graphregion(color(white)) plotregion(color(white)) scheme(s2manual) ///
					  title("") ring(2) `ycom' `xcom' `gphoptionscommon' 
	}
	
		
	qui erase "__scpi__stata_plot.dta"
	qui erase "__scpi__stata_plot.csv"
	
	if !mi("`gphsave'") graph export "`gphsave'.png", replace 
	if !mi("`savedata'") save "`savedata'.dta", replace

	restore

end


version 17.0
python:
import pickle, numpy, pandas, urllib, luddite
from scpi_pkg import version as lver
from scpi_pkg.scplotMulti import scplotMulti
from sfi import Macro
from copy import deepcopy

def scplot_loader(last_object, ptype, joint, e_out, e_method):
	filename    = last_object
	filehandler = open(filename, 'rb') 
	result      = pickle.load(filehandler)
	class_input = result.__class__.__name__

	if e_out == "False":
		e_out_bool = False
	else:
		e_out_bool = True
	
	if joint == "False":
		joint_bool = False
	else:
		joint_bool = True
	
	p = scplotMulti(result=result, ptype=ptype, e_out=e_out_bool, joint=joint_bool, save_data="__scpi__stata_plot", verbose=False)
	df = pandas.read_csv("__scpi__stata_plot.csv")
	df.drop(columns='Unnamed: 0', inplace=True)
	df.to_stata("__scpi__stata_plot.dta", write_index = False)
	
	Macro.setLocal("effect", result.effect)


def version_checker():
	# try to connect to pypi and get the latest version of scpi_pkg
	try:
		local_version = str(lver.__version__)
		pypi_version = luddite.get_version_pypi("scpi_pkg")
		if local_version == pypi_version:
			alert_version = "n"
		else:
			alert_version = "y"
	except urllib.error.URLError:
		alert_version = "n"
		pypi_version = "none"

	Macro.setLocal("alert_version", alert_version)
	Macro.setLocal("python_local_version", local_version)
	Macro.setLocal("python_pypi_version", pypi_version)
	
end
