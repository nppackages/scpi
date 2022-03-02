*! Date        : 1 Mar 2022
*! Version     : 0.2.1
*! Authors     : Filippo Palomba
*! Email       : fpalomba@princeton.edu
*! Description : Plot Synthetic Control


capture program drop scplot
program define scplot, eclass         
version 17.0           

	syntax , [scest uncertainty(string) dots_tr_col(string) dots_tr_symb(string) dots_tr_size(string)     ///
										dots_sc_col(string) dots_sc_symb(string) dots_sc_size(string)     ///
										line_tr_col(string) line_tr_patt(string) line_tr_width(string)    ///
										line_sc_col(string) line_sc_patt(string) line_sc_width(string)    ///
										spike_sc_col(string) spike_sc_patt(string) spike_sc_width(string) ///
										gphoptions(string) gphsave(string) savedata(string) pypinocheck]

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
			di as error "{err}The option 'uncertainty' should be of 'insample', 'gaussian', 'ls', or 'qreg'!"
			exit 198
		}
	}
	
	if mi("`dots_tr_col'") {
		local dots_tr_col "gray"
	}
	if mi("`dots_tr_size'") {
		local dots_tr_size "small"
	}
	if mi("`dots_tr_symb'") {
		local dots_tr_symb "Dh"
	}
	if mi("`dots_sc_col'") {
		local dots_sc_col "blue"
	}
	if mi("`dots_sc_size'") {
		local dots_sc_size "small"
	}
	if mi("`dots_sc_symb'") {
		local dots_sc_symb "Dh"
	}
	
	if mi("`line_tr_col'") {
		local line_tr_col "gray"
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
		local line_sc_patt "dash"
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
	
	python: scplot_loader("`last_object'")
	
	if mi("`uncertainty'"){
		local uncertainty "`e_method'"
	}
	
	preserve
	use "__scpi__stata_plot.dta", clear
	
	if "`last_object'" == "__scest__output.obj" {
		la var y_act "Treated Unit"
		la var y_sc "Synthetic Control Unit"
		la var time "Time index"

		local tline = Tdate[1]
		
		twoway (scatter y_act time, msymbol(`dots_tr_symb') msize(`dots_tr_size') mcolor(`dots_tr_col'))  ///
			   (scatter y_sc time,  msymbol(`dots_sc_symb') msize(`dots_sc_size') mcolor(`dots_sc_col'))  ///
			   (line y_act time, lpattern(`line_tr_patt') lwidth(`line_tr_width') lcolor(`line_tr_col'))  ///
			   (line y_sc time,  lpattern(`line_sc_patt') lwidth(`line_tr_width') lcolor(`line_sc_col')), ///
			   ytitle("Outcome Variable") xtitle("Time") ylabel(,nogrid)    							  ///
			   xline(`tline', lcolor(black) lpattern(dash))                 							  ///
			   legend(order(3 4) lab(3 "Treated") lab(4 "Synthetic Control") region(style(none)) nobox)   ///
			   graphregion(color(white)) plotregion(color(white)) scheme(s2manual) 						  ///
			   title("Synthetic Control Prediction") `gphoptions'

	}
	
	if "`last_object'" == "__scpi__output.obj" {
		la var y_act "Treated Unit"
		la var y_sc "Synthetic Control Unit"
		la var time "Time index"
		la var Tdate "Treatment Date"
		la var lb0 "lower - Only In-Sample"
		la var ub0 "upper - Only In-Sample"
		la var lb1 "lower - In + Out (gaussian)"
		la var ub1 "upper - In + Out (gaussian)"
		la var lb2 "lower - In + Out (location-scale)"
		la var ub2 "upper - In + Out (location-scale)"
		la var lb3 "lower - In + Out (quantile)"
		la var ub3 "upper - In + Out (quantile)"
		la var Tdate "Treatment Date"			
		
		
		
		if "`uncertainty'" == "insample" {
			local tline = Tdate[1]
		
			twoway (rcap lb0 ub0 time, lpattern(`spike_sc_patt') lwidth(`spike_sc_width') lcolor(`spike_sc_col')) ///
				   (scatter y_act time, msymbol(`dots_tr_symb') msize(`dots_tr_size') mcolor(`dots_tr_col'))  		///
				   (scatter y_sc time,  msymbol(`dots_sc_symb') msize(`dots_sc_size') mcolor(`dots_sc_col'))  		///
			       (line y_act time, lpattern(`line_tr_patt') lwidth(`line_tr_width') lcolor(`line_tr_col'))  		///
			       (line y_sc time,  lpattern(`line_sc_patt') lwidth(`line_tr_width') lcolor(`line_sc_col')), 		///
				   ytitle("Outcome Variable") xtitle("Time") ylabel(,nogrid)    									///
				   xline(`tline', lcolor(black) lpattern(dash))                 									///
				   legend(order(4 5) lab(4 "Treated") lab(5 "Synthetic Control") region(style(none)) nobox) 		///
				   graphregion(color(white)) plotregion(color(white)) scheme(s2manual) 								///
				   note("In-sample Uncertainty")  `gphoptions'

		} 
		
		if "`uncertainty'" == "gaussian" | "`uncertainty'" == "all" {
			local tline = Tdate[1]
		
			twoway (rcap lb1 ub1 time, lpattern(`spike_sc_patt') lwidth(`spike_sc_width') lcolor(`spike_sc_col')) ///
				   (scatter y_act time, msymbol(`dots_tr_symb') msize(`dots_tr_size') mcolor(`dots_tr_col'))  		///
				   (scatter y_sc time,  msymbol(`dots_sc_symb') msize(`dots_sc_size') mcolor(`dots_sc_col'))  		///
			       (line y_act time, lpattern(`line_tr_patt') lwidth(`line_tr_width') lcolor(`line_tr_col'))  		///
			       (line y_sc time,  lpattern(`line_sc_patt') lwidth(`line_tr_width') lcolor(`line_sc_col')), 		///
				   ytitle("Outcome Variable") xtitle("Time") ylabel(,nogrid)    									///
				   xline(`tline', lcolor(black) lpattern(dash))                 									///
				   legend(order(4 5) lab(4 "Treated") lab(5 "Synthetic Control") region(style(none)) nobox) 		///
				   graphregion(color(white)) plotregion(color(white)) scheme(s2manual) 								///
				   note("In and Out of Sample Uncertainty - Subgaussian Bounds") `gphoptions'

		} 
		
		if "`uncertainty'" == "ls" | "`uncertainty'" == "all" {
			local tline = Tdate[1]
		
			twoway (rcap lb2 ub2 time, lpattern(`spike_sc_patt') lwidth(`spike_sc_width') lcolor(`spike_sc_col'))  ///
				   (scatter y_act time, msymbol(`dots_tr_symb') msize(`dots_tr_size') mcolor(`dots_tr_col'))  		 ///
				   (scatter y_sc time,  msymbol(`dots_sc_symb') msize(`dots_sc_size') mcolor(`dots_sc_col')) 		 ///
			       (line y_act time, lpattern(`line_tr_patt') lwidth(`line_tr_width') lcolor(`line_tr_col'))  		 ///
			       (line y_sc time,  lpattern(`line_sc_patt') lwidth(`line_tr_width') lcolor(`line_sc_col')), 		 ///
				   ytitle("Outcome Variable") xtitle("Time") ylabel(,nogrid)    									 ///
				   xline(`tline', lcolor(black) lpattern(dash))                 									 ///
				   legend(order(4 5) lab(4 "Treated") lab(5 "Synthetic Control") region(style(none)) nobox) 		 ///
				   graphregion(color(white)) plotregion(color(white)) scheme(s2manual) 								 ///
				   note("In and Out of Sample Uncertainty - Location-scale Model")  `gphoptions'

		} 
		
		if "`uncertainty'" == "qreg" | "`uncertainty'" == "all" {
			local tline = Tdate[1]
		
			twoway (rcap lb3 ub3 time, lpattern(`spike_sc_patt') lwidth(`spike_sc_width') lcolor(`spike_sc_col')) ///
				   (scatter y_act time, msymbol(`dots_tr_symb') msize(`dots_tr_size') mcolor(`dots_tr_col'))        ///
				   (scatter y_sc time,  msymbol(`dots_sc_symb') msize(`dots_sc_size') mcolor(`dots_sc_col'))        ///
			       (line y_act time, lpattern(`line_tr_patt') lwidth(`line_tr_width') lcolor(`line_tr_col'))        ///
			       (line y_sc time,  lpattern(`line_sc_patt') lwidth(`line_tr_width') lcolor(`line_sc_col')),       ///
				   ytitle("Outcome Variable") xtitle("Time") ylabel(,nogrid)    									///
				   xline(`tline', lcolor(black) lpattern(dash))                 									///
				   legend(order(4 5) lab(4 "Treated") lab(5 "Synthetic Control") region(style(none)) nobox) 	 	///
				   graphregion(color(white)) plotregion(color(white)) scheme(s2manual) 								///
				   note("In and Out of Sample Uncertainty - Quantile Regression")  `gphoptions'

		} 
		

	}		
	
	if !mi("`gphsave'") graph export "`gphsave'.png", replace 
	if !mi("`savedata'") save "`savedata'", replace
	
	restore
	
	erase "__scpi__stata_plot.dta"

end


version 17.0
python:
import pickle, numpy, pandas, urllib, luddite
from sfi import Macro

def scplot_loader(last_object):
	filename    = last_object
	filehandler = open(filename, 'rb') 
	result      = pickle.load(filehandler)
	class_input = result.__class__.__name__
	
	period_pre  = result.period_pre
	period_post = result.period_post
	
	time    = numpy.concatenate([period_pre, period_post])
	T0      = period_pre[len(period_pre)-1]
	Tdate   = pandas.DataFrame(numpy.array([T0]*len(time)))
	y_act   = pandas.concat([result.Y_pre, result. Y_post]).to_numpy().flatten()
	y_sc_df = pandas.concat([result.Y_pre_fit, result.Y_post_fit])
	y_sc_na = pandas.DataFrame(numpy.array([numpy.nan]*len(time)))
	
	if class_input == "scest_output":
		not_miss_plot = [t in y_sc_df.index.tolist() for t in time]
		y_sc_na.loc[not_miss_plot,] = y_sc_df.iloc[:,[0]].to_numpy()
		data_points_act     = pandas.DataFrame({'time' : time, 'y_act': y_act})
		data_points         = pandas.concat([data_points_act, y_sc_na, Tdate], axis = 1)
		data_points.columns = ['time', 'y_act', 'y_sc', 'Tdate']
		
	if class_input == "scpi_output":
		e_method = result.e_method
		sc_l_0 = result.CI_in_sample.iloc[:,[0]].to_numpy()
		sc_r_0 = result.CI_in_sample.iloc[:,[1]].to_numpy()
		
		if e_method in ["all", "gaussian"]:
			sc_l_1 = result.CI_all_gaussian.iloc[:,[0]].to_numpy()
			sc_r_1 = result.CI_all_gaussian.iloc[:,[1]].to_numpy()
		else:
			sc_l_1 = sc_r_1 = numpy.empty((len(period_post), 1))

		if e_method in ["all", "ls"]:
			sc_l_2 = result.CI_all_ls.iloc[:,[0]].to_numpy()
			sc_r_2 = result.CI_all_ls.iloc[:,[1]].to_numpy()
		else:
			sc_l_2 = sc_r_2 = numpy.empty((len(period_post), 1))

		if e_method in ["all", "qreg"]:
			sc_l_3 = result.CI_all_qreg.iloc[:,[0]].to_numpy()
			sc_r_3 = result.CI_all_qreg.iloc[:,[1]].to_numpy()
		else:
			sc_l_3 = sc_r_3 = numpy.empty((len(period_post), 1))
			
			
		sc_l_0_na = pandas.DataFrame(numpy.array([numpy.nan]*len(time)))
		sc_r_0_na = pandas.DataFrame(numpy.array([numpy.nan]*len(time)))
		sc_l_1_na = pandas.DataFrame(numpy.array([numpy.nan]*len(time)))
		sc_r_1_na = pandas.DataFrame(numpy.array([numpy.nan]*len(time)))
		sc_l_2_na = pandas.DataFrame(numpy.array([numpy.nan]*len(time)))
		sc_r_2_na = pandas.DataFrame(numpy.array([numpy.nan]*len(time)))
		sc_l_3_na = pandas.DataFrame(numpy.array([numpy.nan]*len(time)))
		sc_r_3_na = pandas.DataFrame(numpy.array([numpy.nan]*len(time)))
		
		not_miss_plot = [t in y_sc_df.index.tolist() for t in time]
		not_miss_ci   = [t in result.CI_in_sample.index.tolist() for t in time]
		
		
		y_sc_na.loc[not_miss_plot, ] = y_sc_df.iloc[:,[0]].to_numpy()
		sc_l_0_na.loc[not_miss_ci, ] = sc_l_0
		sc_r_0_na.loc[not_miss_ci, ] = sc_r_0
		sc_l_1_na.loc[not_miss_ci, ] = sc_l_1
		sc_r_1_na.loc[not_miss_ci, ] = sc_r_1
		sc_l_2_na.loc[not_miss_ci, ] = sc_l_2
		sc_r_2_na.loc[not_miss_ci, ] = sc_r_2
		sc_l_3_na.loc[not_miss_ci, ] = sc_l_3.astype(numpy.float64)
		sc_r_3_na.loc[not_miss_ci, ] = sc_r_3.astype(numpy.float64)

				
		data_points_act = pandas.DataFrame({'time': time, 'y_act': y_act})		
		data_points = pandas.concat([data_points_act, y_sc_na, 
									sc_l_0_na, sc_r_0_na,
									sc_l_1_na, sc_r_1_na,
									sc_l_2_na, sc_r_2_na,
									sc_l_3_na, sc_r_3_na,
									Tdate], axis = 1)
									
		data_points.columns = ['time', 'y_act', 'y_sc', 'lb0', 'ub0', 'lb1', 'ub1', 'lb2', 'ub2', 'lb3', 'ub3', 'Tdate']
		Macro.setLocal("e_method", e_method)

	data_points.to_stata("__scpi__stata_plot.dta", write_index = False)

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