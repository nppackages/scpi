*! Date        : 1 Mar 2022
*! Version     : 0.2.1
*! Authors     : Filippo Palomba
*! Email       : fpalomba@princeton.edu
*! Description : Synthetic control inference

/*
To do:
- implement V
- implement P
*/

capture program drop scpi
program define scpi, eclass         
version 17.0           

	syntax , dfname(string) [p(integer 1) direc(string) q(real -11.92) lb(string) name(string) u_missp u_sigma(string) u_order(integer 1) u_lags(integer 0) u_alpha(real 0.05) sims(integer 200) ///
			 e_method(string) e_order(integer 1) e_lags(integer 0) e_alpha(real 0.05) rho(real -11) rho_max(real -11) cores(integer 1) opt_est(string) opt_inf(string) pypinocheck]

	if mi("`pypinocheck'") & mi("$scpi_version_checked") {
		python: version_checker()
		if "`alert_version'" == "y" {
			di as error "The current version of scpi_pkg in Python is `python_local_version', but version `python_pypi_version' needed! Please update the package in Python and restart Stata!"
			exit 198
		}
		global scpi_version_checked "yes"
	}			 
			 
	if mi("`name'") {
		local name "None"
	}
	
	if mi("`direc'") {
		if `p' > 0 {
			local direc "=="
		} 
		if `p' == 0 {
			local direc "None"
		}
	}
	
	if mi("`opt_est'") {
		local opt_est "None"
	}

	if mi("`opt_inf'") {
		local opt_inf "None"
	}	
	
	if !mi("`p'") {
		if `p' == 0 {
			local p_str "no norm"
		} 
		if `p' == 1 {
			local p_str "L1"
		}
		if `p' == 2 {
			local p_str "L2"
		}
	}
	
	if mi("`lb'") {
		local lb "0"
	} 
	
	if "`lb'" != "0" & "`lb'" != "-inf" {
		di as error "The option lb should be either 'zero' or '-inf'."
		exit 198
	}
	
	if mi("`u_missp'") {
		local u_missp "False"
	}
	
	if mi("`u_sigma'") {
		local u_sigma "HC1"
	}
	
	if mi("`e_method'") {
		local e_method "all"
	}
	
	python: executionTime(`cores', `sims', "`dfname'")
	
	di ""
	di "Linking to Python - Results will be printed soon!"
	di "Estimating weights and Quantifying Uncertainty"
	di "`toprint'"

	
	sleep 500

	python: scpi_wrapper("`p_str'", "`direc'", `q', "`lb'", "`name'", "`u_missp'", "`u_sigma'", `u_order', `u_lags', `u_alpha', "`e_method'", `e_order', `e_lags', ///
						 `e_alpha', `sims', `rho', `rho_max', `cores', "`opt_est'", "`opt_inf'", "`dfname'")
	
	
	ereturn local e_alpha           = "`e_alpha'"
	ereturn local e_lags            = "`e_lags'"
	ereturn local e_order           = "`e_order'"
	ereturn local e_sigma           = "`e_sigma'"
	ereturn local e_method          = "`e_method'"
	ereturn local u_alpha           = "`u_alpha'"
	ereturn local u_sigma           = "`u_sigma'"
	ereturn local u_lags            = "`u_lags'"
	ereturn local u_order           = "`u_order'"
	ereturn local u_missp           = "`u_missp'"
	ereturn local name              = "`name'"
	ereturn local dir               = "`dire'"
	ereturn local p                 = "`p_str'"
	ereturn local cointegrated_data = "`cointegrated_data'"
	ereturn local constant          = "`constant'"
	ereturn local outcomevar        = "`outcomevar'"
	ereturn local features          = "`features'"
	
	ereturn local rho = scalar(rho)
	capture ereturn local q   = scalar(q_est)  
	ereturn scalar M  = scalar(M)
	ereturn scalar KM = scalar(KM)
	ereturn scalar J  = scalar(J)
	ereturn scalar T1 = scalar(T1)
	
	ereturn matrix failed_sims     = failed_sims
	capture ereturn matrix e_var = e_var
	capture ereturn matrix e_mean = e_mean
	capture ereturn matrix u_var = u_var
	capture ereturn matrix u_mean = u_mean
	capture ereturn matrix CI_all_qreg     = CI_all_qreg
	capture ereturn matrix CI_all_ls       = CI_all_ls
	capture ereturn matrix CI_all_gaussian = CI_all_gaussian
	ereturn matrix CI_in_sample    = CI_in_sample
	ereturn matrix Y_pre_fit       = Y_pre_fit
	ereturn matrix Y_pre           = Y_pre
	ereturn matrix Y_post_fit      = Y_post_fit
	ereturn matrix Y_post          = Y_post
	ereturn matrix beta            = bmat
	capture ereturn matrix r = rmat
	ereturn matrix w    = wmat
	ereturn matrix res  = res
	ereturn matrix pred = pred
	capture ereturn matrix C = C
	ereturn matrix B  = B
	ereturn matrix A  = A	
	ereturn matrix T0 = T0_features
	
	
end


version 17.0
python:
import pickle, numpy, urllib, luddite
from scpi_pkg.scpi import scpi
from sfi import Scalar, Matrix, Macro
from math import ceil


def scpi_wrapper(p, dir, q, lb, name, u_missp, u_sigma, u_order, u_lags, u_alpha, e_method, e_order, e_lags, e_alpha, sims, rho, rho_max, cores, opt_est, opt_inf, dfname):

	filename = dfname + '.obj'
	filehandler = open(filename, 'rb') 
	df = pickle.load(filehandler)
	
	if lb == "0":
		lb = 0
	else:
		lb = -numpy.inf

	if dir == "None":
		dire = None
	else:
		dire = dir		
		
	if q == -11.92:
		Q = None
	else:
		Q = q
	
	if name == "None":
		if Q is None:
			w_constr = {'p': p, 'dir': dir, 'Q': 1, 'lb': lb}
		else:
			w_constr = {'p': p, 'dir': dir, 'Q': Q, 'lb': lb}
	else:
		if Q is None:
			w_constr = {'name': str(name)}
		else:
			w_constr = {'name': str(name), 'Q': Q}
				
	if opt_est == "None":
		opt_est = {}
		
	if opt_inf == "None":
		opt_inf = {}
		
	if rho == -11:
		rho = None
	if rho_max == -11:
		rho_max = None
		
	# Parse u_missp
	if u_missp == "False":
		u_missp_bool = False
	else:
		u_missp_bool = True
		
	res_pi = scpi(df, w_constr, None, None, u_missp_bool, u_sigma, u_order, u_lags, None, u_alpha, str(e_method), e_order, e_lags, None, e_alpha, sims, rho, rho_max, cores, False, None, None, None, opt_est, opt_inf)
	print(res_pi)
	
	filename = '__scpi__output.obj'
	file     = open(filename, 'wb')
	pickle.dump(res_pi, file)
	
	
	Scalar.setValue("T1",    res_pi.T1_outcome)
	Scalar.setValue("J",     res_pi.J)
	Scalar.setValue("KM",    res_pi.KM)
	Scalar.setValue("M",     res_pi.M)
	Scalar.setValue("rho",   res_pi.rho)
	if res_pi.w_constr['Q'] is not None:
		Scalar.setValue("q_est", res_pi.w_constr['Q'])
	
	Matrix.create("wmat", res_pi.J, 1, 0)
	Matrix.store("wmat", res_pi.w.values)
	names = [str(row) for row in res_pi.w.index.tolist()]
	Matrix.setRowNames("wmat", names)
	
	if res_pi.KM > 0:
		Matrix.create("rmat", res_pi.KM, 1, 0)
		Matrix.store("rmat", res_pi.r.values)
		names = [str(row) for row in res_pi.r.index.tolist()]
		Matrix.setRowNames("rmat", names)
	
	Matrix.create("bmat", len(res_pi.b), 1, 0)
	Matrix.store("bmat", res_pi.b.values)
	names = [str(row) for row in res_pi.b.index.tolist()]
	Matrix.setRowNames("bmat", names)
	
	Matrix.create("Y_post", len(res_pi.Y_post), 1, 0)
	Matrix.store("Y_post", res_pi.Y_post.values)
	names = [str(row) for row in res_pi.Y_post.index.tolist()]
	Matrix.setRowNames("Y_post", names)

	Matrix.create("Y_post_fit", len(res_pi.Y_post_fit), 1, 0)
	Matrix.store("Y_post_fit", res_pi.Y_post_fit.values)
	names = [str(row) for row in res_pi.Y_post_fit.index.tolist()]
	Matrix.setRowNames("Y_post_fit", names)

	Matrix.create("Y_pre", len(res_pi.Y_pre), 1, 0)
	Matrix.store("Y_pre", res_pi.Y_pre.values)
	names = [str(row) for row in res_pi.Y_pre.index.tolist()]
	Matrix.setRowNames("Y_pre", names)
	
	Matrix.create("Y_pre_fit", len(res_pi.Y_pre_fit), 1, 0)
	Matrix.store("Y_pre_fit", res_pi.Y_pre_fit.values)
	names = [str(row) for row in res_pi.Y_pre_fit.index.tolist()]
	Matrix.setRowNames("Y_pre_fit", names)

	Matrix.create("res", len(res_pi.res), 1, 0)
	Matrix.store("res", res_pi.res.values)
	names = [str(row) for row in res_pi.res.index.tolist()]
	Matrix.setRowNames("res", names)

	Matrix.create("pred", len(res_pi.A_hat), 1, 0)
	Matrix.store("pred", res_pi.A_hat.values)
	names = [str(row) for row in res_pi.A_hat.index.tolist()]
	Matrix.setRowNames("pred", names)

	Matrix.create("A", len(res_pi.A), 1, 0)
	Matrix.store("A", res_pi.A.values)
	names = [str(row) for row in res_pi.A.index.tolist()]
	Matrix.setRowNames("A", names)
	
	Matrix.create("B", len(res_pi.B), res_pi.J, 0)
	Matrix.store("B", res_pi.B.values)
	names = [str(row) for row in res_pi.B.index.tolist()]
	Matrix.setRowNames("B", names)	
	names = [str(col) for col in res_pi.B.columns.tolist()]
	Matrix.setColNames("B", names)	

	if res_pi.KM > 0:
		Matrix.create("C", len(res_pi.C), res_pi.KM, 0)
		Matrix.store("C", res_pi.C.values)
		names = [str(row) for row in res_pi.C.index.tolist()]
		Matrix.setRowNames("C", names)	
		names = [str(col) for col in res_pi.C.columns.tolist()]
		Matrix.setColNames("C", names)	
	
	Matrix.create("CI_in_sample", len(res_pi.CI_in_sample), len(res_pi.CI_in_sample.columns), 0)
	Matrix.store("CI_in_sample", res_pi.CI_in_sample.values)
	names = [str(row) for row in res_pi.CI_in_sample.index.tolist()]
	Matrix.setRowNames("CI_in_sample", names)	
	names = [str(col) for col in res_pi.CI_in_sample.columns.tolist()]
	Matrix.setColNames("CI_in_sample", names)	
	
	if e_method == "gaussian" or e_method == "all":
		Matrix.create("CI_all_gaussian", len(res_pi.CI_all_gaussian), len(res_pi.CI_all_gaussian.columns), 0)
		Matrix.store("CI_all_gaussian", res_pi.CI_all_gaussian.values)
		names = [str(row) for row in res_pi.CI_all_gaussian.index.tolist()]
		Matrix.setRowNames("CI_all_gaussian", names)	
		names = [str(col) for col in res_pi.CI_all_gaussian.columns.tolist()]
		Matrix.setColNames("CI_all_gaussian", names)	

		Matrix.create("e_mean", len(res_pi.CI_all_gaussian), 1, 0)
		Matrix.store("e_mean", res_pi.e_mean)
		Matrix.setRowNames("e_mean", names)

		Matrix.create("e_var", len(res_pi.CI_all_gaussian), 1, 0)
		Matrix.store("e_var", res_pi.e_var)
		Matrix.setRowNames("e_var", names)

	if e_method == "ls" or e_method == "all":
		Matrix.create("CI_all_ls", len(res_pi.CI_all_ls), len(res_pi.CI_all_ls.columns), 0)
		Matrix.store("CI_all_ls", res_pi.CI_all_ls.values)
		names = [str(row) for row in res_pi.CI_all_ls.index.tolist()]
		Matrix.setRowNames("CI_all_ls", names)	
		names = [str(col) for col in res_pi.CI_all_ls.columns.tolist()]
		Matrix.setColNames("CI_all_ls", names)	

		Matrix.create("e_mean", len(res_pi.CI_all_gaussian), 1, 0)
		Matrix.store("e_mean", res_pi.e_mean)
		Matrix.setRowNames("e_mean", names)

		Matrix.create("e_var", len(res_pi.CI_all_gaussian), 1, 0)
		Matrix.store("e_var", res_pi.e_var)
		Matrix.setRowNames("e_var", names)
		
	if e_method == "qreg" or e_method == "all":
		Matrix.create("CI_all_qreg", len(res_pi.CI_all_qreg), len(res_pi.CI_all_qreg.columns), 0)
		Matrix.store("CI_all_qreg", res_pi.CI_all_qreg.values)
		names = [str(row) for row in res_pi.CI_all_qreg.index.tolist()]
		Matrix.setRowNames("CI_all_qreg", names)	
		names = [str(col) for col in res_pi.CI_all_qreg.columns.tolist()]
		Matrix.setColNames("CI_all_qreg", names)	

	if u_missp_bool:
		Matrix.create("u_mean", res_pi.u_mean.shape[0], res_pi.u_mean.shape[1], 0)
		Matrix.store("u_mean", res_pi.u_mean)

	Matrix.create("u_var", res_pi.u_var.shape[0], res_pi.u_var.shape[1], 0)
	Matrix.store("u_var", res_pi.u_var)

	Matrix.create("failed_sims",res_pi.failed_sims.shape[0], res_pi.failed_sims.shape[1], 0)
	Matrix.store("failed_sims", res_pi.failed_sims.values)
	Matrix.setRowNames("failed_sims", ['lb','ub'])	

	Matrix.create("T0_features", len(res_pi.features), 1, 0)
	Matrix.store("T0_features", numpy.array([c for c in res_pi.T0_features.values()]))
	Matrix.setRowNames("T0_features", res_pi.features)
		
	Macro.setLocal("cointegrated_data", str(res_pi.cointegrated_data))
	Macro.setLocal("constant", str(res_pi.glob_cons))
	Macro.setLocal("outcomevar", res_pi.outcome_var)
	Macro.setLocal("features", ' '.join([str(elem) for elem in res_pi.features]))
	Macro.setLocal("u_alpha", str(res_pi.u_alpha))
	Macro.setLocal("u_order", str(res_pi.u_order))
	Macro.setLocal("u_lags", str(res_pi.u_lags))
	Macro.setLocal("u_missp", str(res_pi.u_missp))
	Macro.setLocal("u_sigma", res_pi.u_sigma)
	Macro.setLocal("e_alpha", str(res_pi.e_alpha))
	Macro.setLocal("e_order", str(res_pi.e_order))
	Macro.setLocal("e_lags", str(res_pi.e_lags))
	Macro.setLocal("e_method", res_pi.e_method)
	Macro.setLocal("dire", str(res_pi.w_constr['dir']))

	
	
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

	
def executionTime(cores, sims, dfname):
	filename = dfname + '.obj'
	filehandler = open(filename, 'rb') 
	df = pickle.load(filehandler)

	T0 = df.T0_features
	T1 = df.T1_outcome
	J = df.J
	Ttot = sum(T0.values())
	tincr = Ttot / 1000

	coefsJ = numpy.array([-0.54755616, 0.09985644])

	time = numpy.array([1, J]) @ coefsJ
	time = time * sims / 10
	time = time / cores
	time = time * tincr
	time = time * T1
	time = time / 60
	time = ceil(time)
	time = 2 * time

	if time < 1:
		toprint = "Maximum expected execution time: less than a minute."
	elif time == 1:
		toprint = "Maximum expected execution time: " + str(time) + " minute."
	else:
		toprint = "Maximum expected execution time: " + str(time) + " minutes."

	Macro.setLocal("toprint", toprint)
	
end