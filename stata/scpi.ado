*! Date        : 20 Jan 2023
*! Version     : 2.1.1
*! Authors     : Filippo Palomba
*! Email       : fpalomba@princeton.edu
*! Description : Synthetic control inference

capture program drop scpi
program define scpi, eclass         
version 16.0           

	syntax , dfname(string) [p(integer 1) direc(string) q(real -11.92) lb(string) V(string) name(string) u_missp u_sigma(string) u_order(integer 1) u_lags(integer 0) u_alpha(real 0.05) sims(integer 200) ///
			 e_method(string) e_order(integer 1) e_lags(integer 0) e_alpha(real 0.05) lgapp(string) rho(real -11) rho_max(real -11) cores(integer 1) pypinocheck]

	if mi("`pypinocheck'") & mi("$scpi_version_checked") {
		python: version_checker()
		if "`alert_version'" == "y" {
			di as error "The current version of scpi_pkg in Python is `python_local_version', but version `python_pypi_version' needed! Please update the package in Python and restart Stata!"
			exit 198
		}
		global scpi_version_checked "yes"
	}			 

	if mi("`name'") & mi("`p'") & mi("`direc'") & mi("`q'") & mi("`lb'") {
		di as error "One between 'name' and the set of options 'p', 'direc', 'lb', and 'q' has to be specified!"
		exit 198
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

	if mi("`V'") {
		local V = "separate"
	} 
	else {
		if !inlist("`V'", "separate", "pooled") {
			di as error "The option V should be either 'separate' or 'pooled'!"
			exit 198
		}
	}	

	if mi("`lgapp'") {
		local lgapp = "generalized"
	} 
	else {
		if !inlist("`lgapp'", "generalized", "linear") {
			di as error "The option lgapp should be either 'generalized' or 'linear'!"
			exit 198
		}
	}		
	
	if mi("`u_missp'") {
		local u_missp "False"
	}
	
	if mi("`u_sigma'") {
		local u_sigma "HC1"
	}
	
	if mi("`e_method'") {
		local e_method "gaussian"
	}
	
	python: executionTime(`cores', `sims', "`dfname'")
	
	di ""
	di "Linking to Python - Results will be printed soon!"
	di "Estimating weights and Quantifying Uncertainty"
	di "`toprint'"

	
	sleep 500

	python: scpi_wrapper("`p_str'", "`direc'", `q', "`lb'", "`name'", "`V'", "`u_missp'", "`u_sigma'", `u_order', `u_lags', `u_alpha', "`e_method'", `e_order', `e_lags', ///
						 `e_alpha', `sims', "`lgapp'", `rho', `rho_max', `cores', "`dfname'")
	
	ereturn clear
	
	ereturn local lgapp             = "`lgapp'"
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
	ereturn local donors            = "`donors'"
	ereturn local anticipation      = "`anticipation'"
	ereturn local constant          = "`constant'"
	ereturn local features          = "`features'"
	ereturn local outcomevar        = "`outcomevar'"
	
	ereturn scalar KMI = scalar(KMI)
	ereturn scalar I = 1
	
	ereturn matrix failed_sims     = failed_sims
	capture ereturn matrix e_var = e_var
	capture ereturn matrix e_mean = e_mean
	capture ereturn matrix u_var = u_var
	capture ereturn matrix u_mean = u_mean
	capture ereturn matrix CI_all_qreg     = CI_all_qreg
	capture ereturn matrix CI_all_ls       = CI_all_ls
	capture ereturn matrix CI_all_gaussian = CI_all_gaussian
	ereturn matrix CI_in_sample    = CI_in_sample

	ereturn matrix rho = rho
	capture ereturn matrix Qstar = Qstar
	ereturn matrix M = M
	ereturn matrix KM = KM
	ereturn matrix J = J
	ereturn matrix T1 = T1
	ereturn matrix T0 = T0_features
	ereturn matrix Y_post_fit = Y_post_fit
	ereturn matrix Y_pre_fit = Y_pre_fit
	
	if "`class_type'" == "scpi_output" {
		ereturn matrix Y_post = Y_post		
		ereturn matrix Y_pre = Y_pre		
	}

	ereturn matrix beta = bmat
	capture ereturn matrix r = rmat
	ereturn matrix w = wmat
	ereturn matrix res = res
	ereturn matrix pred = pred
	ereturn matrix P = P
	capture ereturn matrix C = C
	ereturn matrix B = B
	ereturn matrix A = A	
	
end


version 16.0
python:
import pickle, numpy, urllib, luddite
from scpi_pkg.scpi import scpi
from scpi_pkg import version as lver
from sfi import Scalar, Matrix, Macro
from math import ceil, floor


def scpi_wrapper(p, dir, q, lb, name, V, u_missp, u_sigma, u_order, u_lags, u_alpha, e_method, e_order, e_lags, e_alpha, sims, lgapp, rho, rho_max, cores, dfname):

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
				
	if rho == -11:
		rho = None
	if rho_max == -11:
		rho_max = None
		
	# Parse u_missp
	if u_missp == "False":
		u_missp_bool = False
	else:
		u_missp_bool = True

	res_pi = scpi(df, w_constr, V, None, u_missp_bool, u_sigma, u_order, u_lags, None, u_alpha, str(e_method), e_order, e_lags, None, e_alpha, sims, rho, rho_max, lgapp, cores, False, None, None, True, True)
	class_type = res_pi.__class__.__name__
	Macro.setLocal("class_type", class_type)
	
	if class_type == "scpi_output":
		res_pi.iota = 1
	
	filename = '__scpi__output.obj'
	file     = open(filename, 'wb')
	pickle.dump(res_pi, file)	
	
	if res_pi.KMI > 0:
		Matrix.create("rmat", len(res_pi.r), 1, 0)
		Matrix.store("rmat", res_pi.r.values)
		names = [ix2rn(row) for row in res_pi.r.index.tolist()]
		Matrix.setRowNames("rmat", names)
	
	Matrix.create("bmat", len(res_pi.b), 1, 0)
	Matrix.store("bmat", res_pi.b.values)
	names = [ix2rn(row) for row in res_pi.b.index.tolist()]
	Matrix.setRowNames("bmat", names)
	
	if class_type == "scpi_output":
		Matrix.create("Y_post", len(res_pi.Y_post), 1, 0)
		Matrix.store("Y_post", res_pi.Y_post.values)
		names = [ix2rn(row) for row in res_pi.Y_post.index.tolist()]
		Matrix.setRowNames("Y_post", names)

		Matrix.create("Y_pre", len(res_pi.Y_pre), 1, 0)
		Matrix.store("Y_pre", res_pi.Y_pre.values)
		names = [ix2rn(row) for row in res_pi.Y_pre.index.tolist()]
		Matrix.setRowNames("Y_pre", names)
		
	Matrix.create("Y_post_fit", len(res_pi.Y_post_fit), 1, 0)
	Matrix.store("Y_post_fit", res_pi.Y_post_fit.values)
	names = [ix2rn(row) for row in res_pi.Y_post_fit.index.tolist()]
	Matrix.setRowNames("Y_post_fit", names)
	
	Matrix.create("Y_pre_fit", len(res_pi.Y_pre_fit), 1, 0)
	Matrix.store("Y_pre_fit", res_pi.Y_pre_fit.values)
	names = [ix2rn(row) for row in res_pi.Y_pre_fit.index.tolist()]
	Matrix.setRowNames("Y_pre_fit", names)

	Matrix.create("res", len(res_pi.res), 1, 0)
	Matrix.store("res", res_pi.res.values)
	names = [ix2rn(row) for row in res_pi.res.index.tolist()]
	Matrix.setRowNames("res", names)

	Matrix.create("pred", len(res_pi.A_hat), 1, 0)
	Matrix.store("pred", res_pi.A_hat.values)
	names = [ix2rn(row) for row in res_pi.A_hat.index.tolist()]
	Matrix.setRowNames("pred", names)

	Matrix.create("A", len(res_pi.A), 1, 0)
	Matrix.store("A", res_pi.A.values)
	names = [ix2rn(row) for row in res_pi.A.index.tolist()]
	Matrix.setRowNames("A", names)
	
	Matrix.create("B", len(res_pi.B), len(res_pi.B.columns), 0)
	Matrix.store("B", res_pi.B.values)
	names = [ix2rn(row) for row in res_pi.B.index.tolist()]
	Matrix.setRowNames("B", names)	
	names = [ix2rn(col) for col in res_pi.B.columns.tolist()]
	Matrix.setColNames("B", names)	

	if res_pi.KMI > 0:
		Matrix.create("C", len(res_pi.C), len(res_pi.C.columns), 0)
		Matrix.store("C", res_pi.C.values)
		names = [ix2rn(row) for row in res_pi.C.index.tolist()]
		Matrix.setRowNames("C", names)	
		names = [ix2rn(col) for col in res_pi.C.columns.tolist()]
		Matrix.setColNames("C", names)	

	Matrix.create("P", len(res_pi.P), len(res_pi.P.columns), 0)
	Matrix.store("P", res_pi.P.values)
	names = [ix2rn(row) for row in res_pi.P.index.tolist()]
	Matrix.setRowNames("P", names)	
	names = [ix2rn(col) for col in res_pi.P.columns.tolist()]
	Matrix.setColNames("P", names)				

	wc_est = res_pi.w_constr[res_pi.treated_units[0]]	
	names = res_pi.treated_units
	if wc_est['Q'] is not None:
		Matrix.create("Qstar", res_pi.iota, 1, 0)
		Matrix.store("Qstar", [w['Q'] for w in res_pi.w_constr.values()])
		Matrix.setRowNames("Qstar", names)

	Matrix.create("rho", res_pi.iota, 1, 0)
	Matrix.store("rho", [r for r in res_pi.rho.values()])
	Matrix.setRowNames("rho", names)	
	
	Matrix.create("CI_in_sample", len(res_pi.CI_in_sample), len(res_pi.CI_in_sample.columns), 0)
	Matrix.store("CI_in_sample", res_pi.CI_in_sample.values)
	names = [ix2rn(row) for row in res_pi.CI_in_sample.index.tolist()]
	Matrix.setRowNames("CI_in_sample", names)	
	names = [ix2rn(col) for col in res_pi.CI_in_sample.columns.tolist()]
	Matrix.setColNames("CI_in_sample", names)	
	
	if e_method == "gaussian" or e_method == "all":
		Matrix.create("CI_all_gaussian", len(res_pi.CI_all_gaussian), len(res_pi.CI_all_gaussian.columns), 0)
		Matrix.store("CI_all_gaussian", res_pi.CI_all_gaussian.values)
		names = [ix2rn(col) for col in res_pi.CI_all_gaussian.columns.tolist()]
		Matrix.setColNames("CI_all_gaussian", names)
		names = [ix2rn(row) for row in res_pi.CI_all_gaussian.index.tolist()]
		Matrix.setRowNames("CI_all_gaussian", names)	
		
		Matrix.create("e_mean", len(res_pi.CI_all_gaussian), 1, 0)
		Matrix.store("e_mean", res_pi.e_mean)
		Matrix.setRowNames("e_mean", names)
		
		Matrix.create("e_var", len(res_pi.CI_all_gaussian), 1, 0)
		Matrix.store("e_var", res_pi.e_var)
		Matrix.setRowNames("e_var", names)
		
	if e_method == "ls" or e_method == "all":
		Matrix.create("CI_all_ls", len(res_pi.CI_all_ls), len(res_pi.CI_all_ls.columns), 0)
		Matrix.store("CI_all_ls", res_pi.CI_all_ls.values)
		names = [ix2rn(col) for col in res_pi.CI_all_ls.columns.tolist()]
		Matrix.setColNames("CI_all_ls", names)	
		names = [ix2rn(row) for row in res_pi.CI_all_ls.index.tolist()]
		Matrix.setRowNames("CI_all_ls", names)	

		Matrix.create("e_mean", len(res_pi.CI_all_gaussian), 1, 0)
		Matrix.store("e_mean", res_pi.e_mean)
		Matrix.setRowNames("e_mean", names)

		Matrix.create("e_var", len(res_pi.CI_all_gaussian), 1, 0)
		Matrix.store("e_var", res_pi.e_var)
		Matrix.setRowNames("e_var", names)
		
	if e_method == "qreg" or e_method == "all":
		Matrix.create("CI_all_qreg", len(res_pi.CI_all_qreg), len(res_pi.CI_all_qreg.columns), 0)
		Matrix.store("CI_all_qreg", res_pi.CI_all_qreg.values)
		names = [ix2rn(row) for row in res_pi.CI_all_qreg.index.tolist()]
		Matrix.setRowNames("CI_all_qreg", names)	
		names = [ix2rn(col) for col in res_pi.CI_all_qreg.columns.tolist()]
		Matrix.setColNames("CI_all_qreg", names)	

	if u_missp_bool:
		Matrix.create("u_mean", res_pi.u_mean.shape[0], res_pi.u_mean.shape[1], 0)
		Matrix.store("u_mean", res_pi.u_mean)

	Matrix.create("u_var", res_pi.u_var.shape[0], res_pi.u_var.shape[1], 0)
	Matrix.store("u_var", res_pi.u_var)
	Matrix.create("failed_sims",res_pi.failed_sims.shape[0], res_pi.failed_sims.shape[1], 0)
	Matrix.store("failed_sims", res_pi.failed_sims.values)
	Matrix.setRowNames("failed_sims", ['lb','ub'])	

	names = res_pi.treated_units
	
	if class_type == "scpi_output":
		Matrix.create("wmat", len(res_pi.w), 1, 0)
		Matrix.store("wmat", res_pi.w.values)
		nms = [str(row) for row in res_pi.w.index.tolist()]
		Matrix.setRowNames("wmat", nms)
	
		Matrix.create("T1", 1, 1, 0)
		Matrix.store("T1", res_pi.T1_outcome)
		Matrix.setRowNames("T1", names)

		Matrix.create("J", 1, 1, 0)
		Matrix.store("J", res_pi.J)
		Matrix.setRowNames("J", names)

		Matrix.create("KM", 1, 1, 0)
		Matrix.store("KM", res_pi.KM)
		Matrix.setRowNames("KM", names)

		Matrix.create("M", 1, 1, 0)
		Matrix.store("M", res_pi.M)
		Matrix.setRowNames("M", names)
		
		Matrix.create("T0_features", len(res_pi.features), 1, 0)
		Matrix.store("T0_features", numpy.array([c for c in res_pi.T0_features.values()]))
		Matrix.setRowNames("T0_features", res_pi.features)		
			
		Macro.setLocal("cointegrated_data", str(res_pi.cointegrated_data))
		Macro.setLocal("constant", str(res_pi.glob_cons))
		Macro.setLocal("features", ' '.join([str(elem) for elem in res_pi.features]))
		Macro.setLocal("anticipation", str(res_pi.anticipation))

	
	if class_type == "scpi_multi_output":
		Matrix.create("wmat", len(res_pi.w), 1, 0)
		Matrix.store("wmat", res_pi.w.values)
		names = [ix2rn(row) for row in res_pi.w.index.tolist()]
		Matrix.setRowNames("wmat", names)
	
		Matrix.create("T1", len(res_pi.T1_outcome), 1, 0)
		names = [str(k) for k in res_pi.T1_outcome.keys()]
		Matrix.store("T1", [v for v in res_pi.T1_outcome.values()])
		Matrix.setRowNames("T1", names)

		Matrix.create("J", len(res_pi.J), 1, 0)
		names = [str(k) for k in res_pi.J.keys()]
		Matrix.store("J", [v for v in res_pi.J.values()])
		Matrix.setRowNames("J", names)

		Matrix.create("KM", len(res_pi.KM), 1, 0)
		names = [str(k) for k in res_pi.KM.keys()]
		Matrix.store("KM", [v for v in res_pi.KM.values()])
		Matrix.setRowNames("KM", names)

		Matrix.create("M", len(res_pi.M), 1, 0)
		names = [str(k) for k in res_pi.M.keys()]
		Matrix.store("M", [v for v in res_pi.M.values()])
		Matrix.setRowNames("M", names)

		unique_feats = []
		for v in res_pi.T0_features.values():
			for k in v.keys():
				if not k in unique_feats:
					unique_feats.append(k)
		aux = numpy.zeros((res_pi.iota, len(unique_feats)))
		for i in range(res_pi.iota):
			for j in range(len(unique_feats)):
				aux[i,j] = res_pi.T0_features[res_pi.treated_units[i]][unique_feats[j]]

		Matrix.create("T0_features", res_pi.iota, len(unique_feats), 0)
		Matrix.store("T0_features", aux)
		Matrix.setRowNames("T0_features", res_pi.treated_units)				
		Matrix.setColNames("T0_features", unique_feats)				
		
		Macro.setLocal("cointegrated_data", dict2glob(res_pi.cointegrated_data.keys(), res_pi.cointegrated_data.values()))
		Macro.setLocal("constant", dict2glob(res_pi.glob_cons.keys(), res_pi.glob_cons.values()))
		Macro.setLocal("features", dict2glob(res_pi.features.keys(), res_pi.features.values()))
		Macro.setLocal("anticipation", dict2glob(res_pi.anticipation.keys(), res_pi.anticipation.values()))
		
	Macro.setLocal("u_alpha", str(res_pi.u_alpha))
	Macro.setLocal("u_order", str(res_pi.u_order))
	Macro.setLocal("u_lags", str(res_pi.u_lags))
	Macro.setLocal("u_missp", str(res_pi.u_missp))
	Macro.setLocal("u_sigma", res_pi.u_sigma)
	Macro.setLocal("e_alpha", str(res_pi.e_alpha))
	Macro.setLocal("e_order", str(res_pi.e_order))
	Macro.setLocal("e_lags", str(res_pi.e_lags))
	Macro.setLocal("e_method", res_pi.e_method)
	Macro.setLocal("dire", str(wc_est['dir']))
	Macro.setLocal("outcomevar", res_pi.outcome_var)
	Macro.setLocal("donors", dict2glob(res_pi.donors_dict.keys(), res_pi.donors_dict.values()))
	Scalar.setValue("KMI", res_pi.KMI)
	Scalar.setValue("I", res_pi.iota)
	
	
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

def ix2rn(s):
	return str(s).replace('(','').replace(')','').replace("'",'').replace(", ","_")

def dict2glob(keys, vals):
    klist = [k for k in keys]
    vlist = [v for v in vals]
    for i in range(len(vlist)):
        if isinstance(vlist[i], (list)):
            vlist[i] = ' '.join([str(v) for v in vlist[i]])
        if i==0:
            dct = str(klist[i]) + ": " + str(vlist[i]) 
        else:
            dct = dct + " || " + str(klist[i]) + ": " + str(vlist[i])
    return dct
	
def executionTime(cores, sims, dfname):

	filename = dfname + '.obj'
	filehandler = open(filename, 'rb') 
	df = pickle.load(filehandler)
	class_type = df.__class__.__name__

	if class_type == "scdata_output":
		J = df.J
		Ttot = sum(df.T0_features.values())
		T1 = df.T1_outcome
		I = 1
	else:
		T1 = sum(df.T1_outcome.values())
		J = sum(df.J.values())
		T0_M = {}
		for n, v in df.T0_features.items():
			T0_M[n] = sum(v.values())
		Ttot = sum(T0_M.values())
		I = df.iota

	tincr = Ttot / 1000

	coefsJ = numpy.array([-0.54755616, 0.09985644])

	time = numpy.array([1, J]) @ coefsJ
	time = ceil(time) * sims / 10
	time = time / cores
	time = time * tincr
	time = time * T1
	time = time * I
	time = time / 60
	time = ceil(time)
	time = 2 * time

	if time < 60:
		if time < 1:
			toprint = "Maximum expected execution time: less than a minute."
		elif time == 1:
			toprint = "Maximum expected execution time: " + str(time) + " minute."
		else:
			toprint = "Maximum expected execution time: " + str(time) + " minutes."
	else:
		hours = floor(time/60)
		if hours == 1: 
			toprint = "Maximum expected execution time: " + str(hours) + " hour."
		else:
			toprint = "Maximum expected execution time: " + str(hours) + " hours."
			
	Macro.setLocal("toprint", toprint)
	
end
