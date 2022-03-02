*! Date        : 1 Mar 2022
*! Version     : 0.2.1
*! Authors     : Filippo Palomba
*! Email       : fpalomba@princeton.edu
*! Description : Synthetic control estimation

/*
To do:
- implement V
*/

capture program drop scest
program define scest, eclass         
version 17.0           

	syntax , dfname(string) [name(string) p(integer 1) direc(string) q(real -11.92) lb(string) opt(string) pypinocheck]
	
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
	
	if mi("`opt'") {
		local opt "None"
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
	
	if `q' == -11.92 {
		local q "None"
	}
	
	python: scest_wrapper("`p_str'", "`direc'", `q', "`lb'", "`name'", "`opt'", "`dfname'")
	
	ereturn local name              = "`name'"
	ereturn local dir               = "`dire'"
	ereturn local p                 = "`p_str'"
	ereturn local cointegrated_data = "`cointegrated_data'"
	ereturn local constant          = "`constant'"
	ereturn local outcomevar        = "`outcomevar'"
	ereturn local features          = "`features'"
	
	capture ereturn scalar q  = scalar(q_est)
	ereturn scalar M  = scalar(M)
	ereturn scalar KM = scalar(KM)
	ereturn scalar J  = scalar(J)
	ereturn scalar T1 = scalar(T1)
	
	ereturn matrix Y_pre_fit  = Y_pre_fit
	ereturn matrix Y_pre      = Y_pre
	ereturn matrix Y_post_fit = Y_post_fit
	ereturn matrix Y_post     = Y_post
	ereturn matrix beta       = bmat
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
from scpi_pkg.scest import scest
from sfi import Scalar, Matrix, Macro

def scest_wrapper(p, dir, Q, lb, name, opt, dfname):
	
	filename = dfname + '.obj'
	filehandler = open(filename, 'rb') 
	df = pickle.load(filehandler)
	
	if dir == "None":
		dire = None
	else:
		dire = dir
		
	if lb == "0":
		lb = 0
	else:
		lb = -numpy.inf
		
	if name == "None":
		w_constr = {'p': p, 'dir': dire, 'Q': Q, 'lb': lb}
	else:
		if Q is None:
			w_constr = {'name': str(name)}
		else:
			w_constr = {'name': str(name), 'Q': Q}
				
	if opt == "None":
		opt = {}

	res_est = scest(df, w_constr = w_constr, V = None, opt_dict = opt)
	print(res_est)
	
	filename = '__scest__output.obj'
	file     = open(filename, 'wb')
	pickle.dump(res_est, file)	
	
	Scalar.setValue("T1",    res_est.T1_outcome)
	Scalar.setValue("J",     res_est.J)
	Scalar.setValue("KM",    res_est.KM)
	Scalar.setValue("M",     res_est.M)
	if res_est.w_constr['Q'] is not None:
		Scalar.setValue("q_est", res_est.w_constr['Q'])
	
	Matrix.create("wmat", res_est.J, 1, 0)
	Matrix.store("wmat", res_est.w.values)
	names = [str(row) for row in res_est.w.index.tolist()]
	Matrix.setRowNames("wmat", names)
	
	if res_est.KM > 0:
		Matrix.create("rmat", res_est.KM, 1, 0)
		Matrix.store("rmat", res_est.r.values)
		names = [str(row) for row in res_est.r.index.tolist()]
		Matrix.setRowNames("rmat", names)
	
	Matrix.create("bmat", len(res_est.b), 1, 0)
	Matrix.store("bmat", res_est.b.values)
	names = [str(row) for row in res_est.b.index.tolist()]
	Matrix.setRowNames("bmat", names)
	
	Matrix.create("Y_post", len(res_est.Y_post), 1, 0)
	Matrix.store("Y_post", res_est.Y_post.values)
	names = [str(row) for row in res_est.Y_post.index.tolist()]
	Matrix.setRowNames("Y_post", names)

	Matrix.create("Y_post_fit", len(res_est.Y_post_fit), 1, 0)
	Matrix.store("Y_post_fit", res_est.Y_post_fit.values)
	names = [str(row) for row in res_est.Y_post_fit.index.tolist()]
	Matrix.setRowNames("Y_post_fit", names)

	Matrix.create("Y_pre", len(res_est.Y_pre), 1, 0)
	Matrix.store("Y_pre", res_est.Y_pre.values)
	names = [str(row) for row in res_est.Y_pre.index.tolist()]
	Matrix.setRowNames("Y_pre", names)
	
	Matrix.create("Y_pre_fit", len(res_est.Y_pre_fit), 1, 0)
	Matrix.store("Y_pre_fit", res_est.Y_pre_fit.values)
	names = [str(row) for row in res_est.Y_pre_fit.index.tolist()]
	Matrix.setRowNames("Y_pre_fit", names)

	Matrix.create("res", len(res_est.res), 1, 0)
	Matrix.store("res", res_est.res.values)
	names = [str(row) for row in res_est.res.index.tolist()]
	Matrix.setRowNames("res", names)

	Matrix.create("pred", len(res_est.A_hat), 1, 0)
	Matrix.store("pred", res_est.A_hat.values)
	names = [str(row) for row in res_est.A_hat.index.tolist()]
	Matrix.setRowNames("pred", names)

	Matrix.create("A", len(res_est.A), 1, 0)
	Matrix.store("A", res_est.A.values)
	names = [str(row) for row in res_est.A.index.tolist()]
	Matrix.setRowNames("A", names)
	
	Matrix.create("B", len(res_est.B), res_est.J, 0)
	Matrix.store("B", res_est.B.values)
	names = [str(row) for row in res_est.B.index.tolist()]
	Matrix.setRowNames("B", names)	
	names = [str(col) for col in res_est.B.columns.tolist()]
	Matrix.setColNames("B", names)	

	if res_est.KM > 0:
		Matrix.create("C", len(res_est.C), res_est.KM, 0)
		Matrix.store("C", res_est.C.values)
		names = [str(row) for row in res_est.C.index.tolist()]
		Matrix.setRowNames("C", names)	
		names = [str(col) for col in res_est.C.columns.tolist()]
		Matrix.setColNames("C", names)	

	Matrix.create("T0_features", len(res_est.features), 1, 0)
	Matrix.store("T0_features", numpy.array([c for c in res_est.T0_features.values()]))
	Matrix.setRowNames("T0_features", res_est.features)		
		
	Macro.setLocal("cointegrated_data", str(res_est.cointegrated_data))
	Macro.setLocal("constant", str(res_est.glob_cons))
	Macro.setLocal("outcomevar", res_est.outcome_var)
	Macro.setLocal("features", ' '.join([str(elem) for elem in res_est.features]))
	Macro.setLocal("dire", str(res_est.w_constr['dir']))
	
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

























