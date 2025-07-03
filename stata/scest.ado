*! Date        : 3 Jul 2025
*! Version     : 3.0.1
*! Authors     : Filippo Palomba
*! Email       : fpalomba@princeton.edu
*! Description : Synthetic control estimation

capture program drop scest
program define scest, eclass         
version 16.0           

	syntax , dfname(string) [name(string) p(integer 1) direc(string) q(real -11.92) lb(string) V(string) pypinocheck]
	
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
	
	if `q' == -11.92 {
		local q "None"
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
	
	python: scest_wrapper("`p_str'", "`direc'", `q', "`lb'", "`name'", "`V'", "`dfname'")

	ereturn clear

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
	ereturn scalar I = scalar(I)
	
	ereturn matrix Qstar  = Qstar
	ereturn matrix M  = M
	ereturn matrix KM = KM
	ereturn matrix J  = J
	ereturn matrix T1 = T1
	ereturn matrix T0 = T0_features
	ereturn matrix Y_post_fit = Y_post_fit
	ereturn matrix Y_pre_fit  = Y_pre_fit

	if "`class_type'" == "scest_output" {
		ereturn matrix Y_post      = Y_post		
		ereturn matrix Y_pre      = Y_pre		
	}
	
	ereturn matrix beta = bmat
	capture ereturn matrix r = rmat
	ereturn matrix w    = wmat
	ereturn matrix res  = res
	ereturn matrix pred = pred
	ereturn matrix P  = P
	capture ereturn matrix C = C
	ereturn matrix B = B
	ereturn matrix A = A	
	
end



version 16.0
python:

import pickle, numpy, urllib, luddite
from collections import Counter
from scpi_pkg.scest import scest
from scpi_pkg import version as lver
from sfi import Scalar, Matrix, Macro

def scest_wrapper(p, dir, QQ, lb, name, V, dfname):
	
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
		w_constr = {'p': p, 'dir': dire, 'Q': QQ, 'lb': lb}
	else:
		if QQ is None:
			w_constr = {'name': str(name)}
		else:
			w_constr = {'name': str(name), 'Q': QQ}

	res_est = scest(df, w_constr = w_constr, Vmat = None, V = V)

	class_type = res_est.__class__.__name__
	Macro.setLocal("class_type", class_type)

	if class_type == "scest_output":
		print(res_est)
	
	filename = '__scest__output.obj'
	file     = open(filename, 'wb')
	pickle.dump(res_est, file)	
			
	Matrix.create("wmat", len(res_est.w), 1, 0)
	Matrix.store("wmat", res_est.w.values)
	names = [str(row) for row in res_est.w.index.tolist()]
	Matrix.setRowNames("wmat", names)
	
	if res_est.KMI > 0:
		Matrix.create("rmat", len(res_est.r), 1, 0)
		Matrix.store("rmat", res_est.r.values)
		names = [ix2rn(row) for row in res_est.r.index.tolist()]
		Matrix.setRowNames("rmat", names)
	
	Matrix.create("bmat", len(res_est.b), 1, 0)
	Matrix.store("bmat", res_est.b.values)
	names = [ix2rn(row) for row in res_est.b.index.tolist()]
	Matrix.setRowNames("bmat", names)

	if class_type == "scest_output":
		Matrix.create("Y_post", len(res_est.Y_post), 1, 0)
		Matrix.store("Y_post", res_est.Y_post.values)
		names = [ix2rn(row) for row in res_est.Y_post.index.tolist()]
		Matrix.setRowNames("Y_post", names)

		Matrix.create("Y_pre", len(res_est.Y_pre), 1, 0)
		Matrix.store("Y_pre", res_est.Y_pre.values)
		names = [ix2rn(row) for row in res_est.Y_pre.index.tolist()]
		Matrix.setRowNames("Y_pre", names)
		
	Matrix.create("Y_post_fit", len(res_est.Y_post_fit), 1, 0)
	Matrix.store("Y_post_fit", res_est.Y_post_fit.values)
	names = [ix2rn(row) for row in res_est.Y_post_fit.index.tolist()]
	Matrix.setRowNames("Y_post_fit", names)
	
	Matrix.create("Y_pre_fit", len(res_est.Y_pre_fit), 1, 0)
	Matrix.store("Y_pre_fit", res_est.Y_pre_fit.values)
	names = [ix2rn(row) for row in res_est.Y_pre_fit.index.tolist()]
	Matrix.setRowNames("Y_pre_fit", names)

	Matrix.create("res", len(res_est.res), 1, 0)
	Matrix.store("res", res_est.res.values)
	names = [ix2rn(row) for row in res_est.res.index.tolist()]
	Matrix.setRowNames("res", names)

	Matrix.create("pred", len(res_est.A_hat), 1, 0)
	Matrix.store("pred", res_est.A_hat.values)
	names = [ix2rn(row) for row in res_est.A_hat.index.tolist()]
	Matrix.setRowNames("pred", names)

	Matrix.create("A", len(res_est.A), 1, 0)
	Matrix.store("A", res_est.A.values)
	names = [ix2rn(row) for row in res_est.A.index.tolist()]
	Matrix.setRowNames("A", names)
	
	Matrix.create("B", len(res_est.B), len(res_est.B.columns), 0)
	Matrix.store("B", res_est.B.values)
	names = [ix2rn(row) for row in res_est.B.index.tolist()]
	Matrix.setRowNames("B", names)	
	names = [ix2rn(col) for col in res_est.B.columns.tolist()]
	Matrix.setColNames("B", names)	

	if res_est.KMI > 0:
		Matrix.create("C", len(res_est.C), len(res_est.C.columns), 0)
		Matrix.store("C", res_est.C.values)
		names = [ix2rn(row) for row in res_est.C.index.tolist()]
		Matrix.setRowNames("C", names)	
		names = [ix2rn(col) for col in res_est.C.columns.tolist()]
		Matrix.setColNames("C", names)	

	Matrix.create("P", len(res_est.P), len(res_est.P.columns), 0)
	Matrix.store("P", res_est.P.values)
	names = [ix2rn(row) for row in res_est.P.index.tolist()]
	Matrix.setRowNames("P", names)	
	names = [ix2rn(col) for col in res_est.P.columns.tolist()]
	Matrix.setColNames("P", names)				

	wc_est = res_est.w_constr[res_est.treated_units[0]]
	Macro.setLocal("dire", str(wc_est['dir']))
	Macro.setLocal("outcomevar", res_est.outcome_var)
	Scalar.setValue("KMI", res_est.KMI)
	Scalar.setValue("I", res_est.iota)

	names = res_est.treated_units
	Matrix.create("Qstar", res_est.iota, 1, 0)
	if p == "no norm" or name == "ols":
		Matrix.store("Qstar", [numpy.nan] * res_est.iota)
	else:
		Matrix.store("Qstar", [w['Q'] for w in res_est.w_constr.values()])
	Matrix.setRowNames("Qstar", names)

	Macro.setLocal("donors", dict2glob(res_est.donors_dict.keys(), res_est.donors_dict.values()))
	Macro.setLocal("anticipation", dict2glob(res_est.anticipation.keys(), res_est.anticipation.values()))
	
	if class_type == "scest_output":
	
		Matrix.create("T1", 1, 1, 0)
		Matrix.store("T1", res_est.T1_outcome)
		Matrix.setRowNames("T1", names)

		Matrix.create("J", 1, 1, 0)
		Matrix.store("J", res_est.J)
		Matrix.setRowNames("J", names)

		Matrix.create("KM", 1, 1, 0)
		Matrix.store("KM", res_est.KM)
		Matrix.setRowNames("KM", names)

		Matrix.create("M", 1, 1, 0)
		Matrix.store("M", res_est.M)
		Matrix.setRowNames("M", names)
		
		Matrix.create("T0_features", len(res_est.features), 1, 0)
		Matrix.store("T0_features", numpy.array([c for c in res_est.T0_features.values()]))
		Matrix.setRowNames("T0_features", res_est.features)		
			
		Macro.setLocal("cointegrated_data", str(res_est.cointegrated_data))
		Macro.setLocal("constant", str(res_est.glob_cons))
		Macro.setLocal("features", ' '.join([str(elem) for elem in res_est.features]))
		Macro.setLocal("anticipation", str(res_est.anticipation))


	if class_type == "scest_multi_output":
	
		Matrix.create("T1", len(res_est.T1_outcome), 1, 0)
		names = [str(k) for k in res_est.T1_outcome.keys()]
		Matrix.store("T1", [v for v in res_est.T1_outcome.values()])
		Matrix.setRowNames("T1", names)

		Matrix.create("J", len(res_est.J), 1, 0)
		names = [str(k) for k in res_est.J.keys()]
		Matrix.store("J", [v for v in res_est.J.values()])
		Matrix.setRowNames("J", names)

		Matrix.create("KM", len(res_est.KM), 1, 0)
		names = [str(k) for k in res_est.KM.keys()]
		Matrix.store("KM", [v for v in res_est.KM.values()])
		Matrix.setRowNames("KM", names)

		Matrix.create("M", len(res_est.M), 1, 0)
		names = [str(k) for k in res_est.M.keys()]
		Matrix.store("M", [v for v in res_est.M.values()])
		Matrix.setRowNames("M", names)

		unique_feats = []
		for v in res_est.T0_features.values():
			for k in v.keys():
				if not k in unique_feats:
					unique_feats.append(k)
		aux = numpy.zeros((res_est.iota, len(unique_feats)))
		for i in range(res_est.iota):
			for j in range(len(unique_feats)):
				aux[i,j] = res_est.T0_features[res_est.treated_units[i]][unique_feats[j]]

		Matrix.create("T0_features", res_est.iota, len(unique_feats), 0)
		Matrix.store("T0_features", aux)
		Matrix.setRowNames("T0_features", res_est.treated_units)				
		Matrix.setColNames("T0_features", unique_feats)				
		
		Macro.setLocal("cointegrated_data", dict2glob(res_est.cointegrated_data.keys(), res_est.cointegrated_data.values()))
		Macro.setLocal("constant", dict2glob(res_est.glob_cons.keys(), res_est.glob_cons.values()))
		Macro.setLocal("features", dict2glob(res_est.features.keys(), res_est.features.values()))
		Macro.setLocal("anticipation", dict2glob(res_est.anticipation.keys(), res_est.anticipation.values()))


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
end

























