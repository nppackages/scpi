*! Date        : 14 Mar 2023
*! Version     : 2.2.1
*! Authors     : Filippo Palomba
*! Email       : fpalomba@princeton.edu
*! Description : Data preparation for scest or scpi

capture program drop scdata
program define scdata, eclass         
version 16.0           
		
	syntax varlist [if] [in], id(varname) time(varname) outcome(varname) treatment(varname) dfname(string) ///
							  [covadj(string) anticipation(integer 0) cointegrated pypinocheck constant]

	if mi("`pypinocheck'") & mi("$scpi_version_checked") {
		python: version_checker()
		if "`alert_version'" == "y" {
			di as error "The current version of scpi_pkg in Python is `python_local_version', but version `python_pypi_version' needed! Please update the package in Python and restart Stata!"
			exit 198
		}
		global scpi_version_checked "yes"
	}
		

	local features "`varlist'"	  

	if mi("`covadj'") {
		local covadj "None"
	}
	
	if mi("`cointegrated'") {
		local cointegrated "False"
	} 
	
	if mi("`constant'") {
		local constant "False"
	} 
	
	qui export delimited using "__scpi__data_to_python.csv", replace
	
	python: scdata_wrapper("`features'", "`id'", "`time'", "`outcome'", "`covadj'", `anticipation', "`cointegrated'", "`constant'", "False", "`treatment'", "`dfname'")
	
	erase "__scpi__data_to_python.csv"
	
	ereturn clear

	ereturn scalar KM = scalar(KM)
	ereturn scalar J = scalar(J)
	
	ereturn local effect = "unit-time"
	ereturn local cointegrated_data = "`cointegrated_data'"
	ereturn local constant          = "`constant'"
	ereturn local outcomevar        = "`outcomevar'"
	ereturn local features          = "`features'"
	
	ereturn matrix P  = P
	capture ereturn matrix C = C
	ereturn matrix B  = B
	ereturn matrix A  = A	  
end





version 16.0
python:
import pandas, pickle, numpy, urllib, luddite
from scpi_pkg.scdata import scdata
from scpi_pkg import version as lver
from sfi import Scalar, Matrix, Macro

def ix2rn(s):
	return str(s).replace('(','').replace(')','').replace("'",'').replace(", ","_")

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

def scdata_wrapper(features, id_var, time_var, outcome_var, covadj, anticipation, cointegrated, constant, reportmissing, treatment, dfname):
	
	# Create dataframe
	df = pandas.read_csv('__scpi__data_to_python.csv')

	# Create treatment and control groups
	NperiodsT = df[[treatment, id_var]].groupby([id_var]).sum()
	unit_tr   = NperiodsT.loc[NperiodsT[treatment] > 0,].index.tolist()
	unit_co   = NperiodsT.loc[NperiodsT[treatment] == 0,].index.tolist()

	# Create period_pre and period_post
	NunitsT  = df[[treatment, time_var]].groupby([time_var]).sum()
	aux_post = NunitsT.loc[NunitsT[treatment] > 0,].index.tolist()
	aux_pre  = NunitsT.loc[NunitsT[treatment] == 0,].index.tolist()	
	period_post = numpy.array(aux_post)
	period_pre  = numpy.array(aux_pre)
	
	# Parse features
	aux = features.split()
	f_list = [f.strip() for f in aux]
	
	# Parse cointegrated
	if cointegrated == "False":
		cointegrated_bool = False
	else:
		cointegrated_bool = True
		
	# Parse constant
	if constant == "False":
		constant_bool = False
	else:
		constant_bool = True
	
	# Parse covadj
	if covadj == "None":
		cov_adj_list = None
	else:
		lists = covadj.split(";")
		if len(lists) == 1:
			covs = lists[0].split(",")
			cov_adj_list = []
			for cov in covs:
				cov_adj_list.append(cov.strip())
		else:
			cov_adj_list = []
			for l in lists:
				covs = l.split(",")
				cov_list = []
				for cov in covs:
					cov_list.append(cov.strip())
				cov_adj_list.append(cov_list)
		
		
	data_prep = scdata(df, id_var, time_var, outcome_var, period_pre, period_post, unit_tr[0], unit_co, f_list, cov_adj_list, cointegrated_bool, anticipation, constant_bool)
	
	filename = dfname + '.obj'
	file     = open(filename, 'wb')
	pickle.dump(data_prep, file, protocol = pickle.HIGHEST_PROTOCOL)
	
	Matrix.create("A", len(data_prep.A), 1, 0)
	Matrix.store("A", data_prep.A.values)
	names = [ix2rn(row) for row in data_prep.A.index.tolist()]
	Matrix.setRowNames("A", names)
	
	Matrix.create("B", len(data_prep.B), data_prep.J, 0)
	Matrix.store("B", data_prep.B.values)
	names = [ix2rn(row) for row in data_prep.B.index.tolist()]
	Matrix.setRowNames("B", names)	
	names = [str(col) for col in data_prep.B.columns.tolist()]
	Matrix.setColNames("B", names)	

	if data_prep.KM > 0:
		Matrix.create("C", len(data_prep.C), data_prep.KM, 0)
		Matrix.store("C", data_prep.C.values)
		names = [ix2rn(row) for row in data_prep.C.index.tolist()]
		Matrix.setRowNames("C", names)	
		names = [str(col) for col in data_prep.C.columns.tolist()]
		Matrix.setColNames("C", names)	

	Matrix.create("P", len(data_prep.P), len(data_prep.P.columns), 0)
	Matrix.store("P", data_prep.P.values)
	names = [ix2rn(row) for row in data_prep.P.index.tolist()]
	Matrix.setRowNames("P", names)	
	names = [str(col) for col in data_prep.P.columns.tolist()]
	Matrix.setColNames("P", names)	
	
	Scalar.setValue("KM", data_prep.KM)
	Scalar.setValue("J", data_prep.J)
	Macro.setLocal("cointegrated_data", str(data_prep.cointegrated_data))
	Macro.setLocal("constant", str(data_prep.glob_cons))
	Macro.setLocal("outcomevar", data_prep.outcome_var)
	Macro.setLocal("features", ' '.join([str(elem) for elem in data_prep.features]))
	
end

