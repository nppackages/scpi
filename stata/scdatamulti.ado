*! Date        : 28 Jul 2022
*! Version     : 1.0
*! Authors     : Filippo Palomba
*! Email       : fpalomba@princeton.edu
*! Description : Data preparation for scest or scpi

capture program drop scdatamulti
program define scdatamulti, eclass         
version 17.0           
		
	syntax anything [if] [in], id(varname) time(varname) outcome(varname) treatment(varname) dfname(string) ///
							  [covadj(string) cointegrated(string) constant(string) anticipation(string) effect(string) pypinocheck]

	local features "`anything'"  // backup copy

	if mi("`covadj'") {
		local covadj "None"
	}
	
	if mi("`cointegrated'") {
		local cointegrated "False"
	} 
	
	if mi("`constant'") {
		local constant "False"
	} 
	
	if mi("`anticipation'") {
		local anticipation "0"
	}
	
	if mi("`effect'") {
		local effect = "unit-time"
	}
	else {
		if !inlist("`effect'", "unit-time", "unit") {
			di as error "The option 'effect' should be either 'unit' or 'unit-time'!"
		}
	}
	
	if mi("`pypinocheck'") & mi("$scpi_version_checked") {
		python: version_checker()
		if "`alert_version'" == "y" {
			di as error "The current version of scpi_pkg in Python is `python_local_version', but version `python_pypi_version' needed! Please update the package in Python and restart Stata!"
			exit 198
		}
		global scpi_version_checked "yes"
	}
	
	qui export delimited using "__scpi__data_to_python.csv", replace

	python: scdatamulti_wrapper("`features'", "`id'", "`time'", "`outcome'", "`treatment'", "`covadj'", "`anticipation'", "`cointegrated'", "`constant'", "`dfname'", "`effect'")
	
	erase "__scpi__data_to_python.csv"
	
	ereturn clear

	ereturn scalar I = I
	ereturn scalar KMI = KMI
	ereturn matrix KM = KM
	ereturn matrix J = J
	
	ereturn local cointegrated_data = "`cointegrated_data'"
	ereturn local constant          = "`constant'"
	ereturn local outcomevar        = "`outcomevar'"
	ereturn local features          = "`features'"
	
	ereturn matrix P = P
	capture ereturn matrix C = C
	ereturn matrix B = B
	ereturn matrix A = A	  
	
end



version 17.0
python:
import pandas, pickle, numpy, urllib, luddite, re
from scpi_pkg.scdataMulti import scdataMulti
from scpi_pkg import version as lver
from sfi import Scalar, Matrix, Macro


def scdatamulti_wrapper(features, id_var, time_var, outcome_var, treatment, covadj, anticipation, cointegrated, constant, dfname, effect):
	
	# Create dataframe
	df = pandas.read_csv('__scpi__data_to_python.csv')

	# Compute number of treatment periods for each unit to distinguish between single and multiple treated units
	NperiodsT = df[[treatment, id_var]].groupby([id_var]).sum()
	unit_tr   = NperiodsT.loc[NperiodsT[treatment] > 0,].index.tolist()
	
	char_rem = "['()]"
	
	## Parse features	
	if features[0] != "(":  # common spec
		aux = features.split()
		f_list = [f.strip() for f in aux]
		f_dict = {'features': f_list}
		
	else: # unit by unit specification
		
		parsed_features = re.findall('\(.*?\)', features)
		
		i = 1
		for feat in parsed_features:
			aux = re.sub(char_rem,"",feat).split(":")
			# name of treated unit
			tr_name = aux[0].strip()
			
			# list of features of treated unit
			tr_feat = aux[1].strip().split()
			
			if i == 1:
				f_dict = {tr_name: tr_feat}
			else:
				f_dict[tr_name] = tr_feat
			i += 1


	## Parse covariates for adjustment
	if covadj == "None":
		cov_dict = None
	else:
		if covadj[0] != "(":  # common spec
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
						covv = cov.strip()
						if covv == "None":
							break
						cov_list.append(cov.strip())
					cov_adj_list.append(cov_list)
			cov_dict = {'covs': cov_adj_list}
		
		else: 
			parsed_covs = re.findall('\(.*?\)', covadj)
			
			i = 1
			for cov in parsed_covs:
				aux = re.sub(char_rem,"",cov).split(":")
				
				# name of treated unit
				tr_name = aux[0].strip()
				
				# list of covariates of treated unit
				lists = aux[1].split(";")
				
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
							covv = cov.strip()
							if covv == "None":
								break
							cov_list.append(cov.strip())
						cov_adj_list.append(cov_list)

				if i == 1:
					cov_dict = {tr_name: cov_adj_list}
				else:
					cov_dict[tr_name] = cov_adj_list
				i += 1

				
	## Parse constant (if boolean scdataMulti handle it automatically)
	if constant[0] == "(":  			
		parsed_constant = re.findall('\(.*?\)', constant)
		
		i = 1
		for const in parsed_constant:
			aux = re.sub(char_rem,"",const).split(":")
			
			# name of treated unit
			tr_name = aux[0].strip()
			
			# list of features of treated unit
			tr_cons = aux[1].strip()
			
			if tr_cons == "True":
				tr_cons = True
			else:
				tr_cons = False
			
			if i == 1:
				cons_dict = {tr_name: tr_cons}
			else:
				cons_dict[tr_name] = tr_cons
			i += 1
	else:
		if constant == "True":
			cons_dict = True
		else:
			cons_dict = False

	## Parse cointegrated_data (if boolean scdataMulti handle it automatically)
	if cointegrated[0] == "(":  			
		parsed_cointegration = re.findall('\(.*?\)', cointegrated)
		
		i = 1
		for coint in parsed_cointegration:
			aux = re.sub(char_rem,"",coint).split(":")
			
			# name of treated unit
			tr_name = aux[0].strip()
			
			# list of features of treated unit
			tr_coint = aux[1].strip()
			
			if tr_coint == "True":
				tr_coint = True
			else:
				tr_coint = False
			
			if i == 1:
				coint_dict = {tr_name: tr_cons}
			else:
				coint_dict[tr_name] = tr_cons
			i += 1
	else:
		if cointegrated == "True":
			coint_dict = True
		else:
			coint_dict = False
			
			
	## Parse anticipation (if integer scdataMulti handle it automatically)
	if anticipation[0] == "(":  			
		parsed_anticipation = re.findall('\(.*?\)', anticipation)
		
		i = 1
		for antic in parsed_anticipation:
			aux = re.sub(char_rem,"",antic).split(":")
			
			# name of treated unit
			tr_name = aux[0].strip()
			
			# list of features of treated unit
			tr_antic = int(aux[1].strip())
							
			if i == 1:
				antic_dict = {tr_name: tr_cons}
			else:
				antic_dict[tr_name] = tr_cons
			i += 1
	else:
		antic_dict = int(anticipation)
	

	data_prep = scdataMulti(df=df, id_var=id_var, time_var=time_var, outcome_var=outcome_var, treatment_var=treatment,
							features=f_dict, cov_adj=cov_dict, constant=cons_dict, anticipation=antic_dict, cointegrated_data=coint_dict, effect=effect)


	# Save data and store matrices in stata
	
	filename = dfname + '.obj'
	file     = open(filename, 'wb')
	pickle.dump(data_prep, file, protocol = pickle.HIGHEST_PROTOCOL)
	
	Matrix.create("A", len(data_prep.A), 1, 0)
	Matrix.store("A", data_prep.A.values)
	names = [re.sub("['()]","",str(row)).replace(" ","").replace(",","_") for row in data_prep.A.index.tolist()]
	Matrix.setRowNames("A", names)
	
	Matrix.create("B", len(data_prep.B), len(data_prep.B.columns), 0)
	Matrix.store("B", data_prep.B.values)
	names = [re.sub("['()]","",str(row)).replace(" ","").replace(",","_") for row in data_prep.B.index.tolist()]
	Matrix.setRowNames("B", names)	
	names = [str(col) for col in data_prep.B.columns.tolist()]
	Matrix.setColNames("B", names)	
	
	if data_prep.KMI > 0:
		Matrix.create("C", len(data_prep.C), len(data_prep.C.columns), 0)
		Matrix.store("C", data_prep.C.values)
		names = [re.sub("['()]","",str(row)).replace(" ","").replace(",","_") for row in data_prep.C.index.tolist()]
		Matrix.setRowNames("C", names)	
		names = [str(col) for col in data_prep.C.columns.tolist()]
		Matrix.setColNames("C", names)	

	Matrix.create("P", len(data_prep.P), len(data_prep.P.columns), 0)
	Matrix.store("P", data_prep.P.values)
	names = [str(row) for row in data_prep.P.index.tolist()]
	Matrix.setRowNames("P", names)	
	names = [str(col) for col in data_prep.P.columns.tolist()]
	Matrix.setColNames("P", names)	
	
	Scalar.setValue("I", data_prep.iota)
	Scalar.setValue("KMI", data_prep.KMI)
	
	klist = []; vlist = []
	for k,v in data_prep.KM.items():
		klist.append(k) 
		vlist.append(v)
	mat = numpy.array(vlist)
	Matrix.create("KM", 1, len(data_prep.KM), 0)
	Matrix.store("KM", mat)
	Matrix.setColNames("KM", klist)
	
	klist = []; vlist = []
	for k,v in data_prep.KM.items():
		klist.append(k) 
		vlist.append(v)
	mat = numpy.array(vlist)
	Matrix.create("J", 1, len(data_prep.J), 0)
	Matrix.store("J", mat)
	Matrix.setColNames("J", klist)
	
	Macro.setLocal("cointegrated_data", str(data_prep.cointegrated_data))
	Macro.setLocal("constant", str(data_prep.glob_cons))
	Macro.setLocal("outcomevar", data_prep.outcome_var)
	Macro.setLocal("features", ' '.join([str(elem) for elem in data_prep.features]))						   
									

		
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
