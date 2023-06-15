{smcl}
{* *!version 2.2.1 2023-03-14}{...}
{viewerjumpto "Syntax" "scdatamulti##syntax"}{...}
{viewerjumpto "Description" "scdatamulti##description"}{...}
{viewerjumpto "Options" "scdatamulti##options"}{...}
{viewerjumpto "Examples" "scdatamulti##examples"}{...}
{viewerjumpto "Details" "scdatamulti##details"}{...}
{viewerjumpto "Stored results" "scdatamulti##stored_results"}{...}
{viewerjumpto "References" "scdatamulti##references"}{...}
{viewerjumpto "Authors" "scdatamulti##authors"}{...}

{title:Title}

{p 4 8}{cmd:scdatamulti} {hline 2} Data Preparation for Synthetic Control Methods with Staggered Adoption.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:scdatamulti } {it:features} {ifin} 
{cmd:,} 
{cmd:id(}{it:idvar}{cmd:)} 
{cmd:time(}{it:timevar}{cmd:)}
{cmd:outcome(}{it:outcomevar}{cmd:)}
{cmd:treatment(}{it:treatmentvar}{cmd:)}  
{cmd:dfname(}{it:string}{cmd:)}
[{cmd:covadj(}{it:string}{cmd:)}
{cmd:cointegrated(}{it:string}{cmd:)}
{cmd:constant(}{it:string}{cmd:)}
{cmd:anticipation(}{it:string}{cmd:)}
{cmd:effect(}{it:string}{cmd:)}
{cmd:post_est(}{it:string}{cmd:)}
{cmd:units_est(}{it:string}{cmd:)}
{cmd:donors_est(}{it:string}{cmd:)}
{cmd:pypinocheck}]{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:scdatamulti} prepares the data to be used by {help scest:scest} or {help scpi:scpi} to implement estimation and inference procedures for Synthetic Control (SC) methods
in the general case of multiple treated units and staggered adoption. It allows the user to specify for each treated unit the features 
to be matched, covariate-adjustment feature by feature, anticipation effects, and presence of cointegration. The command follows the terminology proposed in 
{browse "https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA.pdf":Cattaneo, Feng, and Titiunik (2021)}. The command is a wrapper of 
the companion Python package. As such, the user needs to have a running version of Python with the package installed. A tutorial on how to install Python and link it to Stata can be found {browse "https://nppackages.github.io/scpi/":here}.{p_end}

{p 8 8} Companion {browse "www.r-project.org":R} and {browse "https://www.python.org/":Python} packages are described in 
{browse "https://arxiv.org/abs/2202.05984":Cattaneo, Feng, Palomba and Titiunik (2022)}.{p_end}

{p 8 8} Companion commands are: {help scdata:scdata} for data preparation in the single treated unit case, 
{help scest:scest} for point estimation, {help scpi:scpi} for inference procedures, and {help scplot:scplot} for SC plots.{p_end}

{p 4 8}Related Stata, R, and Python packages useful for inference in SC designs are described in the following website:{p_end}

{p 8 8}{browse "https://nppackages.github.io/scpi/":https://nppackages.github.io/scpi/}{p_end}

{p 4 8}For an introduction to synthetic control methods, see {browse "https://economics.mit.edu/files/17847":Abadie (2021)} and 
references therein.{p_end}

{marker options}{...}
{title:Options}

{dlgtab:Variables}

{p 4 8}{cmd:id(}{it:idvar}{cmd:)} specifies the variable containing the identifier for each unit.{p_end}

{p 4 8}{cmd:time(}{it:timevar}{cmd:)} specifies the variable containing the time period of each observation.{p_end}

{p 4 8}{cmd:outcome(}{it:outcomevar}{cmd:)} specifies the outcome variable of interest. Note that {it:outcomevar} may not be among the {it:features} specified.{p_end}

{p 4 8}{cmd:treatment(}{it:treatmentvar}{cmd:)} specifies the treatment indicator.{p_end}

{dlgtab:Estimator}

{p 4 8}{cmd:covadj(}{it:string}{cmd:)} specifies the variable to be used for adjustment for each features for each treated unit. 
    If the user wants to specify the same set of covariates for all features, a string should be provided according to the following format: {opt covadj("cov1, cov2")}. 
    If instead a different set of covariates per feature has to be specified, then the following format should be used {opt covadj("cov1, cov2; cov1, cov3")}. 
    Note that in this latter case the number of sub-lists delimited by ";" must be equal to the number of {it:features}. Moreover, the order of the 
    sub-lists matters, in the sense that the first sub-list is interpreted as the set of covariates used for adjustment for the first 
    feature, and so on. Finally, the user can specify 'constant' and 'trend' as covariates even if they are not 
    present in the loaded dataset, the former includes a constant, whilst the latter a linear deterministic trend. 
    See Details section for more.{p_end}


{p 4 8}{cmd:cointegrated(}{it:string}{cmd:)} a logical value (the input should be either True or False) that 
specifies the presence of a cointegrating relationship between the features of the 
treated unit(s) and the the features of the donors. Default is {cmd:cointegrated("False")}. It can be specified for each treated unit.
See Details section for more.{p_end}

{p 4 8}{cmd:constant(}{it:string}{cmd:)} a logical value (the input should be either True or False) that  includes a common constant 
term across features. Default is {cmd:constant("False"}}. It can be specified for each treated unit.
See Details section for more.{p_end}

{p 4 8}{cmd:anticipation(}{it:string}{cmd:)} specifies the number of periods of potential anticipation effects. Default is no anticipation.
Note that it has to be a string, e.g. {cmd: anticipation("1")}. It can be specified for each treated unit.
See Details section for more.{p_end}

{p 4 8}{cmd:effect(}{it:string}{cmd:)} a string indicating the type of treatment effect to be estimated. Options are: 'unit-time', which estimates
        treatment effects for each treated unit-time combination; 'unit', which estimates the treatment effect for each unit by averaging post-treatment features over time;
        'time', which estimates the average treatment effect on the treated at various horizons.{p_end}

{p 4 8}{cmd:post_est(}{it:string}{cmd:)} a string specifying the number of post-treatment periods for which treatment effects have to be estimated for each treated unit. If
        effect = "unit" it indicates the number of periods over which the average post-treatment effect is computed. Note that it has to be a string, e.g. {cmd: post_est("1")}.{p_end}

{p 4 8}{cmd:units_est(}{it:string}{cmd:)} a string specifying the treated units for which treatment effects have to be estimated. Treated 
        units must be separated by commas, e.g. {cmd:units_est("unit1, unit2, unit3")}.{p_end}

{p 4 8}{cmd:donors_est(}{it:string}{cmd:)} a string specifying the donors units to be used. Note that all treated units share the same
        potential donors. If this is not desired, the donor pool can be separately specified for each treated unit. See Details section for more.{p_end}

{dlgtab:Others}

{p 4 8}{cmd:dfname(}{it:string}{cmd:)} specifies the name of the Python object that is saved and that will be passed to {help scest:scest} or {help scpi:scpi}.{p_end}

{p 4 8}{cmd:pypinocheck)} if specified avoids to check that the version of scpi_pkg in Python is the one required by {cmd:scdata} in Stata. When not specified performs the check and stores a macro called to avoid checking it multiple times.{p_end}

    {hline}

{marker details}{...}
{title:Details}

{p 4 8}This section describes how to use {cmd:scdatamulti} in two cases: first, when the user wants a common specification across treated units;
second, when the user wants to tailor her specification for each treated unit. {p_end}

{dlgtab:Common Specification}

{p 4 8}Let's start first with the simple case of common specification across treated units. Suppose, for the sake of the example, that 
there are just two treated units and two features to be matched on. The command would simply be{p_end}

{p 8 8}{cmd:scdatamulti feature1 feature2, id(idvar) outcome(feature1) treatment(trvar) time(timevar)}{p_end}

{p 4 8}If covariate adjustment, cointegration, anticipation effects, and a global constant need to be specified for each treated unit, then{p_end}

{p 8 8}{cmd:scdatamulti feature1 feature2, id(idvar) outcome(feature1) treatment(trvar) time(timevar) ///}{p_end}
{p 8 8}{cmd:constant(True) cointegrated(True) anticipation(1) covadj("constant, trend")}{p_end}


{dlgtab:Heterogeneous Specification}

{p 4 8}Again, suppose there are two treated units and an individual specification is desired. In particular, we would like to 
match one feature of unit one and two features of the second unit. Then{p_end}

{p 8 8}{cmd:scdatamulti (unit1: feature1) (unit2: feature1 feature2), id(idvar) outcome(feature1) treatment(trvar) time(timevar) ///}{p_end}
{p 8 8}{cmd:            constant($cons_spec) cointegrated($coint_spec) anticipation($ant_spec) covadj($cov_spec)}{p_end}

{p 4 8}Where the globals are defined as follows{p_end}

{p 6 8}First, we specify covariate adjustment just for the first feature of both treated units adding a linear trend for the first unit and a constant
term for the second unit.{p_end}

{p 8 8}{cmd:global cov_spec = "(unit1: trend) (unit2: constant; None)"}

{p 6 8}Second, we add a global constant for both treated units. There are two equivalent ways to do it:{p_end} 

{p 8 8}{cmd:global cons_spec = "True"}{p_end}
{p 8 8}{cmd:global cons_spec = "(unit1: True) (unit2: True)"}{p_end}

{p 6 8} Similarly,{p_end}

{p 8 8}{cmd:global coint_spec = "(unit1: True) (unit2: True)"}{p_end}
{p 8 8}{cmd:global ant_spec = "(unit1: 0) (unit2: 1)"}{p_end}
{p 8 8}{cmd:global donors_spec = "(unit1: donor1 donor2) (unit2: donor2 donor3)"}{p_end}

{marker examples}{...}
{title:Example: Germany Data}

{p 4 8}Setup{p_end}
{p 8 8}{cmd:. use scpi_germany.dta}{p_end}

{p 4 8}Prepare data{p_end}
{p 8 8}{cmd:. scdata gdp, dfname("python_scdata") id(country) outcome(gdp) time(year) treatment(status) cointegrated}{p_end}


{marker stored_results}{...}
{title:Stored results}

{p 4 8}{cmd:scdata} stores the following in {cmd:e()}:

{synoptset 25 tabbed}{...}

{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(I)}}number of treated units{p_end}
{synopt:{cmd:e(KMI)}}total number of covariates used for adjustment{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(features)}}name of features{p_end}
{synopt:{cmd:e(outcomevar)}}name of outcome variable{p_end}
{synopt:{cmd:e(constant)}}logical indicating the presence of a common constant across features{p_end}
{synopt:{cmd:e(cointegrated)}}logical indicating cointegration{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(A)}}pre-treatment features of the treated unit{p_end}
{synopt:{cmd:e(B)}}pre-treatment features of the control units{p_end}
{synopt:{cmd:e(C)}}covariates used for adjustment{p_end}
{synopt:{cmd:e(P)}}predictor matrix{p_end}
{synopt:{cmd:e(J)}}number of donors for each treated unit{p_end}
{synopt:{cmd:e(KM)}}total number of covariates used for adjustment for each treated unit{p_end}

{marker references}{...}
{title:References}

{p 4 8}Abadie, A. 2021. 
{browse "https://economics.mit.edu/files/17847":Using synthetic controls: Feasibility, data requirements, and methodological aspects.} 
{it:Journal of Economic Literature}, 59(2), 391-425.{p_end}

{p 4 8}Cattaneo, M. D., Feng, Y., and Titiunik, R. 2021. 
{browse "https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA.pdf":Prediction intervals for synthetic control methods}. 
{it:Journal of the American Statistical Association}, 116(536), 1865-1880.{p_end}

{p 4 8}Cattaneo, M. D., Feng, Y., Palomba F., and Titiunik, R. 2022. 
{browse "https://arxiv.org/abs/2202.05984":scpi: Uncertainty Quantification for Synthetic Control Estimators}, {it:arXiv}:2202.05984.{p_end}

{p 4 8}Cattaneo, M. D., Feng, Y., Palomba F., and Titiunik, R. 2022. 
{browse "https://arxiv.org/abs/2210.05026":Uncertainty Quantification in Synthetic Controls with Staggered Treatment Adoption}, {it:arXiv}:2210.05026. {p_end}

{marker authors}{...}
{title:Authors}

{p 4 8}Matias D. Cattaneo, Princeton University, Princeton, NJ.
{browse "mailto:cattaneo@princeton.edu":cattaneo@princeton.edu}.{p_end}

{p 4 8}Yingjie Feng, Tsinghua University, Beijing, China.
{browse "mailto:fengyj@sem.tsinghua.edu.cn":fengyj@sem.tsinghua.edu.cn}.

{p 4 8}Filippo Palomba, Princeton University, Princeton, NJ.
{browse "mailto:fpalomba@princeton.edu":fpalomba@princeton.edu}.

{p 4 8}Rocio Titiunik, Princeton University, Princeton, NJ.
{browse "mailto:titiunik@princeton.edu":titiunik@princeton.edu}.{p_end}