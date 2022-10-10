{smcl}
{* *!version 0.3 2022-03-07}{...}
{viewerjumpto "Syntax" "scest##syntax"}{...}
{viewerjumpto "Description" "scest##description"}{...}
{viewerjumpto "Options" "scest##options"}{...}
{viewerjumpto "Examples" "scest##examples"}{...}
{viewerjumpto "Stored results" "scest##stored_results"}{...}
{viewerjumpto "References" "scest##references"}{...}
{viewerjumpto "Authors" "scest##authors"}{...}

{title:Title}

{p 4 8}{cmd:scest} {hline 2} Estimation for Synthetic Control Methods.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:scest } 
{cmd:,} 
{cmd:dfname(}{it:string}{cmd:)} 
[{cmd:p(}{it:#}{cmd:)}
{cmd:direc(}{it:string}{cmd:)}
{cmd:Q(}{it:#}{cmd:)}  
{cmd:lb(}{it:#}{cmd:)}  
{cmd:name(}{it:string}{cmd:)}
{cmd:V(}{it:string}{cmd:)}
{cmd:opt(}{it:string}{cmd:)}
{cmd:pypinocheck}]{p_end}

{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:scest} implements estimation procedures for Synthetic Control (SC) methods using least squares, lasso, ridge, or simplex-type constraints according to
{browse "https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf":Cattaneo, Feng, and Titiunik (2021)}. The command is a wrapper of the companion Python package. 
As such, the user needs to have a running version of Python with the package installed. A tutorial on how to install Python and link it to Stata
can be found {browse "https://nppackages.github.io/scpi/":here}.{p_end}

{p 8 8} Companion {browse "www.r-project.org":R} and {browse "https://www.python.org/":Python} packages are described in 
{browse "https://arxiv.org/abs/2202.05984":Cattaneo, Feng, Palomba and Titiunik (2022)}.{p_end}

{p 8 8} Companion commands are: {help scdata:scdata} for data preparation, {help scpi:scpi} for inference procedures, and {help scplot:scplot} for SC plots.{p_end}

{p 4 8}Related Stata, R, and Python packages useful for inference in SC designs are described in the following website:{p_end}

{p 8 8}{browse "https://nppackages.github.io/scpi/":https://nppackages.github.io/scpi/}{p_end}

{p 4 8}For an introduction to synthetic control methods, see {browse "https://economics.mit.edu/files/17847":Abadie (2021)} and 
references therein.{p_end}

{marker options}{...}
{title:Options}

{p 4 8}{cmd:dfname(}{it:string}{cmd:)} specifies the name of the Python object containing the processed data created with {help scdata:scdata}.{p_end}

{dlgtab:Loss Function and Constraints}

{p 2 4} These options let the user specify the type of constraint to be imposed to estimate the SC weights and the loss function. The user controls the lower bound on the weights (option {opt lb}),
the norm of the weights to be constrained (option {opt p}), the direction of the constraint on the norm (option {opt dir}), the size of the constraint on the norm (option {opt q}), and
the shape of the wieghting matrix in the loss function (option {opt V}). Alternatively,
some popular constraints can be selected through the option {opt name}. A detailed description of the popular constraints implemented can be found in 
{browse "https://nppackages.github.io/references/Cattaneo-Feng-Palomba-Titiunik_2022_scpi.pdf":Cattaneo, Feng, Palomba and Titiunik (2022)}. {p_end}

{p 4 8}{cmd:lb(}{it:#}{cmd:)} specifies the lower bound on the weights. The default is {cmd:lb(0)}. {p_end}

{p 4 8}{cmd:p(}{it:#}{cmd:)} sets the type of norm to be constrained. Options are:{p_end}
{p 8 12} {opt 0} no constraint on the norm of the weights is imposed. {p_end}
{p 8 12} {opt 1} a constraint is imposed on the L1 norm of the weights (the default). {p_end}
{p 8 12} {opt 2} a constraint is imposed on the L2 norm of the weights. {p_end}

{p 4 8}{cmd:direc(}{it:string}{cmd:)} specifies the direction of the constraint on the norm of the weights. Options are:{p_end}
{p 8 12} {opt "<="} the constraint on the norm of the weights is an inequality constraint. {p_end}
{p 8 12} {opt "=="} the constraint on the norm of the weights is an equality constraint (the default). {p_end}

{p 4 8}{cmd:Q(}{it:#}{cmd:)} specifies the size of the constraint on the norm of the weights. {p_end}

{p 4 8}{cmd:name(}{it:string}{cmd:)} specifies the name of the constraint to be used. Options are:{p_end}
{p 8 12}{opt simplex} classic synthetic control estimator where the weights are constrained to be non-negative and their L1 norm must be equal to 1.{p_end}
{p 8 12}{opt lasso} weights are estimated using a Lasso-type penalization{p_end}
{p 8 12}{opt ridge} weights are estimated using a Ridge-type penalization.{p_end}
{p 8 12}{opt L1-L2} weights are estimated using a Simplex-type constraint and a Ridge-type penalization.{p_end}
{p 8 12}{opt ols} weights are estimated without constraints using least squares{p_end}

{p 4 8}{cmd:V(}{it:string}{cmd:)} specifies the weighting matrix to be used in the loss function. The default is the identity matrix (option {cmd:V("separate")}), so equal weight is given to all observations. The other possibility is to 
	specify {cmd:V("pooled")} for the pooled fit.{p_end}

{dlgtab:Others}

{p 4 8}{cmd:opt(}{it:string}{cmd:)} a string specifying the stopping criteria used by the underling optimizer ({browse "https://nlopt.readthedocs.io/en/latest/NLopt_Python_Reference/":nlopt}) 
    for point estimation. The default is a sequential quadratic programming (SQP) algorithm for nonlinearly constrained gradient-based optimization ('SLSQP'). 
    The default value is {cmd:opt("'maxeval' = 5000, 'xtol_rel' = 1e-8, 'xtol_abs' = 1e-8, 'ftol_rel' = 1e-12, 'ftol_abs' = 1e-12, 'tol_eq' = 1e-8, 'tol_ineq' = 1e-8")}.
    In case a lasso-type constraint is implemented, a different optimizer ({browse "https://www.cvxpy.org/":cvxpy}) is used and stopping criteria cannot be changed. {p_end}

{p 4 8}{cmd:pypinocheck)} if specified avoids to check that the version of scpi_pkg in Python is the one required by {cmd:scest} in Stata. When not specified performs the check and stores a macro called to avoid checking it multiple times.{p_end}

    {hline}


{marker examples}{...}
{title:Example: Germany Data}

{p 4 8}Setup{p_end}
{p 8 8}{cmd:. use scpi_germany.dta}{p_end}

{p 4 8}Prepare data{p_end}
{p 8 8}{cmd:. scdata gdp, dfname("python_scdata") id(country) outcome(gdp) time(year) treatment(status) cointegrated}{p_end}

{p 4 8}Estimate Synthetic Control with a simplex constraint{p_end}
{p 8 8}{cmd:. scest, dfname("python_scdata") name(simplex)}{p_end}


{marker stored_results}{...}
{title:Stored results}

{p 4 8}{cmd:scest} stores the following in {cmd:e()}:

{synoptset 25 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(KMI)}}number of covariates used for adjustment{p_end}
{synopt:{cmd:e(I)}}number of treated units{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(outcomevar)}}name of outcome variable{p_end}
{synopt:{cmd:e(features)}}name of features{p_end}
{synopt:{cmd:e(constant)}}logical indicating the presence of a common constant across features{p_end}
{synopt:{cmd:e(anticipation)}}logical indicating the extent of anticipation effects{p_end}
{synopt:{cmd:e(donors)}}donor units for each treated unit{p_end}
{synopt:{cmd:e(cointegrated_data)}}logical indicating cointegration{p_end}
{synopt:{cmd:e(p)}}type of norm of the weights used in constrained estimation{p_end}
{synopt:{cmd:e(dir)}}direction of the constraint on the norm of the weights{p_end}
{synopt:{cmd:e(name)}}name of constraint used in estimation{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(A)}}pre-treatment features of the treated unit{p_end}
{synopt:{cmd:e(B)}}pre-treatment features of the control units{p_end}
{synopt:{cmd:e(C)}}covariates used for adjustment{p_end}
{synopt:{cmd:e(P)}}covariates used to predict the out-of-sample series for the synthetic unit{p_end}
{synopt:{cmd:e(pred)}}predicted values of the features of the treated unit{p_end}
{synopt:{cmd:e(res)}}residuals {cmd:e(A)} - {cmd:e(pred)}{p_end}
{synopt:{cmd:e(w)}}weights of the controls{p_end}
{synopt:{cmd:e(r)}}coefficients of the covariates used for adjustment{p_end}
{synopt:{cmd:e(beta)}}stacked version of {cmd:e(w)} and {cmd:e(r)}{p_end}
{synopt:{cmd:e(Y_pre)}}pre-treatment outcome of the treated unit (only returned if one treated unit){p_end}
{synopt:{cmd:e(Y_post)}}post-treatment outcome of the treated unit (only returned if one treated unit){p_end}
{synopt:{cmd:e(Y_pre_fit)}}estimate pre-treatment outcome of the treated unit{p_end}
{synopt:{cmd:e(Y_post_fit)}}estimated post-treatment outcome of the treated unit{p_end}
{synopt:{cmd:e(T0)}}number of pre-treatment periods per feature for each treated unit{p_end}
{synopt:{cmd:e(M)}}number of features for each treated unit{p_end}
{synopt:{cmd:e(KM)}}number of covariates used for adjustment for each treated unit{p_end}
{synopt:{cmd:e(J)}}number of donors for each treated unit{p_end}
{synopt:{cmd:e(T1)}}number of post-treatment periods for each treated unit{p_end}
{synopt:{cmd:e(Qstar)}}regularized constraint on the norm{p_end}

{marker references}{...}
{title:References}

{p 4 8}Abadie, A. 2021. 
{browse "https://economics.mit.edu/files/17847":Using synthetic controls: Feasibility, data requirements, and methodological aspects.} 
{it:Journal of Economic Literature}, 59(2), 391-425.{p_end}

{p 4 8}Cattaneo, M. D., Feng, Y., and Titiunik, R. 2021. 
{browse "https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf":Prediction Intervals for Synthetic Sontrol Methods.}
{it:Journal of the American Statistical Association}, 116(536), 1865-1880.{p_end}

{p 4 8}Cattaneo, M. D., Feng, Y., Palomba F., and Titiunik, R. 2022. 
{browse "https://arxiv.org/abs/2202.05984":scpi: Uncertainty Quantification for Synthetic Control Estimators, {it:arXiv}:2202.05984.}.{p_end}

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
