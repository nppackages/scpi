{smcl}
{* *!version 0.1 2022-01-25}{...}
{viewerjumpto "Syntax" "scpi##syntax"}{...}
{viewerjumpto "Description" "scpi##description"}{...}
{viewerjumpto "Options" "scpi##options"}{...}
{viewerjumpto "Examples" "scpi##examples"}{...}
{viewerjumpto "Stored results" "scpi##stored_results"}{...}
{viewerjumpto "References" "scpi##references"}{...}
{viewerjumpto "Authors" "scpi##authors"}{...}

{title:Title}

{p 4 8}{cmd:scpi} {hline 2} Estimation and Inference for Synthetic Control Methods.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:scpi } 
{cmd:,} 
{cmd:dfname(}{it:string}{cmd:)} 
[{cmd:p(}{it:#}{cmd:)}
{cmd:direc(}{it:string}{cmd:)}
{cmd:Q(}{it:#}{cmd:)}  
{cmd:lb(}{it:#}{cmd:)}  
{cmd:name(}{it:string}{cmd:)}
{cmd:u_missp)}
{cmd:u_sigma(}{it:string}{cmd:)}
{cmd:u_order(}{it:#}{cmd:)}
{cmd:u_lags(}{it:#}{cmd:)}
{cmd:u_alpha(}{it:#}{cmd:)}
{cmd:sims(}{it:#}{cmd:)}
{cmd:e_method(}{it:string}{cmd:)}
{cmd:e_order(}{it:#}{cmd:)}
{cmd:e_lags(}{it:#}{cmd:)}
{cmd:e_alpha(}{it:#}{cmd:)}
{cmd:rho(}{it:#}{cmd:)}
{cmd:rho_max(}{it:#}{cmd:)}
{cmd:opt_est(}{it:string}{cmd:)}]
{cmd:opt_inf(}{it:string}{cmd:)}]
{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:scpi} implements estimation and inference procedures for Synthetic Control (SC) methods using least squares, lasso, ridge, or simplex-type constraints according to
{browse "https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf":Cattaneo, Feng, and Titiunik (2021)}. The command is a wrapper of the companion Python package. 
As such, the user needs to have a running version of Python with the package installed. A tutorial on how to install Python and link it to Stata
can be found {browse "https://nppackages.github.io/scpi/":here}.{p_end}

{p 8 8} Companion {browse "www.r-project.org":R} and {browse "https://www.python.org/":Python} packages are described in 
{browse "https://arxiv.org/abs/2202.05984":Cattaneo, Feng, Palomba and Titiunik (2022)}.{p_end}

{p 8 8} Companion commands are: {help scdata:scdata} for data preparation, {help scest:scest} for estimation procedures, and {help scplot:scplot} for SC plots.{p_end}

{p 4 8}Related Stata, R, and Python packages useful for inference in SC designs are described in the following website:{p_end}

{p 8 8}{browse "https://nppackages.github.io/scpi/":https://nppackages.github.io/scpi/}{p_end}

{p 4 8}For an introduction to synthetic control methods, see {browse "https://economics.mit.edu/files/17847":Abadie (2021)} and 
references therein.{p_end}

{marker options}{...}
{title:Options}

{p 4 8}{cmd:dfname(}{it:string}{cmd:)} specifies the name of the Python object containing the processed data created with {help scdata:scdata}.{p_end}

{dlgtab:Constraint}

{p 2 4} These options let the user specify the type of constraint to be imposed to estimate the SC weights. The user controls the lower bound on the weights (option {opt lb}), the norm of the weights to
be constrained (option {opt p}), the direction of the constraint on the norm (option {opt dir}), and the size of the constraint on the norm (option {opt q}). Alternatively,
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
{p 8 12}{opt ols} weights are estimated without constraints using least squares{p_end}

{dlgtab:In-sample Uncertainty}

{p 2 4} This set of options allows the user to specify her preferred approch to model in-sample uncertainty, that is uncertainty that stems from the estimation the weights.{p_end}

{p 4 8}{cmd:u_missp} if specified indicates that model misspecification should be taken into account.{p_end}
{p 4 8}{cmd:u_sigma(}{it:string}{cmd:)} specifies the type of variance-covariance estimator to be used when estimating the conditional variance of the pseudo-residuals. 
Options are: {opt HC0}, {opt HC1} (default), {opt HC2}, and {opt HC3}.{p_end}
{p 4 8}{cmd:u_order(}{it:#}{cmd:)} specifies the order of the polynomial in the predictors used to estimate conditional moments of the pseudo-residuals. Default is {cmd:u_order(1)}.{p_end}
{p 4 8}{cmd:u_lags(}{it:#}{cmd:)} specifies the lags of the predictors used to estimate conditional moments of the pseudo-residuals. Default is {cmd:u_lags(0)}.{p_end}
{p 4 8}{cmd:u_alpha(}{it:#}{cmd:)} specifies the confidence level for in-sample uncertainty. Default is {cmd:u_alpha(0.05)}.{p_end}
{p 4 8}{cmd:sims(}{it:#}{cmd:)} specifies the number of simulations to be used in quantifying in-sample uncertainty. Default is {cmd:sims(200)}.{p_end}

{dlgtab:Out-of-sample Uncertainty}

{p 2 4} This set of options allows the user to specify her preferred approch to model out-of-sample uncertainty.{p_end}

{p 4 8}{cmd:e_method(}{it:#}{cmd:)} specifies the method to be used to quantify out-of-sample uncertainty. Options are:{p_end}
{p 8 12} {opt gaussian} conditional subgaussian bounds. {p_end}
{p 8 12} {opt ls} location-scale model. {p_end}
{p 8 12} {opt qreg} quantile regression. {p_end}
{p 8 12} {opt all} all of the above (the default). {p_end}

{p 4 8}{cmd:e_order(}{it:#}{cmd:)} specifies the order of the polynomial in the predictors used to estimate conditional moments of the out-of-sample error. Default is {cmd:e_order(1)}.{p_end}
{p 4 8}{cmd:e_lags(}{it:#}{cmd:)} specifies the lags of the predictors used to estimate conditional moments of the out-of-sample error. Default is {cmd:e_lags(0)}.{p_end}
{p 4 8}{cmd:e_alpha(}{it:#}{cmd:)} specifies the confidence level for out-fo-sample uncertainty. Default is {cmd:e_alpha(0.05)}.{p_end}

{dlgtab:Regularization}

{p 4 8}{cmd:rho(}{it:#}{cmd:)} specifies the regularizing parameter that imposes sparsity on the estimated vector of weights. If unspecified, the tuning parameter is computed 
based on optimization inequalities.{p_end}
{p 4 8}{cmd:rho_max(}{it:#}{cmd:)} specifies the maximum value attainable by the tuning parameter {opt rho}. {p_end}

{dlgtab:Others}

{p 4 8}{cmd:opt_est(}{it:string}{cmd:)} a string specifying the stopping criteria used by the underling optimizer (nlopt) for point estimation. 
    The default is a sequential quadratic programming (SQP) algorithm for nonlinearly constrained gradient-based optimization ('SLSQP'). 
    In case a lasso-type constraint is implemented, the method of moving asymptotes (MMA) is used. The default value is 
    {cmd:opt("'maxeval' = 5000, 'xtol_rel' = 1e-8, 'xtol_abs' = 1e-8, 'ftol_rel' = 1e-12, 'ftol_abs' = 1e-12, 'tol_eq' = 1e-8, 'tol_ineq' = 1e-8")}.{p_end}

{p 4 8}{cmd:opt_inf(}{it:string}{cmd:)} a string specifying the stopping criteria used by the underling optimizer (nlopt) for point estimation. 
    The default is a sequential quadratic programming (SQP) algorithm for nonlinearly constrained gradient-based optimization ('SLSQP'). 
    In case a lasso-type constraint is implemented, the method of moving asymptotes (MMA) is used. The default value is 
    {cmd:opt("'maxeval' = 5000, 'xtol_rel' = 1e-8, 'xtol_abs' = 1e-8, 'ftol_rel' = 1e-4, 'ftol_abs' = 1e-4, 'tol_eq' = 1e-8, 'tol_ineq' = 1e-8")}.{p_end}

    {hline}


{marker examples}{...}
{title:Example: Germany Data}

{p 4 8}Setup{p_end}
{p 8 8}{cmd:. use scpi_germany.dta}{p_end}

{p 4 8}Prepare data{p_end}
{p 8 8}{cmd:. scdata gdp, dfname("python_scdata") id(country) outcome(gdp) time(year) treatment(status) cointegrated}{p_end}

{p 4 8}Estimate Synthetic Control with a simplex constraint and quantify uncertainty{p_end}
{p 8 8}{cmd:. scpi, dfname("python_scdata") name(simplex) u_missp}{p_end}

marker stored_results}{...}
{title:Stored results}

{p 4 8}{cmd:scpi} stores the following in {cmd:e()}:

{synoptset 25 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(M)}}number of features{p_end}
{synopt:{cmd:e(KM)}}number of covariates used for adjustment{p_end}
{synopt:{cmd:e(J)}}number of donors{p_end}
{synopt:{cmd:e(T1)}}number of post-treatment periods{p_end}
{synopt:{cmd:e(q)}}size of the constraint on the norm{p_end}
{synopt:{cmd:e(rho)}}post-estimation regularization parameter{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(features)}}name of features{p_end}
{synopt:{cmd:e(outcomevar)}}name of outcome variable{p_end}
{synopt:{cmd:e(constant)}}logical indicating the presence of a common constant across features{p_end}
{synopt:{cmd:e(cointegrated_data)}}logical indicating cointegration{p_end}
{synopt:{cmd:e(p)}}type of norm of the weights used in constrained estimation{p_end}
{synopt:{cmd:e(dir)}}direction of the constraint on the norm of the weights{p_end}
{synopt:{cmd:e(name)}}name of constraint used in estimation{p_end}
{synopt:{cmd:e(u_missp)}}a logical indicating whether the model has been treated as misspecified or not{p_end}
{synopt:{cmd:e(u_order)}}order of the polynomial in the predictors used to estimate conditional moments of the pseudo-residuals{p_end}
{synopt:{cmd:e(u_lags)}}lags of the predictors used to estimate conditional moments of the pseudo-residuals{p_end}
{synopt:{cmd:e(u_sigma)}}estimator of the conditional variance-covariance matrix of the pseudo-residuals{p_end}
{synopt:{cmd:e(u_alpha)}}confidence level for in-sample uncertainty{p_end}
{synopt:{cmd:e(e_method)}}method used to quantify out-of-sample uncertainty{p_end}
{synopt:{cmd:e(e_order)}}order of the polynomial in the predictors used to estimate conditional moments of the out-of-sample error{p_end}
{synopt:{cmd:e(e_lags)}}order of the polynomial in the predictors used to estimate conditional moments of the out-of-sample error{p_end}
{synopt:{cmd:e(e_alpha)}}confidence level for out-of-sample uncertainty{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(T0)}}number of pre-treatment periods per feature{p_end}
{synopt:{cmd:e(A)}}pre-treatment features of the treated unit{p_end}
{synopt:{cmd:e(B)}}pre-treatment features of the control units{p_end}
{synopt:{cmd:e(C)}}covariates used for adjustment{p_end}
{synopt:{cmd:e(pred)}}predicted values of the features of the treated unit{p_end}
{synopt:{cmd:e(res)}}residuals {cmd:e(A)} - {cmd:e(pred)}{p_end}
{synopt:{cmd:e(w)}}weights of the controls{p_end}
{synopt:{cmd:e(r)}}coefficients of the covariates used for adjustment{p_end}
{synopt:{cmd:e(beta)}}stacked version of {cmd:e(w)} and {cmd:e(r)}{p_end}
{synopt:{cmd:e(Y_post)}}post-treatment outcome of the treated unit{p_end}
{synopt:{cmd:e(Y_post_fit)}}estimated post-treatment outcome of the treated unit{p_end}
{synopt:{cmd:e(Y_pre)}}pre-treatment outcome of the treated unit{p_end}
{synopt:{cmd:e(Y_pre_fit)}}estimate pre-treatment outcome of the treated unit{p_end}
{synopt:{cmd:e(CI_in_sample)}}prediction intervals taking only in-sample uncertainty into account{p_end}
{synopt:{cmd:e(CI_all_gaussian)}}prediction intervals taking in- and out-of-sample uncertainty into account{p_end}
{synopt:{cmd:e(CI_all_ls)}}prediction intervals taking in- and out-of-sample uncertainty into account{p_end}
{synopt:{cmd:e(CI_all_qreg)}}prediction intervals taking in- and out-of-sample uncertainty into account{p_end}
{synopt:{cmd:e(u_mean)}}estimated conditional mean of the pseudo-residuals{p_end}
{synopt:{cmd:e(u_var)}}estimated conditional variance-covariance of the pseudo-residuals{p_end}
{synopt:{cmd:e(e_mean)}}estimated conditional mean of the out-of-sample error{p_end}
{synopt:{cmd:e(e_var)}}estimated conditional variance of the out-of-sample error{p_end}
{synopt:{cmd:e(failed_sims)}}percentage of failed simulations per post-treatment period to estimate lower and upper bounds.{p_end}

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
