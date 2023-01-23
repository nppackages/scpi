{smcl}
{* *!version 2.2.0 2023-01-23}{...}
{viewerjumpto "Syntax" "scplot##syntax"}{...}
{viewerjumpto "Description" "scplot##description"}{...}
{viewerjumpto "Options" "scplot##options"}{...}
{viewerjumpto "Examples" "scplot##examples"}{...}
{viewerjumpto "Stored results" "scplot##stored_results"}{...}
{viewerjumpto "References" "scplot##references"}{...}
{viewerjumpto "Authors" "scplot##authors"}{...}

{title:Title}

{p 4 8}{cmd:scplotmulti} {hline 2} Synthetic Control Methods with Multiple Treated Units Plots.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:scplot } 
{cmd:,} 
[{cmd:scest}
{cmd:uncertainty(}{it:string}{cmd:)}
{cmd:uncertainty(}{it:string}{cmd:)}
{cmd:joint}
{cmd:yscalefree}
{cmd:xscalefree}
{cmd:dots_tr_col(}{it:{help colorstyle:colorstyle}}{cmd:)}
{cmd:dots_tr_symb(}{it:{help symbolstyle:symbolstyle}}{cmd:)}
{cmd:dots_tr_size(}{it:{help markersizestyle:markersizestyle}}{cmd:)}
{cmd:dots_sc_col(}{it:{help colorstyle:colorstyle}}{cmd:)}
{cmd:dots_sc_symb(}{it:{help symbolstyle:symbolstyle}}{cmd:)}
{cmd:dots_sc_size(}{it:{help markersizestyle:markersizestyle}}{cmd:)}
{cmd:line_tr_col(}{it:{help colorstyle:colorstyle}}{cmd:)}
{cmd:line_tr_patt(}{it:{help linepatternstyle:linepatternstyle}}{cmd:)}
{cmd:line_tr_width(}{it:{help linewidthstyle:linewidthstyle}}{cmd:)}
{cmd:line_sc_col(}{it:{help colorstyle:colorstyle}}{cmd:)}
{cmd:line_sc_patt(}{it:{help linepatternstyle:linepatternstyle}}{cmd:)}
{cmd:line_sc_width(}{it:{help linewidthstyle:linewidthstyle}}{cmd:)}
{cmd:spike_sc_col(}{it:{help colorstyle:colorstyle}}{cmd:)}
{cmd:spike_sc_patt(}{it:{help linepatternstyle:linepatternstyle}}{cmd:)}
{cmd:spike_sc_width(}{it:{help linewidthstyle:linewidthstyle}}{cmd:)}
{cmd:gphoptions(}{it:string}{cmd:)}
{cmd:gphcombineoptions(}{it:{help graph combine:graph combine}}{cmd:)}
{cmd:gphsave(}{it:string}{cmd:)}
{cmd:savedata(}{it:dta_name}{cmd:)}
{cmd:keepsingleplots}
{cmd:pypinocheck}]{p_end}
{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:scplot} implements several Synthetic Control (SC) plots even in the presence of multiple treated units and staggered adoption.
The command is designed te be called after {help scest:scest} or {help scpi:scpi} which implement  
estimation and inference procedures for SC methods using least squares, lasso, ridge, or simplex-type constraints according to
{browse "https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf":Cattaneo, Feng, and Titiunik (2021)}. The command is a wrapper of the companion Python package. 
As such, the user needs to have a running version of Python with the package installed. A tutorial on how to install Python and link it to Stata
can be found {browse "https://nppackages.github.io/scpi/":here}.{p_end}

{p 8 8} Companion {browse "www.r-project.org":R} and {browse "https://www.python.org/":Python} packages are described in 
{browse "https://arxiv.org/abs/2202.05984":Cattaneo, Feng, Palomba and Titiunik (2022)}.{p_end}

{p 8 8} Companion commands are: {help scdata:scdata} for data preparation, {help scest:scest} for estimation procedures, and {help scpi:scpi} for inference procedures.{p_end}

{p 4 8}Related Stata, R, and Python packages useful for inference in SC designs are described in the following website:{p_end}

{p 8 8}{browse "https://nppackages.github.io/scpi/":https://nppackages.github.io/scpi/}{p_end}

{p 4 8}For an introduction to synthetic control methods, see {browse "https://economics.mit.edu/files/17847":Abadie (2021)} and 
references therein.{p_end}

{marker options}{...}
{title:Options}

{dlgtab:Type of Plot}

{p 4 8}{cmd:scest} if specified {cmd:scplot} must be called after {help scest:scest}. Otherwise, it is presumed that {cmd:scplot} is called after {help scpi:scpi}.{p_end}

{p 4 8}{cmd:uncertainty(}{it:string}{cmd:)} specifies which prediction intervals are plotted. It does not affect the plot if {opt scest} is specified. Options are:{p_end}
{p 8 12} {opt insample} prediction intervals quantify only in-sample uncertainty. {p_end}
{p 8 12} {opt gaussian} prediction intervals quantify in-sample and out-of-sample uncertainty using conditional subgaussian bounds. {p_end}
{p 8 12} {opt ls} prediction intervals quantify in-sample and out-of-sample uncertainty imposing a location-scale model. {p_end}
{p 8 12} {opt qreg} prediction intervals quantify in-sample and out-of-sample uncertainty using quantile regressions. {p_end}

{p 4 8}{cmd:ptype(}{it:string}{cmd:)} specifies the type of plot to be produced. If set to 'treatment', then treatment effects are
        plotted. If set to 'series' (default), the actual and synthetic time series are reported.{p_end}

{p 4 8}{cmd:joint} if specified simultaneous prediction intervals are included in the plot(s).{p_end}

{dlgtab:Scale Options}

{p 4 8}{cmd:yscalefree} if specified each graph has its own scale for the y axis.{p_end}
{p 4 8}{cmd:xscalefree} if specified each graph has its own scale for the x axis.{p_end}


{dlgtab:Marker Options}

{p 2 4} These options let the user specify color, size, and form of the markers in the plot.{p_end}

{p 4 8} {cmd:dots_tr_col(}{it:{help colorstyle:colorstyle}}{cmd:)} specifies the color of the markers for the treated unit.{p_end}
{p 4 8} {cmd:dots_tr_symb(}{it:{help symbolstyle:symbolstyle}}{cmd:)} specifies the form of the markers for the treated unit.{p_end}
{p 4 8} {cmd:dots_tr_size(}{it:{help markersizestyle:markersizestyle}}{cmd:)} specifies the size of the markers for the treated unit.{p_end}
{p 4 8} {cmd:dots_sc_col(}{it:{help colorstyle:colorstyle}}{cmd:)} specifies the color of the markers for the SC unit.{p_end}
{p 4 8} {cmd:dots_sc_symb(}{it:{help symbolstyle:symbolstyle}}{cmd:)} specifies the form of the markers for the SC unit.{p_end}
{p 4 8} {cmd:dots_sc_size(}{it:{help markersizestyle:markersizestyle}}{cmd:)} specifies the size of the markers for the SC unit.{p_end}

{dlgtab:Line Options}

{p 2 4} These options let the user specify color, pattern, and width of the lines in the plot.{p_end}

{p 4 8} {cmd:line_tr_col(}{it:{help colorstyle:colorstyle}}{cmd:)} specifies the color of the line for the treated unit.{p_end}
{p 4 8} {cmd:line_tr_patt(}{it:{help linepatternstyle:linepatternstyle}}{cmd:)} specifies the pattern of the line for the treated unit.{p_end}
{p 4 8} {cmd:line_tr_width(}{it:{help linewidthstyle:linewidthstyle}}{cmd:)} specifies the width of the line for the treated unit.{p_end}
{p 4 8} {cmd:line_sc_col(}{it:{help colorstyle:colorstyle}}{cmd:)} specifies the color of the line for the SC unit.{p_end}
{p 4 8} {cmd:line_sc_patt(}{it:{help linepatternstyle:linepatternstyle}}{cmd:)} specifies the pattern of the line for the SC unit.{p_end}
{p 4 8} {cmd:line_sc_width(}{it:{help linewidthstyle:linewidthstyle}}{cmd:)} specifies the width of the line for the SC unit.{p_end}

{dlgtab:Bar Options}

{p 2 4} These options let the user specify color, pattern, and width of the bar (spikes) in the plot. These options do not have effect if {opt scest} is specified.{p_end}

{p 4 8} {cmd:spike_sc_col(}{it:{help colorstyle:colorstyle}}{cmd:)} specifies the color of the bars for the SC unit.{p_end}
{p 4 8} {cmd:spike_sc_patt(}{it:{help linepatternstyle:linepatternstyle}}{cmd:)} specifies the pattern of the bars for the SC unit.{p_end}
{p 4 8} {cmd:spike_sc_width(}{it:{help linewidthstyle:linewidthstyle}}{cmd:)} specifies the width of the bars for the SC unit.{p_end}

{dlgtab:Others}

{p 4 8}{cmd:gphoptions(}{it:string}{cmd:)} specifies additional options to modify individual plots.{p_end}
{p 4 8}{cmd:gphcombineoptions(}{it:{help graph combine:graph combine}}{cmd:)} specifies additional options to modify the final combined plot.{p_end}
{p 4 8}{cmd:gphsave(}{it:string}{cmd:)} specifies the path and the name of the {it:.gph} file that is saved by the command.{p_end}
{p 4 8}{cmd:savedata(}{it:dta_name}{cmd:)} saves a {it:dta_name.dta} file containing the processed data used to produce the plot.{p_end}
{p 4 8}{cmd:keepsingleplots)} if specified saves the individual plots in .gph format in the current directory.{p_end}

{p 4 8}{cmd:pypinocheck)} if specified avoids to check that the version of scpi_pkg in Python is the one required by {cmd:scplot} in Stata. When not specified performs the check and stores a macro called to avoid checking it multiple times.{p_end}

    {hline}


{marker examples}{...}
{title:Example: Germany Data}

{p 4 8}Setup{p_end}
{p 8 8}{cmd:. use scpi_germany.dta}{p_end}

{p 4 8}Prepare data{p_end}
{p 8 8}{cmd:. scdata gdp, dfname("python_scdata") id(country) outcome(gdp) time(year) treatment(status) cointegrated}{p_end}

{p 4 8}Estimate Synthetic Control with a simplex constraint and quantify uncertainty{p_end}
{p 8 8}{cmd:. scpi, dfname("python_scdata") name(simplex) u_missp}{p_end}

{p 4 8}Plot Synthetic Control Estimate with Prediction Intervals{p_end}
{p 8 8}{cmd:. scplot, gphsave("plot_scpi")}{p_end}

{marker references}{...}
{title:References}

{p 4 8}Abadie, A. 2021. 
{browse "https://economics.mit.edu/files/17847":Using synthetic controls: Feasibility, data requirements, and methodological aspects.} 
{it:Journal of Economic Literature}, 59(2), 391-425.{p_end}

{p 4 8}Cattaneo, M. D., Feng, Y., and Titiunik, R. 2021. 
{browse "https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf":Prediction Intervals for Synthetic Sontrol Methods}. 
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
