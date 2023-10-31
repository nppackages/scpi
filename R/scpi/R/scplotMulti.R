#' @title Plot Synthetic Control Point Estimates and Prediction Interval With Multiple Treated units and Staggered Adoption
#'
#' @description The command produces a wide range of plots of Synthetic Control estimates and corresponding prediction intervals. The command allows form multiple treated units and staggered adoption.
#' Prediction intervals can take into account either in-sample uncertainty only or in-sample and
#' out-of-sample uncertainty using the techniques developed in \href{https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, Feng, and Titiunik (2021)}. \code{\link{scpi}}.
#' The input object should come from the command \code{\link{scest}} or from the command \code{\link{scpi}}.
#'
#' Companion \href{https://www.stata.com/}{Stata} and \href{https://www.python.org/}{Python} packages are described in \href{https://arxiv.org/abs/2202.05984}{Cattaneo, Feng, Palomba, and Titiunik (2022)}.
#'
#' Companion commands are:  \link{scdata} and \link{scdataMulti} for data preparation in the single and multiple treated unit(s) cases, respectively,
#' \link{scest} for point estimation, \link{scpi} for inference procedures, and \link{scplotMulti} for plots with multiple treated units.
#'
#' Related Stata, R, and Python packages useful for inference in SC designs are described in the following website:
#'
#' \href{ https://nppackages.github.io/scpi/}{ https://nppackages.github.io/scpi/}
#'
#' For an introduction to synthetic control methods, see \href{https://www.aeaweb.org/articles?id=10.1257/jel.20191450}{Abadie (2021)} and references therein.
#'
#' @param result a class 'scest' object, obtained by calling \code{\link{scest}}, or a class
#' 'scpi' object, obtained by calling \code{\link{scpi}}. The data object given as input to this command has to be
#' processed with \code{\link{scdataMulti}}.
#' @param type a character that specifies the type of plot to be produced. If set to 'treatment' then treatment effects are plotted. 
#' If set to 'series' (default), the actual and synthetic time series are reported.
#' @param e.out a logical specifying whether out-of-sample uncertainty should be included in the plot(s).
#' @param joint a logical specifying whether simultaneous prediction intervals should be included in the plot(s). It requires \code{e.out = TRUE}.
#' @param col.treated a string specifying the color for the treated unit series. Find the full list at \url{http://sape.inf.usi.ch/quick-reference/ggplot2/colour}.
#' @param col.synth a string specifying the color for the synthetic unit series. Find the full list at \url{http://sape.inf.usi.ch/quick-reference/ggplot2/colour}. 
#' @param scales should axes scales be fixed ("fixed", the default), free ("free"), or free in one dimension ("free_x", "free_y")?
#' @param point.size a scalar controlling the size of points in the scatter plot. Default is 1.5.
#' @param ncols an integer controlling the number of columns in the plot.
#' @param save.data a character specifying the name and the path of the saved dataframe containing the processed data used to produce the plot. 
#' @param verbose if \code{TRUE} prints additional information in the console.
#'
#' @return
#' \item{plots}{a list containing standard ggplot object(s) that can be used for further customization.}
#'
#' @author
#' Matias Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Yingjie Feng, Tsinghua University. \email{fengyj@sem.tsinghua.edu.cn}.
#'
#' Filippo Palomba, Princeton University (maintainer). \email{fpalomba@princeton.edu}.
#'
#' Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu}.
#'
#' @references
#' \itemize{
#' \item{\href{https://www.aeaweb.org/articles?id=10.1257/jel.20191450}{Abadie, A. (2021)}. Using synthetic controls: Feasibility, data requirements, and methodological aspects.
#' \emph{Journal of Economic Literature}, 59(2), 391-425.}
#' \item{\href{https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, M. D., Feng, Y., and Titiunik, R. 
#' (2021)}. Prediction intervals for synthetic control methods. \emph{Journal of the American Statistical Association}, 116(536), 1865-1880.} 
#' \item{\href{https://arxiv.org/abs/2202.05984}{Cattaneo, M. D., Feng, Y., Palomba F., and Titiunik, R. (2022)},
#' scpi: Uncertainty Quantification for Synthetic Control Methods, \emph{arXiv}:2202.05984.}
#' \item{\href{https://arxiv.org/abs/2210.05026}{Cattaneo, M. D., Feng, Y., Palomba F., and Titiunik, R. (2022).}
#' Uncertainty Quantification in Synthetic Controls with Staggered Treatment Adoption, \emph{arXiv}:2210.05026.}
#' }
#'
#' @seealso \code{\link{scdata}}, \code{\link{scdataMulti}}, \code{\link{scest}}, \code{\link{scpi}}, \code{\link{scplotMulti}}
#'
#' @examples
#'
#' datager <- scpi_germany
#'
#' datager$tr_id <- 0
#' datager$tr_id[(datager$country == "West Germany" & datager$year > 1990)] <- 1
#' datager$tr_id[(datager$country == "Italy" & datager$year > 1992)] <- 0
#'
#' outcome.var <- "gdp"
#' id.var <- "country"
#' treatment.var <- "tr_id"
#' time.var <- "year"
#' df.unit <- scdataMulti(datager, id.var = id.var, outcome.var = outcome.var,
#'                        treatment.var = treatment.var,
#'                        time.var = time.var, features = list(c("gdp", "trade")),
#'                		    cointegrated.data = TRUE, constant = TRUE)
#'
#' res.unit <- scpi(df.unit, sims = 10, cores = 1)
#' scplotMulti(res.unit, joint = TRUE)
#'
#' @export

scplotMulti <- function(result, type = "series", e.out = TRUE, joint = FALSE,
                        col.treated = "black", col.synth = "mediumblue", scales = "fixed",
                        point.size = 1.5, ncols = 3, save.data = NULL, verbose = TRUE) {

  Treatment <- ID <- Time <- Tdate <- Type <- tstd <- NULL

  if (!(result$data$specs$class.type %in% c("scpi_scest_multi", "scpi_scpi_multi"))) {
    stop("data should be the object returned by running scest or scpi when the data have been processed through scdataMulti!")
  }

  if (!(type %in% c("series", "treatment"))) {
    stop("'type' should be either 'series' or 'treatment'!")
  }

  plot.type <- result$data$specs$effect
  I <- result$data$specs$I
  if (plot.type == "time") {
    I <- 1
  }
  
  if (I > 20 && (plot.type != "unit" || type != "treatment") && verbose) {
    warning(paste0(I, " treated units detected, therefore some graphs might be too crowded! ",
                   "We suggest saving the data with the option save.data, consult our replication ",
                   "files at https://nppackages.github.io/scpi/, and reproduce the same graph for just",
                   " a fraction of the sample at a time!"), call. = FALSE, immediate. = TRUE)
  }

  class.type      <- result$data$specs$class.type
  sparse.matrices <- result$data$specs$sparse.matrices
  treated.units   <- result$data$specs$treated.units
  anticipation    <- result$data$specs$anticipation
  period.post     <- result$data$specs$period.post
  units.est       <- result$data$specs$units.est
  T1.outcome      <- result$data$specs$T1.outcome
  Y.df            <- result$data$Y.df
  Y.pre.fit       <- result$est.results$Y.pre.fit # matrix
  Y.post.fit      <- result$est.results$Y.post.fit # matrix

  # create to plot object
  toplot.obj <- outcomeGet(Y.pre.fit=Y.pre.fit, Y.post.fit=Y.post.fit, Y.df=Y.df,
                           units.est=units.est, treated.units=treated.units, plot.type=plot.type,
                           anticipation=anticipation, period.post=period.post,
                           sparse.matrices = sparse.matrices)

  toplot <- toplot.obj$toplot
  treated.reception <- toplot.obj$treated.reception
  plot.type <- toplot.obj$plot.type

  toplot$Effect <- toplot$Actual - toplot$Synthetic # compute treatment effect

  toplot <- merge(toplot, treated.reception, by = "ID") # compute periods since treated
  
  if (plot.type == 'unit-time' && type == "series") { # plot series
    toplot <- toplot[c("ID","Time","Actual", "Synthetic")]
    toplot <- reshape2::melt(toplot, id=c("ID","Time"))
    names(toplot) <- c("ID", "Time", "Type", "Y")

    toplot <- merge(toplot, treated.reception, by="ID")

    plot <- ggplot(toplot) + xlab("Date") + ylab("Outcome") +
            geom_line(aes(x=.data$Time, y=.data$Y, colour=.data$Type)) + 
            geom_point(aes(x=.data$Time, y=.data$Y, colour=.data$Type), size=point.size) + 
            geom_vline(data = treated.reception, aes(xintercept=.data$Tdate)) +
            facet_wrap(~.data$ID, ncol = ncols, scales = scales) + theme(legend.position="bottom") +
            scale_color_manual(name = "", values = c(col.treated, col.synth),
                               labels = c("Treated", "Synthetic Control"))

  } else if (plot.type == "unit-time" && type == "treatment") {  # plot treatment effects
    plot <- ggplot(toplot) + xlab("Date") + ylab("Effect") +
            geom_line(aes(x=.data$Time, y=.data$Effect), colour = col.synth) +
            geom_point(aes(x=.data$Time, y=.data$Effect), colour = col.synth, size=point.size) +
            geom_vline(data = treated.reception, aes(xintercept=.data$Tdate)) +
            geom_hline(aes(yintercept= 0), linetype = "dashed") +
            facet_wrap(~.data$ID, ncol = ncols, scales = scales) + theme(legend.position="bottom")

  } else if (plot.type == "unit" && type == "series") {

    # add number of post-treatment periods that have been averaged
    auxdf <- data.frame(ID = names(T1.outcome), T1 = unlist(T1.outcome))
    treated.reception <- merge(treated.reception, auxdf, by = "ID")

    toplot <- toplot[c("ID", "Time", "Actual", "Synthetic")]
    toplot <- reshape2::melt(toplot, id=c("ID", "Time"))
    names(toplot) <- c("ID", "Time", "Type", "Y")

    toplot <- merge(toplot, treated.reception, by="ID")
    plot <- ggplot(toplot) + xlab("Date") + ylab("Outcome") +
            geom_line(data=subset(toplot, Time < Tdate), aes(x=.data$Time, y=.data$Y, colour=.data$Type)) +
            geom_point(aes(x=.data$Time, y=.data$Y, colour=.data$Type), size=point.size) +
            geom_text(aes(x = Inf, y = -Inf, label = .data$T1), vjust = "inward", hjust = "inward") +
            facet_wrap(~.data$ID, ncol = ncols, scales = scales) + theme(legend.position = "bottom") +
            scale_color_manual(name = "", values = c(col.treated, col.synth),
                               labels = c("Treated", "Synthetic Control"))

  } else if (plot.type == "unit" && type == "treatment") {
    auxdf <- data.frame(ID = names(T1.outcome), Periods = unlist(T1.outcome))
    toplot <- merge(toplot, auxdf, by = "ID")
    plot <- ggplot(subset(toplot, Treatment == 1)) + xlab("Unit") + ylab("Effect") +
            geom_point(aes(x=.data$ID, y=.data$Effect, size=.data$Periods), colour = col.synth, shape = 19) +
            geom_hline(aes(yintercept= 0), linetype = "dashed") + coord_flip()

  }


  #############################################################################
  ## ADD UNCERTAINTY
  #############################################################################

  if (class.type == 'scpi_scpi_multi') {

    e.method <- result$inference.results$e.method

    if (type == "series") {
      toplot <- merge(toplot, ci2df(result$inference.results$CI.in.sample, type = "insample"),
                      by = c("ID", "Time"), all = TRUE)

      if (e.method %in% c("all", "gaussian")) {
        toplot <- merge(toplot, ci2df(result$inference.results$CI.all.gaussian, type = "gaussian"),
                        by = c("ID", "Time"), all = TRUE)
      }

      if (e.method %in% c("all", "ls")) {
        toplot <- merge(toplot, ci2df(result$inference.results$CI.all.ls, type = "ls"),
                        by = c("ID", "Time"), all = TRUE)
      }

      if (e.method %in% c("all", "qreg")) {
        toplot <- merge(toplot, ci2df(result$inference.results$CI.all.qreg, type = "qreg"),
                        by = c("ID", "Time"), all = TRUE)
      }

      if (joint == TRUE) {
        toplot <- merge(toplot, ci2df(result$inference.results$bounds$joint, type = "bd"),
                        by = c("ID", "Time"), all = TRUE)
        toplot$lb.joint <- toplot$Y + toplot$lb.bd
        toplot$ub.joint <- toplot$Y + toplot$ub.bd
        toplot[c("ub.bd", "lb.bd")] <- list(NULL)
      }

    } else if (type == "treatment") {
      toplot <- merge(toplot, ci2df(result$inference.results$bounds$insample, type = "bd"),
                      by = c("ID", "Time"), all = TRUE)
      toplot$lb.insample <- toplot$Effect - toplot$ub.bd
      toplot$ub.insample <- toplot$Effect - toplot$lb.bd
      toplot[c("ub.bd", "lb.bd")] <- list(NULL)

      if (e.method %in% c("all", "gaussian")) {
        toplot <- merge(toplot, ci2df(result$inference.results$bounds$subgaussian, type = "bd"),
                        by = c("ID", "Time"), all = TRUE)
        toplot$lb.gaussian <- toplot$Effect - toplot$ub.bd
        toplot$ub.gaussian <- toplot$Effect - toplot$lb.bd
        toplot[c("ub.bd", "lb.bd")] <- list(NULL)
      }

      if (e.method %in% c("all", "ls")) {
        toplot <- merge(toplot, ci2df(result$inference.results$bounds$ls, type = "bd"),
                        by = c("ID", "Time"), all = TRUE)
        toplot$lb.ls <- toplot$Effect - toplot$ub.bd
        toplot$ub.ls <- toplot$Effect - toplot$lb.bd
        toplot[c("ub.bd", "lb.bd")] <- list(NULL)
      }

      if (e.method %in% c("all", "qreg")) {
        toplot <- merge(toplot, ci2df(result$inference.results$bounds$qreg, type = "bd"),
                        by = c("ID", "Time"), all = TRUE)
        toplot$lb.qreg <- toplot$Effect - toplot$ub.bd
        toplot$ub.qreg <- toplot$Effect - toplot$lb.bd
        toplot[c("ub.bd", "lb.bd")] <- list(NULL)
      }

      if (joint == TRUE) {
        toplot <- merge(toplot, ci2df(result$inference.results$bounds$joint, type = "bd"),
                        by = c("ID", "Time"), all = TRUE)
        toplot$lb.joint <- toplot$Effect - toplot$ub.bd
        toplot$ub.joint <- toplot$Effect - toplot$lb.bd
        toplot[c("ub.bd", "lb.bd")] <- list(NULL)
      }

    }
    if (e.out == FALSE) {   # Only in-sample uncertainty

      if (type == "treatment" && plot.type == "unit") {
        plot.w <- plot +
                  geom_errorbar(data = toplot,
                                aes(x = .data$ID, ymin = .data$lb.insample, ymax = .data$ub.insample),
                                colour = "col.synth", width = 0.25, linetype = "solid") + 
                  ggtitle("In-sample Uncertainty")

      } else {
        
        plot.w <- plot + geom_errorbar(data = toplot,
                                       aes(x = .data$Time, ymin = .data$lb.insample, ymax = .data$ub.insample), colour = col.synth,
                                       width = 0.5, linetype = 'solid')  +  ggtitle("In-sample Uncertainty")
      }

    }

    if (e.out == TRUE && e.method %in% c("all", "gaussian")) {
      if (type == "treatment" && plot.type == "unit") {
        plot.w1 <- plot + geom_errorbar(data = toplot,
                                        aes(x = .data$ID, ymin = .data$lb.gaussian, ymax = .data$ub.gaussian), colour = col.synth,
                                        width = 0.25, linetype = 'solid') + ggtitle("In and Out of Sample Uncertainty - Subgaussian Bounds")
      } else {
        plot.w1 <- plot + geom_errorbar(data = toplot,
                                        aes(x = .data$Time, ymin = .data$lb.gaussian, ymax = .data$ub.gaussian), colour = col.synth,
                                        width = 0.5, linetype = 1) + ggtitle("In and Out of Sample Uncertainty - Subgaussian Bounds")
      }

      if (joint == TRUE && plot.type == "unit-time") {
        if (type == "treatment") {
          timebd <- min(toplot[toplot$Treatment == 1,]$Time) - anticipation
          plotdf <- subset(toplot, Time >= timebd)
        } else if (type == "series") {
          plotdf <- subset(toplot, Type == "Synthetic")
        }

        plot.w1 <- plot.w1 + geom_ribbon(data=plotdf, 
                                         aes(x=.data$Time, ymin=.data$lb.joint, ymax=.data$ub.joint), fill=col.synth, alpha=0.1)

      } else if (joint == TRUE && plot.type == "unit") {
        if (type == "treatment") {
          plotdf <- subset(toplot, Treatment == 1)
          plot.w1 <- plot.w1 + geom_errorbar(data=plotdf, 
                                             aes(x=.data$ID, ymin=.data$lb.joint, ymax=.data$ub.joint), 
                                             colour="darkred", width = 0.25)
          plot.w1$layers <- c(plot.w1$layers[[2]], plot.w1$layers[[4]], plot.w1$layers[[3]],plot.w1$layers[[1]])

        } else if (type == "series") {
          plotdf <- subset(toplot, Type == "Synthetic")
          plot.w1 <- plot.w1 + geom_errorbar(data=plotdf, aes(x=.data$Time, ymin=.data$lb.joint, ymax=.data$ub.joint), 
                                             colour="darkred", width = 0.6)
          plot.w1$layers <- c(plot.w1$layers[[3]], plot.w1$layers[[5]], plot.w1$layers[[4]], plot.w1$layers[[2]],plot.w1$layers[[1]])
        }
      }
    }

    if (e.out == TRUE && e.method %in% c("all", "ls")){
      if (type == "treatment" && plot.type == "unit") {
        plot.w2 <- plot + geom_errorbar(data = toplot,
                                        aes(x = .data$ID, ymin = .data$lb.ls, ymax = .data$ub.ls), colour = col.synth,
                                        width = 0.25, linetype = 1) + ggtitle("In and Out of Sample Uncertainty - Location-scale Model")
      } else {
        plot.w2 <- plot + geom_errorbar(data = toplot,
                                        aes(x = .data$Time, ymin = .data$lb.ls, ymax = .data$ub.ls), colour = col.synth,
                                        width = 0.5, linetype = 1) + ggtitle("In and Out of Sample Uncertainty - Location-scale Model")
      }

      if (joint == TRUE && plot.type == "unit-time") {
        if (type == "treatment") {
          timebd <- min(toplot[toplot$Treatment==1,]$Time) - anticipation
          plotdf <- subset(toplot, Time >= timebd)
        }
        else if (type == "series") plotdf <- subset(toplot, Type == "Synthetic")

        plot.w2 <- plot.w2 + geom_ribbon(data=plotdf, 
                                         aes(x=.data$Time, ymin=.data$lb.joint, ymax=.data$ub.joint), fill=col.synth, alpha=0.1)
      } else if (joint == TRUE && plot.type == "unit") {
        if (type == "treatment") {
          plotdf <- subset(toplot, Treatment == 1)

          plot.w2 <- plot.w2 + geom_errorbar(data=plotdf, 
                                             aes(x=.data$ID, ymin=.data$lb.joint, ymax=.data$ub.joint), 
                                             colour="darkred", width = 0.25)
          plot.w2$layers <- c(plot.w2$layers[[2]], plot.w2$layers[[4]], plot.w2$layers[[3]],plot.w2$layers[[1]])

        } else if (type == "series") {
          plotdf <- subset(toplot, Type == "Synthetic")
          plot.w2 <- plot.w2 + geom_errorbar(data=plotdf, aes(x=.data$Time, ymin=.data$lb.joint, ymax=.data$ub.joint), 
                                             colour="darkred", width = 0.6)
        }
      }
    }

    if (e.out == TRUE & e.method %in% c("all", "qreg")){

      if (type == "treatment" & plot.type == "unit") {
        plot.w3 <- plot + geom_errorbar(data = toplot,
                                      aes(x = .data$ID, ymin = .data$lb.qreg, ymax = .data$ub.qreg), colour = col.synth,
                                      width = 0.25, linetype = 1) + ggtitle("In and Out of Sample Uncertainty - Quantile Regression")
      } else {
        plot.w3 <- plot + geom_errorbar(data = toplot,
                                        aes(x = .data$Time, ymin = .data$lb.qreg, ymax = .data$ub.qreg), colour = col.synth,
                                        width = 0.5, linetype = 1) + ggtitle("In and Out of Sample Uncertainty - Quantile Regression")
      }

      if (joint == TRUE & plot.type == "unit-time") {
        if (type == "treatment") {
          timebd <- min(toplot[toplot$Treatment==1,]$Time) - anticipation
          plotdf <- subset(toplot, Time >= timebd)
        }
        else if (type == "series") plotdf <- subset(toplot, Type == "Synthetic")

        plot.w3 <- plot.w3 + geom_ribbon(data=plotdf, 
                                         aes(x=.data$Time, ymin=.data$lb.joint, ymax=.data$ub.joint), fill=col.synth, alpha=0.1)
      } else if (joint == TRUE & plot.type == "unit") {
        if (type == "treatment") {
          plotdf <- subset(toplot, Treatment == 1)
          plot.w3 <- plot.w3 + geom_errorbar(data=plotdf, 
                                             aes(x=.data$ID, ymin=.data$lb.joint, ymax=.data$ub.joint), 
                                             colour="darkred", width = 0.25)
          plot.w3$layers <- c(plot.w3$layers[[2]], plot.w3$layers[[4]], plot.w3$layers[[3]],plot.w3$layers[[1]])

        } else if (type == "series") {
          plotdf <- subset(toplot, Type == "Synthetic")
          plot.w3 <- plot.w3 + geom_errorbar(data=plotdf, aes(x=.data$Time, ymin=.data$lb.joint, ymax=.data$ub.joint), 
                                             colour="darkred", width = 0.6)
        }
      }
    }
  }

  if (class.type == 'scpi_scpi_multi') {
    if (e.out == FALSE) {plots <- list('plot_in'= plot.w)}
    else if (e.method == 'gaussian') {plots <- list('plot_out'= plot.w1)}
    else if (e.method == 'ls') {plots <- list('plot_out'= plot.w2)}    
    else if (e.method == 'qreg') {plots <- list('plot_out'= plot.w3)}   
    else if (e.method == 'all') {plots <- list('plot_out_gau'=plot.w1,
                                               'plot_out_ls'=plot.w2,
                                               'plot_out_qr'=plot.w3)}
  } else {
    plots <- list('plot_out' = plot)
  }

  ## Save data to reproduce plot
  if (!is.null(save.data)) {
    save(toplot, file = paste0(save.data, ".RData"))
  }

  return(plots)
}
