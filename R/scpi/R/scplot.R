#' @title Plot Synthetic Control Point Estimates and Prediction Interval
#'
#' @description The command plots the actual pre-treatment and post-treatment series of the treated
#' unit and the estimated counterfactual synthetic control unit with corresponding prediction intervals.
#' Prediction intervals can take into account either in-sample uncertainty only or in-sample and
#' out-of-sample uncertainty using the techniques developed in \href{https://mdcattaneo.github.io/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, Feng, and Titiunik (2021)}. \code{\link{scpi}}.
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
#' 'scpi' object, obtained by calling \code{\link{scpi}}.
#' @param fig.path a string indicating the path where the plot(s) should be saved.
#' @param fig.name a string indicating the name of the plot(s). If multiple plots will be saved the command automatically
#' generates a numeric suffix to avoid overwriting them.
#' @param fig.format a string indicating the format in which the plot(s) should be saved.
#' @param e.out a logical specifying whether out-of-sample uncertainty should be included in the plot(s).
#' @param joint a logical specifying whether simultaneous prediction intervals should be included in the plot(s). It requires \code{e.out = TRUE}.
#' @param col.treated a string specifying the color for the treated unit series. Find the full list at \url{http://sape.inf.usi.ch/quick-reference/ggplot2/colour}.
#' @param col.synth a string specifying the color for the synthetic unit series. Find the full list at \url{http://sape.inf.usi.ch/quick-reference/ggplot2/colour}. 
#' @param label.xy a character list with two elements indicating the name of the axes
#' (eg. label.xy = list(x.lab = "Year", y.lab = "GDP growth (%)")).
#' @param plot.range a numeric array indicating the time range of the plot(s).
#' @param x.ticks a numeric list containing the location of the ticks on the x axis.
#' @param event.label a list containing a character object ('lab') indicating the label of the event and
#' a numeric object indicating the height of the label in the plot.
#' @param plot.specs a list containing some specifics to be passed to ggsave (eg. img.width, img.height, dpi)
#' @param save.data a character specifying the name and the path of the saved dataframe containing the processed data used
#' to produce the plot.
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
#' \item{\href{https://www.aeaweb.org/articles?id=10.1257/jel.20191450}{Abadie, A. (2021)}. Using synthetic controls:
#' Feasibility, data requirements, and methodological aspects.
#' \emph{Journal of Economic Literature}, 59(2), 391-425.}
#' \item{\href{https://mdcattaneo.github.io/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, M. D., Feng, Y., and Titiunik, R. 
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
#' data <- scpi_germany
#'
#' df <- scdata(df = data, id.var = "country", time.var = "year",
#'              outcome.var = "gdp", period.pre = (1960:1990),
#'              period.post = (1991:2003), unit.tr = "West Germany",
#'              unit.co = setdiff(unique(data$country), "West Germany"),
#'              constant = TRUE, cointegrated.data = TRUE)
#'
#' result <- scest(df, w.constr = list(name = "simplex", Q = 1))
#'
#' scplot(result)
#'
#' @export
#'
#'

scplot  <- function(result, fig.path = NULL, fig.name = NULL, fig.format = "png", e.out = TRUE, joint = FALSE,  
                    col.treated = "black", col.synth = "mediumblue", label.xy = NULL, plot.range = NULL, 
                    x.ticks = NULL, event.label = NULL, plot.specs = NULL, save.data = NULL) {

  if (!is.null(fig.path)) {
    if (is.character(fig.path) == FALSE) stop("The object 'fig.path' should be a character!")
  }

  if (!is.null(fig.name)) {
    if (is.character(fig.name) == FALSE) stop("The object 'fig.name' should be a character!")
  }

  if (is.character(fig.format) == FALSE) stop("The object 'fig.format' should be a character!")

  if (!is.null(save.data)) {
    if (is.character(save.data) == FALSE) stop("The object 'save.data' should be a character!")
  }

  if (is.null(label.xy) == FALSE){
    if (is.list(label.xy) == FALSE) {
      stop("label.xy should be a list of two elements! (eg. label.xy = list(x.lab = 'Year', y.lab = 'GDP')")
    }

    if (!all(names(label.xy) %in% c("x.lab", "y.lab"))) {
      stop("label.xy should be a list of two elements named 'x.lab' ad 'y.lab'!")
    }

    if (!all(sapply(c(label.xy$x.lab, label.xy$y.lab), is.character))) {
      stop("label.xy should be a list of two character elements!")
    }
  }

  if (is.null(label.xy) == TRUE) {
    x.lab <- "Time"
    y.lab <- "Outcome Variable"
  } else {
    x.lab <- label.xy$x.lab
    y.lab <- label.xy$y.lab
  }

  if (is.null(event.label) == FALSE){
    if (is.list(event.label) == FALSE) {
      stop("event.label should be a list of two elements! (eg. event.label = list(lab = 'Event', y.lab = 10)")
    }

    if (!all(names(event.label) %in% c("lab", "height"))) {
      stop("event.label should be a list of two elements named 'lab' ad 'height'!")
    }

    if (!all(is.character(event.label$lab), is.numeric(event.label$height))) {
      stop("event.label should be a list of two elements: one numeric ('height'), one character ('lab')!")
    }
  } 

  if (is.null(plot.specs) == FALSE){
    if (is.list(plot.specs) == FALSE) {
      stop("plot.specs should be a list of elements! (eg. plot.specs = list(img.height = 4.5, img.width = 6, dpi = 1000)")
    }

    if (!all(names(plot.specs) %in% c("img.height", "img.width", "dpi"))) {
      stop("plot.specs should be a list of elements named 'img.height', 'img.width', 'dpi'!")
    }

    img.width  <- plot.specs$img.width
    img.height <- plot.specs$img.height
    dpi        <- plot.specs$dpi
  } else {
    img.width  <- 6
    img.height <- 4.5
    dpi        <- 1000

  }

  # Some plot specifics used if the plot is saved
  if (any(is.null(fig.path), is.null(fig.name))) {
      save.plot  <- FALSE
  } else {
      save.plot  <- TRUE
  }

  if (result$data$specs$class.type == "scpi_scest") {    # Result comes from scest

    if (save.plot == TRUE) {
      if (!fig.format %in% c("eps","ps","tex","pdf","jpeg","tiff","png","bmp","svg","wmf")) {
        stop("The specified format is not valid. See ?ggsave to check valid formats.")
      } else {
        plot.fitted <- paste(fig.path,"/",fig.name,".",fig.format, sep = "")
      }
    }

    y.fit   <- rbind(result$est.results$Y.pre.fit, result$est.results$Y.post.fit)
    y.act   <- rbind(result$data$Y.pre, result$data$Y.post)


    period.pre  <- result$data$specs$period.pre
    period.post <- result$data$specs$period.post
    T0 <- period.pre[length(period.pre)] # intercept

    if (is.null(plot.range)) {
      plot.range <- c(period.pre, period.post)
    }


    # Actual data
    dat    <- data.frame(t     = c(period.pre, period.post),
                         Y.act = c(y.act),
                         sname = "Treated")

    # Fill with NAs Y.fit and confidence bounds where missing
    Y.fit.na  <- matrix(NA, nrow = length(c(period.pre, period.post)))
    names <- strsplit(rownames(y.fit), "\\.")
    rnames <- unlist(lapply(names, "[[", 2))
    not.missing.plot <- as.character(c(period.pre,period.post)) %in% rnames

    Y.fit.na[not.missing.plot, 1]  <- y.fit

    # Synthetic unit data
    dat.sc <- data.frame(t        = c(period.pre, period.post),
                         Y.sc     = Y.fit.na,
                         sname    = "SC Unit")

    if (is.null(x.ticks)) {
      x.ticks <- c(seq(plot.range[1], plot.range[length(plot.range)], length.out = 5), T0)
      x.ticks <- round(unique(x.ticks))
    }

    if (is.null(event.label)) {
      event.lab <- ""
      event.lab.height <- dat.sc$Y.sc[dat.sc$t == T0]
    } else {
      event.lab <- paste("\n", event.label$lab, sep = "")
      event.lab.height <- event.label$height
    }


    dat.plot    <- dat[dat[,'t'] %in% plot.range, ]
    dat.sc.plot <- dat.sc[dat.sc[,'t'] %in% plot.range, ]

    plotdf <- dplyr::left_join(dat.plot,dat.sc.plot, by = 't')

    ## Plot specs
    plot <- ggplot() + theme_bw() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  axis.text.x = element_text(vjust = 0.5, hjust=1)) +
            labs(x = x.lab, y = y.lab) +
            theme(legend.position = "bottom", legend.box = "horizontal", legend.title = element_blank(),
                  legend.background = element_rect(fill = "white", color = "black"))

    ## Series to plot
    plot  <- plot + geom_line( data = plotdf, aes(x = t, y = .data$Y.act, colour = .data$sname.x), linetype = 'solid') +
                    geom_point(data = plotdf, aes(x = t, y = .data$Y.act, colour = .data$sname.x), shape = 1) +
                    geom_line( data = plotdf, aes(x = t, y = .data$Y.sc,  colour = .data$sname.y), linetype = 'dashed') +
                    geom_point(data = plotdf, aes(x = t, y = .data$Y.sc,  colour = .data$sname.y), shape = 19) +
                    geom_vline(xintercept = T0, linetype = "dashed") +
                    geom_text(aes(x = T0, label = event.lab, y = event.lab.height), angle = 90, size = 4) +
                    scale_x_continuous(breaks = x.ticks) +
                    scale_color_manual(name = "", values = c(col.synth, col.treated),
                                       labels = c("Synthetic Control", "Treated"),
                                       guide = guide_legend(override.aes = list(
                                         linetype = c('dashed','solid'), shape = c(19, 1)))) +
                    ggtitle("Synthetic Control Prediction")

    if (save.plot == TRUE) {
      suppressWarnings(ggsave(filename = plot.fitted, plot = plot, width = img.width, height = img.height, dpi = dpi))
      cat("Plot saved at '", fig.path,"'\n", sep = "")
      cat("File name: ", plot.fitted,"\n", sep = "")
    }


    ## Return plot object to be modified by the user if needed
    plots <- list('plot_out' = plot)

  } else if (result$data$specs$class.type == "scpi_scpi") {    # Result comes from scpi

    e.method <- result$inference.results$e.method
    period.post <- result$data$specs$period.post
    
    if (save.plot == TRUE) {
      if (!fig.format %in% c("eps","ps","tex","pdf","jpeg","tiff","png","bmp","svg","wmf")) {
        stop("The specified format is not valid. See ?ggsave to check valid formats.")
      } else {
        plot.in.sample     <- paste(fig.path,"/",fig.name,".",fig.format, sep = "")

        if (e.method == "all") {
          plot.in.out.sample <- c(paste(fig.path,"/",fig.name,"_1.",fig.format, sep = ""),
                                  paste(fig.path,"/",fig.name,"_2.",fig.format, sep = ""),
                                  paste(fig.path,"/",fig.name,"_3.",fig.format, sep = ""))
        } else {
          plot.in.out.sample <- rep(paste(fig.path,"/",fig.name,".",fig.format, sep = ""), 3)
        }
      }
    }

    y.fit <- rbind(result$est.results$Y.pre.fit, result$est.results$Y.post.fit)

    sc.l.0  <- result$inference.results$CI.in.sample[, 1, drop = FALSE]
    sc.r.0  <- result$inference.results$CI.in.sample[, 2, drop = FALSE]

    if (e.method %in% c("all","gaussian")) {
      sc.l.1  <- result$inference.results$CI.all.gaussian[, 1, drop = FALSE]
      sc.r.1  <- result$inference.results$CI.all.gaussian[, 2, drop = FALSE]
    } else {
      sc.l.1 <- sc.r.1 <- rep(NA, length(period.post))
    }

    if (e.method %in% c("all","ls")) {
      sc.l.2  <- result$inference.results$CI.all.ls[, 1, drop = FALSE]
      sc.r.2  <- result$inference.results$CI.all.ls[, 2, drop = FALSE]
    } else {
      sc.l.2 <- sc.r.2 <- rep(NA, length(period.post))
    }

    if (e.method %in% c("all","qreg")) {
      sc.l.3  <- result$inference.results$CI.all.qreg[, 1, drop = FALSE]
      sc.r.3  <- result$inference.results$CI.all.qreg[, 2, drop = FALSE]
    } else {
      sc.l.3 <- sc.r.3 <- rep(NA, length(period.post))
    }

    y.act   <- rbind(result$data$Y.pre, result$data$Y.post)

    period.pre  <- result$data$specs$period.pre
    period.post <- result$data$specs$period.post

    T0 <- period.pre[length(period.pre)] # intercept

    if (is.null(plot.range)) {
      plot.range <- c(period.pre, period.post)
    }

    # Actual data
    dat    <- data.frame(t     = c(period.pre, period.post),
                         Y.act = c(y.act),
                         sname = "Treated")

    # Fill with NAs Y.fit and confidence bounds where missing
    Y.fit.na  <- matrix(NA, nrow = length(c(period.pre, period.post)))
    sc.l.0.na <- matrix(NA, nrow = length(c(period.pre, period.post)))
    sc.r.0.na <- matrix(NA, nrow = length(c(period.pre, period.post)))
    sc.l.1.na <- matrix(NA, nrow = length(c(period.pre, period.post)))
    sc.r.1.na <- matrix(NA, nrow = length(c(period.pre, period.post)))
    sc.l.2.na <- matrix(NA, nrow = length(c(period.pre, period.post)))
    sc.r.2.na <- matrix(NA, nrow = length(c(period.pre, period.post)))
    sc.l.3.na <- matrix(NA, nrow = length(c(period.pre, period.post)))
    sc.r.3.na <- matrix(NA, nrow = length(c(period.pre, period.post)))
    sc.l.j.na <- matrix(NA, nrow = length(c(period.pre, period.post)))
    sc.r.j.na <- matrix(NA, nrow = length(c(period.pre, period.post)))

    names <- strsplit(rownames(y.fit), "\\.")
    not.missing.plot <- c(period.pre,period.post) %in% unlist(lapply(names, "[[", 2))
    names <- strsplit(rownames(sc.l.0), "\\.")
    not.missing.ci   <- c(period.pre,period.post) %in% unlist(lapply(names, "[[", 2))

    Y.fit.na[not.missing.plot, 1]  <- y.fit
    sc.l.0.na[not.missing.ci, 1]   <- sc.l.0
    sc.r.0.na[not.missing.ci, 1]   <- sc.r.0
    sc.l.1.na[not.missing.ci, 1]   <- sc.l.1
    sc.r.1.na[not.missing.ci, 1]   <- sc.r.1
    sc.l.2.na[not.missing.ci, 1]   <- sc.l.2
    sc.r.2.na[not.missing.ci, 1]   <- sc.r.2
    sc.l.3.na[not.missing.ci, 1]   <- sc.l.3
    sc.r.3.na[not.missing.ci, 1]   <- sc.r.3
    sc.l.j.na[not.missing.ci, 1]   <- sc.l.0 + result$inference.results$bounds$joint[, 1]
    sc.r.j.na[not.missing.ci, 1]   <- sc.r.0 + result$inference.results$bounds$joint[, 2]

    # Synthetic unit data
    dat.sc <- data.frame(t        = c(period.pre, period.post),
                         Y.sc     = Y.fit.na,
                         lb0      = c(sc.l.0.na), ub0 = c(sc.r.0.na),
                         lb1      = c(sc.l.1.na), ub1 = c(sc.r.1.na),
                         lb2      = c(sc.l.2.na), ub2 = c(sc.r.2.na),
                         lb3      = c(sc.l.3.na), ub3 = c(sc.r.3.na),
                         lbj      = c(sc.l.j.na), ubj = c(sc.r.j.na),
                         sname    = "SC Unit")

    if (is.null(x.ticks)) {
      x.ticks <- c(seq(plot.range[1], plot.range[length(plot.range)], length.out = 5), T0)
      x.ticks <- round(unique(x.ticks))
    }

    if (is.null(event.label)) {
      event.lab <- ""
      event.lab.height <- dat.sc$Y.sc[dat.sc$t == T0]
    } else {
      event.lab <- paste("\n", event.label$lab, sep = "")
      event.lab.height <- event.label$height
    }

    dat.plot    <- subset(dat,    t %in% plot.range)
    dat.sc.plot <- subset(dat.sc, t %in% plot.range)

    plotdf <- dplyr::left_join(dat.plot, dat.sc.plot, by = 't')

    ## Plot specs
    plot <- ggplot() + theme_bw() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  axis.text.x = element_text(vjust = 0.5, hjust=1)) +
            labs(x = x.lab, y = y.lab) +
            theme(legend.position = "bottom", legend.box = "horizontal", legend.title = element_blank(),
                  legend.background = element_rect(fill = "white", color = "black"))

    ## Series to plot
    plot <- plot + 
            geom_line( data = plotdf, aes(x = t, y = .data$Y.act, colour = .data$sname.x), linetype = 'solid') +
            geom_point(data = plotdf, aes(x = t, y = .data$Y.act, colour = .data$sname.x), shape = 1) +
            geom_line( data = plotdf, aes(x = t, y = .data$Y.sc,  colour = .data$sname.y), linetype = 'dashed') +
            geom_point(data = plotdf, aes(x = t, y = .data$Y.sc,  colour = .data$sname.y), shape = 19) +
            geom_vline(xintercept = T0, linetype = "dashed") +
            geom_text(aes(x = T0, label = event.lab, y = event.lab.height), angle = 90, size = 4) +
            scale_x_continuous(breaks = x.ticks) + 
            scale_color_manual(name = "", values = c(col.synth, col.treated),
                               labels = c("Synthetic Control", "Treated"),
                               guide = guide_legend(override.aes = list(
                                 linetype = c('dashed','solid'), shape = c(19, 1))))
          
    
    if (e.out == FALSE) {   # Only in-sample uncertainty
      plot.w <- plot + geom_errorbar(data = plotdf,
                                     aes(x = .data$t, ymin = .data$lb0, ymax = .data$ub0, colour = .data$sname.y),
                                     width = 0.5, linetype = 'solid')  +  ggtitle("In-sample Uncertainty")
      
      if (save.plot == TRUE) {
        suppressWarnings(ggsave(filename = plot.in.sample, plot = plot.w, width = img.width, height = img.height, dpi = dpi))
        cat("Plot saved at '", fig.path,"'\n", sep = "")
        cat("File name: ", plot.in.sample,"\n", sep = "")
        }
    }

    if (e.out == TRUE & e.method %in% c("all", "gaussian")){
      plot.w1 <- plot + geom_errorbar(data = plotdf,
                                     aes(x = .data$t, ymin = .data$lb1, ymax = .data$ub1, colour = .data$sname.y),
                                     width = 0.5, linetype = 1) + ggtitle("In and Out of Sample Uncertainty - Subgaussian Bounds")
      if (joint == TRUE) {
        plot.w1 <- plot.w1 + geom_ribbon(data=plotdf, 
                                         aes(x=.data$t, ymin=.data$lbj, ymax=.data$ubj), fill=col.synth, alpha=0.1)
      }

      if (save.plot == TRUE) {
        suppressWarnings(ggsave(filename = plot.in.out.sample[1], plot = plot.w1, width = img.width, height = img.height, dpi = dpi))
        cat("Plot saved at '", fig.path,"'\n", sep = "")
        cat("File name: ", plot.in.out.sample[1],"\n", sep = "")
        }
    }



    if (e.out == TRUE & e.method %in% c("all", "ls")){
      plot.w2 <- plot + geom_errorbar(data = plotdf,
                                      aes(x = .data$t, ymin = .data$lb2, ymax = .data$ub2, colour = .data$sname.y),
                                      width = 0.5, linetype = 1) + ggtitle("In and Out of Sample Uncertainty - Location-scale Model")
      
      if (joint == TRUE) {
        plot.w2 <- plot.w2 + geom_ribbon(data=plotdf, 
                                         aes(x=.data$t, ymin=.data$lbj, ymax=.data$ubj), fill=col.synth, alpha=0.1)
      }

      if (save.plot == TRUE) {
        suppressWarnings(ggsave(filename = plot.in.out.sample[2], plot = plot.w2, width = img.width, height = img.height, dpi = dpi))
        cat("Plot saved at '", fig.path,"'\n", sep = "")
        cat("File name: ", plot.in.out.sample[2],"\n", sep = "")
      }
    }



    if (e.out == TRUE && e.method %in% c("all", "qreg")){
      plot.w3 <- plot + geom_errorbar(data = plotdf,
                                      aes(x = .data$t, ymin = .data$lb3, ymax = .data$ub3, colour = .data$sname.y),
                                      width = 0.5, linetype = 1) + ggtitle("In and Out of Sample Uncertainty - Quantile Regression")

      if (joint == TRUE) {
        plot.w3 <- plot.w3 + geom_ribbon(data=plotdf, 
                                         aes(x=.data$t, ymin=.data$lbj, ymax=.data$ubj), fill=col.synth, alpha=0.1)
      }
      
      if (save.plot == TRUE) {
        suppressWarnings(ggsave(filename = plot.in.out.sample[3], plot = plot.w3, width = img.width, height = img.height, dpi = dpi))
        cat("Plot saved at '", fig.path,"'\n", sep = "")
        cat("File name: ", plot.in.out.sample[3],"\n", sep = "")
        }
    }
  
    
    if (e.out == FALSE) plots = list('plot_in'= plot.w)
    else if (e.method == 'gaussian') plots = list('plot_out'= plot.w1)
    else if (e.method == 'ls') plots = list('plot_out'= plot.w2)
    else if (e.method == 'qreg') plots = list('plot_out'= plot.w3)
    else if (e.method == 'all') {
      plots = list('plot_out_gau'=plot.w1,
                   'plot_out_ls'=plot.w2,
                   'plot_out_qr'=plot.w3)
    }
  } else {
    stop("The object 'result' should be the output of scest or scpi!")
  }

  ## Save data to reproduce plot
  if (!is.null(save.data)) {
    save(plotdf, file = paste0(save.data, ".RData"))
  }

  return(plots)
}
