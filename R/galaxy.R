#' @import plotrix 
#' @import mvmeta 
#' @title Galaxy Plot: A New Visualization Tool of Bivariate Meta-Analysis Studies
#' @description  A new visualization method that simultaneously presents the effect sizes of bivariate outcomes and their standard errors in a two-dimensional space.
#' @author Chuan Hong, Chongliang Luo, Yong Chen
#' @usage galaxy(data, y1, s1, y2, s2, scale1, scale2, scale.adj, 
#'               corr, group, study.label, annotate, xlab, ylab, main, legend.pos)
#' 
#' @param data dataset with at least 4 columns for the effect sizes of the two outcomes and their standard errors
#' @param y1 column name for outcome 1, default is 'y1'
#' @param s1 column name for standard error of \code{y1}, default is 's1'
#' @param y2 column name for outcome 2, default is 'y2'
#' @param s2 column name for standard error of \code{y2}, default is 's2' 
#' @param scale1  parameter for the length of the cross hair: the ellipse width is scale1 / s1 * scale.adj
#' @param scale2  parameter for the length of the cross hair: the ellipse height is scale2 / s2 * scale.adj 
#' @param scale.adj a pre-specified parameter to adjust for \code{scale1} and \code{scale2}
#' @param corr column name for within-study correlation 
#' @param group column name for study group
#' @param study.label  column name for study label
#' @param annotate  logical specifying whether study label should be added to the plot, default is FALSE.
#' @param xlab   x axis label, default \code{y1}
#' @param ylab   y axis label, default \code{y2}
#' @param main   main title 
#' @param legend.pos  The position of the legend for study groups if \code{group} is specified, see \code{legend}, default is 'bottomright'.
#' @details This function returns the galaxy plot to visualize bivariate meta-analysis data, 
#' which faithfully retains the information in two separate funnel plots, while providing 
#' useful insights into outcome correlations, between-study heterogeneity and joint asymmetry. 
#' Galaxy plot: a new visualization tool of bivariate meta-analysis studies.
#' Funnel plots have been widely used to detect small study effects in the results of 
#' univariate meta-analyses. However, there is no existing visualization tool that is 
#' the counterpart of the funnel plot in the multivariate setting. We propose a new 
#' visualization method, the galaxy plot, which can simultaneously present the effect sizes 
#' of bivariate outcomes and their standard errors in a two-dimensional space. 
#' The galaxy plot is an intuitive visualization tool that can aid in interpretation 
#' of results of multivariate meta-analysis. It preserves all of the information presented 
#' by separate funnel plots for each outcome while elucidating more complex features that 
#' may only be revealed by examining the joint distribution of the bivariate outcomes.
#' @return NULL
#'          
#' @references Hong, C., Duan, R., Zeng, L., Hubbard, R., Lumley, T., Riley, R., Chu, H., Kimmel, S., and Chen, Y. (2020) 
#' Galaxy Plot: A New Visualization Tool of Bivariate Meta-Analysis Studies, American Journal of Epidemiology, https://doi.org/10.1093/aje/kwz286.
#' @examples 
#' data(sim_dat)
#' galaxy(data=sim_dat, scale.adj = 0.9, corr = 'corr', group = 'subgroup', 
#'         study.label = 'study.id', annotate = TRUE, main = 'galaxy plot')
#' 
#' @export
galaxy <- function(data, y1='y1', s1='s1', y2='y2', s2='s2', scale1, scale2, scale.adj=1, 
                   corr=NULL, group=NULL, study.label=NULL, annotate=F, 
                   xlab, ylab, main, legend.pos='bottomright'){
  if(!y1 %in% names(data)) stop('column y1 = ', y1, ' not found, please specify it!')
  if(!s1 %in% names(data)) stop('column s1 = ', s1, ' not found, please specify it!')
  if(!y2 %in% names(data)) stop('column y2 = ', y2, ' not found, please specify it!')
  if(!s2 %in% names(data)) stop('column s2 = ', s2, ' not found, please specify it!')
  names(data)[names(data)==y1] <- 'y1' 
  names(data)[names(data)==s1] <- 's1' 
  names(data)[names(data)==y2] <- 'y2' 
  names(data)[names(data)==s2] <- 's2' 
  n.remove <- sum(is.na(data$y1) | is.na(data$s1) | is.na(data$y2) | is.na(data$s2) | data$s1==0 | data$s2==0)
  if(n.remove > 0) warning(n.remove, ' studies removed due to NA or zero se!' )
  data <- data[!(is.na(data$y1) | is.na(data$s1) | is.na(data$y2) | is.na(data$s2) | data$s1==0 | data$s2==0),]
  n <- nrow(data)
  
  ## within-study correlation
  if(!is.null(corr)){
    names(data)[names(data) == corr] <- 'corr'
    angle1 <- asin(data$corr)
  }
  ## group variable
  if(!is.null(group)){
    names(data)[names(data) == group] <- 'group'
    data$group <- factor(data$group)
  }
  ## study label
  if(!is.null(study.label)){
    names(data)[names(data) == study.label] <- 'study.label' 
  }
  ## scale of the ellipses for y1 and y2
  if(missing(scale1)){
    scale1 <- (max(data$y1) - min(data$y1) ) / n * quantile(data$s1, 0.2)
  } 
  if(missing(scale2)){
    scale2 <- (max(data$y2) - min(data$y2) ) / n * quantile(data$s2, 0.2)
  }
  scale1 <- scale1 * scale.adj
  scale2 <- scale2 * scale.adj
  message('scale of the ellipses are: scale1 = ', round(scale1,4), ', scale2 = ', round(scale2,4))
  
  ## galaxy center
  junk1 = mvmeta(data$y1, data$s1^2, method = "fixed")
  junk2 = mvmeta(data$y2, data$s2^2, method = "fixed")
  b = c(junk1$coefficients, junk2$coefficients)
   
  ## endpoints of the hairs 
  ax0 <- data$y1 - scale1 * 1.2/data$s1
  ax1 <- data$y1 + scale1 * 1.2/data$s1
  ay0 <- data$y2 - scale2 * 1.2/data$s2
  ay1 <- data$y2 + scale2 * 1.2/data$s2
  ## span of x and y axis
  dx <- max(ax1) - min(ax0)
  dy <- max(ay1) - min(ay0)
  
  res = plot(# c(min(data$y1) - 1, max(data$y1) + 1), c(min(data$y2) - 1, max(data$y2) + 1), 
             c(min(ax0), max(ax1)) + 0.05 * c(-1, 1) * dx, c(min(ay0), max(ay1)) + 0.05 * c(-1, 1) * dy, 
             type = "n", xlab = ifelse(missing(xlab), y1, xlab), 
             ylab = ifelse(missing(ylab), y2, ylab), 
             cex.main = 2, main = ifelse(missing(main), '', main)) 
  draw.ellipse(data$y1, data$y2, scale1/data$s1, scale2/data$s2, col=as.numeric(data$group)+1)
  if(!is.null(group)){ # color fill ellipse for group of studies
    legend(legend.pos, pch=19, col=seq_along(levels(data$group))+1, legend = levels(data$group))
  }
  if(!is.null(corr)){  # arrow for within-study correlation
    ll <- sqrt((scale1/data$s1)^2+(scale2/data$s2)^2) / sqrt(dx^2+dy^2) * 0.7
    arrows(data$y1, data$y2, data$y1+ll*cos(angle1)*dx, data$y2+ll*sin(angle1)*dy, length = 0.05)
  }
  points(x = b[1], y = b[2], col = "red", pch = 8, cex = 2, lwd = 3)
  segments(data$y1, ay0, data$y1, ay1 )
  segments(ax0, data$y2, ax1, data$y2 )
  # segments(data$y1, data$y2, data$y1, (data$y2 + scale2 * 1.2/data$s2) )
  # segments(data$y1, data$y2, data$y1, (data$y2 - scale2 * 1.2/data$s2) )
  # segments(data$y1, data$y2, (data$y1 + scale1 * 1.2/data$s1), data$y2 )
  # segments(data$y1, data$y2, (data$y1 - scale1 * 1.2/data$s1), data$y2 )
  if(annotate==T) text(data$y1, data$y2, data$study.label, adj = c(0,1))  # annotate study id's
  
  return(res)
}
