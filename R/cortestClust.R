####################################
## Clustered correlation function
####################################

#' Test for Marginal Association Between Paired Clustered Data
#'
#' Test for association between paired samples in clustered data with potentially
#' informative cluster size.
#'
#' @param x,y numeric vectors of data values.
#' @param id a vector or factor object which identifies the clusters. The length of \code{id} should be the
#' same as the number of observations.
#' @param method a character string indicating which correlation coefficient is to be used for the test.
#' One of "\code{pearson}", "\code{kendall}", or "\code{spearman}". Can be abbreviated.
#' @param alternative indicates the alternative hypothesis and must be one of "\code{two.sided}", "\code{greater}",
#' or "\code{less}".You can specify just the initial letter. "\code{greater}" corresponds to positive association,
#' "\code{less}" to negative association.
#' @param conf.level confidence level for the returned confidence interval.
#' @param formula a formula of the form ~ u + v, where each of \code{u} and \code{v} are numeric variables
#' giving the data values for one sample. The samples must be of the sample length.
#' @param data an optional matrix or data frame containing variables in the formula \code{formula}. By default the
#' variables are taken from \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param na.action a function which indicates what should happen when data contain \code{NA}s. Defaults to
#' \code{getOption("na.action")}.
#' @param ... further arguments to be passed to or from methods.
#' @return A list with class "\code{htest}" containing the following compoments:
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{estimate}{the estimated measure of marginal association, with name "\code{cluster-weighted cor}",
#' "\code{cluster-weighted tau}", or "\code{cluster-weighted rho}" corresponding to the method employed.}
#' \item{null.value}{the value of the assocation measure under the null hypothesis, always 0.}
#' \item{conf.int}{a confidence interval for the measure of association.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating how the association was measured.}
#' \item{data.name}{a character string giving the name(s) of the data and the total number of clusters.}
#' \item{M}{the number of clusters.}
#' @references
#' Lorenz, D., Datta, S., Harkema, S. (2011) Marginal association measures for clustered data.
#' \emph{Statistics in Medicine}, \bold{30}, 3181--3191.
#'
#' Lorenz, D., Levy, S., Datta, S. (2018) Inferring marginal association with paired and unpaired
#' clustered data. \emph{Stat. Methods Med. Res.}, \bold{27}, 1806--1817.
#' @details The three methods each estimate the marginal association between paired observatons from
#' clustered data and compute a test of the value being zero.
#'
#' If \code{method} is "\code{pearson}" ("\code{kendall}"), the test statistic is based on the
#' Pearson product-moment (Kendall concordance coefficient) analog of Lorenz \emph{et al.} (2011).
#'
#' If \code{method} is "\code{spearman}", the test statistic
#' is based on the Spearman coefficient analog of Lorenz \emph{et al.} (2018) modified for paired data.
#'
#' @examples
#' data(screen8)
#' cortestClust(screen8$read, screen8$math, screen8$sch.id)
#'
#' ## Formula interface.
#' cortestClust(~ math + read, sch.id, data=screen8, method="kendall")
#'
#' @export
cortestClust <- function(x, ...) {
  UseMethod("cortestClust")
}


#' @rdname cortestClust
#' @export
cortestClust.default <- function(x, y, id, method=c("pearson", "kendall", "spearman"),
                           alternative = c("two.sided", "less", "greater"), conf.level=.95,...) {
  alternative <- match.arg(alternative)
  M <- length(unique(id))
  method <- match.arg(method)
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  if (length(x) != length(y))
    stop("'x' and 'y' must have the same length")
  if (!is.numeric(x))
    stop("'x' must be a numeric vector")
  if (!is.numeric(y))
    stop("'y' must be a numeric vector")
  OK <- stats::complete.cases(x, y, id)
  x <- x[OK]
  y <- y[OK]
  id <- id[OK]
  id <- factor(as.factor(id), levels=unique(as.factor(id)))
  id.list <- unique(id)
  cl.size <- as.numeric(tapply(id, id, length))
  NVAL <- 0

  if (method == "pearson") {
    METHOD <- "Cluster-weighted Pearson's product-moment correlation"
    names(NVAL) <- "correlation"
    clmean.x <- tapply(x, id, mean)
    clmean.y <- tapply(y, id, mean)
    clmean.xy <- tapply(x*y, id, mean)
    clmean.xx <- tapply(x^2, id, mean)
    clmean.yy <- tapply(y^2, id, mean)
    est <- (mean(clmean.xy) - mean(clmean.x)*mean(clmean.y))/
      sqrt((mean(clmean.xx)-mean(clmean.x)^2)*(mean(clmean.yy)-mean(clmean.y)^2))
    moment.var <- stats::var(cbind(clmean.x, clmean.y, clmean.xy, clmean.xx, clmean.yy))
    grad <- c(est*mean(clmean.x)/(mean(clmean.xx)-mean(clmean.x)^2)-
                mean(clmean.y)/sqrt((mean(clmean.xx)-mean(clmean.x)^2)*(mean(clmean.yy)-mean(clmean.y)^2)),
              est*mean(clmean.y)/(mean(clmean.yy)-mean(clmean.y)^2)-
                mean(clmean.x)/sqrt((mean(clmean.xx)-mean(clmean.x)^2)*(mean(clmean.yy)-mean(clmean.y)^2)),
              1/sqrt((mean(clmean.xx)-mean(clmean.x)^2)*(mean(clmean.yy)-mean(clmean.y)^2)),
              -est/(2*(mean(clmean.xx)-mean(clmean.x)^2)), -est/(2*(mean(clmean.yy)-mean(clmean.y)^2)))
    var.est <- (1/M) * t(grad) %*% moment.var %*% grad
  }
  else {
    if (method == "kendall") {
      METHOD <- "Cluster-weighted Kendall's rank correlation tau"

      ## Required funtions for Kendall
      ecdf.clust <- function(x, q, cluster, cl.size) {
        temp <- outer(x, q, "<=")*1
        sum.x <- rowsum(temp, cluster, reorder=F)
        ecdf.x <- sum.x/cl.size
        colSums(ecdf.x)/(nrow(ecdf.x)+1)
      }
      hajek.var2 <- function(x, y, cluster, cl.size, M) {
        if (M<=500) {
          temp.less.x <- outer(x, x, "<")*1
          temp.greater.x <- 1-temp.less.x
          temp.less.y <- outer(y, y, "<")*1
          temp.greater.y <- 1-temp.less.y
          temp <- temp.less.x*temp.less.y + temp.greater.x*temp.greater.y
          temp.sum <- colSums(temp)/length(cluster)
        }
        else{
          half <- ceiling(M/2)
          temp.less.x1 <- outer(x, x[cluster%in%id.list[1:half]], "<")*1
          temp.greater.x1 <- 1-temp.less.x1
          temp.less.y1 <- outer(y, y[cluster%in%id.list[1:half]], "<")*1
          temp.greater.y1 <- 1-temp.less.y1
          temp1 <- temp.less.x1*temp.less.y1 + temp.greater.x1*temp.greater.y1
          rm(temp.less.x1, temp.greater.x1, temp.less.y1, temp.greater.y1)

          temp.less.x2 <- outer(x, x[cluster%in%id.list[(half+1):M]], "<")*1
          temp.greater.x2 <- 1-temp.less.x2
          temp.less.y2 <- outer(y, y[cluster%in%id.list[(half+1):M]], "<")*1
          temp.greater.y2 <- 1-temp.less.y2
          temp2 <- temp.less.x2*temp.less.y2 + temp.greater.x2*temp.greater.y2
          rm(temp.less.x2, temp.greater.x2, temp.less.y2, temp.greater.y2)
          temp.sum <- c(colSums(temp1), colSums(temp2))/length(cluster)
        }
        tapply(temp.sum, cluster, sum)/cl.size
      }
      ni.cl <- outer(1/cl.size, 1/cl.size)


      if (M<=500) {
        ecdf.x <- ecdf.clust(x, x, id, cl.size)
        ecdf.y <- ecdf.clust(y, y, id, cl.size)

        temp.x <- outer(ecdf.x, ecdf.x, function(x1, x2) sign(x1 - x2))
        temp.y <- outer(ecdf.y, ecdf.y, function(y1, y2) sign(y1 - y2))

        temp <- rowsum(temp.x*temp.y, id, reorder=F)
      }
      else {
        half <- ceiling(M/2)
        ecdf.x1 <- ecdf.clust(x, x[id%in%id.list[1:half]], id, cl.size)
        ecdf.x2 <- ecdf.clust(x, x[id%in%id.list[(half+1):M]], id, cl.size)
        ecdf.x <- c(ecdf.x1, ecdf.x2)
        rm(ecdf.x1, ecdf.x2)

        ## perform ecdf.clust for y, split into two
        ecdf.y1 <- ecdf.clust(y, y[id%in%id.list[1:half]], id, cl.size)
        ecdf.y2 <- ecdf.clust(y, y[id%in%id.list[(half+1):M]], id, cl.size)
        ecdf.y <- c(ecdf.y1, ecdf.y2)
        rm(ecdf.y1, ecdf.y2)

        ## outer of ecdf.x by two
        temp.x1 <- outer(ecdf.x, ecdf.x[id%in%id.list[1:half]], function(x1, x2) sign(x1-x2))
        temp.x2 <- outer(ecdf.x, ecdf.x[id%in%id.list[(half+1):M]], function(x1, x2) sign(x1-x2))

        temp.y1 <- outer(ecdf.y, ecdf.y[id%in%id.list[1:half]], function(x1, x2) sign(x1-x2))
        temp.y2 <- outer(ecdf.y, ecdf.y[id%in%id.list[(half+1):M]], function(x1, x2) sign(x1-x2))

        temp1 <- rowsum(temp.x1*temp.y1, id, reorder=F)
        temp2 <- rowsum(temp.x2*temp.y2, id, reorder=F)
        rm(temp.x1, temp.x2, temp.y1, temp.y2)
        temp <- cbind(temp1, temp2)
      }
      sum.cl <- t(rowsum(t(temp), id, reorder=F))
      diag(sum.cl) <- 0

      est <- sum(ni.cl*sum.cl)/(M*(M-1))
      var.est <- (16/M)*stats::var(hajek.var2(x, y, id, cl.size, M))
    }
    else {
      METHOD <- "Cluster-weighted Spearman's rank correlation rho"
      rank.clust <- function(z, cluster, ni, fnt = c("less", "lesseq"), dat=c("all", "half")) {
        fnt <- match.arg(fnt)
        dat <- match.arg(dat)
        fnt.in <- switch(fnt, less=paste("<"), lesseq=paste("<="))
        if (dat=="all") {
          temp <- 1*outer(z, z, fnt.in)
          return(colMeans(rowsum(temp,cluster)/ni))
        }
        else {
          m <- length(unique(cluster))
          cut <- ceiling(m/2)
          id.list <- unique(cluster)
          temp1 <- 1*outer(z, z[cluster%in%id.list[1:cut]], fnt.in)
          temp2 <- 1*outer(z, z[cluster%in%id.list[(cut+1):m]], fnt.in)
          return(c(colMeans(rowsum(temp1,cluster)/ni), colMeans(rowsum(temp2,cluster)/ni)))
        }
      }

      if (M<=500) {
        F.x <- rank.clust(x, id, cl.size, fnt="lesseq")
        F.xl <- rank.clust(x, id, cl.size, fnt="less")
        F.y <- rank.clust(y, id, cl.size, fnt="lesseq")
        F.yl <- rank.clust(y, id, cl.size, fnt="less")
      }
      else {
        F.x <- rank.clust(x, id, cl.size, fnt="lesseq", dat="half")
        F.xl <- rank.clust(x, id, cl.size, fnt="less", dat="half")
        F.y <- rank.clust(y, id, cl.size, fnt="lesseq", dat="half")
        F.yl <- rank.clust(y, id, cl.size, fnt="less", dat="half")
      }
      r.xij <- (F.x+F.xl)/2
      r.yij <- (F.y+F.yl)/2

      clmean.x <- tapply(r.xij, id, mean)
      clmean.y <- tapply(r.yij, id, mean)
      clmean.xy <- tapply(r.xij*r.yij, id, mean)
      clmean.xx <- tapply(r.xij^2, id, mean)
      clmean.yy <- tapply(r.yij^2, id, mean)
      est <- (mean(clmean.xy) - mean(clmean.x)*mean(clmean.y))/
        sqrt((mean(clmean.xx)-mean(clmean.x)^2)*(mean(clmean.yy)-mean(clmean.y)^2))
      moment.var <- stats::var(cbind(clmean.x, clmean.y, clmean.xy, clmean.xx, clmean.yy))
      grad <- c(est*mean(clmean.x)/(mean(clmean.xx)-mean(clmean.x)^2)-
                  mean(clmean.y)/sqrt((mean(clmean.xx)-mean(clmean.x)^2)*(mean(clmean.yy)-mean(clmean.y)^2)),
                est*mean(clmean.y)/(mean(clmean.yy)-mean(clmean.y)^2)-
                  mean(clmean.x)/sqrt((mean(clmean.xx)-mean(clmean.x)^2)*(mean(clmean.yy)-mean(clmean.y)^2)),
                1/sqrt((mean(clmean.xx)-mean(clmean.x)^2)*(mean(clmean.yy)-mean(clmean.y)^2)),
                -est/(2*(mean(clmean.xx)-mean(clmean.x)^2)), -est/(2*(mean(clmean.yy)-mean(clmean.y)^2)))
      var.est <- (1/M) * t(grad) %*% moment.var %*% grad

    }
  }

  z <- est/sqrt(var.est)
  if (!missing(conf.level) && (length(conf.level) !=
                               1 || !is.finite(conf.level) || conf.level < 0 ||
                               conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  alpha <- 1-conf.level
  cint <- switch(alternative, less = c(-1, est+stats::qnorm(conf.level)*sqrt(var.est)),
                 greater = c(est - stats::qnorm(conf.level)*sqrt(var.est), 1),
                 two.sided = cint <- c(est - stats::qnorm(1-alpha/2)*sqrt(var.est),
                                       est + stats::qnorm(1-alpha/2)*sqrt(var.est)))
  pval <- switch(alternative, less = stats::pnorm(z), greater = stats::pnorm(z, lower.tail = FALSE),
                 two.sided = 2*(1-stats::pnorm(abs(z))))

  names(est) <- switch(method, pearson="cluster-weighted cor", kendall="cluster-weighted tau",
                       spearman="cluster-weighted rho")
  names(z) <- "z"
  attr(cint, "conf.level") <- conf.level
  DNAME <- paste0(paste0(DNAME, ", M = "), as.character(M))

  RVAL <- list(statistic = z,
               p.value = pval, estimate = est, null.value = NVAL,
               conf.int = cint, alternative = alternative, method = METHOD,
               data.name = DNAME, M = M)
  class(RVAL) <- "htest"
  if (M < 30)
    warning('Number of clusters < 30. Normal approximation may be incorrect')
  return(RVAL)
}

#' @rdname cortestClust
#' @export
cortestClust.formula <- function (formula, id, data, subset, na.action, ...)
{
  if (missing(formula) || !inherits(formula, "formula") ||
      length(formula) != 2L)
    stop("'formula' missing or invalid")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, environment(formula))
  #if (length(mf) != 2L)
  #  stop("invalid formula")
  DNAME <- paste(names(mf)[1:2], collapse = " and ")
  names(mf) <- c("x", "y", "id")
  y <- do.call("cortestClust", c(mf, list(...)))
  y$data.name <- paste0(paste0(DNAME, ", M = "), as.character(length(unique(mf$id))))
  y
}
