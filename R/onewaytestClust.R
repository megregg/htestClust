##########################################################
## REWEIGHTED ANOVA TEST FOR CLUSTERED DATA
##
## From 'ANOVA CC.R'
##########################################################

#' Test for Equal Marignal Means in Clustered Data
#'
#' Test whether two or more intra-cluster groups have the same marginal means in clustered data. Reweighted to
#' correct for potential cluster- or group size informativeness.
#'
#' @param x  a two-dimensional matrix or data frame containing the within-cluster group means, where rows are the clusters
#' and columns are the group means.
#' @param formula a formula of the form \code{lhs} ~ \code{rhs}, where \code{lhs} is a numberic variable
#' giving the data values and \code{rhs} a numeric or factor with two or more levels giving the correponding groups.
#' @param id a vector or factor object denoting cluster membership.
#' @param data an optional matrix or data frame containing variables in the formula \code{formula} and \code{id}.
#' By default the variables are taken from \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param ... further arguments to be passed to or from methods.
#' @return A list with class "\code{htest}" containing the following compoments:
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{estimate}{the estimated marginal group means.}
#' \item{parameter}{the degrees of freedom of the chi square distribution.}
#' \item{method}{a character string indicating the test performed.}
#' \item{data.name}{a character string giving the name of the data and the total number of clusters.}
#' \item{M}{the number of clusters.}
#' @references
#' Gregg, M., Marginal methods and software for clustered data with cluster- and group-size informativeness.
#' PhD dissertation, University of Louisville, 2020.
#'
#'@details The null hypothesis is that all levels of \code{group} have equal marginal means.
#'
#'If \code{x} is a matrix or data frame, the dimension of \code{x} should be MxK, where M is the
#'number of clusters and K is the number of groups. Each row of \code{x} correponds to a cluster, where
#'the column values contain the respective group means from that cluster. Clusters which do not contain
#'observations in a particular group should have \code{NA} in the correponding column.
#'
#' @examples
#' data(screen8)
#' ## do average reading scores differ across after-school activities?
#' ## test using a table
#' read.tab <- tapply(screen8$read, list(screen8$sch.id, screen8$activity), mean)
#' onewaytestClust(read.tab)
#'
#' ## test using formula method
#' onewaytestClust(read~activity, id=sch.id, data=screen8)
#'
#' @export
onewaytestClust <- function(x, ...) {
  UseMethod("onewaytestClust")
}

#' @rdname onewaytestClust
#' @export
onewaytestClust.default <- function(x, ...) {
  DNAME <- deparse(substitute(x))
  if (!is.matrix(x) & !is.data.frame(x))
    stop("'x' must be a matrix or data frame")
  k <- ncol(x)
  if(k < 2L)
    stop("not enough groups")
  ## create contrast matrix
  cfun <- function(pos) {
    1*(seq(1:k)==pos) - 1*(seq(1:k)==(pos+1))
  }
  cmat <- matrix(unlist(lapply(1:(k-1), cfun)), ncol=k, byrow=T)

  ## check if matrix has rows with all columns NA and remove
  rrm <- as.numeric(apply(x, 1, function(y) sum(1*is.na(y))))
  keep <- (rrm!=k)
  x <- x[keep,]
  M <- nrow(x)

  ## Internal functions
  wtgrp.fun <- function(dat){
    ## Get matrix of indicators for group
    grp.mat <- 1-1*is.na(dat)
    ## multiply respective group indicator column (e.g, column 1 for group 1) by the rowSums
    cl.typ.mat <- apply(grp.mat, 2, function(x) x*rowSums(grp.mat))
    ## Caclulate denominator values for each group
    D.grp <- apply(cl.typ.mat, 2, function(x) sum(table(x[x!=0])/as.numeric(names(table(x[x!=0])))))

    ## function to get the interior weighted sums
    ## lapply across the groups and divide by D.grp to get weighted group means
    sum.fun <- function(coln){
      ok <- stats::complete.cases(dat[,coln])
      sum(dat[ok,coln]/cl.typ.mat[ok,coln])
    }
    wt.means <- as.numeric(unlist(lapply(1:ncol(dat), sum.fun))/D.grp)
    return(wt.means)
  }
  ## Jackknife functions
  wtgrp.fun.single <- function(dat, grp){
    wtgrp.fun(dat)[grp]
  }

  jk.fun <- function(x, xdat, grp) {
    tdat <- xdat[x,]
    wtgrp.fun.single(tdat, grp)
  }

  ## Group means (appropriately weighted)
  theta.hat <- wtgrp.fun(x)

  ## Jackknifed variance
  jk.apply <- function(grp) {
    bootstrap::jackknife(1:M, jk.fun, x, grp)$jack.values
  }
  theta.hat.i <- matrix(unlist(lapply(1:k, jk.apply)), ncol=k, byrow=F)
  theta.bar <- colMeans(theta.hat.i)
  t.dmat <- t(apply(theta.hat.i, 1, function(x) x - theta.bar))
  t.dmat2 <- apply(t.dmat, 1, function(x) x%*%t(x))
  JK.var <- (M/(M - k))*(M-1)*matrix(apply(t.dmat2, 1, mean), nrow=k)

  ## Rest of the stuff
  STATISTIC <- t(cmat%*%t(t(theta.hat)))%*%MASS::ginv(cmat%*%JK.var%*%t(cmat))%*%(cmat%*%t(t(theta.hat)))
  PARAMETER <- k-1
  PVAL <- stats::pchisq(STATISTIC, df=PARAMETER, lower.tail=F)

  names(STATISTIC) <- "X-squared"
  names(PARAMETER) <- "df"
  ESTIMATE <- theta.hat
  names(ESTIMATE) <- colnames(x)
  METHOD <- "Reweighted one-way analysis of means for clustered data"
  DNAME <- paste0(paste0(DNAME, ", M = "), as.character(M))
  RVAL <- list(statistic = STATISTIC, parameter = PARAMETER, estimate=ESTIMATE,
               p.value = PVAL, method = METHOD, data.name = DNAME, M=M)
  class(RVAL) <- "htest"
  if (M < 30)
    warning('Number of clusters < 30. Normal approximation may be incorrect')
  return(RVAL)
}


#' @rdname onewaytestClust
#' @export
onewaytestClust.formula <- function (formula, id, data, subset,...)
{
  if (missing(formula) || (length(formula) != 3L) || (length(attr(stats::terms(formula[-2L]),
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  dp <- as.character(formula)
  if (length(dp) != 3L)
    stop("a two-sided formula is required")
  DNAME <- paste(dp[[2L]], "and", dp[[3L]])
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  mf <- eval(m, parent.frame())
  names(mf) <- c("r", "group", "id")
  mf$g <- factor(mf$group)
  k <- nlevels(mf$g)
  if (k < 2L)
    stop("not enough groups")
  df <- tapply(mf$r, list(mf$id, mf$g), mean)
  y <- do.call("onewaytestClust", args=list(df))
  y$data.name <- paste0(paste0(DNAME, ", M = "), as.character(length(unique(mf$id))))
  y
}
