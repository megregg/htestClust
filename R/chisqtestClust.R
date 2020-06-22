####################################
## Chisq Tests for clustered data
##
## From 'Chisq CC.R'
####################################


#' Chi-squared Test for Clustered Count Data
#'
#' \code{chisqtestClust} performs chi-squared contingency table tests and goodness-of-fit tests for
#' clustered data with potentially informative cluster size.
#'
#' @usage chisqtestClust(x, y = NULL, id, p = NULL,
#'              variance = c("MoM", "sand.null", "sand.est", "emp"))
#'
#' @param x  a numeric vector or factor. Can also be a table or data frame.
#' @param y a numeric vector or factor of the same length as \code{x}. Ignored if \code{x} is a table or data frame.
#' @param id a numeric vector or factor which identifies the clusters; ignored if \code{x} is a table or data frame.
#' The length of \code{id} must be the same as the length of \code{x}.
#' @param p a vector of probabilities with length equal to the number of unique categories of \code{x}
#' if \code{x} is a vector, or equal to the number of columns of \code{x} if \code{x} is a table or data frame.
#' @param variance character string specifying the method of variance estimation. Must be one of "\code{sand.null}",
#' "\code{sand.est}", "\code{emp}", or "\code{MoM}".
#' @return A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{parameter}{the degrees of freedom of the approximate chi-squared distrubtion of the test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{method}{a character string indicating the test performed, and which variance estimation method was used.}
#' \item{data.name}{a character string giving the name(s) of the data and the total number of clusters.}
#' \item{M}{the number of clusters.}
#' \item{observed}{the observed reweighted proportions.}
#' \item{expected}{the expected proportions under the null hypothesis.}
#' @references
#' Gregg, M., Datta, S., Lorenz, D. (2020) Variance estimation in tests of clustered categorical data
#' with informative cluster size.
#' \emph{Statistical Methods in Medical Research}, doi:10.1177/0962280220928572.
#'
#' @details If \code{x} is 2-dimensional table or data frame, or if \code{x} is a vector or factor and \code{y}
#' is not given, then the cluster-weighted \emph{goodness-of-fit test} is performed. When \code{x} is a table
#' or data frame, the rows of \code{x} must give the aggregate category counts across the clusters. In this case,
#' the hypothesis tested is whether the marginal population probabilities equal those in \code{p}, or are all
#' equal if \code{p} is not given.
#'
#' When \code{x}, \code{y}, and \code{id} are all given as vectors or factors, the cluster-weighted
#' \emph{chi-squared test of independence} is performed. The lengths of \code{x}, \code{y}, and \code{id}
#' must be equal. In this case, the hypothesis tested is that the joint probabilities of \code{x} and \code{y} are
#' equal to the product of the marginal probabilities.
#'
#' @examples
#' data(screen8)
#' ## is the marginal extracurricular activity participation evenly distributed across categories?
#' ## Goodness of Fit test using vectors.
#' chisqtestClust(x=screen8$activity, id=screen8$sch.id)
#'
#' ## Goodness of Fit test using table.
#' act.table <- table(screen8$sch.id, screen8$activity)
#' chisqtestClust(act.table)
#'
#' ## test if extracurricular activity participation and gender are independent
#' chisqtestClust(screen8$gender, screen8$activity, screen8$sch.id)
#'
#' @export
chisqtestClust <- function(x, y=NULL, id, p = NULL,
                             variance = c("MoM", "sand.null", "sand.est", "emp")) {
  variance <- match.arg(variance)
  if (is.null(y)) { ## GOF TEST
    DNAME <- deparse(substitute(x))
    ##
    ## x is a table
    if (is.table(x)|is.data.frame(x)) {
      if (ncol(x) < 2L)
        stop("'x' must have at least 2 columns")
      x <- x[stats::complete.cases(x),]
      if (!(all(rowSums(x)>0)))
        stop("all clusters must have counts > 0")
      if (!(all(x)>=0))
        stop("elements of x must be nonnegative ")
      if (is.null(p))
        p <- rep(1/ncol(x), ncol(x))
      if (ncol(x)!=length(p))
        stop("length of p must equal the number of columns of x")
      tmpdat <- x
    }
    ##
    ## x is a vector
    else {
      if (length(x) != length(id))
        stop("'id' and 'x' must have the same length")
      OK <- stats::complete.cases(x, id)
      x <- factor(x[OK])
      id <- factor(id[OK])
      if (length(unique(x)) == 1L)
        stop("'x' must at least have 2 categories")
      if (is.null(p))
        p = rep(1/length(unique(x)), length(unique(x)))
      if (length(unique(x)) != length(p))
        stop("length of p must must equal number of categories of x")

      tmpdat <- table(id, x)
    }
    if (any(p < 0))
      stop("probabilities must be non-negative.")
    if (abs(sum(p) - 1) > sqrt(.Machine$double.eps))
      stop("probabilities must sum to 1.")
    M <- as.integer(nrow(tmpdat))
    cl.n <- apply(tmpdat, 1, sum)
    phat <- tmpdat/cl.n
    pbar <- apply(phat,2,mean)
    dmat <- phat - matrix(rep(p,M), nrow=M, byrow=T)
    dhat <- colMeans(dmat)

    if (variance=="emp") {
      var.hat <- stats::var(dmat)
    }
    else {
      U <- switch(variance,
                  MoM = t(t(dmat)/dhat),
                  sand.null = t(t(phat)/p),
                  sand.est = t(t(phat)/pbar))
      Vtmp <- apply(U, 1, function(x) x%*%t(x))
      V <- matrix(apply(Vtmp, 1, mean), nrow=length(p))
      Hinv <- switch(variance,
                     MoM = MASS::ginv(diag(-1/dhat)),
                     sand.null = MASS::ginv(diag(-pbar/p^2)),
                     sand.est = MASS::ginv(diag(-1/pbar)))
      var.hat <- Hinv%*%V%*%Hinv
    }
    OBSERVED <- pbar
    EXPECTED <- p
    names(EXPECTED) <- names(pbar)
    PARAMETER <- ncol(tmpdat) - 1
    STATISTIC <- M*t(dhat)%*%MASS::ginv(var.hat)%*%(dhat)
    PVAL <- stats::pchisq(STATISTIC, df=PARAMETER, lower.tail=F)
    DNAME <- paste0(paste0(DNAME, ", M = "), as.character(M))

    names(STATISTIC) <- "X-squared"
    names(PARAMETER) <- "df"
    METHOD <- paste0("Cluster-weighted chi-squared test for given probabilities with variance est: ", variance)

    if (any(M*p < 5) && is.finite(PARAMETER))
      warning("M*p < 5: Chi-squared approximation may be incorrect")

    RVAL <- list(statistic = STATISTIC, parameter = PARAMETER,
                 p.value = PVAL, method = METHOD, data.name = DNAME, observed = OBSERVED,
                 expected = EXPECTED, M=M)
  }
  ##
  ## INDEPENDENCE TEST
  else {
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    if (length(x) != length(y))
      stop("'x' and 'y' must have the same length")
    if (length(x) != length(id))
      stop("'id' must have the same length as x and y")
    OK <- stats::complete.cases(id,x,y)
    METHOD <- paste0("Cluster-weighted Chi-squared test of independence with variance est: ", variance)
    x <- factor(x[OK])
    y <- factor(y[OK])
    id <- factor(id[OK])
    M <- length(unique(id))
    if ((nlevels(x) < 2L) || (nlevels(y) < 2L))
      stop("'x' and 'y' must have at least 2 levels")

    tmp <- table(x, y, id)
    o <- prop.table(tmp,3)
    e <- apply(tmp,3,function(x) {
      outer(rowSums(x)/sum(x), colSums(x)/sum(x))
    })
    dim(e) <- dim(o)
    o.mat <- matrix(as.vector(o), nrow=M, byrow=T)
    e.mat <- matrix(as.vector(e), nrow=M, byrow=T)
    d.mat <- o.mat-e.mat
    o.hat <- colMeans(o.mat)
    e.hat <- colMeans(e.mat)
    d.hat <- colMeans(d.mat)

    OBSERVED <- matrix(o.hat, nrow=nlevels(x), byrow=F)
    EXPECTED <- matrix(e.hat, nrow=nlevels(x), byrow=F)
    rownames(OBSERVED) <- rownames(EXPECTED) <- levels(x)
    colnames(OBSERVED) <- colnames(EXPECTED) <- levels(y)

    if (variance=="emp") {
      var.hat <- stats::var(d.mat)
    } else {
      U <- switch(variance,
                  MoM = t(t(d.mat)/d.hat),
                  sand.null =  t(t(o.mat)/e.hat),
                  sand.est = t(t(o.mat)/o.hat))

      Vtmp <- apply(U, 1, function(x) x%*%t(x))
      V <- matrix(apply(Vtmp, 1, mean), nrow=(nlevels(y)*nlevels(x)))
      Hinv <- switch(variance,
                     MoM = MASS::ginv(diag(-1/d.hat)),
                     sand.null = MASS::ginv(diag(-o.hat/e.hat^2)),
                     sand.est = MASS::ginv(diag(-1/o.hat)))

      var.hat <- Hinv%*%V%*%Hinv
    }
    STATISTIC <- M*(t(d.hat) %*% MASS::ginv(var.hat) %*% d.hat)
    PARAMETER <- (nlevels(x)-1)*(nlevels(y)-1)
    PVAL <- stats::pchisq(STATISTIC, df=PARAMETER, lower=F)
    DNAME <- paste0(paste0(DNAME, ", M = "), as.character(M))
    names(STATISTIC) <- "X-squared"
    names(PARAMETER) <- "df"
    RVAL <- list(statistic = STATISTIC, parameter = PARAMETER,
                 p.value = PVAL, method = METHOD, data.name = DNAME, M = M, observed = OBSERVED, expected = EXPECTED)
  }
  DNAME <- paste0(paste0(DNAME, "M = "), as.character(M))
  class(RVAL) <- "htest"
  if (M < 30)
    warning('Number of clusters < 30. Chi-squared approximation may be incorrect')
  return(RVAL)
}
