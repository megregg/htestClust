##########################################################
## REWEIGHTED F TEST ANALOG FOR CLUSTERED DATA
##
## From 'Var test CC.R'
##########################################################

#' Reweighted Test to Compare Two Variances in Clustered Data
#'
#' Performs a reweighted test to compare marginal variances of intra-cluster groups in clustered data.
#' Reweighted to correct for potential cluster- or group-size informativeness.
#'
#' @param x,y numeric vectors of data values.
#' @param idx vector or factor object denoting cluster membership for \code{x} observations.
#' Length must be equal to length of \code{x}.
#' @param idy vector or factor object denoting cluster membership for \code{y} observations. Length must be equal
#' to length of \code{y}
#' @param difference the hypothesized difference of the marginal population variances of \code{x} and \code{y}.
#' @param alternative indicates the alternative hypothesis and must be one of "\code{two.sided}", "\code{greater}",
#' or "\code{less}".You can specify just the initial letter.
#' @param conf.level confidence level of the interval.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} is a numeric variable giving the data values and
#' \code{rhs} a factor with two levels giving the corresponding groups.
#' @param id a vector or factor giving the corresponding cluster membership.
#' @param data an optional matrix or data frame containing variables in the formula \code{formula} and \code{id}.
#' By default the variables are taken from \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param na.action a function which indicates what should happen when data contain \code{NA}s. Defaults to
#' \code{getOption("na.action")}.
#' @param ... further arguments to be passed to or from methods.
#' @return A list with class "\code{htest}" containing the following compoments:
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{conf.int}{a confidence interval for the difference of the population marginal variances.}
#' \item{estimate}{the difference in reweighted sample variances of \code{x} and \code{y}.}
#' \item{null.value}{the difference of population marginal variances under the null.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating the test performed.}
#' \item{data.name}{a character string giving the name of the data and the total number of clusters.}
#' \item{M}{the number of clusters.}
#' @references
#' Gregg, M., Marginal methods and software for clustered data with cluster- and group-size informativeness.
#' PhD dissertation, University of Louisville, 2020.
#'
#' @details The null hypothesis is that the difference of the marginal variances of the populations of
#' intra-cluster groups from which \code{x} and \code{y} were drawn is equal to \code{difference}.
#'
#' Using the default method, \code{difference} is the difference of the reweighted sample variances of \code{x}
#' and \code{y}. When using the formula method, the order of the difference is determined by the order of the
#' factor levels of \code{rhs}.
#'
#'
#' @examples
#' data(screen8)
#' boys <- subset(screen8, gender=='M')
#' girls <- subset(screen8, gender=='F')
#'
#' ## Do boys and girls have the same variability in math scores?
#' ## Test using vectors
#' vartestClust(x=boys$math, y=girls$math, idx=boys$sch.id, idy=girls$sch.id)
#'
#' ## Test using formula method.
#' vartestClust(math~gender, id=sch.id, data=screen8)
#'
#' ## Note that in this example, the sign of the estimate returned when using the formula
#' ## method is opposite to that when the test was performed using vectors. This is due to
#' ## the order of the gender factor levels
#'
#' @export
vartestClust <- function(x, ...) {
  UseMethod("vartestClust")
}

#' @rdname vartestClust
#' @export
vartestClust.default <- function(x, y, idx, idy, difference = 0, alternative = c("two.sided", "less", "greater"),
                                   conf.level = 0.95, ...) {
  METHOD <- "Reweighted test to compare two intra-cluster group variances"
  alternative <- match.arg(alternative)
  alpha <- 1-conf.level
  if (!missing(difference) && (length(difference) != 1 ||  is.na(difference)))
    stop("'difference' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  if (!is.numeric(x) || !is.numeric(y))
    stop("'x' and 'y' must be numeric")
  if(length(x)!=length(idx))
    stop("'x' and 'idx' must be of equal length")
  if(length(y)!=length(idy))
    stop("'y' and 'idy' must be of equal length")

  ## INTERNAL FUNCTIONS
  svec.fun <- function(dat.c){

    dat <- dat.c[,1:2]
    dat2 <- dat.c[,3:4]
    ## Get matrix of indicators for group
    grp.mat <- 1-1*is.na(dat)
    ## multiply respective group indicator column (e.g, column 1 for group 1) by the rowSums
    cl.typ.mat <- apply(grp.mat, 2, function(x) x*rowSums(grp.mat))

    ## Caclulate denominator values for each group
    D.grp <- apply(cl.typ.mat, 2, function(x) sum(table(x[x!=0])/as.numeric(names(table(x[x!=0])))))

    cl.typ.mat2 <- cbind(cl.typ.mat, cl.typ.mat)

    ## Calculate interior sum of group averages, average squared values
    sum.fun3 <- function(coln){
      ok <- stats::complete.cases(dat.c[,coln])
      sum(dat.c[ok,coln]/cl.typ.mat2[ok,coln])
    }
    return(as.numeric(unlist(lapply(1:ncol(dat.c), sum.fun3))/rep(D.grp,ncol(dat2))))
  }
  svec.fun.single <- function(dat, grp){
    svec.fun(dat)[grp]
  }
  jk.fun <- function(x, xdat, grp) {
    tdat <- xdat[x,]
    svec.fun.single(tdat, grp)
  }

  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  ## COLLECT DATA INTO DATA FRAME
  xok <- stats::complete.cases(x, idx)
  xdat <- data.frame(cbind(idx[xok], x[xok]))
  colnames(xdat) <- c("id", "t")
  xdat$group <- 0

  yok <- stats::complete.cases(y, idy)
  ydat <- data.frame(cbind(idy[yok], y[yok]))
  colnames(ydat) <- c("id", "t")
  ydat$group <- 1

  cdat <- rbind(xdat, ydat)
  cdat <- cdat[order(cdat$id),]
  M <- length(unique(cdat$id))

  cw.mat <- tapply(cdat$t, list(cdat$id, cdat$group), mean)
  cw.mat2 <- tapply(cdat$t, list(cdat$id, cdat$group), function(x) mean(x^2))
  s.hat <- svec.fun(cbind(cw.mat, cw.mat2))

  ## Jackknife variance
  jk.apply <- function(grp) {
    bootstrap::jackknife(1:M, jk.fun, cbind(cw.mat, cw.mat2), grp)$jack.values
  }
  theta.i <- matrix(unlist(lapply(1:4, jk.apply)), ncol=4, byrow=F)
  ivec.fun <- function(svec) {
    return((svec[3] - svec[1]^2) - (svec[4] - svec[2]^2))
  }
  Ivec.i <- apply(theta.i, 1, ivec.fun)
  Ivec.jk.bar <- mean(Ivec.i)
  var.F.CWd <- (M-1)*mean((Ivec.i - Ivec.jk.bar)^2)
  ESTIMATE <- ivec.fun(s.hat)
  STATISTIC <- (ESTIMATE - difference)/sqrt(var.F.CWd)
  names(STATISTIC) <- "z"
  CINT <- switch(alternative, less = c(-Inf, ESTIMATE + stats::qnorm(1-alpha)*sqrt(var.F.CWd)),
                 greater = c(ESTIMATE - stats::qnorm(1-alpha)*sqrt(var.F.CWd), Inf),
                 two.sided = c(ESTIMATE - stats::qnorm(1-alpha/2)*sqrt(var.F.CWd),
                               ESTIMATE + stats::qnorm(1-alpha/2)*sqrt(var.F.CWd)))
  attr(CINT, "conf.level") <- conf.level
  names(ESTIMATE) <- names(difference) <- "difference of variances"
  PVAL <- switch(alternative, less = stats::pnorm(STATISTIC), greater = stats::pnorm(STATISTIC, lower.tail = FALSE),
                 two.sided = 2*(1-stats::pnorm(abs(STATISTIC))))
  DNAME <- paste0(paste0(DNAME, ", M = "), as.character(M))
  RVAL <- list(statistic = STATISTIC, p.value = PVAL,
               conf.int = CINT, estimate = ESTIMATE , null.value = difference,
               alternative = alternative, method = METHOD, data.name = DNAME, M=M)
  class(RVAL) <- "htest"
  if (M < 30)
    warning('Number of clusters < 30. Normal approximation may be incorrect')
  return(RVAL)
}

#' @rdname vartestClust
#' @export
vartestClust.formula <- function (formula, id, data, subset, na.action, ...)
{
  if (missing(formula) || (length(formula) != 3L) || (length(attr(stats::terms(formula[-2L]),
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf)[1:2], collapse = " by ")
  names(mf) <- c("r", "group", "id")
  g <- factor(mf$group)
  if (nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  DATA <- stats::setNames(split(mf[,-2], g), c("x", "y"))
  y <- do.call("vartestClust", args=c(list(DATA$x$r, DATA$y$r, DATA$x$id, DATA$y$id,...)))
  y$data.name <- paste0(paste0(DNAME, ", M = "), as.character(length(unique(mf$id))))
  y
}
