##########################################################
## REWEIGHTED T-TEST FOR CLUSTERED DATA
##
## From 't test CC.R'
##########################################################

#' Test of Marginal Means in Clustered Data
#'
#' Performs one and two sample tests of marignal means in clustered data, reweighted to correct
#' for potential cluster- or group-size informativeness.
#'
#' @param x,y numeric vectors of data values.
#' @param idx vector or factor object denoting cluster membership for \code{x} observations (or cluster
#' membership for paired observations when \code{paired} is \code{TRUE}). Length must be equal
#' to length of \code{x}.
#' @param idy vector or factor object denoting cluster membership for \code{y} observations. Length must be equal
#' to length of \code{y}
#' @param alternative indicates the alternative hypothesis and must be one of "\code{two.sided}", "\code{greater}",
#' or "\code{less}".You can specify just the initial letter.
#' @param mu a number specifying an optional parameter used to form the null hypothesis.
#' @param paired a logical indicating whether \code{x} and \code{y} are paired.
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
#' \item{conf.int}{a confidence interval for the mean appropriate to the specified alternative hypothesis}
#' \item{estimate}{the estimated mean or difference in means, depending on whether it was a one-sample or two-sample test.}
#' \item{null.value}{the specified hypothesized value of the mean or mean difference.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of reweighted test of means was performed.}
#' \item{data.name}{a character string giving the name of the data and the total number of clusters.}
#' \item{M}{the number of clusters.}
#' @references
#' Gregg, M., Marginal methods and software for clustered data with cluster- and group-size informativeness.
#' PhD dissertation, University of Louisville, 2020.
#'
#' @details The formula interface is only applicable for the 2-sample tests.
#'
#' If \code{paired} is \code{TRUE} then \code{x}, \code{y}, and \code{idx} must be given and be of the same length.
#' \code{idy} is ignored.
#'
#' @examples
#' data(screen8)
#' ## One sample test
#' ## Test if marginal math scores are equal to 70
#' ttestClust(x=screen8$math, idx=screen8$sch.id, mu = 70)
#'
#' ## paired test
#' ## Test is marginal math scores have equal mean to marginal reading scores
#' ttestClust(x=screen8$math, y=screen8$read, idx=screen8$sch.id, paired=TRUE)
#'
#' ## unpaired test
#' ## Test if boys and girls have equal marginal math scores
#' boys <- subset(screen8, gender=='M')
#' girls <- subset(screen8, gender=='F')
#' ttestClust(x=boys$math, y=girls$math, idx=boys$sch.id, idy=girls$sch.id)
#'
#' ## unpaired test using formula method
#' ttestClust(math~gender, id=sch.id, data=screen8)
#'
#' @export ttestClust
ttestClust <- function(x, ...) {
  UseMethod("ttestClust")
}

#' @rdname ttestClust
#' @export
ttestClust.default <- function(x, y=NULL, idx, idy=NULL, alternative = c("two.sided", "less", "greater"),
                                 mu = 0, paired = FALSE, conf.level = 0.95, ...) {
  alternative <- match.arg(alternative)
  alpha <- 1-conf.level
  if (!missing(mu) && (length(mu) != 1 ||  is.na(mu)))
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  if (!is.null(y)) {
    if (!is.numeric(y))
      stop("'y' must be numeric")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

    if (paired) { ## IF PAIRED, PREP DATA FOR ONE-SAMPLE TEST
      if (length(x) != length(y) | length(x) != length(idx))
        stop("'x', 'y', and 'idx' must have the same length")
      METHOD <- "Paired cluster-weighted test of means"
      OK <- stats::complete.cases(x, y, idx)
      x <- x[OK] - y[OK]
      idx <- idx[OK]
      y <- idy <- NULL
      M <- length(unique(idx))
    }
    else { ## TWO-SAMPLE TEST: USING ANALYTIC JK VARIANCE
      if(is.null(idy))
        stop("'idx' and 'idy' must both be provided for two sample test")
      METHOD <- "Two sample group-weighted test of means"
      ## Function to compute two-sample test
      jk.fun <- function(datset, M) {
        id.x <- unique(datset$id[datset$group==0])
        id.y <- unique(datset$id[datset$group==1])
        id.c <- intersect(id.x, id.y)
        ## clusters are all complete
        if (M == length(id.c)) {
          xbar.i <- stats::aggregate(datset$t[datset$group==0], list(datset$id[datset$group==0]), mean)[,2]
          ybar.i <- stats::aggregate(datset$t[datset$group==1], list(datset$id[datset$group==1]), mean)[,2]

          mx <- mean(xbar.i)
          my <- mean(ybar.i)
          V.T <- (M/(M-1)^2)*stats::var(xbar.i - ybar.i)
        }
        else {
          ## incomplete clusters
          dat.c <- datset[(datset$id)%in%id.c,]

          ## x only clusters
          idx.d <- setdiff(id.x, id.y)
          datx.d <- datset[(datset$id)%in%idx.d,]

          ## y only clusters
          idy.d <- setdiff(id.y, id.x)
          daty.d <- datset[(datset$id)%in%idy.d,]

          M.c <- length(id.c)
          M.x <- length(idx.d)
          M.y <- length(idy.d)

          ## Get xbar.i, ybar.i values
          xbar.i.c <- stats::aggregate(dat.c$t[dat.c$group==0], list(dat.c$id[dat.c$group==0]), mean)[,2]
          ybar.i.c <- stats::aggregate(dat.c$t[dat.c$group==1], list(dat.c$id[dat.c$group==1]), mean)[,2]

          xbar.i.d <- if(M.x!=0) {stats::aggregate(datx.d$t, list(datx.d$id), mean)[,2]} else {0}
          ybar.i.d <- if(M.y!=0) {stats::aggregate(daty.d$t, list(daty.d$id), mean)[,2]} else {0}

          xbar.c <- mean(xbar.i.c)
          ybar.c <- mean(ybar.i.c)
          xbar.d <- mean(xbar.i.d)
          ybar.d <- mean(ybar.i.d)

          ## Denominators
          D1 <- ((M.c-1)+2*M.x)
          D2 <- ((M.c-1)+2*M.y)
          D.x1 <- (M.c+2*(M.x-1))
          D.x2 <- (M.c+2*M.y)
          D.y1 <- (M.c + 2*M.x)
          D.y2 <- (M.c+2*(M.y-1))

          sum.comp <- (M.c/M)*((2*M.x*mean(xbar.i.d))/(D1) - (2*M.y*mean(ybar.i.d))/(D2) + (M.c-1)*(mean(xbar.i.c)/D1 - mean(ybar.i.c)/D2))
          sum.xonly <- (M.x/M)*(M.c*(mean(xbar.i.c)/D.x1 - mean(ybar.i.c)/D.x2) + 2*((mean(xbar.i.d)*(M.x-1))/D.x1 - (M.y*mean(ybar.i.d))/D.x2))
          sum.yonly <- (M.y/M)*(M.c*(mean(xbar.i.c)/D.y1 - mean(ybar.i.c)/D.y2) + 2*(M.x*mean(xbar.i.d)/D.y1 - (mean(ybar.i.d)*(M.y-1))/D.y2))

          ## K values
          K.c <- 2*((M.x*mean(xbar.i.d))/D1 - (M.y*mean(ybar.i.d))/D2)*((M.c-M)/M) +
            M.c*(mean(xbar.i.c)/D1 - mean(ybar.i.c)/D2)*((M.c-M-1)/M) + sum.xonly + sum.yonly

          K.x <- ((M.x-M)/M)*(M.c*(xbar.c/D.x1 - ybar.c/D.x2)-(2*M.y*ybar.d)/D.x2) +
            ((2*xbar.d*M.x)/D.x1)*((M.x-M-1)/M) + sum.comp + sum.yonly

          K.y <- ((M.y-M)/M)*(M.c*(xbar.c/D.y1 - ybar.c/D.y2) + (2*M.x*xbar.d)/D.y1) -
            ((M.y-M-1)/M)*((2*M.y*ybar.d)/D.y2) + sum.comp + sum.xonly

          ## Pieces of V(T)
          V.comp <- sum((xbar.i.c/D1 - ybar.i.c/D2)^2) + 2*K.c*((M.c*xbar.c)/D1 - (M.c*ybar.c)/D2) + M.c*(K.c^2)
          V.x <- 4*sum((xbar.i.d/D.x1)^2) + ((4*K.x)/D.x1)*(M.x*xbar.d) + M.x*(K.x^2)
          V.y <- 4*sum((ybar.i.d/D.y2)^2) - ((4*K.y)/D.y2)*(M.y*ybar.d) + M.y*(K.y^2)
          V.T <- (M/(M-1))*(V.comp + V.x + V.y)

          ## statistic
          mx <- (1/(M.c + 2*M.x))*(M.c*mean(xbar.i.c) + 2*M.x*mean(xbar.i.d))
          my <- (1/(M.c + 2*M.y))*(M.c*mean(ybar.i.c) + 2*M.y*mean(ybar.i.d))
        }
        stat <- mx - my
        est <- c(mx, my)
        RVAL <- list(statistic=stat, variance=V.T, estimate=est)
        return(RVAL)
      }

      ## Get data in correct format to use 'jk.fun'
      ## (must have column names "id", "t", "group", and group must be 0/1)
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

      ## Apply function that returns statistic and variance
      cw.vals <- jk.fun(cdat, M=M)
      DEST <- cw.vals$statistic
      var.hat <- cw.vals$variance
      STATISTIC <- (DEST - mu)/sqrt(var.hat)
      CINT <- switch(alternative, less = c(-Inf, DEST + stats::qnorm(1-alpha)*sqrt(var.hat)),
                     greater = c(DEST - stats::qnorm(1-alpha)*sqrt(var.hat)),
                     two.sided = c(DEST - stats::qnorm(1-alpha/2)*sqrt(var.hat),
                                   DEST + stats::qnorm(1-alpha/2)*sqrt(var.hat)))
      ESTIMATE <- cw.vals$est
      names(ESTIMATE) <- c("weighted mean of x", "weighted mean of y")
    }
  }
  else { ## CLEAN ONE-SAMPLE NON-PAIRED DATA
    METHOD <- "One sample cluster-weighted test of means"
    DNAME <- deparse(substitute(x))
    if (paired)
      stop("'y' is missing for paired test")
    ok <- stats::complete.cases(x, idx)
    x <- x[ok]
    idx <- idx[ok]
    M <- length(unique(idx))
  }
  if (is.null(y)) { ## APPLY ONE-SAMPLE/PAIRED TEST
    xbar.i <- stats::aggregate(x, list(idx), mean)[,2]
    mx <- mean(xbar.i)
    var.hat <- mean((xbar.i-mu)^2)/M
    STATISTIC <- (mx-mu)/sqrt(var.hat)
    CINT <- switch(alternative, less = c(-Inf, mx + stats::qnorm(1-alpha)*sqrt(var.hat)),
                   greater = c(mx - stats::qnorm(1-alpha)*sqrt(var.hat), Inf),
                   two.sided = c(mx - stats::qnorm(1-alpha/2)*sqrt(var.hat),
                                 mx + stats::qnorm(1-alpha/2)*sqrt(var.hat)))
    ESTIMATE <- stats::setNames(mx, if (paired)
      "cluster-weighted mean of the differences"
      else "cluster-weighted mean of x")
  }

  PVAL <- switch(alternative, less = stats::pnorm(STATISTIC), greater = stats::pnorm(STATISTIC, lower.tail = FALSE),
                 two.sided = 2*(1-stats::pnorm(abs(STATISTIC))))
  names(mu) <- if (paired || !is.null(y))
    "difference in means"
  else "mean"
  attr(CINT, "conf.level") <- conf.level
  names(STATISTIC) <- "z"
  DNAME <- paste0(paste0(DNAME, ", M = "), as.character(M))
  RVAL <- list(statistic = STATISTIC, p.value = PVAL,
               conf.int = CINT, estimate = ESTIMATE , null.value = mu,
               alternative = alternative, method = METHOD, data.name = DNAME, M=M)
  class(RVAL) <- "htest"
  if (M < 30)
    warning('Number of clusters < 30. Normal approximation may be incorrect')
  return(RVAL)
}



#' @rdname ttestClust
#' @export
ttestClust.formula <- function (formula, id, data, subset, na.action, ...)
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
  y <- do.call("ttestClust", args=c(list(DATA$x$r, DATA$y$r, DATA$x$id, DATA$y$id,...)))
  #y$data.name <- DNAME
  y$data.name <- paste0(paste0(DNAME, ", M = "), as.character(length(unique(mf$id))))
  names(y$estimate) <- paste("weighted mean in group", levels(g))
  y
}

