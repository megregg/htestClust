############################
## McNemar test for ICS
##
## From 'McNemar cleaned.R'
############################

#' Test of Marginal Homogeneity for Clustered Data
#'
#' Performs a test of marginal homogeneity of paired clustered data with potentially informative
#' cluster size.
#' @param x,y  vector or factor objects of equal length.
#' @param id a vector or factor object which identifies the clusters. The length of \code{id} must be the same
#' as the length of \code{x}.
#' @param variance character string specifying the method of variance estimation. Must be one of "\code{MoM}"
#' or "\code{emp}".
#' @return A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{parameter}{the degrees of freedom of the approxiate chi-squared distribution of the test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{method}{a character string indicating the test performed and which variance estimation method was used.}
#' \item{data.name}{a character string giving the name(s) of the data and the total number of clusters.}
#' \item{M}{the number of clusters.}
#' @references
#' Durkalski, V., Palesch, Y., Lipsitz, S., Rust, P. (2003). Analysis of clustered matched pair data.
#' \emph{Statistics in Medicine}, \bold{22}, 2417--2428.
#'
#' Gregg, M., Marginal methods and software for clustered data with cluster- and group-size informativeness.
#' PhD dissertation, University of Louisville, 2020.
#'
#' @details The null is that the marginal probabilities of being classified into cells \code{[i,j]} and
#' \code{[j,i]} are equal.
#'
#' Arguments \code{x}, \code{y}, and \code{id} must be vectors or factors of the same length.
#' Incomplete cases are removed.
#'
#' When \code{variance} is \code{MoM}, a method of moments variance estimate evaluated under the null is used.
#' This is equivalent to the test by Durkaslski \emph{et al.} (2003). When \code{variance} is \code{emp},
#' an emperical variance estimate is used. See Gregg (2020) for details.
#'
#' @examples
#' data(screen8)
#' ## Is marginal proportion of students in lowest fitness category
#' ## at the end of year equal to the beginning of year?
#' screen8$low.start <- 1*(screen8$qfit.s=='Q1')
#' screen8$low.end <- 1*(screen8$qfit=='Q1')
#' mcnemartestClust(screen8$low.start, screen8$low.end, screen8$sch.id)
#'
#' @export

mcnemartestClust <- function(x, y, id, variance = c("MoM", "emp")) {
  if (length(id) != length(x) | length(x) != length(y))
    stop("'id', 'x' and 'y' must have the same length")
  variance <- match.arg(variance)
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  OK <- stats::complete.cases(id, x, y)
  id <- as.factor(id[OK])
  x <- as.factor(x[OK])
  y <- as.factor(y[OK])
  r <- nlevels(x)
  m <- length(unique(id))
  DNAME <- paste(paste(DNAME,",", sep=""), paste("M =", as.character(m)))
  if ((r != 2) || (nlevels(y) != r))
    stop("'x' and 'y' must have 2 levels")
  if (!setequal(levels(x), levels(y)))
    stop("'x' and 'y' values must be equivocally named")

  cat1 <- levels(x)[1]
  cat2 <- levels(x)[2]
  ## Durkalski test
  if (variance=="MoM") {
    nk <- tapply(id, id, length)
    bk <- tapply(x==cat2 & y==cat1, id, sum)
    ck <- tapply(x==cat1 & y==cat2, id, sum)
    weighted.chisq <- sum((bk-ck)/nk)^2/sum(((bk-ck)/nk)^2)
    METHOD <- "Cluster-weighted test of marginal homogeneity with variance est: MoM"
  }  else {   ##My test
    I21.bar <- tapply((x==cat1 & y==cat2), id, mean)
    I12.bar <- tapply((x==cat2 & y==cat1), id, mean)
    z <- sum(I21.bar) - sum(I12.bar)
    #var.z <- length(unique(id))*(var(I21.bar) + var(I12.bar) - 2*cov(I21.bar, I12.bar))
    var.z <- m*stats::var(I21.bar-I12.bar)
    weighted.chisq <- (z^2)/var.z
    METHOD <- "Cluster-weighted test of marginal homogeneity with variance est: emp"
  }

  PARAMETER <- 1
  names(PARAMETER) <- "df"
  names(weighted.chisq) <- "X-squared"
  pval <- stats::pchisq(weighted.chisq, 1, lower.tail=F)
  RVAL <- list(statistic = weighted.chisq, parameter = PARAMETER,
               p.value = pval, method = METHOD, data.name = DNAME, M = m)
  class(RVAL) <- "htest"
  if (m < 30) {
    warning('Number of clusters < 30. Chi-squared approximation may be incorrect')
    return(RVAL)
  }
  return(RVAL)
}
