% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ttestClust.R
\name{ttestClust}
\alias{ttestClust}
\alias{ttestClust.default}
\alias{ttestClust.formula}
\title{Test of Marginal Means in Clustered Data}
\usage{
ttestClust(x, ...)

\method{ttestClust}{default}(
  x,
  y = NULL,
  idx,
  idy = NULL,
  alternative = c("two.sided", "less", "greater"),
  mu = 0,
  paired = FALSE,
  conf.level = 0.95,
  ...
)

\method{ttestClust}{formula}(formula, id, data, subset, na.action, ...)
}
\arguments{
\item{x, y}{numeric vectors of data values.}

\item{...}{further arguments to be passed to or from methods.}

\item{idx}{vector or factor object denoting cluster membership for \code{x} observations (or cluster
membership for paired observations when \code{paired} is \code{TRUE}). Length must be equal
to length of \code{x}.}

\item{idy}{vector or factor object denoting cluster membership for \code{y} observations. Length must be equal
to length of \code{y}}

\item{alternative}{indicates the alternative hypothesis and must be one of "\code{two.sided}", "\code{greater}",
or "\code{less}".You can specify just the initial letter.}

\item{mu}{a number specifying an optional parameter used to form the null hypothesis.}

\item{paired}{a logical indicating whether \code{x} and \code{y} are paired.}

\item{conf.level}{confidence level of the interval.}

\item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} is a numeric variable giving the data values and
\code{rhs} a factor with two levels giving the corresponding groups.}

\item{id}{a vector or factor giving the corresponding cluster membership.}

\item{data}{an optional matrix or data frame containing variables in the formula \code{formula} and \code{id}.
By default the variables are taken from \code{environment(formula)}.}

\item{subset}{an optional vector specifying a subset of observations to be used.}

\item{na.action}{a function which indicates what should happen when data contain \code{NA}s. Defaults to
\code{getOption("na.action")}.}
}
\value{
A list with class "\code{htest}" containing the following compoments:
\item{statistic}{the value of the test statistic.}
\item{p.value}{the p-value of the test.}
\item{conf.int}{a confidence interval for the mean appropriate to the specified alternative hypothesis}
\item{estimate}{the estimated mean or difference in means, depending on whether it was a one-sample or two-sample test.}
\item{null.value}{the specified hypothesized value of the mean or mean difference.}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{a character string indicating what type of reweighted test of means was performed.}
\item{data.name}{a character string giving the name of the data and the total number of clusters.}
\item{M}{the number of clusters.}
}
\description{
Performs one and two sample tests of marignal means in clustered data, reweighted to correct
for potential cluster- or group-size informativeness.
}
\details{
The formula interface is only applicable for the 2-sample tests.

If \code{paired} is \code{TRUE} then \code{x}, \code{y}, and \code{idx} must be given and be of the same length.
\code{idy} is ignored.
}
\examples{
data(screen8)
## One sample test
## Test if marginal math scores are equal to 70
ttestClust(x=screen8$math, idx=screen8$sch.id, mu = 70)

## paired test
## Test is marginal math scores have equal mean to marginal reading scores
ttestClust(x=screen8$math, y=screen8$read, idx=screen8$sch.id, paired=TRUE)

## unpaired test
## Test if boys and girls have equal marginal math scores
boys <- subset(screen8, gender=='M')
girls <- subset(screen8, gender=='F')
ttestClust(x=boys$math, y=girls$math, idx=boys$sch.id, idy=girls$sch.id)

## unpaired test using formula method
ttestClust(math~gender, id=sch.id, data=screen8)

}
\references{
Gregg, M., Marginal methods and software for clustered data with cluster- and group-size informativeness.
PhD dissertation, University of Louisville, 2020.
}
