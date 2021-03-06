% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chisqtestClust.R
\name{chisqtestClust}
\alias{chisqtestClust}
\title{Chi-squared Test for Clustered Count Data}
\usage{
chisqtestClust(x, y = NULL, id, p = NULL,
             variance = c("MoM", "sand.null", "sand.est", "emp"))
}
\arguments{
\item{x}{a numeric vector or factor. Can also be a table or data frame.}

\item{y}{a numeric vector or factor of the same length as \code{x}. Ignored if \code{x} is a table or data frame.}

\item{id}{a numeric vector or factor which identifies the clusters; ignored if \code{x} is a table or data frame.
The length of \code{id} must be the same as the length of \code{x}.}

\item{p}{a vector of probabilities with length equal to the number of unique categories of \code{x}
if \code{x} is a vector, or equal to the number of columns of \code{x} if \code{x} is a table or data frame.}

\item{variance}{character string specifying the method of variance estimation. Must be one of "\code{sand.null}",
"\code{sand.est}", "\code{emp}", or "\code{MoM}".}
}
\value{
A list with class "\code{htest}" containing the following components:
\item{statistic}{the value of the test statistic.}
\item{parameter}{the degrees of freedom of the approximate chi-squared distrubtion of the test statistic.}
\item{p.value}{the p-value of the test.}
\item{method}{a character string indicating the test performed, and which variance estimation method was used.}
\item{data.name}{a character string giving the name(s) of the data and the total number of clusters.}
\item{M}{the number of clusters.}
\item{observed}{the observed reweighted proportions.}
\item{expected}{the expected proportions under the null hypothesis.}
}
\description{
\code{chisqtestClust} performs chi-squared contingency table tests and goodness-of-fit tests for
clustered data with potentially informative cluster size.
}
\details{
If \code{x} is 2-dimensional table or data frame, or if \code{x} is a vector or factor and \code{y}
is not given, then the cluster-weighted \emph{goodness-of-fit test} is performed. When \code{x} is a table
or data frame, the rows of \code{x} must give the aggregate category counts across the clusters. In this case,
the hypothesis tested is whether the marginal population probabilities equal those in \code{p}, or are all
equal if \code{p} is not given.

When \code{x}, \code{y}, and \code{id} are all given as vectors or factors, the cluster-weighted
\emph{chi-squared test of independence} is performed. The lengths of \code{x}, \code{y}, and \code{id}
must be equal. In this case, the hypothesis tested is that the joint probabilities of \code{x} and \code{y} are
equal to the product of the marginal probabilities.
}
\examples{
data(screen8)
## is the marginal extracurricular activity participation evenly distributed across categories?
## Goodness of Fit test using vectors.
chisqtestClust(x=screen8$activity, id=screen8$sch.id)

## Goodness of Fit test using table.
act.table <- table(screen8$sch.id, screen8$activity)
chisqtestClust(act.table)

## test if extracurricular activity participation and gender are independent
chisqtestClust(screen8$gender, screen8$activity, screen8$sch.id)

}
\references{
Gregg, M., Datta, S., Lorenz, D. (2020) Variance estimation in tests of clustered categorical data
with informative cluster size.
\emph{Statistical Methods in Medical Research}, doi:10.1177/0962280220928572.
}
