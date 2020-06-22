#' Example data for informative cluster size
#'
#' Simulated hypothetical clustered data created for illustration of functions in the \code{htest.clust} package.
#'
#' @name screen8
#' @docType data
#' @author Mary Gregg
#' @usage data(screen8)
#'
#' @format A data frame with 2224 rows and 12 columns:
#' \describe{
#'   \item{sch.id}{identification variable for school (clusters).}
#'   \item{stud.id}{identification variable for students within schools (observations within clusters).}
#'   \item{age}{student age in years.}
#'   \item{gender}{binary student gender.}
#'   \item{height}{student height in inches.}
#'   \item{weight}{student weight in lbs.}
#'   \item{math}{score from standardized math test.}
#'   \item{read}{score from standardized reading test.}
#'   \item{phq2}{ordinal (0-6) score from a mental health screening; higher scores correspond to higher levels of depression.}
#'   \item{qfit}{age-adjusted fitness quartile from physical health assessment at end of school year.}
#'   \item{qfit.s}{age-adjusted fitness quartile from physical health assessment at beginning of school year.}
#'   \item{activity}{student's primary after-school activity.}
#' }
#'
#' @details
#' Hypothetical data simulated for the following scenario.
#' An urban school district has collected demographic, biometric, and academic performance data
#' from graduating 8th grade students. \code{screen8} contains a sample of this data from 2224 students across
#' 73 schools. Student-level observations are clustered within schools.
#' The school district has implemented an incentive program in which schools with higher participation rates are
#' prioritized for classroom and technology upgrades. Cluster size could be informative in this data, as resource-poor
#' schools might have higher participation rates (larger cluster size), but also tend to have worse health metrics and
#' lower standardized test scores.
#'
#' @examples
#' data(screen8)
#' head(screen8)
#'
#' ## plot average math scores by cluster size
#' cl.size <- as.numeric(table(screen8$sch.id))
#' ave.math <- tapply(screen8$math, list(screen8$sch.id), mean)
#' plot(cl.size, ave.math)
#'
#' @keywords datasets
#'
"screen8"
