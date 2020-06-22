############################
## Test for ICS
##
## From 'testICS clean.R'
############################

#' Test for Informative Cluster Size
#'
#' Performs a test for informative cluster size.
#' @param x  a vector of numeric responses. Can also be a data frame.
#' @param id a vector or factor object which identifies the clusters; ignored if \code{x} is a data frame.
#' The length of \code{id} must be the same as the length of \code{x}.
#' @param test.method character string specifying the method of construction for the test statistic.
#' Must be one of "\code{TF}" or "\code{TCM}".
#' @param B the number of bootstrap iterations.
#' @param print.it a logical indicating whether to print the progression of bootstrap iterations.
#' @return A list with class "\code{htest}" containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{method}{a character string indicating the test performed and the method of construction.}
#' \item{data.name}{a character string giving the name(s) of the data.}
#' @references
#' Nevalainen, J., Oja, H., Datta, S. (2017) Tests for informative cluster size using a novel balanced bootstrap scheme.
#' \emph{Statistics in Medicine}, \bold{36}, 2630--2640.
#'
#' @details The null is that the marginal distributions of the responses are independent of the cluster sizes.
#' A small p-value is evidence for the presence of informative cluster size.
#'
#' When \code{test.method = "TF"}, the test statistic is constructed based on differences between the null and
#' alternative distribution functions. "\code{TF}" is the suggested method when there are a large
#' number of unique cluster sizes and the number of clusters of each size is small. When \code{test.method = "TCM"},
#' the test statistic is a multisample Cramer von Mises-based test. This method is recommended
#' when there are a small number of possible cluster sizes. See Nevalainen \emph{et al.} (2017) for more details.
#'
#' When \code{x} is a data frame, the first column should contain values denoting cluster membership and
#' the second column the responses.
#'
#' This test is computationally intensive and can take significant time to execute. \code{print.it} defaults to
#' \code{TRUE} to identify the bootstrap progression.
#'
#' @examples
#' \dontrun{
#'   data(screen8)
#'   ## using vectors
#'   ## test if cluster size is related to math scores
#'   icstestClust(screen8$math, screen8$sch.id, B=100)
#'
#'   ## same test, but using a data frame and supressing iterations
#'   tdat <- data.frame(screen8$sch.id, screen8$math)
#'   icstestClust(tdat, B=100, print.it = F)
#' }
#' @export

icstestClust <- function(x, id, test.method=c('TF', 'TCM'), B=1000,
                     print.it=TRUE) {
  bi <- B
  if (bi<=0)
    stop("Number of bootstrap iterations must be > 0")
  if (is.data.frame(x)) {
    if (ncol(x) != 2L)
      stop("'x' must have 2 columns")
    DNAME <- deparse(substitute(x))
    x <- x[stats::complete.cases(x),]
    colnames(x) <- c("id", "response")
    x$id <- factor(as.factor(x$id), levels=unique(as.factor(x$id)))
    if (is.numeric(x$response)==FALSE)
      stop("'response' must be numeric")
    dat <- x[order(x$id),]
  }
  else {
    if ((l <- length(x)) != length(id))
      stop("'x' and 'id' must have the same length")
    if (is.numeric(x)==FALSE)
      stop("'response' must be numeric")
    DNAME <- deparse(substitute(x))
    OK <- stats::complete.cases(x, id)
    x <- x[OK]
    id <- id[OK]
    id <- factor(as.factor(id), levels=unique(as.factor(id)))
    dat <- data.frame(cbind(id,x))
    colnames(dat)<- c("id", "response")
    dat <- dat[order(dat$id),]
  }
  TEST.METH <- match.arg(test.method)

  ##
  ## FUNCTIONS FROM NEVALAINEN TEST OF ICS CODE
  ## ('testICS' modified so it only computes TF or TCF)
  Fhat<-function(Y,k,weights=FALSE)
  {
    Y0<-matrix(Y[Y[,2]==k,],ncol=4)
    w0<-rep(1,nrow(Y0))
    if (weights==TRUE) w0<-Y0[,4]/sum(Y0[,4])
    points<-matrix(c(sort(Y[,1]),rep(0,nrow(Y))),ncol=2)
    for (i in 1:nrow(points))
    {
      points[i,2] <- stats::weighted.mean((Y0[,1]<=points[i,1]),w0)
    }
    points
  }

  F1mF2<-function(Y1)
  {
    Y1[,2]<-999
    # a cluster size weighted cdf
    F1<-Fhat(Y1,999,weights=TRUE)
    # the usual cdf
    F2<-Fhat(Y1,999)
    d<-abs(F1[,2]-F2[,2])
    max(d)
  }


  cvmF1mF2<-function(Y)
  {
    Y1<-Y
    # mean cdf
    Y1[,2]<-999
    Fbar<-Fhat(Y1,999,weights=TRUE)
    K<-max(Y[,2])
    Tk<-matrix(nrow=1,ncol=K)
    for (i in 1:K)
    {
      Mk<-length(unique(Y[Y[,2]==i,3]))
      if (Mk==0) Tk[,i]<-0
      else
      {
        Fk<-Fhat(Y,i)
        d<-c(Fk[,1],max(Fk[,1]))-c(min(Fk[,1]),Fk[,1])
        d<-d[2:length(d)]
        Tk[,i]<-sum(i*Mk*(d*(Fk[,2]-Fbar[,2])**2))
      }
    }
    sum(Tk)
  }

  testICS2<-function(Y,Bn=1000, meth=c('TF', 'TCM'))
  {
    METH <- match.arg(meth)
    ni<-data.frame(table(Y$id))[,2]
    Ymat<-matrix(ncol=4,nrow=nrow(Y))
    Ymat[,1]<-Y$response
    Ymat[,2]<-rep(ni,ni)
    Ymat[,3]<-Y$id
    Ymat[,4]<-1/rep(ni,ni)
    colnames(Ymat)<-c("response","ni","id","w")
    M<-length(ni)

    # which function are we using
    fun <- switch(METH, TF = F1mF2, TCM = cvmF1mF2)

    #computation of the test statistic
    stat <- fun(Ymat)

    # resample the Vi's in a balanced way
    bs<-NULL
    for (b in 1:Bn)
    {
      if(print.it) print(b)
      # First permute the Y's within each cluster, whole sample
      Yp<-matrix(nrow=0,ncol=4)
      for (l in 1:M)
      {
        z1<-matrix(Ymat[Ymat[,3]==l,],nrow=ni[l],ncol=4)
        z2<-sample(z1[,1],z1[1,2])
        Yp<-rbind(Yp,cbind(matrix(z2,nrow=ni[l],ncol=1),matrix(z1[,2:4],nrow=ni[l],ncol=3)))
      }
      # Bootstrap sample
      Yb<-NULL
      for (q in 1:M)
      {
        # Cluster size
        nib<-ni[q]
        Vis<-matrix(Yp,ncol=4)
        # Sample one cluster index
        i<-sample(unique(c(Vis[,3])),1)
        # First sampled cluster
        Vib1<-matrix(Vis[Vis[,3]==i,],ncol=4)
        # If nib smaller than in the sampled cluster, take nib first elements
        if (nib<=Vib1[1,2])
        {
          Vib2<-matrix(Vib1[1:nib,],ncol=4)
        }
        # If desired nib is larger than the size of the sampled cluster,
        # attach the remaining elements from the most similar cluster
        if (nib>Vib1[1,2])
        {
          Vib2<-Vib1
          # Sufficiently large clusters and their indices
          Vipass<-matrix(Ymat[Ymat[,2]>=nib],ncol=4)
          indpass<-unique(Vipass[,3])
          # Compute distances
          D<-matrix(ncol=2,nrow=length(indpass))
          dind<-0
          for (c in indpass)
          {
            dind<-dind+1
            A<-Vipass[Vipass[,3]==c,]
            D[dind,2]<-sum((Vib2[,1]-A[1:nrow(Vib2),1])**2)
            D[dind,1]<-c
          }
          # Replacement cluster id minimizing the distance, in case of ties, sample one
          r<-D[D[,2]==min(D[,2]),1]
          if (length(r)>1) r<-sample(r,1)
          R<-Vipass[Vipass[,3]==r,]
          Vib2<-rbind(Vib2,R[(nrow(Vib2)+1):nib,])
        }
        Vib2[,2]<-nib
        Vib2[,4]<-1/nib
        Vib2[,3]<-q
        Yb<-rbind(Yb,Vib2)
      }
      bs[b] <- fun(Yb)
      #bs[b,1]<-F1mF2(Yb)
      #bs[b,2]<-cvmF1mF2(Yb)
    }
    #pvalues
    p <- mean(stat<=bs)
    results<-list(statistic = stat, p.value=p, method=METH)
    return(results)
  }

  ret <- testICS2(Y=dat, Bn=bi, meth=TEST.METH)
  pval <- ret$p.value
  STAT <- ret$statistic
  names(STAT) <- ret$method

  METHOD <- paste("Test of informative cluster size", paste('(',ret$method,')',sep=""))

  RVAL <- list(statistic = STAT, p.value = pval,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}
