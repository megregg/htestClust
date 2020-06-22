####################################
## Clustered Wilcoxon function
####################################

#' Rank Sum and Signed Rank Tests for Clustered Data
#'
#' Performs a one-sample or paired cluster-weighted signed rank test, or a cluster- or
#' group-weighted rank sum test. These tests are appropriate for clustered data with
#' potentially informative cluster size.
#'
#' @param x,y numeric vectors of data values.
#' @param idx vector or factor object denoting cluster membership for \code{x} observations (or cluster
#' membership for paired observations when \code{paired} is \code{TRUE}). Length must be equal
#' to length of \code{x}.
#' @param idy vector or factor object denoting cluster membership for \code{y} observations. Length must be equal
#' to length of \code{y}
#' @param alternative indicates the alternative hypothesis and must be one of "\code{two.sided}", "\code{greater}",
#' or "\code{less}".You can specify just the initial letter.
#' @param mu a number specifying an optional parameter used to form the null hypothesis. Ignored when
#' performing a rank-sum test. See 'Details'.
#' @param paired a logical indicating whether \code{x} and \code{y} are paired. When \code{TRUE}, the
#' cluster-weighted signed rank test is performed.
#' @param method a character string specifying the method of rank sum test to be performed. See 'Details'.
#' @param formula a formula of the form \code{lhs} ~ \code{rhs}, where \code{lhs} is a numberic variable
#' giving the data values and \code{rhs} a numeric or factor with two levels giving the correponding groups.
#' @param id a vector or factor object denoting cluster membership.
#' @param data an optional matrix or data frame containing variables in the formula \code{formula} and \code{id}.
#' By default the variables are taken from \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param na.action a function which indicates what should happen when data contain \code{NA}s. Defaults to
#' \code{getOption("na.action")}.
#' @param ... further arguments to be passed to or from methods.
#' @return A list with class "\code{htest}" containing the following compoments:
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{null.value}{the location parameter \code{mu}. Always 0 for rank sum test.}
#' \item{data.name}{a character string giving the name(s) of the data and the total number of clusters.}
#' \item{method}{a character string indicating the test performed and method of construction.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{M}{the number of clusters.}
#' @references
#' Datta, S., Satten, G. (2005) Rank-sum tests for clustered data.
#' \emph{J. Am. Stat. Assoc.}, \bold{100}, 908--915.
#'
#' Datta, S., Satten, G. (2008) A signed-rank test for clustered data.
#' \emph{Biometrics}, \bold{64}, 501--507.
#'
#' Dutta, S., Datta, S. (2015) A rank-sum test for clustered data when the number of subjects in a
#' group within a cluster is informative. \emph{Biometrics}, \bold{72}, 432--440.
#'
#' @details The formula interface is only applicable for the 2-sample rank-sum tests.
#'
#' If only \code{x} and \code{idx} are given, a cluster-weighted signed rank test of the null that
#' the distribution of \code{x} is symmetric about \code{mu} is performed.
#'
#' If \code{x} and \code{y} are given and \code{paired} is \code{TRUE}, only \code{idx} is necessary (\code{idy}
#' is ignored). In this case, a cluster-weighted signed-rank test of the null that the distribution of \code{x - y}
#' is symmetric about \code{mu} is performed.
#'
#' When \code{method} is \code{cluster}, the cluster-weighted rank sum test of Datta and Satten (2005) is performed.
#' The data must have complete intra-cluster group distribution (i.e., all clusters must contain observations
#' belonging to both groups) for this test to be performed.
#'
#' When \code{method} is \code{group}, the group-weighted rank-sum test of Dutta and Datta (2015) is performed.
#' This test is appropriate for clustered data with potentially informative intra-cluster group size. Incomplete
#' intra-cluster group distribution is permitted.
#'
#' For the rank sum tests, the null is that the two groups follow the same marginal distribution. \code{mu} is
#' ignored when performing these tests.
#'
#' The tests performed by this function involve computation of reweighted empirical CDFs. This is computationally
#' intensive and can result in lengthy execution time for large data sets.
#'
#' @examples
#' \dontrun{
#'   data(screen8)
#'   ## One-sample signed rank test
#'   wilcoxtestClust(x=screen8$math, idx=screen8$sch.id, mu=70)
#'
#'   ## Paired signed rank test
#'   wilcoxtestClust(x=screen8$math, y=screen8$read, idx=screen8$sch.id, paired=TRUE, mu=10)
#'
#'   ## Cluster-weighted rank sum test
#'   wilcoxtestClust(math~gender, id=sch.id, data=screen8)
#'
#'
#'   ## Group-weighted rank sum test
#'   boys <- subset(screen8, gender=='M')
#'   girls <- subset(screen8, gender=='F')
#'   wilcoxtestClust(x=boys$math, y=girls$math, idx=boys$sch.id, idy=girls$sch.id, method="group")
#'
#'   ## Group-weighted rank sum test using formula method
#'   wilcoxtestClust(math~gender, id=sch.id, data=screen8, method="group")
#'}
#' @export

wilcoxtestClust <- function(x, ...) {
  UseMethod("wilcoxtestClust")
}

#' @rdname wilcoxtestClust
#' @export
wilcoxtestClust.default <- function (x, y = NULL, idx, idy=NULL,
                                       alternative = c("two.sided", "less", "greater"),
                                       mu = 0, paired = FALSE, method = c("cluster", "group"), ...)
{
  meth <- match.arg(method)
  alternative <- match.arg(alternative)
  if (!missing(mu) && ((length(mu) > 1L) || !is.finite(mu)))
    stop("'mu' must be a single number")
  if (!is.numeric(x))
    stop("'x' must be numeric")
  if (!is.null(y)) {
    if (!is.numeric(y))
      stop("'y' must be numeric")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

    if (paired) { ## SIGNED RANK: PAIRED
      if (length(x) != length(y))
        stop("'x' and 'y' must have the same length")
      OK <- stats::complete.cases(x, y, idx)
      x <- x[OK] - y[OK]
      idx <- idx[OK]
      y <- idy <- NULL
      m <- length(unique(idx))
      METHOD <- "Paired cluster-weighted signed rank test"
    }
    else { ## RANK SUM: TWO SAMPLE CASE
      ## Set up data for rank sum tests
      ## Remove missing values, set up in data frame with ID, X, Grp (required by Dutta/Datta code)
      xok <- stats::complete.cases(x, idx)
      yok <- stats::complete.cases(y, idy)
      x <- x[xok]
      idx <- factor(as.factor(idx[xok]), levels=unique(as.factor(idx[xok])))
      y <- y[yok]
      idy <- factor(as.factor(idy[yok]), levels=unique(as.factor(idy[yok])))
      xdat <- data.frame(idx, x, rep(0,length(x)))  ## numeric group variable of 0/1 required by Dutta/Datta code
      ydat <- data.frame(idy, y, rep(1,length(y)))
      colnames(xdat) <- c("ID", "X", "grp")
      colnames(ydat) <- c("ID", "X", "grp")
      dat <- rbind(xdat, ydat)
      dat[,1] <- as.numeric(dat[,1])  ## Dutta/Datta code uses numeric IDs
      dat <- dat[with(dat, order(ID)),]
    }
  }
  else { ## SIGNED RANK: ONE SAMPLE CASE
    DNAME <- deparse(substitute(x))
    if (paired)
      stop("'y' is missing for paired test")
    xok <- stats::complete.cases(x, idx)
    yok <- NULL
    x <- x[xok]
    idx <- idx[xok]
    m <- length(unique(idx))
    METHOD <- "One sample cluster-weighted signed rank test"
  }

  if (length(x) < 1L)
    stop("not enough (finite) 'x' observations")
  ## PERFORM SIGNED RANK TEST (ONE SAMPLE OR PAIRED)
  if (is.null(y)) {
    ## The 'Zeros' code came from R code from iid signed rank. But Datta/Satten don't assume
    ## F is continuous, and have ties correction
    x <- x - mu
    #ZEROES <- any(x == 0)
    #if (ZEROES)  {
    #  x <- x[x != 0]
    #  idx <- idx[x!=0]
    #}
    ## Ensuring data is order by cluster ID and cluster ID is a factor
    datx <- cbind.data.frame(idx, x)
    datx <- datx[with(datx, order(idx)),]
    x <- datx$x
    idx <- factor(as.factor(datx$idx), levels=unique(as.factor(datx$idx)))

    ## INSERTING DATTA CODE HERE
    ##
    ## Functions required for Signed Rank code (Datta)
    ##
    Fi <- function(x,i) { Xi <- Xij[(cni[i]+1):(cni[i+1])];
    (sum(abs(Xi)<=x)+sum(abs(Xi)<x))/(2*ni[i])}

    Ftot <- function(x) { st <- 0;
    for (i in 1:g) st <- st + Fi(x,i);
    return(st)}

    Fcom <- function(x) { st <- 0;
    for (i in 1:g) st <- st + Fi(x,i)*ni[i];
    return(st/n)}

    Xij <- x
    ni <- as.vector(table(idx))
    g <- length(ni)
    n <- sum(ni)
    cni <- c(0, cumsum(ni))

    TS <- VTS <- 0
    for (i in 1:g) {

      Xi <- Xij[(cni[i]+1):(cni[i+1])]
      first <- (sum(Xi>0)-sum(Xi<0))/length(Xi)
      second <- 0
      third <- 0
      for (x in Xi) { second <- second + sign(x)*(Ftot(abs(x))-Fi(abs(x),i));
      third <- third + sign(x)*Fcom(abs(x))}

      TS <- TS + first+second/length(Xi)
      VTS <- VTS + (first+ (g-1)*third/length(Xi))^2
    }
    Z <- TS/sqrt(VTS)
    PVAL <- switch(alternative, less = stats::pnorm(Z), greater = stats::pnorm(Z,
                                                                 lower.tail = FALSE), two.sided = 2*(1-stats::pnorm(abs(Z))))
  }
  ## PERFORM RANK SUM TEST (CLUSTER OR GROUP WEIGHTED)
  else {
    ## only tests if distribution of x and y are equal (can't test if the shift is different than zero)
    mu <- 0
    m <- length(unique(dat[,1]))
    ## CLUSTER WEIGHTED (DATTA/SATTEN CODE)
    if (meth=="cluster") {
      METHOD <- "Cluster-weighted rank sum test"

      clus.rank.sum<-function(Cluster,X,grp) {
        ## group-weighted test doesn't have incomplete cluster modification
        if (length(unique(Cluster[grp==1]))!=length(unique(Cluster[grp==0])))
          stop ("Incomplete intra-cluster group structure: can not apply cluster-weighted rank sum test")

        #####calculate quantity 2 (using the pooled estimate of F)
        n<-length(X)
        F.hat<-numeric(n)
        for (i in 1:n){
          F.hat[i]<-(sum(X<=X[i])+sum(X<X[i]))/(2*n)
        }
        #####calculate quantity 1 (using ECD-F for each cluster)
        #### M is No. of clusters, n is No. of observations
        M<-length(unique(Cluster))
        n.i<-table(Cluster)
        F.prop<-numeric(n)
        for(ii in 1:n){
          F.j<-numeric(M)
          for (i in 1:M){
            F.j[i]<-(sum(X[Cluster==i]<X[ii])+0.5*sum(X[Cluster==i]==X[ii]))/(n.i[i])
          }
          F.prop[ii]<-sum(F.j[-Cluster[ii]])
        }

        ###########calculate S=E(W*|X,g)
        a<-numeric(M)
        b<-1+F.prop
        for (i in 1:M){
          a[i]<-sum((grp[Cluster==i]*b[Cluster==i])/(n.i[i]))
        }
        c<-1/(M+1)
        S<-c*sum(a)
        ########note: for m groups maybe can use grp[Cluster==i&grp=m]

        #########Calculate E(S)=E(W*)
        #n.i1<-table(Cluster[grp==1])
        n.i1 <- table(Cluster, grp)[,2]
        d<-n.i1/n.i
        E.S<-(1/2)*sum(d)

        #######Calculate estimate of variance of S
        W.hat<-numeric(M)        #####first calculate W.hat for each cluster
        a<-n.i1/n.i
        for (i in 1:M){
          b<-1/(n.i[i]*(M+1))
          c<-(grp[Cluster==i])*(M-1)
          d<-sum(a[-i])
          W.hat[i]<-b*sum((c-d)*F.hat[Cluster==i])
        }
        a<-n.i1/n.i
        E.W<-(M/(2*(M+1)))*(a-sum(a)/M)    ##second, calculate E(W)

        var.s<-sum((W.hat-E.W)^2) #calculate var(s)
        stat<-(S-E.S)/sqrt(var.s)   #calculate the test statistic
        #p.value<-2*pnorm(abs(stat),lower.tail=F)
        #list(S=S,E.S=E.S,Var.S=var.s,z.stat=stat,p.value=p.value)
        stat
      }
      Z <- clus.rank.sum(dat$ID, dat$X, dat$grp)
      PVAL <- switch(alternative, less = stats::pnorm(Z), greater = stats::pnorm(Z, lower.tail = FALSE),
                     two.sided = 2*(1-stats::pnorm(abs(Z))))

    }
    else {
      ## GROUP WEIGHTED RANK SUM (DUTTA/DATTA)
      METHOD <- "Group-weighted rank sum test"
      ## Required functions
      rn<-function(dv){
        ik=dv[1]
        x=dv[2]
        ds1=dat[dat[,3]==1,]
        vs1=(kh==2)*(ds1[,2]<x)+(kh==1)*(ds1[,2]<=x)
        ic1 <- subset(unique(dat[,1]), !(unique(dat[,1]) %in% unique(ds1[,1])))
        if (length(ic1)==0) {
          sl1=stats::aggregate(vs1,list(ds1[,1]),mean)[,2]
        }
        else {
          cmp1 <- stats::aggregate(vs1,list(ds1[,1]),mean)
          incp1 <- data.frame(cbind(ic1, rep(0, length(ic1))))
          colnames(incp1) <- colnames(cmp1)
          tmp1 <- rbind(cmp1, incp1)
          sl1 <- tmp1[with(tmp1, order(tmp1[,1])), 2]
        }
        ds2=dat[dat[,3]==0,]
        vs2=(kh==2)*(ds2[,2]<x)+(kh==1)*(ds2[,2]<=x)
        ic2 <- subset(unique(dat[,1]), !(unique(dat[,1]) %in% unique(ds2[,1])))
        if (length(ic2)==0) {
          sl2=stats::aggregate(vs2,list(ds2[,1]),mean)[,2]
        }
        else {
          cmp2 <- stats::aggregate(vs2,list(ds2[,1]),mean)
          incp2 <- data.frame(cbind(ic2, rep(0, length(ic2))))
          colnames(incp2) <- colnames(cmp2)
          tmp2 <- rbind(cmp2, incp2)
          sl2 <- tmp2[with(tmp2, order(tmp2[,1])), 2]
        }
        Fwt <- 1*unique(dat[,1])%in%c(ic1,ic2) + 2*(1-1*unique(dat[,1])%in%c(ic1,ic2))

        fg=(sl1+sl2)/Fwt
        fg[ik]=0
        return(fg)
      }

      rst <- function(il){
        ly=sum(mat[-which(dw[,1]==il),-il])
        return(ly)
      }

      m <- g <- length(unique(dat[,1]))
      dw <- dat[(dat[,3]==1),]
      ns <- (dw[,1])
      incc1 <- subset(unique(dat[,1]), !(unique(dat[,1]) %in% unique(dat[(dat[,3]==0),][,1])))
      incc0 <- subset(unique(dat[,1]), !(unique(dat[,1]) %in% unique(dat[(dat[,3]==1),][,1])))
      nv <- 4*(1-1*(ns%in%incc1))*as.vector(table(ns)[match(ns,names(table(ns)))]) + 2*(1*(ns%in%incc1))

      kh <- 1
      mat <- t(apply(cbind(dw[,1:2]),1,rn))/nv
      vf1 <- apply(cbind(seq(1,m)),1,rst)
      sFs1 <- sum(mat)

      kh <- 2
      mat <- t(apply(cbind(dw[,1:2]),1,rn))/nv
      vf2 <- apply(cbind(seq(1,m)),1,rst)
      sFs2 <- sum(mat)

      I <- sum(1*unique(dat[,1])%in%c(incc1, incc0))
      C <- sum(1-1*unique(dat[,1])%in%c(incc1, incc0))
      v1=(sFs1+sFs2)+(C+2*I)/2
      vd= (vf1+vf2)+((C+2*I)-1)/2

      h=1
      TS <- v1
      E.T<- 0.25*m*(m+1)

      test=(m/m^h)*v1-((m-1)/(m-1)^h)*vd
      v.test=stats::var(test)
      v_hat=(((m^h)^2)/(m-1))*v.test
      v.hat=ifelse(v_hat==0,0.00000001,v_hat)
      Z <- (TS-E.T)/sqrt(v.hat)
      PVAL <- switch(alternative, less = stats::pnorm(Z), greater = stats::pnorm(Z, lower.tail = FALSE),
                     two.sided = 2*(1-stats::pnorm(abs(Z))))
    }
  }
  names(Z) <- "z"
  #names(g) <- "m"
  DNAME <- paste0(paste0(DNAME, ", M = "), as.character(m))
  names(mu) <- if (paired || !is.null(y))
    "location shift"
  else "location"

  rval <- list(statistic=Z, p.value=PVAL, null.value = mu, data.name=DNAME, method=METHOD,
               alternative=alternative, M = m)
  class(rval) <- "htest"
  if (m < 30)
    warning('Number of clusters < 30. Normal approximation may be incorrect')
  return(rval)
}

#' @rdname wilcoxtestClust
#' @export
wilcoxtestClust.formula <- function (formula, id, data, subset, na.action, ...)
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
  y <- do.call("wilcoxtestClust", args=c(list(DATA$x$r, DATA$y$r, DATA$x$id, DATA$y$id,...)))
  #y$data.name <- DNAME
  y$data.name <- paste0(paste0(DNAME, ", M = "), as.character(length(unique(mf$id))))
  y
}
