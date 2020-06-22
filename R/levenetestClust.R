##########################################################
## REWEIGHTED F TEST ANALOG FOR CLUSTERED DATA
##
## From 'Var test CC.R'
##########################################################

#' Reweighted Levene's Test for Homogeneity of Variance in Clustered Data
#'
#' Performs a reweighted test for homogeneity of marginal variances across intra-cluster groups in clustered data.
#' Reweighted to correct for potential cluster- or group size informativness.
#'
#' @param y vector of numeric responses.
#' @param group vector or factor object defining groups.
#' @param id vector or factor object denoting cluster membership for \code{y} reponses.
#' @param center The name of a function to compute the center of each group. If \code{mean}, the reweighted
#' group means will be used. The default \code{median} is the suggested measure of center, as it provides a more
#' robust test.
#' @param trim optional numeric argument taking values \[0, 0.5\] to specify the percentage trimmed mean.
#' Ignored if \code{center = median}.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} is a numeric variable giving the data values and
#' \code{rhs} a factor with two or more levels giving the corresponding groups.
#' @param data an optional matrix or data frame containing variables in the formula \code{formula} and \code{id}.
#' By default the variables are taken from \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param na.action a function which indicates what should happen when data contain \code{NA}s. Defaults to
#' \code{getOption("na.action")}.
#' @param ... further arguments to be passed to or from methods.
#' @return A list with class "\code{htest}" containing the following compoments:
#' \item{statistic}{the value of the test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{parameter}{the degrees of freedom of the chi square distribution.}
#' \item{method}{a character string indicating the test performed.}
#' \item{data.name}{a character string giving the name of the data and the total number of clusters.}
#' \item{M}{the number of clusters.}
#' @references
#' Gregg, M., Marginal methods and software for clustered data with cluster- and group-size informativeness.
#' PhD dissertation, University of Louisville, 2020.
#'
#' @details The null hypothesis is that all levels of \code{group} have equal marginal variances.
#'
#'
#' @examples
#' data(screen8)
#'
#' ## Do boys and girls have the same variability in math scores?
#' ## Test using vectors
#' levenetestClust(y=screen8$math, group=screen8$gender, id=screen8$sch.id)
#'
#' ## Test using formula method
#' levenetestClust(math~gender, id=sch.id, data=screen8)
#'
#' ## Using 10% trimmed mean
#' levenetestClust(math~gender, id=sch.id, data=screen8, center="mean", trim=.1)
#'
#' @export levenetestClust

levenetestClust <- function(y, ...) {
  UseMethod("levenetestClust")
}

#' @rdname levenetestClust
#' @export
levenetestClust.default <- function(y, group, id, center=c("median", "mean"), trim=NA, ...){
  if (!is.numeric(y))
    stop(deparse(substitute(y)), " is not a numeric variable")
  if (is.numeric(trim) && (trim < 0 || trim > 0.5))
    stop("'trim' must be between [0, 0.5]")
  DNAME <- paste(deparse(substitute(y)), "by", deparse(substitute(group)))
  OK <- stats::complete.cases(y, group, id)
  y <- y[OK]
  group <- as.numeric(group[OK])
  id <- id[OK]
  center.meth <- match.arg(center)
  if (is.numeric(trim) && center=="median") {
    trim <- 0
  }

  ##
  ## FUNCTIONS
  ##
  ## Function to compute reweighted group means/trimmed means
  grpwt.center.fun <- function(t.in, id.in, grp.in, trim.in) {
    ## calculate which groups are present in each cluster
    tmp <- tapply(t.in, list(id.in, grp.in), mean)
    grp.typ <- 1-1*is.na(tmp)
    cl.typ.mat <- apply(grp.typ, 2, function(x) x*rowSums(grp.typ))
    grpsize.mat <- table(id.in, grp.in)
    wt.mat <- cl.typ.mat*grpsize.mat
    M.tilde <- apply(cl.typ.mat, 2, function(x) sum(table(x[x!=0])/as.numeric(names(table(x[x!=0])))))
    dat.in <- data.frame(cbind(id.in, grp.in, t.in))

    if (trim.in > 0) {
      ntrim <- trim.in
      ##
      ## Get reweighted rank of each observation, by group
      ecdf.grp <- function(x, q, cluster, wt) {
        temp <- outer(x, q, "<=")*1
        sum.x <- rowsum(temp, cluster, reorder=F)
        ecdf.x <- sum.x/wt
        colSums(ecdf.x)
      }
      rank.fun <- function(grp){
        ecdf.grp(t.in[grp.in==grp], t.in[grp.in==grp], id.in[grp.in==grp], wt.mat[,grp][wt.mat[,grp]!=0])/M.tilde[grp]
      }
      Rlist <- lapply(1:ncol(grp.typ), rank.fun)  # list of ranks for each group
      ##
      ## Function for trimmed reweighted mean
      grpwt.mean <- function(grp) {
        dat.g <- dat.in[dat.in$grp.in==grp,]
        t.tmp <- t.in[grp.in==grp]
        r.tmp <- unlist(Rlist[grp])
        ## data frame with rank
        dat.gr <- data.frame(cbind(dat.g, r.tmp))
        ## sort by rank
        dat.ord <- dat.gr[order(r.tmp),]
        n.tmp <- nrow(dat.ord)
        lo <- floor(ntrim*n.tmp)
        trimdat <- dat.ord[(lo+1):(n.tmp-lo),]

        trimdat.all <- rbind(trimdat[,1:3], dat.in[dat.in$grp.in!=grp,])

        ## redo the weights after removing the observations (only the number of observations in the cluster changes; the
        ## cl.typ.mat values stay the same (b/c/ we're only dealing with one group at a time))
        t.trim <- trimdat.all$t.in
        id.trim <- trimdat.all$id.in
        group.trim <- trimdat.all$grp.in

        tmp.trim <- tapply(t.trim, list(id.trim, group.trim), mean)
        grp.typ.trim <- 1-1*is.na(tmp.trim)
        cl.typ.mat.trim <- apply(grp.typ.trim, 2, function(x) x*rowSums(grp.typ.trim))
        wt.trim <- (cl.typ.mat.trim[,grp][cl.typ.mat.trim[,grp]!=0])*table(trimdat$id.in)

        grpmean.tmp <- tapply(trimdat$t.in, list(trimdat$id.in), sum)/wt.trim

        M.tilde.trim <- apply(cl.typ.mat.trim, 2, function(x) sum(table(x[x!=0])/as.numeric(names(table(x[x!=0])))))[grp]
        trimmean <- sum(grpmean.tmp)/M.tilde.trim
        return(trimmean)
      }
      centerval <- unlist(lapply(1:ncol(grp.typ), grpwt.mean))
    }
    else { ## Reweighted mean, no trim
      ## function to get the interior weighted sums
      ## Takes tmp and c.typ.mat, and loops through each group to get sum
      sum.fun <- function(dat2, typ) {
        grp.sum <- NULL
        for (i in 1:ncol(dat2)) {
          ok <- stats::complete.cases(dat2[,i])
          grp.sum[i] <- sum(dat2[ok,i]/typ[ok,i])
        }
        return(grp.sum)
      }
      ## Apply the function to get the interior sums for each group
      sum.grp <- sum.fun(tmp, cl.typ.mat)
      ## Divide interior sums by respective group denominator value
      centerval <- sum.grp/M.tilde
    }
    return(centerval)
  }

  ##
  ## Function to compute reweighted ANOVA test (will be applied to transformed values)
  anova.fun <- function(cwmat.in, cmat.in) {
    wtgrp.fun <- function(dat){
      ## Get matrix of indicators for group
      grp.mat <- 1-1*is.na(dat)
      ## multiply respective group indicator column (e.g, column 1 for group 1) by the rowSums
      cl.typ.mat <- apply(grp.mat, 2, function(x) x*rowSums(grp.mat))

      ## Caclulate denominator values for each group
      D.grp <- apply(cl.typ.mat, 2, function(x) sum(table(x[x!=0])/as.numeric(names(table(x[x!=0])))))

      ## function to get the interior weighted sums
      ## Takes cw.mat and c.typ.mat, and loops through each group to get sum
      sum.fun <- function(dat2, typ) {
        grp.sum <- NULL
        for (i in 1:ncol(dat2)) {
          ok <- stats::complete.cases(dat2[,i])
          grp.sum[i] <- sum(dat2[ok,i]/typ[ok,i])
        }
        return(grp.sum)
      }

      ## Apply the function to get the interior sums for each group
      sum.grp <- sum.fun(dat, cl.typ.mat)
      ## Divide interior sums by respective group denominator value
      wt.means <- sum.grp/D.grp
      return(wt.means)
    }

    wtgrp.fun.single <- function(dat, grp){
      wtgrp.fun(dat)[grp]
    }

    jk.fun <- function(x, xdat, grp) {
      tdat <- xdat[x,]
      wtgrp.fun.single(tdat, grp)
    }
    Min <- nrow(cwmat.in)
    Kin <- ncol(cwmat.in)

    theta.hat <- wtgrp.fun(cwmat.in)
    jk.apply <- function(grp) {
      bootstrap::jackknife(1:Min, jk.fun,cwmat.in, grp)$jack.values
    }
    theta.hat.i <- matrix(unlist(lapply(1:Kin, jk.apply)), ncol=Kin, byrow=F)
    theta.bar <- colMeans(theta.hat.i)
    t.dmat <- t(apply(theta.hat.i, 1, function(x) x - theta.bar))
    t.dmat2 <- apply(t.dmat, 1, function(x) x%*%t(x))
    ## corrected variance, statisic, pval
    JK.var <- (Min/(Min - Kin))*Min*(Min-1)*matrix(apply(t.dmat2, 1, mean), nrow=Kin)
    stat <- Min*t(cmat.in%*%t(t(theta.hat)))%*%MASS::ginv(cmat.in%*%JK.var%*%t(cmat.in))%*%(cmat.in%*%t(t(theta.hat)))
    pval <- stats::pchisq(stat, df=(Kin-1), lower.tail=F)
    return(list(statistic=stat, p.value=pval))
  }

  ##
  ## TRANSFORM THE DATA AND PERFORM THE TEST
  trim.tmp <- ifelse(is.numeric(trim) && trim > 0, trim, 0)
  if(center.meth=="mean") {
    meds <- grpwt.center.fun(t.in = y, id.in = id, grp.in = group, trim.in=trim.tmp)
  } else{
    meds <- tapply(y, list(group), stats::median)
  }
  y.trans <- abs(y - meds[group])
  ytrans.tab <- tapply(y.trans, list(id, group), mean)
  ## create contrast matrix
  k <- ncol(ytrans.tab)
  cfun <- function(pos) {
    1*(seq(1:k)==pos) - 1*(seq(1:k)==(pos+1))
  }
  cmat <- matrix(unlist(lapply(1:(k-1), cfun)), ncol=k, byrow=T)
  CW <- anova.fun(ytrans.tab, cmat)
  dots <- deparse(substitute(...))
  method.tmp <- ifelse(trim.tmp > 0, paste0(center.meth,': ', trim.tmp), center.meth)
  METHOD <- paste0("Reweighted Levene's Test for Homogeneity of Variance in Clustered Data (center = ",
                   method.tmp, ')')
  STATISTIC <- CW$statistic
  PARAMETER <- k-1
  M <- nrow(ytrans.tab)
  PVAL <- CW$p.value
  names(STATISTIC) <- "X-squared"
  names(PARAMETER) <- "df"
  DNAME <- paste0(paste0(DNAME, ", M = "), as.character(M))
  RVAL <- list(statistic = STATISTIC, p.value = PVAL,
               parameter = PARAMETER, method = METHOD, data.name = DNAME, M = M)
  class(RVAL) <- "htest"
  if (M < 30)
    warning('Number of clusters < 30. Normal approximation may be incorrect')
  return(RVAL)
}


#' @rdname levenetestClust
#' @export
levenetestClust.formula <- function (formula, id, data, subset, na.action, ...)
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
  z <- do.call("levenetestClust", args=c(list(mf$r, mf$group, mf$id, ...)))
  z$data.name <- paste0(paste0(DNAME, ", M = "), as.character(length(unique(mf$id))))
  z
}
