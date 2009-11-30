#' Holes index.
#'
#' Calculates the holes index. See Cook and Swayne (2007)
#' Interactive and Dynamic Graphics for Data Analysis for equations.
#'
#' @param mat matrix being used
#' @keywords hplot
holes <- function(mat) {
  n <- nrow(mat)
  d <- ncol(mat)

  num <- 1 - 1/n * sum(exp(-0.5 * rowSums(mat ^ 2)))
  den <- 1 - exp(-d / 2)
  
  num/den
}



#' Central mass index.
#'
#' Calculates the central mass index.  See Cook and Swayne (2007)
#' Interactive and Dynamic Graphics for Data Analysis for equations.
#'
#' @param mat matrix being used
#' @keywords hplot
cm <- function(mat)
  1 - holes(mat)



#' LDA projection pursuit index.
#'
#' Calculate the LDA projection pursuit index.  See Cook and Swayne (2007)
#' Interactive and Dynamic Graphics for Data Analysis for equations.
#' 
#' @param mat matrix being used
#' @param cl class to be used.  Such as "color"
#' @keywords hplot
#'
lda_pp <- function(cl) {
    if (length(unique(cl)) == 0)
      stop("ERROR: You need to select the class variable!")
    if (length(unique(cl)) == 1)
      stop("LDA index needs at least two classes!")

  function(mat) {
    fit <- manova(mat ~ cl)                      

    1 - summary(fit,test="Wilks")$stats[[3]]    
  }
}

#' PDA projection pursuit index.
#'
#' Calculate the PDA projection pursuit index.  See Lee and Cook (2009)
#' A Projection Pursuit Index for Large p, Small n Data
#' 
#' @param mat matrix being used
#' @param cl class to be used.  Such as "color"
#' @keywords hplot
#'
pda_pp <- function(cl, lambda=0.2)
{
  if (is.null(lambda)) 
    stop("ERROR : You need to use parameter lambda!")
  if (length(unique(cl)) == 0)
    stop("ERROR: You need to select the class variable!")
  if (length(unique(cl)) < 2)
    stop("PDA index needs at least two classes!")

  function(mat) {
    mat <- as.matrix(mat);
    cl <- as.matrix(cl);
    if (ncol(cl)!=1) cl<-t(cl)

    p<-ncol(mat);n<-nrow(mat);	  # n: no. of obs, p: no. of variables
    ngroup<-table(cl)             # different classes
    groups<-length(ngroup)	  # no. of classes
    gn<-names(table(cl))	  # names of classes

    cl.i<-matrix(rep(0,n),n)
    for (i in 1:n) {
      for (j in 1:groups) if (cl[i, 1] == gn[j]) 
        cl.i[i, 1] <- j
    }

    cl.sort<-sort.list(cl.i)
# the class label after sorted, num 1,2,3... rep. different classes
    cl<-cl.i[cl.sort]
    mat<-mat[cl.sort,]   # sort mat according to class label 		  
    ngroup<-table(cl.i)  # the obs. no. in each class, stored in table	  
    groups<-length(ngroup)  # no. of classes
    gname<-names(ngroup)	  # names of classes, now is integer 1,2,3...

    CalIndex(n, p, groups, mat, cl, gname, as.integer(ngroup), lambda)
# n: no. of obs,
# p: no. of variables,
# groups: no. of classes,
# mat: sorted mat,
# cl: sorted class label of all mat 
# gname: names of classes, now is integer 1,2,3...;
#    as.integer(ngroup): the obs. no. in each class
  }
  #return(pda_index)
}   

# fval <==> data; groupraw <==> class: sorted class label of all data;
# group: class label of all data; 
# ngroup: the obs. no. in each class
CalIndex<-function(n, p, groups, fvals, groupraw, gname, ngroup, lambda)
{ 
  # int i, j, k, right, left
  g = groups
  mean = matrix(rep(0,g*p),g,p)
  ovmean = matrix(rep(0,p),p)
  group = matrix(rep(0,n),n)  # the class label of each obs

  right = n-1
  left = 0
 
  group = groupraw;
  
  val = 0

# Calculate mean for within class and the overall mean
  for (i in 1:n) {
    for (j in 1:p) {
      mean[group[i],j] = mean[group[i],j] + fvals[i,j]/ngroup[group[i]]
      ovmean[j] = ovmean[j] + fvals[i,j]/n
    }
  }

  cov = matrix(rep(0,p*p),p,p)
  tempcov = matrix(rep(0,p*p),p,p)

  for (i in 1:n) {
    for (j in 1:p) {
      for (k in 1:j) {
        cov[k,j] = cov[k,j] + (1-lambda)*(fvals[i,j] -
             mean[group[i],j]) * (fvals[i,k]-mean[group[i],k])
        cov[j,k] = cov[k,j]
      }
    }
  }
  
  for (j in 1:p)
    cov[j,j] = cov[j,j] + lambda*n

  tempcov = cov

# =========================================================
  pivot = matrix(rep(0,p),p)
  det = 0

  val = ludcomp(tempcov, p, pivot)
# ===========================================================

  for (j in 1:p)
    for (k in 1:p)
      for (i in 1:g)
	cov[j,k] = cov[j,k] + (1-lambda) * ngroup[i] *
          (mean[i,j] - ovmean[j]) * (mean[i,k] - ovmean[k])

  tempcov = cov
  tempval = ludcomp(tempcov, p, pivot)
  
  if (tempval < 0.00000001) {
    val = 0
    print ("ZERO VARIANCE!")
  }
  else
    val = 1 - val/tempval
  
  return(val)
}

#==============================================================================
ludcomp <- function(a, n, pivot)
{
  det = 1;
  temp = 0;
  s = matrix(rep(0,n),n)
  for (i in 1:n) {
    s[i]=a[i,1+1]
    for (j in 2:n)
      if (s[i]<a[i,j]) 
        s[i] = a[i,j]
  }
  
  for( k in 1:(n-1)) {
    for ( i in k:n )
    {
      temp = abs(a[i,k]/s[i])
      if (i==k) {
        c = temp
        pivot[k] = i
      }
      else if (c<temp)  {
	c = temp
	pivot[k] = i
      }
    }

    if (pivot[k] != k)
    {
      det = det * (-1)
      for (j in k:n) {
        temp = a[k,j]
        a[k,j] = a[pivot[k],j]
        a[pivot[k],j] = temp
      }
      temp = s[k]
      s[k] = s[pivot[k]]
      s[pivot[k]] = temp
    }
	
    for ( i in (k+1):n) {
      temp = a[i,k]/a[k,k]
      a[i,k] = temp
      for (j in (k+1):n)
        a[i,j] = a[i,j] - temp*a[k,j]
      }
      det = det * a[k,k]
  }
  
  det = det * a[n,n]
  
  return(det)

}
