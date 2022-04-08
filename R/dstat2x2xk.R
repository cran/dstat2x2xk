
dstat2x2xk<-function(tab,gamma=1,kappa=NULL,lambda=NULL,rnd=2,warn=TRUE){

  # Check input

  if (!(is.vector(rnd)&(length(rnd)==1)&(rnd>=0))) stop("rnd should be >=0")
  if (!is.array(tab)) stop("tab must be an array")
  if (!(length(dim(tab))==3)) stop("tab must be a 2 x 2 x k array")
  if (!(dim(tab)[1]==2)) stop("tab must be a 2 x 2 x k array")
  if (!(dim(tab)[2]==2)) stop("tab must be a 2 x 2 x k array")
  if (!(dim(tab)[3]>=2)) stop("tab must be a 2 x 2 x k array with k>=2")
  if (!(min(as.vector(tab))>=0)) stop("tab must be a table of nonnegative counts")
  if (!all(as.vector(tab)==round(as.vector(tab)))) stop("tab must be a table of nonnegative integer counts")
  if ((!is.null(kappa))&(!is.null(lambda))){
    stop("Only 1 of kappa and lambda may be nonnull.")
  }
  if (!is.logical(warn)) stop("warn must be TRUE or FALSE")
  if ((is.null(kappa))&(is.null(lambda))){
    kappa<-0
    warning("kappa defaulted to 0 for the usual Mantel-Haenszel test.")
  }
  if (!is.null(lambda)){
    if ((!is.vector(lambda))&(length(lambda)==1)) stop("lambda must be one number")
    if ((lambda<=0)|(lambda>=1)) stop("lambda, if specified, must be between 0 and 1.")
  }
  if (!is.null(kappa)){
    if ((!is.vector(kappa))&(length(kappa)==1)) stop("kappa must be one number")
    if ((kappa<0)) stop("kappa must be nonnegative.")
    if ((kappa>1)) warning("Values of kappa above 1 are not recommended.  Crash possible.")
  }
  initialk<-dim(tab)[3]

  # A 2x2 subtable of the 2x2xk table is degenerate if its row or column
  # totals include a 0 count.  The function check2x2xktable eliminates
  # degenerate subtables, if there are any, and warns that it has done
  # this.
  check2x2xktable <- function(tb) {
    r1 <- tb[1, 1, ] + tb[1, 2, ]
    r2 <- tb[2, 1, ] + tb[2, 2, ]
    c1 <- tb[1, 1, ] + tb[2, 1, ]
    c2 <- tb[1, 2, ] + tb[2, 2, ]
    mc <- pmin(pmin(r1, r2), pmin(c1, c2))==0
    if (sum(mc)==initialk) stop(paste("All",initialk,"tables are degenerate."))
    if (sum(mc)>=1){
      if (warn) warning(paste(sum(mc),"of the",initialk,"subtables are degenerate and were removed."))
      tb<-tb[,,!mc,drop=FALSE]
    }
    tb
  }

  tab<-check2x2xktable(tab) #Check for and remove degenerate subtables
  finalk<-dim(tab)[3]

  # The function one2x2 table computes the null truncated
  # extended hypergeometric distribution for one 2x2 subtable of
  # the full 2x2xk table
  one2x2 <- function(tb, Gamma1 = gamma) {
    # Computations for one of the 2x2 tables
    m1 <- tb[1, 1] + tb[1, 2]
    m2 <- tb[2, 1] + tb[2, 2]
    np <- tb[1, 1] + tb[2, 1]
    mx <- min(np, m1)
    mn <- max(0, np - m2)
    x <- mn:mx
    g <- rep(0, mx + 1)
    pr <- BiasedUrn::dFNCHypergeo(x, m1, m2, np, Gamma1)
    g[(mn + 1):(mx + 1)] <- pr #marginal distribution
    vals<-0:mx
    expect<-sum(g*vals)
    if (!is.null(kappa)) ct<-round(expect*kappa)
    else {
      cusumg<-cumsum(g)
      if (min(cusumg)>lambda) ct<-mn
      else ct<-max(which(cumsum(g)<=lambda))-1
    }
    gc<-(vals>=ct)*g
    if (sum(gc)>0) {
      gc<-gc/sum(gc) #conditional distribution
      cexpect<-sum(vals*gc)
      cvar<-sum(((vals-cexpect)^2)*gc)
    }
    else {
      gc<-NULL
      cexpect<-NULL
      cvar<-NULL
    }
    names(g)<-vals
    list(expect=expect,ct=ct,g=g,gc=gc,cexpect=cexpect,cvar=cvar)
  }

  # Convolution of two probability distributions
  gconv <- function (g1, g2){
    convolve(g1, rev(g2), type = "open")
  }

  o<-matrix(NA,finalk,7)
  colnames(o)<-c("n11","E(n11)","cut","use","E(n11|use=1)","var(n11|use=1)","OR")
  g<-1

  for (i in 1:finalk){
    h<-one2x2(tab[,,i])
    o[i,1]<-tab[1,1,i]
    o[i,2]<-h$expect
    o[i,3]<-h$ct
    o[i,4]<-1*(tab[1,1,i]>=h$ct)
    o[i,5]<-h$cexpect
    o[i,6]<-h$cvar
    o[i,7]<-(tab[1,1,i]*tab[2,2,i])/(tab[1,2,i]*tab[2,1,i])
    if (o[i,4]==1) g<-gconv(g,h$gc)
  }
  o<-as.data.frame(o)
  if (!is.null(dimnames(tab)[3])) rownames(o)<-dimnames(tab)[[3]]
  o[,2]<-round(o[,2],rnd)
  o[,5]<-round(o[,5],rnd)
  o[,6]<-round(o[,6],rnd)
  o[,7]<-round(o[,7],rnd)

  tstat<-sum(o$use*o$n11)
  if (sum(o$use==1)==0){
    list(tstat=tstat,pval=1,detail=o)
  }
  else{
    pval<-sum(g[(tstat+1):length(g)])
    if (pval<(-0.0001)) stop("dstat2x2xk error: negative P-value")
    pval<-max(0,pval)
    list(tstat=tstat,pval=pval,detail=o)
  }
}


