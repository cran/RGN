#******************************************************************
#
# Purpose: Optimize weighted sum-of-squares objective function with Robust Gauss-Newton algorithm
#
# References
# * Qin2018a: Youwei Qin, Kavetski Dmitri, George Kuczera (2018),
#            A Robust Gauss-Newton algorithm for the optimization of hydrological models: From standard Gauss-Newton to Robust Gauss-Newton,
#            Water Resources Research, 54, https://doi.org/10.1029/2017WR022488

# * Qin2018b: Youwei Qin, Kavetski Dmitri, George Kuczera (2018),
#            A Robust Gauss-Newton algorithm for the optimization of hydrological models: Benchmarking against industry-standard algorithms,
#            Water Resources Research, 54, https://doi.org/10.1029/2017WR022489
#
#******************************************************************

# setup rgn variables

EPS=.Machine$double.eps

# Initialize the RGN converge variables
setDefaultRgnConvergeSettings=function(iterMax=NULL, dump=NULL, logFile=NULL){
  #optional arguments iterMax, dump, logFile
  cnvSet=setRgnConvType()
  #---
  #
  cnvSet$iterMax = 100; if(PRESENT(iterMax)) cnvSet$iterMax = iterMax
  cnvSet$noReduction = 4
  cnvSet$noRelChangeFTol = 1.0e-5
  cnvSet$noRelChangeF = 5
  cnvSet$noRelChangeParTol = 1.0e-5
  cnvSet$noRelChangePar = 5
  cnvSet$tolSafe = 1.0e-14
  cnvSet$dumpResults = 0; if(PRESENT(dump)) cnvSet$dumpResults = dump
  cnvSet$logFile = paste0(tempdir(),'rgnLog.txt'); if(PRESENT(logFile)) cnvSet$logFile = logFile
  return(cnvSet)
}

# Initialize the RGN constants
setRgnConstants=function(alpha=NULL, beta=NULL, nls=NULL){
  # NOTE - Side-effect - global variable via <<- operator
  #set = rgnSetType
  rgnConstants = list()
  rgnConstants$gtol = 1.0e-10
  rgnConstants$stol = 1.0e-11
  rgnConstants$ftol = 1.0e-10
  rgnConstants$gtolMin = 0.1
  rgnConstants$tol = 0.1
  rgnConstants$tolFor = 0.1
  rgnConstants$alpha = 0.001; if(PRESENT(alpha)) rgnConstants$alpha = alpha
  rgnConstants$sigma = 1.0
  rgnConstants$c1 = 0.0001
  rgnConstants$rho = 0.6
  rgnConstants$nls = 4; if(PRESENT(nls)) rgnConstants$nls = nls
  rgnConstants$beta = 10.0; if(PRESENT(beta)) rgnConstants$beta = beta
  rgnConstants$hLow = 1.0e-8
  rgnConstants$hHiFrac = 0.5
  rgnConstants$xScale = 10.0
  rgnConstants$fScale = 1.0
  return(rgnConstants)
} # END setRgnConstants


setRgnConvType = function(){
  rgnConvType=list()
  rgnConvType$iterMax=0
  rgnConvType$noReduction=0
  rgnConvType$noRelChangeF=0
  rgnConvType$noRelChangePar=0
  rgnConvType$fail=0
  rgnConvType$noRelChangeFTol=0.0
  rgnConvType$noRelChangeParTol=0.0
  rgnConvType$tolSafe=0.0
  rgnConvType$dumpResults=0
  rgnConvType$logFile=""
  return(rgnConvType)
}

setRgnInfoType = function(){
  rgnInfoType=list()
  rgnInfoType$nIter=0
  rgnInfoType$termFlag=0
  rgnInfoType$nEval=0
  rgnInfoType$f=0.0
  rgnInfoType$cpuTime=0.0
  rgnInfoType$objTime=0.0
  return(rgnInfoType)
}

#############################################
# mimic F90 functionality

PRESENT=function(x=NULL){
  ans=!is.null(x)
  return(ans)
}
SIZE=function(x,k=NULL){
  if(PRESENT(k)){
    return(dim(x)[k]) # matrix
  }else{
    return(length(x)) # vector
  }
}
DOT_PRODUCT=function(u,v){
  return(sum(u*v))
}
MATMUL=function(A,B){
  if(is.matrix(A)&is.matrix(B)){
    return(A%*%B) # both matrices - so proceed with %*%
  }else{
    if(length(A)==1|length(B)==1){
      return(A*B) # a scalar involved, so use regular *
    }else{
      return(A%*%B) # vector/matrix - %*% seems ok
    }
  }
}
TRANSPOSE=function(A){
  return(t(A))
}
MERGE=function(A,val,mask){
  # this implementation manually handles vector or matrix
  if(is.vector(mask)){ # scalar or vector
    n=length(mask)
    val.arr=val
    A.arr=A
    if(length(val.arr)<n) val.arr=rep(val,n) # assume a scalar provided for val and coerce to a vector
    if(length(A.arr)<n) A.arr=rep(A,n) # assume a scalar provided for A and coerce to a vector
  }else{ # assume it is a matrix - we will not support array dim >2
    nr=nrow(mask)
    nc=ncol(mask)
    val.arr=val
    A.arr=A
    if(!is.matrix(val.arr)) val.arr=matrix(val,nr,nc) # assume a scalar provided for val and coerce to a matrix
    if(!is.matrix(A.arr)) A.arr=matrix(A,nr,nc) # assume a scalar provided for A and coerce to a matrix
  }
  #perform the merge
  ind=which(!mask)
  if(length(ind)>0) A.arr[ind]=val.arr[ind]
  return(A.arr)
}
MAXVAL=function(A,mask){max(A[mask])}
MINVAL=function(A,mask){min(A[mask])}
MIN=function(...){
  dots=list(...) # capture the list of inputs
  n=length(dots)
  nindiv=rep(0,n)
  for(i in 1:n)nindiv[i]=length(dots[[i]])
  mxn=max(nindiv) # max length of any inputted item
  if(mxn==1){ # simple
    return(min(...))
  }else{ # we need to perform elementwise
    out=rep(9999999999999,mxn)
    for(j in 1:mxn){ # repeat for all elements
      for(i in 1:n){ # loop over all inputted objects
        if(j<=nindiv[i]) out[j]=min(out[j],dots[[i]][j]) # protect against not all shapes being conformal
      }
    }
    return(out)
  }
}
MAX=function(...){
  dots=list(...) # capture the list of inputs
  n=length(dots)
  nindiv=rep(0,n)
  for(i in 1:n)nindiv[i]=length(dots[[i]])
  mxn=max(nindiv) # max length of any inputted item
  if(mxn==1){ # simple
    return(max(...))
  }else{ # we need to perform elementwise
    out=rep(-9999999999999,mxn)
    for(j in 1:mxn){ # repeat for all elements
      for(i in 1:n){ # loop over all inputted objects
        if(j<=nindiv[i]) out[j]=max(out[j],dots[[i]][j]) # protect against not all shapes being conformal
      }
    }
    return(out)
  }
}
ABS=function(x){abs(x)}

#############################################

goto1=function(procnam){
  return(list(error = 1, message = paste("f-",procnam,"RGN objFunc call failed")))
}

#############################################
# calculate residuals used in objective function

objFuncCall = function(simFunc,x,simTarget,weights,fixParVal=NULL,fixParLoc=NULL,fitParLoc=NULL,...){

  timeObj = vector(length=2)
  timeObj[1] = Sys.time()

  # deal with fixed and fitted pars
  if (!is.null(fixParLoc)){
    xAll = c()
    xAll[fixParLoc] = fixParVal
    xAll[fitParLoc] = x
  } else {
    xAll = x
  }

 sim = simFunc(x=xAll,...)
  # sim = simFunc(x=x,...)
  r = weights*(simTarget-sim)
  # r = r[!is.na(r)]
  f = sum(r^2)
  f = f/2.0

  timeObj[2] = Sys.time()
  timeFunc=timeObj[2]-timeObj[1]
  outObjFunc = list(r=r, f=f, sim=sim, timeFunc=timeFunc) # DM to do: add error number and message to this

  return(outObjFunc)

}

#############################################

#' @title Robust Gauss Newton optimization
#'
#' @description \code{rgn} performs optimization of weighted-sum-of-squares (WSS) objective function using the Robust Gauss Newton algorithm
#' @param simFunc is a function that simulates a (vector) response, with first argument the vector of parameters over which optimization is performed
#' @param simTarget is the target vector that \code{simFunc} is trying to match
#' @param weights is a vector of weights used in the WSS objective function. Defaults to equal weights.
#' @param par is the vector of initial parameters
#' @param lower is the lower bounds on parameters
#' @param upper is the upper bounds on parameters
#' @param control list of RGN settings
#' \itemize{
#' \item{\code{control$n.multi} is number of multi-starts
#'         (i.e. invocations of optimization with different initial parameter estimates). Default is 1.}
#' \item{\code{control$iterMax} is maximum iterations. Default is 100.}
#' \item{\code{control$dump} is level of diagnostic outputs between 0 (none) and 3 (highest). Default is 0.}
#' \item{\code{control$keep.multi} (TRUE/FALSE) controls whether diagnostic output from each multi-start is recorded. Default is FALSE.}
#' \item{\code{control$logFile} is log file name}
#' }
#' @param ... other arguments to \code{simFunc()}
#'
#' @details \code{rgn} minimizes the objective function \code{sum((weights*(simFunc-simTarget)^2))},
#' which is a sum of squared weighted residuals (\code{residuals=weights*(simFunc-simTarget)}).
#' Note \code{simFunc} corresponds to the vector of residuals when default
#' arguments for \code{simTarget} and \code{weights} are used.
#'
#' @return List with
#' \itemize{
#' \item{\code{par}, the optimal parameters}
#' \item{\code{value}, the optimal objective function value}
#' \item{\code{sim}, the simulated vector using optimal parameters}
#' \item{\code{residuals}, the vector of residuals using optimal parameters}
#' \item{\code{counts}, the total number of function calls}
#' \item{\code{convergence}, an integer code indicating reason for completion.
#' \code{1} maximum iterations reached,
#' \code{2} relative reduction in function value small.
#' \code{3} absolute reduction in function value small
#' \code{4} relative change in parameters small}
#' }
#' @examples
#' # Example 1: Rosenbrock
#' simFunc_rosenbrock=function(x) c(1.0-x[1],10.0*(x[2]-x[1]**2))
#' rgnOut = rgn(simFunc=simFunc_rosenbrock,
#'              par=c(-1.0,  0.0), lower=c(-1.5, -1.0), upper=c( 1.5,  3.0),
#'              simTarget=c(0,0))
#' rgnOut$par #optimal parameters
#' rgnOut$value #optimal objective function value
#'
#' # Example 2: Hymod
#' \donttest{
#' data("BassRiver") # load Bass River hydrological data
#' rgnOut = rgn(simFunc=simFunc_hymod,
#'              par=c(400.,0.5,0.1,0.2,0.1),
#'              lower=c(1.,0.1,0.05,0.000001,0.000001),
#'              upper=c(1000.,2.,0.95,0.99999,0.99999),
#'              simTarget=BassRiverData$Runoff.mm.day[365:length(BassRiverData$Date)],
#'              stateVal=c(100.0,30.0,27.0,25.0,30.0,0.0,0.0,0.0), # initial states for hymod
#'              nWarmUp=365,                                       # warmup period
#'              rain=BassRiverData$Rain.mm,                        # precip input
#'              pet=BassRiverData$ET.mm)                           # PET input
#' rgnOut$par #optimal parameters
#' rgnOut$value #optimal objective function value
#' }
#'
#' @useDynLib RGN
#' @export
#'
rgn = function(simFunc, simTarget=0, weights=NULL, par, lower, upper, control=NULL,...){

  if (is.null(control$n.multi)){
    n.multi = 1
  } else {
    n.multi = control$n.multi
  }

  if (is.null(control$keep.multi)){
    keep.multi=F
  } else {
    keep.multi=control$keep.multi
  }

  control.single = control; control.single$n.multi = control.single$keep.multi = NULL

  par.multi=par
  if (!is.null(par.multi)){
    if (is.vector(par.multi)){
      par.multi = matrix(par.multi,nrow=1)
    }
  }

  # number of initial parameter sets
  if (is.null(par.multi)){
    n.par.multi = 0
  } else {
    n.par.multi = dim(par.multi)[1]
  }

  f.best = 9e9; counts = 0
  multistarts = list()
  for (n in 1:n.multi){
    if (n <= n.par.multi){
      par = par.multi[n,]
    } else {
      set.seed(n)
      par = lower + stats::runif(length(lower))*(upper-lower)
    }
    rgn.single.out = rgn.fixPars(simFunc=simFunc, simTarget=simTarget, weights=weights,
                                 par=par,lower=lower,upper=upper,control=control.single,...)
    f.single = (rgn.single.out$value)
    pars.single = rgn.single.out$par
    if (f.single<f.best){
      f.best = f.single
      pars.best = pars.single
      n.best = n
    }
    counts = counts + rgn.single.out$count
    if (keep.multi){
      multistarts[[n]] = rgn.single.out
    }
  }

  if (n.multi==1){
    out = rgn.single.out
  } else {
    out = list(value=f.best,par=pars.best,counts=counts,n.best=n.best)
    if (keep.multi){out$multistarts=multistarts}
  }

  return(out)
}

#############################################
# this function allows fixed parameters to be used (that have same lower and upper bound)
rgn.fixPars = function(simFunc, simTarget, weights=NULL, par, lower, upper, control=NULL, ...){

  fixParLoc = fixParVal = fitParLoc = c()
  for (i in 1:length(par)){
    if (lower[i]==upper[i]){
      fixParLoc = c(fixParLoc,i)
      fixParVal = c(fixParVal,lower[i])
    } else{
      fitParLoc = c(fitParLoc,i)
    }
  }

  tmp = rgn.single(simFunc, simTarget, weights=NULL,
                   par=par[fitParLoc], lower=lower[fitParLoc], upper=upper[fitParLoc],
                   control=control,
                   fixParLoc = fixParLoc, fixParVal = fixParVal, fitParLoc = fitParLoc, ...)

  par = tmp$par
  parAll = c()
  parAll[fixParLoc] = fixParVal
  parAll[fitParLoc] = par

  tmp$par = par

  return(tmp)

}

##########################################################
# running rgn for a single set of initial parameters

rgn.single = function(simFunc, simTarget, weights=NULL, par, lower, upper, control=NULL, ...){

  info = setRgnInfoType()

  cnv = setDefaultRgnConvergeSettings(iterMax=control$iterMax,
                                      dump=control$dump,
                                      logFile=control$logFile)

  x0 = par; xLo = lower; xHi = upper
  p = length(x0)
  n = length(simTarget[!is.na(simTarget)])
  if(is.null(weights)){weights=rep(1,n)}
  nIter=0; i=0; j=0; k=0; m=0; nrls=0; nf=0; iMax=0; nr=0; termCode=0; noReduction=0; noRelChangeF=0; noRelChangePar=0
  forceRelease=FALSE; flag_ls=FALSE; xist=FALSE
  f=0.0; fl=0.0; fh=0.0; fBest=0.0; gMaxFree=0.0; gMaxThawn=0.0; minSingFrac=0.0; sig=0.0; fredExp=0.0; ft=0.0; fls=0.0; cons=0.0; scaledG=0.0; scaledS=0.0; scaledfExp=0.0; scaledfAct=0.0; fredAct=0.0; fOldBest=0.0; maxRelPar=0.0; time=c(0.0,0.0); gradDf=0.0; hessDf=0.0
  #local parameters
  procnam="rgnMain"
  time4fcall=0.0;time4fcallAcc=0.0
  dfm=c("","","","")
  BF=0; BL=1; BFL=2; BH=3; BFH=4;NUL_CON=-1; GRAD_CON=0; SEARCH_CON=1; FRED_CON=2

  #----
  # CONTAINED functions - declared first in R
  updateBest=function(f, x, r,fBest){
    # input real x(:), r(:), f
    update=FALSE;xBest=NULL;rBest=NULL
    if(f < fBest){
      update=TRUE
      fBest = f
      xBest = x
      rBest = r
    }
    return(list(update=update,fBest=fBest,xBest=xBest,rBest=rBest))
  } #END updateBest
  updateHess=function(Hess, k){
    # input integer: k
    # output real Hess(:,:)
    diagK=Hess[k,k]
    Hess[k, ]=zero; Hess[ ,k]=zero
    Hess[k,k]=diagK
    return(Hess)
  } #END updateHess


  error=0                              # Initialize error flag
  message=""                           # Initialize message
  time4fcall=0.0;time4fcallAcc=0.0
  # Allocate work arrays
  h=rep(0.0,p)
  r=rep(0.0,n)
  rBest=rep(0.0,n)
  rl=rep(0.0,n)
  rh=rep(0.0,n)
  xl=rep(0.0,p)
  xh=rep(0.0,p)
  xBest=rep(0.0,p)
  Ja=matrix(0.0,n,p)
  g=rep(0.0,p)
  He=matrix(0.0,p,p)
  as=rep(0,p)
  xScale=rep(0.0,p)
  xp=rep(0.0,p)
  xt=rep(0.0,p)
  delX=rep(0.0,p)
  xls=rep(0.0,p)
  hLo=rep(0.0,p)
  hHi=rep(0.0,p)
  x0ldbest=rep(0.0,p)
  fOptSeries=c()
  delXAct=rep(0.0,p)

  if(cnv$dumpResults >= 1){ # Fortran formats are not relevant in R, need to swap over to fortmat() function
    dfm[1]=paste('(a,', p,'g15.7)',sep="")
    dfm[2]=paste('(a,', p,'i3)',sep="")
    dfm[3]=paste('(33x,', p,'g15.7)',sep="")
    dfm[4]=paste('(a,', p,'(i4,11x))',sep="")
    #write(dfm,file=cnv$logFile)
    write("R output",file=cnv$logFile)
  }

  # Assign constants
  set = setRgnConstants()
  xScale = rep(set$xScale,p)
  hLo = set$hLow; hHi =  set$hHiFrac*(xHi-xLo)
  #
  # Initialize
  as[1:p] = rep(BF,p)  # Assume all initial parameters are free
  nIter = 0
  h = hHi
  forceRelease = FALSE
  noReduction = 0; noRelChangeF = 0; noRelChangePar = 0
  info$termFlag = 0; info$nEval = 0
  time[1]=Sys.time()
  x = x0
  tmp=objFuncCall(simFunc=simFunc,x=x,simTarget=simTarget,weights=weights,...); f=tmp$f;rBest=tmp$r;time4fcall=tmp$timeFunc  #SUB2FUNC conversion

  info$nEval = info$nEval + 1; if(error !=0) return(goto1(procnam))
  fBest = f; xBest = x
  time4fcallAcc=time4fcallAcc+time4fcall
  #CALL userRunTimeMessage ('Starting RGN', -1)
  #-------------
    # RGN iteration loop
  iterLoop=TRUE
  while(iterLoop){
    nIter = nIter + 1
    #
    # Save best result from previous iteration
    fOldBest = fBest; x0ldBest = xBest

    if(cnv$dumpResults >= 1) {
      write('-----------------------------------------',file=cnv$logFile,append=TRUE)
      write(paste('Iteration No.=                 ', nIter),file=cnv$logFile,append=TRUE)
      write(paste('ObjFun Calls=                  ', info$nEval),file=cnv$logFile,append=TRUE)
      write(paste('ObjFun Value f=                 ', f),file=cnv$logFile,append=TRUE)
      write(paste('Parameter set=                  ', paste(x,collapse=" ")),file=cnv$logFile,append=TRUE)
      write(paste('Bound Index=                    ', paste(MERGE(-1, MERGE(0, 1, x[1:p] <= xHi[1:p]), x[1:p] < xLo[1:p]),collapse=" ")),file=cnv$logFile,append=TRUE)
      write(paste('Sampling Scale h=               ', paste(h,collapse=" ")),file=cnv$logFile,append=TRUE)
    }
    #
    # Get Jacobian and update best function result
    xh = x; xl = x; r = rBest

    for(k in 1:p){
      xh[k] = x[k] + h[k]; xh[k] = MIN(xHi[k], xh[k])
      if(cnv$dumpResults >= 2) write(paste('Forward Jacobian sample point:   ', paste(xh,collapse=" ")),file=cnv$logFile,append=TRUE)
      tmp=objFuncCall(simFunc=simFunc,x=xh,simTarget=simTarget,weights=weights,...); fh=tmp$f;rh=tmp$r;time4fcall=tmp$tFunc #SUB2FUNC conversion
      info$nEval = info$nEval + 1; if(error !=0) return(goto1(procnam)); tmp=updateBest(fh, xh, rh,fBest);if(tmp$update){fBest=tmp$fBest;xBest=tmp$xBest;rBest=tmp$rBest} #SUB2FUNC conversion
      xl[k] = x[k] - h[k]; xl[k] = MAX(xLo[k], xl[k])
      time4fcallAcc=time4fcallAcc+time4fcall
      if(cnv$dumpResults >= 2) write(paste('Backward Jacobian sample Point: ', paste(xl,collapse=" ")),file=cnv$logFile,append=TRUE)
      tmp=objFuncCall(simFunc=simFunc,x=xl,simTarget=simTarget,weights=weights,...); fl=tmp$f;rl=tmp$r;time4fcall=tmp$tFunc  #SUB2FUNC conversion
      info$nEval = info$nEval + 1; if(error !=0) return(goto1(procnam)); tmp=updateBest(fl, xl, rl,fBest);if(tmp$update){fBest=tmp$fBest;xBest=tmp$xBest;rBest=tmp$rBest} #SUB2FUNC conversion
      time4fcallAcc=time4fcallAcc+time4fcall
      Ja[ ,k] = (rh-rl)/(xh[k]-xl[k])
      xh[k] = x[k]; xl[k] = x[k]
      if(cnv$dumpResults >= 2) write(paste('Jacobian matrix column:          ', k, fh, fl),file=cnv$logFile,append=TRUE)
    }

    if (sum(abs(Ja))==0){
      info$termFlag = 5; break
    }

    #
    # Calculate gradient and Hessian
    for(i in 1:p){
      g[i] = DOT_PRODUCT(Ja[ ,i],r)
      for(j in 1:i){
        He[i,j] = DOT_PRODUCT(Ja[ ,i],Ja[ ,j])
        He[j,i] = He[i,j]
      }
    }
    #
    # Perform active set update
    for(k in 1:p){
      if(x[k] <= xLo[k] + 10.0*EPS*MAX(ABS(xLo[k]),xScale[k])){
        as[k] = MERGE(BFL, BL, g[k] < 0.0)
        He=updateHess(He, k)
      }else if(x[k] >= xHi[k] - 10.0*EPS*MAX(ABS(xHi[k]),xScale[k])){
        as[k] = MERGE(BFH, BH, g[k] > 0.0)
        He=updateHess(He, k)
      }else{
        as[k] = BF
      }
    }

    if(cnv$dumpResults >= 2) {
      write(paste('Best objective function value f=', fBest),file=cnv$logFile,append=TRUE)
      write(paste('Best parameter x=               ', paste(xBest,collapse=" ")),file=cnv$logFile,append=TRUE)
      write(paste('Gradient g at parameter x=       ', paste(g,collapse=" ")),file=cnv$logFile,append=TRUE)
      write(paste('Hessian at x=                    ', He[1,1]),file=cnv$logFile,append=TRUE)
      if (p>1){
        for(j in 2:p){
          write(He[j,1:j],file=cnv$logFile,append=TRUE)
        }
      }
    }
    #
    # Determine termination code and hence status of forcRelease
    if(nIter > 1){
      termCode = NUL_CON
      scaledG = 0.0
      for(k in 1:p){
        if(as[k] == BF | as[k] == BFL | as[k] == BFH){
          scaledG = MAX (scaledG, ABS(g[k])*MAX(ABS(x[k]),xScale[k])/MAX(f,set$fScale))
        }
      }
      scaledS = 0.0
      for(k in 1:p){
        scaledS = MAX (scaledS, ABS(delXAct[k])/MAX(ABS(x[k]),xScale[k]))
      }
      scaledfExp = ABS(fredExp)/MAX(f,set$fScale)
      scaledfAct = ABS(fredAct)/MAX(f,set$fScale)

      if(scaledG <= set$gtol){
        termCode = GRAD_CON
      }else if(scaledS <= set$stol & scaledG <= set$gtolMin) {
        termCode = SEARCH_CON
      }else if(scaledfExp <= set$ftol & scaledfAct <= set$ftol & scaledG <= set$gtolMin) {
        termCode = FRED_CON
      }

      nf = sum(MERGE(1, 0, as == BF))
      if(nf == 0) {
        forceRelease = TRUE
      }else if(termCode != NUL_CON) {
        forceRelease = TRUE
      }else{
        forceRelease = FALSE
      }
    }

    #
    # Check conditions for releasing parameters
    nrls = sum(MERGE(1, 0, as == BFL | as == BFH))
    nf = sum(MERGE(1, 0, as == BF))
    if(nrls > 0) {
      if(nf > 0) {
        gMaxFree = MAXVAL(ABS(g)*MAX(ABS(x),xScale), mask = as == BF)
        for(k in 1:p){
          if(ABS(g[k])*MAX(ABS(x[k]),xScale[k])>= set$tol*gMaxFree & (as[k] == BFL | as[k] == BFH)) {
            as[k] = BF
          }
        }
      }
      if(forceRelease) {
        iMax = 0
        gMaxThawn = 0.0
        for(k in 1:p){
          if(ABS(g[k]*MAX(ABS(x[k]),xScale[k])) > gMaxThawn & (as[k] == BFL | as[k] == BFH)) {
            gMaxThawn = ABS(g[k]*MAX(ABS(x[k]),xScale[k])); iMax = k
          }
        }
        if(iMax > 0) as[iMax] = BF
        for(k in 1:p){
          if(ABS(g[k]*MAX(ABS(x[k]),xScale[k])) > set$tolFor*gMaxThawn & (as[k] == BFL | as[k] == BFH)) as[k] = BF
        }
      }
    }

    # Solve normal equations after removing non-free parameters
    if(cnv$dumpResults >= 2) write(paste('Active set=                     ', paste(as,collapse=" ")),file=cnv$logFile,append=TRUE)
    nr = sum(MERGE(rep(1,length(as)), 0, as == BF))

    HeRdc=matrix(0.0,nr,nr); delXRdc=rep(0.0,nr); gRdc=rep(0.0,nr); tsv=rep(0.0,nr)

    j = 0
    for(k in 1:p){
      if(as[k] == BF) {
      j = j + 1; gRdc[j] = g[k]
      m = 0
      for(i in 1:p){
        if(as[i] == BF) {
          m = m + 1; HeRdc[m,j] = He[i,k]
        }
      }
      HeRdc[j, ] = HeRdc[ ,j]
      }
    }
    minSingFrac = set$alpha*sqrt(EPS)

    tmp=svdSolve(m=nr, n=nr, A=HeRdc, b=-gRdc, x=delXRdc, minSingFrac=minSingFrac);delXRdc=tmp$x;tsv=tmp$tS;error=tmp$error;message=tmp$message #SUB2FUNC conversion

    if(error !=0) return(goto1(procnam))
    j = 0
    for(k in 1:p){
      if(as[k] == BF) {
      j = j + 1; delX[k] = delXRdc[j]
      }else{
        delX[k] = 0.0
      }
    }
    if(cnv$dumpResults >= 1) write(paste('Truncated SV=                   ', paste(tsv,collapse=" ")),file=cnv$logFile,append=TRUE)
    #
    # Project search step onto box constraints
    if(cnv$dumpResults >= 2) write(paste('SVD delX=                       ', paste(delX,collapse=" ")),file=cnv$logFile,append=TRUE)
    xt = x + delX
    xp = MIN(xHi, MAX(xLo, xt))
    delX = xp - x
    if(cnv$dumpResults >= 2) {
      write(paste('Projected delX=                 ', paste(delX,collapse=" ")),file=cnv$logFile,append=TRUE)
      write(paste('Projected xp=                   ', paste(xp,collapse=" ")),file=cnv$logFile,append=TRUE)
    }
    #
    # Update delXRdc for calculation fredExp
    j = 0
    for(k in 1:p){
      if(as[k] == BF) {
        j = j + 1
        delXRdc[j] = delX[k]
      }
    }
    #
    # Calculate the expected function reduction with constrained delXRdc
    gradDf = DOT_PRODUCT(gRdc,delXRdc)
    hessDf = DOT_PRODUCT(delXRdc,MATMUL(HeRdc,delXRdc))
    fredExp = gradDf + 0.5*hessDf

    #
    # Perform inexact line search
    sig = MIN (set$sigma, 1.0)
    cons = set$c1*MIN(0.0, DOT_PRODUCT(delX,g))
    if(cnv$dumpResults >= 3) write(paste('Cons=                            ', cons),file=cnv$logFile,append=TRUE)
    flag_ls = FALSE
    for(i in 0: set$nls){
      xt = x + sig*delX
      tmp=objFuncCall(simFunc=simFunc,x=xt,simTarget=simTarget,weights=weights,...);rl=tmp$r;ft=tmp$f;time4fcall=tmp$tFunc #SUB2FUNC conversion
      info$nEval = info$nEval + 1; if(error !=0) return(goto1(procnam))
      time4fcallAcc=time4fcallAcc+time4fcall
      if(cnv$dumpResults >= 3) {
        write(paste('xt=                             ', paste(xt,collapse=" ")),file=cnv$logFile,append=TRUE)
        write(paste('ft=                             ', ft),file=cnv$logFile,append=TRUE)
        write(paste( 'ft+sig=                         ',ft + sig*cons),file=cnv$logFile,append=TRUE)
      }
      if(ft < f + sig*cons) {
        xls = xt; fls = ft; flag_ls = TRUE
        if(cnv$dumpResults >= 1) write(paste('Line search successful at iteration', i ,' with sigma=', sig),file=cnv$logFile,append=TRUE)
        break
      }else{
        sig = set$rho*sig
      }
    }
    if(!flag_ls) {
      if(cnv$dumpResults >= 1) write('Line search failed',file=cnv$logFile,append=TRUE)
      fls = f; xls = x
    }
    fredAct = fls - f
    delXAct = xls - x # # Save actual search step for termination calculation
    #
    # Update best variables
    if(fBest < fls) {    # Jacobian evaluation produced better f than line search
      x = xBest; f = fBest
      flag_ls = TRUE
    }else{
      x = xls; f = fls
      #CALL updateBest (fls, xls, rl)
      tmp=updateBest(fls, xls, rl,fBest);if(tmp$update){fBest=tmp$fBest;xBest=tmp$xBest;rBest=tmp$rBest} #SUB2FUNC conversion
    }

    info$nEvalSeries[nIter] = info$nEval

    #
    # Store the value of best objective function for termination check
    fOptSeries[nIter] = fBest
    #
    # Update sampling scale
    if(flag_ls) {
      #         if(flag_ls & !(noReduction > cnv$fail | noRelChangeF > cnv$fail | noRelChangePar > cnv$fail)) {  #GAK enhance
      h = MIN (set$beta*h, hHi)
    }else{
      h = MAX (h/set$beta, hLo)
    }

    #
    # Check for convergence
    if(nIter >= cnv$iterMax) {
      info$termFlag = 1; break
    }
    if(nIter > 1) {
      noRelChangeF = 0
      for(k in MAX(1,nIter - cnv$noRelChangeF + 1):nIter){
        if(ABS((fOptSeries[k]-fOptSeries[nIter])/(fOptSeries[k]+cnv$tolSafe)) <= cnv$noRelChangeFTol) {
          noRelChangeF = noRelChangeF + 1
        }
      }
      if(noRelChangeF >= cnv$noRelChangeF) {
        info$termFlag = 2; break
      }

      noReduction = MERGE (noReduction+1, 0, f >= fOldBest)
      if(noReduction >= cnv$noReduction) {
        info$termFlag = 3; break
      }

      maxRelPar = -hugeRe
      for(k in 1:p){
        if(as[k] == BF) {
          maxRelPar = MAX (maxRelPar, ABS((x0ldBest[k]-x[k])/(x0ldBest[k]+cnv$tolsafe)))
        }
      }
      noRelChangePar = MERGE (noRelChangePar+1, 0, maxRelPar >= 0.0 & maxRelPar < cnv$noRelChangeParTol)
      if(noRelChangePar >= cnv$noRelChangePar) {
        info$termFlag = 4; break
      }
    }

  } #END iterLoop

  #
  # Save optional information
  time[2]=Sys.time(); info$cpuTime = time[2] - time[1];info$objTime=time4fcallAcc
  info$nIter = nIter; info$f = f; info$fOptSeries = fOptSeries
  if(cnv$dumpResults >= 1) {
    write(paste('>>>>> RGN ended with termination code: ', info$termFlag),file=cnv$logFile,append=TRUE)
    write(paste('      f=', info$f),file=cnv$logFile,append=TRUE)
    write(paste('      number of function calls:    ', info$nEval),file=cnv$logFile,append=TRUE)
    write(paste('      cpu time (sec):               ', info$cpuTime),file=cnv$logFile,append=TRUE)
  }

  OFout = objFuncCall(simFunc=simFunc,
                      x=x,
                      simTarget=simTarget,
                      weights=weights,...)

  return(list(par=x,
              value=info$f,
              sim=OFout$sim,
              residuals=OFout$r,
              info=info,
              convergence=info$termFlag,
              counts=info$nEval))
} # END rgn

# ----
svdSolve=function(m, n, A, b, x=NULL, Ainv=NULL, S=NULL, minSingFrac=NULL,minSingVal=NULL,tS=NULL,cn=NULL){
  # Solves Ax=b using SVD decomposition followed setting singular values to zero and then back substitution
  error = 0
  message = 'ok'
  i=0; j=0; k=0; nite=0
  U=matrix(0,m,n)
  SD=matrix(0,n,n)
  W=rep(n,0)
  V=matrix(0,n,n)
  tmp=rep(n,0)
  wMin=0.0

  #----
  # Check consistency of dimension
  if(SIZE(A,1) != m) { error = 1; message = 'm not same as assumed size in A(m,n)' }
  if(SIZE(A,2) != n) { error = 1; message = 'n not same as assumed size in A(m,n)' }
  if(PRESENT(b) & PRESENT(x)) {
    if(SIZE(b) != n)   { error = 1; message = 'n not same as assumed size in b(n)'   }
    if(SIZE(x) != n)   { error = 1; message = 'n not same as assumed size in x(n)'   }
  }
  if(PRESENT(S)) {
    if(SIZE(S) != n)   { error = 1; message = 'n not same as assumed size in S(n)'   }
  }
  if(error!=0) {
    return(list(x=x, Ainv=Ainv, S=S, tS=tS, error, message, error=error,message=message, minSingVal=minSingVal, cn=cn))
  }
    # Perform SVD of A

  # Singular value decomposition
  tmp=svdDecomp (a=A);U=tmp$u;SD=tmp$s;V=tmp$v;nite=tmp$nite #SUB2FUNC conversion

  # Dealing with U and V, in the opposite direction
  U=-U; V=-V
  for(n in 1:SIZE(V,1)) W[n] = SD[n,n]
  if(PRESENT(S)) S = W

  # Zero singular values
  if(PRESENT(minSingFrac)){
    wMin = minSingFrac
  }else{
    wMin = sqrt(EPS)
  }
  wMin = MAXVAL(W)*wMin
  if(PRESENT(minSingVal)) minSingVal = wMin
  W = MERGE(W, 0.0, W > wMin)
  if(PRESENT(tS)) tS = W
  if(PRESENT(cn)) cn = MAXVAL(W)/MINVAL(W)
  #
  # Get x using back substitution
  if(PRESENT(b) & PRESENT(x)) {
    tmp.out=svdBackSub (m=m, n=n, U=U, W=W, V=V, b=b, x=x, error=error, message=message);x=tmp.out$x;error=tmp.out$error;message=tmp.out$message # SUB2FUNC conversion
  }
  #
  # Get inverse
  if(PRESENT(Ainv)){
    if(m == n){
      for(i in 1:n){
        for(j in 1:i){
          for(k in 1:n){
            if(W[k] > 0.0) {
              tmp[k] = V[i,k]/W[k]
            }else{
              tmp[k] = 0.0
            }
          }
          Ainv[i,j] = DOT_PRODUCT(tmp[1:n],U[j,1:n])
          Ainv[j,i] = Ainv[i,j]
        }
      }
    }else{
      error = 1; message = 'cannot get inverse for non-square matrix A'
    }
  }
  #
  return(list(x=x, Ainv=Ainv, S=S, tS=tS, error, message, error=error,message=message, minSingVal=minSingVal, cn=cn))
} #END svdSolve

# Singular value decomposition
svdDecomp=function(a){
  # locals
  u = s = v = q1=matrix(0.0,SIZE(a,1), SIZE(a,2))
  u1=matrix(0.0,SIZE(a,1), SIZE(a,1))
  q = e = matrix(0.0,SIZE(a,2), SIZE(a,2))
  f=rep(0.0,SIZE(a,2))
  err=0.0
  n=0
  # init u,v,u1
  for(n in 1:SIZE(a,1)) u1[n,n] = 1.0
  for(n in 1:SIZE(a,2)) v[n,n] = 1.0

  tmp = .Fortran('Qr_f90',nRow=nrow(a),nCol=ncol(a),a=a,q=q1,r=s); q1=tmp$q; s=tmp$r

  u = MATMUL(u1, q1)
  tmp = .Fortran('Qr_f90',nRow=nrow(TRANSPOSE(s)),nCol=ncol(TRANSPOSE(s)),a=TRANSPOSE(s),q=q1,r=s);q=tmp$q;s=tmp$r
  v = MATMUL(v, q)
  # iterate while converged:
  nite = 1
  #while(TRUE){
  for(i in 1:1000){ # this gives a timeout option when we do not achieve convergence to precision - EPS
    tmp = .Fortran('Qr_f90',nRow=nrow(TRANSPOSE(s)),nCol=ncol(TRANSPOSE(s)),a=TRANSPOSE(s),q=q1,r=s);q=tmp$q;s=tmp$r
    u = MATMUL(u, q)
    tmp = .Fortran('Qr_f90',nRow=nrow(TRANSPOSE(s)),nCol=ncol(TRANSPOSE(s)),a=TRANSPOSE(s),q=q1,r=s);q=tmp$q;s=tmp$r
    v = MATMUL(v, q)
    # check the error:
    e = Triu(s)
    f = Diag_ele(s)
    err = Norm(as.vector(e))/ Norm(f) # RESHAPE conversion, carefully checked instance column from matrix to vector
    nite = nite + 1
    if(err < EPS) break
  }
  if(i>1000) warning(paste("convergence issue, only achieved:",err))
  return(list(u=u, s=s, v=v, nite=nite))
} #END svdDecomp

# L2-norm
Norm=function(x){
  # x is a real vector
  return(sqrt(sum(x**2)))
}
# Diagonal elements
Diag_ele=function(a){
  # A is a real matrix
  # v is a real vector

  i=0;n=0

  n = min(c(SIZE(a,1), SIZE(a,2)))
  v=rep(0.0,n)
  for(i in 1:n){v[i] = a[i,i]}
  return(v)
} #END Diag_ele
# Upper triangular part
Triu=function(a){
  # A is a real matrix
  # au is a real vector

  i=0;j=0;n=0;m=0

  m = SIZE(a,1)
  n = SIZE(a,2)
  au=matrix(0.0,m,n)
  for(i in 1:m){
    for(j in (i+1):n){
      if(i+1 <= n) au[i,j] = a[i,j]
    }
  }
  return(au)
} # END Triu

svdBackSub=function(m, n, U, W, V, b, x, error, message){
  # Solves Ax=b using SVD back substitution
  # Singular value decomposition of A(m,n) = U(m,n) * W(n) *Vtranspose (n,n)
  j=0 # local integer
  tmp=rep(0.0,n) # local real vector
  #----
  # Check consistency of dimension
  error = 0; message = 'ok'

  if(SIZE(U,1) != m) { message = 'm not same as assumed size in U(m,n)' ; error = 1}
  if(SIZE(U,2) != n) { message = 'n not same as assumed size in U(m,n)' ; error = 1}
  if(SIZE(W) != n)   { message = 'n not same as assumed size in W(n)'   ; error = 1}
  if(SIZE(V,1) != n) { message = 'n not same as assumed size in V(n,n1)'; error = 1}
  if(SIZE(V,2) != n) { message = 'n not same as assumed size in V(n1,n)'; error = 1}
  if(SIZE(b) != n)   { message = 'n not same as assumed size in b(n)'   ; error = 1}
  if(SIZE(x) != n)   { message = 'n not same as assumed size in x(n)'   ; error = 1}

  if(error==0){
    # Perform back substitution
    for(j in 1:n){
      if(abs(W[j]) > 0.0000000001){
        tmp[j] = DOT_PRODUCT(U[1:m,j],b[1:m])/W[j]
      }else{
        tmp[j] = 0.0
      }
    }
    for(j in 1:n){
      x[j] = DOT_PRODUCT(V[j,1:n],tmp[1:n])
    }
  }
  return(list(x=x,error=error,message=message))
} #END svdBackSub

