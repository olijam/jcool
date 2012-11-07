#
#

minfo <- function(svar,ivars,count=0,val=1,progress=FALSE,entropy=FALSE) {
  ivars <- as.data.frame(ivars)
  for(j in 1:ncol(ivars)) {
    ivars[,j] <- as.factor(ivars[,j])
  }

  rv <- .Call("R_minfo",
      as.factor(svar),
      ivars,
      as.integer(count),
      as.real(val),
      as.integer(progress),
      as.integer(entropy))
  if(progress & count>0) return(invisible(rv))
  return(rv)
}


minfo.old <- function(svar,ivars,iterate=FALSE,progress=TRUE,count=999999,val=1) {
    rv <- .Call("R_minfo_old",
        as.matrix(ivars),
        as.vector(svar),
        as.integer(iterate),
        as.integer(progress),
        as.integer(count),
        as.real(val))
    if(iterate & progress) return(invisible(rv))
    return(rv)
}

minfo.old.old <- function(svars,ivars) {
    rows         <- nrow(ivars)
    cols         <- ncol(ivars)
    data         <- matrix(nrow=rows,ncol=cols)

    J            <- as.factor(do.call("paste",as.data.frame(svars)))
    state        <- as.integer(as.integer(J) - 1)

    for (j in 1:cols) {
        J                <- factor(ivars[,j],exclude=NULL)
        data[,j]         <- as.integer(as.integer(J) - 1)
    }

    minfo <- .Call("minfo_old_old",data,state);

    names(minfo) <- colnames(ivars);

    return(minfo);
}
