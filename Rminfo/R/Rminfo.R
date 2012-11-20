#
#

minfo_factors <- function(df,args) {
  df <- as.data.frame(df)
  n  <- 1:ncol(df); names(n) <- names(df)
  
  args <- as.list(args)

  if(is.character(args[["svar"]])) args[["svar"]] <- n[args[["svar"]]]
  df[,args[["svar"]]] <- as.character(df[,args[["svar"]]])
  args[["svar"]] <- as.integer(args[["svar"]] -1)

  if(is.character(args[["ivars"]])) args[["ivars"]] <- n[args[["ivars"]]]
  for(v in args[["ivars"]]) df[,v] <- as.character(df[,v])
  args[["ivars"]] <- as.integer(args[["ivars"]] -1)

  args[["var_cnt"]]     <- if(args[["var_cnt"]])     as.integer(args[["var_cnt"]])     else 15
  args[["sample_size"]] <- if(args[["sample_size"]]) as.integer(args[["sample_size"]]) else nrow(df)
  args[["verbose"]]     <- if(args[["verbose"]])     as.integer(args[["verbose"]])     else 1
  args[["threads"]]     <- if(args[["threads"]])     as.integer(args[["threads"]])     else 1
  args[["iterations"]]  <- if(args[["iterations"]])  as.integer(args[["iterations"]])  else 1

#  print(args)
#  for(v in args) print(class(v))

  rv = .Call("R_pminfo_factors",df,args)
  if(args[["verbose"]]) { return(invisible(rv)) }
  else                  { return(rv) }
}

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
