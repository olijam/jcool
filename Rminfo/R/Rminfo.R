#
#

minfo_factors <- function(df,args) {
  df$value <- as.data.frame(df$value)
  n  <- 1:ncol(df$value); names(n) <- names(df$value)
  
  args <- as.list(args)

  if(is.character(args[["svar"]])) args[["svar"]] <- n[args[["svar"]]]
  df$value[,args[["svar"]]] <- as.character(df$value[,args[["svar"]]])
  args[["svar"]] <- as.integer(args[["svar"]] -1)

  if(is.character(args[["ivars"]])) args[["ivars"]] <- n[args[["ivars"]]]
  for(v in args[["ivars"]]) df$value[,v] <- as.character(df$value[,v])
  args[["ivars"]] <- as.integer(args[["ivars"]] -1)

  args[["var_cnt"]]     <- as.integer(if(is.null(args[["var_cnt"]]))     15       else args[["var_cnt"]])
  args[["sample_size"]] <- as.integer(if(is.null(args[["sample_size"]])) nrow(df$value) else args[["sample_size"]])
  args[["verbose"]]     <- as.integer(if(is.null(args[["verbose"]]))     1        else args[["verbose"]])
  args[["threads"]]     <- as.integer(if(is.null(args[["threads"]]))     1        else args[["threads"]])
  args[["iterations"]]  <- as.integer(if(is.null(args[["iterations"]]))  1        else args[["iterations"]])

#  print(args)
#  for(v in args) print(class(v))

  print(".Call R_pminfo_factors")
  rv = .Call("R_pminfo_factors",df$value,args)
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
