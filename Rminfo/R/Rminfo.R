autoload("read.csv.sql",package="sqldf")

Rminfo <- setRefClass("Rminfo",
  fields=list(df="data.frame",summary="integer",names="character",.factors="list"),
  methods=list(
    initialize=function(...) {
      args <- list(...)
      if(!is.null(args$memory.db)) dbname = NULL
      else                         dbname = tempfile()
      if(length(args) > 0) {
        if(is.null(names(args))) {
          df <<- as.data.frame(args,stringsAsFactors=FALSE)
        } else if(!is.null(args$filter)) {
          if(is.null(args$field.types)) {
            n <- read.csv.sql(filter=sprintf("%s | head -1",args$filter))
            f.types <- rep("TEXT",length(n))
            names(f.types) <- names(n)
          } else f.types = NULL
          df <<- read.csv.sql(filter=args$filter,field.types=f.types,dbname=dbname)
        } else if(!is.null(args$file)) {
          if(is.null(args$field.types)) {
            n <- read.csv.sql(file=args$file,nrows=1)
            f.types <- rep("TEXT",length(n))
            names(f.types) <- names(n)
          } else f.types = NULL
          df <<- read.csv.sql(file=args$file,field.types=f.types,dbname=dbname,header=FALSE)
          if(is.null(args$field.types)) names(df) <<- names(n)
        } else if(!is.null(args$X)) {
          df <<- as.data.frame(args$X,stringAsFactors=FALSE)
        }
        for(v in 1:ncol(df))  if(!is.character(df[,v])) df[,v] <<- as.character(df[,v])

        summary <<- dim(df)
        names   <<- names(df)
      }

      gc()
    },
    factors=function(...) {
      args <- list(...)
      n  <- 1:ncol(df); names(n) <- names(df)
      if(is.null(args$ivars)) args$ivars <- args$ivar

      if(is.null(args$svar) | is.null(args$ivars)) stop("svar and ivars fields required")

      if(is.character(args$svar)) args$svar <- n[args$svar]
      args$svar <- as.integer(args$svar -1)

      if(is.character(args$ivars)) args$ivars <- n[args$ivars]
      args$ivars <- as.integer(args$ivars -1)

      args$vars        <- as.integer(if(is.null(args$vars))        15       else args$vars)
      args$samples     <- as.integer(if(is.null(args$samples))     nrow(df) else args$samples)
      args$verbose     <- as.integer(if(is.null(args$verbose))     1        else args$verbose)
      args$threads     <- as.integer(if(is.null(args$threads))     1        else args$threads)
      args$iterations  <- as.integer(if(is.null(args$iterations))  1        else args$iterations)

      .factors <<- .Call("R_pminfo_factors",df,args)
      if(args$verbose) { return(invisible(.factors)) }
      else             { return(.factors) }
    }
  )
)


#setMethod("minfoFactors",signature="Rminfo", function(args,...) {
#  n  <- 1:ncol(df); names(n) <- names(df)
#  
#  args <- as.list(args)
#
#  if(is.character(args[["svar"]])) args[["svar"]] <- n[args[["svar"]]]
#  args[["svar"]] <- as.integer(args[["svar"]] -1)
#
#  if(is.character(args[["ivars"]])) args[["ivars"]] <- n[args[["ivars"]]]
#  args[["ivars"]] <- as.integer(args[["ivars"]] -1)
#
#  args[["var_cnt"]]     <- as.integer(if(is.null(args[["var_cnt"]]))     15       else args[["var_cnt"]])
#  args[["sample_size"]] <- as.integer(if(is.null(args[["sample_size"]])) nrow(df) else args[["sample_size"]])
#  args[["verbose"]]     <- as.integer(if(is.null(args[["verbose"]]))     1        else args[["verbose"]])
#  args[["threads"]]     <- as.integer(if(is.null(args[["threads"]]))     1        else args[["threads"]])
#  args[["iterations"]]  <- as.integer(if(is.null(args[["iterations"]]))  1        else args[["iterations"]])
#
##  print(args)
##  for(v in args) print(class(v))
#
#  print(".Call R_pminfo_factors")
#  factors <<- .Call("R_pminfo_factors",df,args)
#  if(args[["verbose"]]) { return(invisible(factors)) }
#  else                  { return(factors) }
#})

minfo <-  function(svar,ivars,count=0,val=1,progress=FALSE,entropy=FALSE,...) {
  for(j in 1:ncol(ivars)) {
    ivars[,j] <- as.factor(ivars[,j])
  }

  rv <- .Call("R_minfo",
      svar,
      ivars,
      as.integer(count),
      as.real(val),
      as.integer(progress),
      as.integer(entropy))
  if(progress & count>0) return(invisible(rv))
  return(rv)
}


minfo_old <-  function(svar,ivars,iterate=FALSE,progress=TRUE,count=999999,val=1,...) {
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

minfo_old_old <- function(this,svars,ivars,...) {
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
