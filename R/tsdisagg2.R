#' Time Series Disaggregation.
#'
#' @description The "tsdisagg2" function performs temporal disaggregation or interpolation of low frequency to high frequency time series.
#'
#' @param y A data.frame, matrix, list or vector with low frequency data.
#' @param x A data.frame, matrix, list or vector with high frequency data.
#' @param c Constant; If c=1, the model will be estimated with a constant; If c=0, the opposite case. (Default: c=0)
#' @param method Set disaggregation method; Available methods are Boot, Feibes and Lisman (method="bfl1" or method="bfl2"), Chow and Lin (method="cl1" or method="cl2"), Fernandez (method="f") and Litterman (method="l"). Default: method="cl1"
#' @param s Frequency of observations; Available frequencies are 3, 4 or 12. For example, if s=4, we have quaterly observations. (Default: s=4)
#' @param type Type of restriction; Could be "last", "first", "sum" or "average". (Default: type="sum")
#' @param ML Maximum Likelihood (ML=1) or Generalised Least Squares (ML=0) "rho" estimation. (Default: ML=0)
#' @param rho Sets a value for "rho" (Default: rho=0)
#' @param neg If neg=1, will be tested negative for "rho"; If neg=0, only positive values will be tested. (Default: neg=0)
#' @param da First year considered on low frequency data.
#' @param dz Last year considered on low frequency data.
#' @param plots If plots=1, generates the plot of the estimated series and the plot withe Objective function values; If prints=0, the opposite case. (Default: plots=0)
#'
#' @details
#' The function is used  to disaggregate a low frequency to a higher frequency time series, while either the sum, the average, the first or the last value of the resulting high-frequency series is consistent with the low frequency series.
#' Implements the following methods for temporal disaggregation: Boot, Feibes and Lisman (first and second differences), Chow and Lin (independent and AR(1) errors), Fernandez and Litterman.
#' For Boot, Feibes and Lisman methods, the disaggregation can not be performed with help of indicator series. For the remaining methods, desaggregation can be performed with or without the help of one or more indicator series.
#' If the high-frequency indicator(s) cover(s) a longer time span than the low-frequency series, an extrapolation is performed, using the same model as for interpolation.
#'
#'
#' @return
#' The function prints details of the disaggregation (smooth, loglik, ...), estimated parameters (sigma_ols, sigma_gls, model coefficients, ...) and the disaggregated series. The function also returns an invisible list containing all the numeric results.
#'
#' @importFrom graphics grid legend par plot points
#' @importFrom stats pnorm
#'
#' @export
#'
#' @examples
#' anual <- runif( 19, 300, 455 )
#' indicators <- data.frame( runif( 76, 500, 700 ), runif( 76, 800, 980 ) )
#'
#' ### Constant ###
#'
#' tsdisagg2( y=anual, x=indicators, c=1, da=1995, dz=2013, plots=1 )
#' # Estimate model with constant
#'
#' ### Method selection ###
#'
#' tsdisagg2( y=anual, x=indicators, method="f", da=1995, dz=2013, plots=1 )
#' # Use option method
#'
#' ### "rho" value ###
#'
#' tsdisagg2( y=anual, x=indicators, method="cl2", da=1995, dz=2013, plots=1 )
#' # Search for positive optimal "rho" is enabled (if method="cl2" or method="l")
#'
#' tsdisagg2( y=anual, x=indicators, method="cl2", rho=0.35, da=1995, dz=2013, plots=1 )
#' # Set "rho" value manually (the grid search is not performed)
#'
#' ### Interpolation or distribution ###
#'
#' tsdisagg2( y=anual, x=indicators, da=1995, dz=2013, method="f", type="last" )
#' # Performs disaggregation by interpolation with type="last" or type="first"
#'
#' tsdisagg2( y=anual, x=indicators, da=1995, dz=2013, method="f", type="average" )
#' # Performs disaggregation by distribution with type="sum" or type="average"
#'
#' ### Use returned objects ###
#' td <- tsdisagg2( y=anual, x=indicators, da=1995, dz=2013, method="f", type="average" )
#' names(td)
#' td$BETA_ESTIMATION

tsdisagg2 <- function( y, x, c=0, method="cl1", s=4, type="sum", ML=0, rho, neg=0, da, dz, plots=0 ) {

  ########################################## Test and set inputs #################################################

  if( missing(y) ){
    stop(" List of periodic values 'Y' is empty! ")
  }

  mtd <- c( "bfl1", "bfl2", "cl1", "cl2", "f", "l" )

  if ( !(method %in% mtd) ) {
    stop(" Invalid 'method' input! Available methods are 'bfl1', 'bfl2', 'cl1', 'cl2', 'f' or 'l'! ")
  }

  if ( method %in% c("bfl1", "bfl2") & !missing(x) ) {
    stop(" Invalid input! Methods 'bfl1' and 'bfl2' can not use indicators! Check argument 'x' or argument 'method'! ")
  }

  if ( method==mtd[1] ) {
    search <- 0
    rho <- 0
    diff <- 1
    c <- 1
    t <- 0
    tit1 <- "Boot, Feibes and Lisman (first differences)"
  }

  if ( method==mtd[2] ) {
    search <- 0
    rho <- 1
    diff <- 1
    c <- 1
    t <- 1
    tit1 <- "Boot, Feibes and Lisman (second differences)"
  }

  if ( method==mtd[3] ) {
    search <- 0
    rho <- 0
    diff <- 0
    t <- 0
    tit1 <- "Chow and Lin (independent residuals)"
  }

  if ( method==mtd[4] ) {

    if ( missing(rho) ) {
      search <- 1
    }
    if ( !missing(rho) ) {
      search <- 0
    }

    diff <- 0
    t <- 0
    tit1 <- "Chow and Lin (AR1 residuals)"
  }

  if ( method==mtd[5] ) {
    rho <- 0
    search <- 0
    diff <- 1
    t <- 0
    tit1 <- "Fernandez"
  }

  if ( method==mtd[6] ) {
    if ( missing(rho) ) {
      search <- 1
    }
    if ( !missing(rho) ) {
      search <- 0
    }
    diff <- 1
    t <- 0
    tit1 <- "Litterman"
  }

  ct_inp <-c(c,t)
  ct <- as.logical(ct_inp)
  ct_p <- which(ct)

  if ( length(ct_p)==0 & missing(x) ) {
    stop(" Regressor list 'X' is empty! ")
  }

  s_test <- c(3,4,12)

  if ( !(s %in% s_test) ) {
    stop(" Invalid 's' input! This argument should be '3', '4' or '12'! ")
  }

  if ( !missing(rho) ) {
    if ( !is.numeric(rho) | rho < -1 | rho > 1 ) {
      stop(" Invalid 'rho' input! This argument should take a value within '-1' and '1'! ")
    }
  }

  type_test <- c( "sum", "average", "last", "first" )
  if ( !(type %in% type_test) ) {
    stop(" Invalid 'type' input! This argument should be 'sum', 'average', 'last' or 'first'! ")
  }


  if ( type=="sum" | type=="average" ) {
    inter <- 0
  }


  if ( type=="first" | type=="last" ) {
    inter <- 1
  }

  ################################### Y matrix ##################################################


  y.name <- "Y-hat"
  y.file <- as.matrix( data.frame(y) )
  data.y <- y.file

  y <- matrix( data.y[, ncol(y.file)], ncol=1 )
  y1 <- y


  if ( missing(da) | missing(dz) ) {
    if ( ncol(y.file)==1 ) {
      da <- 1
      dz <- nrow(y.file)
    }
  }

  a <- dz-da+1

  ################################### X matrix ##################################################

  if ( missing(x) ) {
    data.x <- c()
    x.names <- c()
    tri <- 0
    du <- dz
  }

  if ( !missing(x) ) {
    x1 <- x

    if ( class(x)!="character" ) {
      x.file <- as.matrix( data.frame(x) )
      data.x <- x.file
      x.names <- c()

      for ( i in 1:ncol(x.file) ) {
        x.names[i] <- c( paste0("X",i) )
      }
    }

    colnames(data.x) <- x.names
    data.x1 <- data.x

    tri <- ( nrow(data.x) - a * s ) %% s
    du <- dz + ( nrow(data.x) - a * s - tri ) / s
  }

  p <- du - dz
  pa <- p + a
  n <- pa * s + tri

  if ( length(ct_p)!=0 ) {
    ct <- cbind( rep(1, n), 1:n )
    colnames(ct) <- c( "c", "t" )
  }

  if ( length(ct_p)==0 ) {
    ct <- c()
  }

  x <- as.matrix( cbind( ct[,ct_p], data.x ) )

  ################################### Auxiliar matrixes ##################################################

  iff <- ML==1
  z <- diff!=1
  inn <- diag(n)
  l1 <- matrix(0,ncol=n,nrow=n)
  di <- diff==1

  e0 <- c()
  for (i in 2:n) {
    e0[i-1] <- (i-2)*n+i
  }

  l1[e0] <- 1

  d <- inn-(l1*di)
  ss <- inn-l1

  ################################### Get "C" matrix ##################################################

  if ( inter==0 & type=="sum" ) {
    v <- rep(1, s)
  }

  if ( inter==0 & type=="average" ) {
    v <- rep(1, s) / s
  }

  if ( inter==1 & type=="first" ) {
    v <- c( 1, rep(0, s-1) )
  }

  if ( inter==1 & type=="last" ) {
    v <- c( rep(0, s-1), 1 )
  }

  petri <- s * p + tri
  C <- kronecker( diag(a), t(v) )

  if ( petri!=0 ) {
    C <- cbind( C, matrix(0, ncol=petri, nrow=a) )
  }

  xc <-  C %*% x

  ################################### Get optimal Rho (input "search") ##################################################

  if ( search == 1 ) {

    ub <- 0.99 * (-1)^( neg == 1 )
    r <- seq( 0.09 * (-1)^( neg == 0 ), ub, by = 0.01* (-1)^( neg == 1 ))
    of <- matrix( NA, nrow = length(r), ncol = 2)

    for ( i in 1:length(r) ) {
      aux <- inn - l1 * r[i]
      aux[1,1] <- sqrt( 1 - z * r[i]^2 )
      omg1 <- t( aux %*% d ) %*% ( aux %*% d )
      detr <- ( C %*% solve( omg1, tol=1e-30 ) %*% t(C) )
      ea <- y - xc %*% solve( t( xc ) %*% solve( detr, tol=1e-30 ) %*% xc, tol=1e-30 ) %*% ( t( xc ) %*% solve( detr, tol=1e-30 ) %*% y )
      of[i,2] <- c( - iff * log( det( detr ) ) - a * log( t( ea ) %*% solve( detr, tol=1e-30 ) %*% ea ) )
      of[i,1] <- c( r[i] )
    }

    of <- data.frame(of)
    colnames(of) <- c( "rho", "value" )
    rho <- r[ which( of[,2]==max( of[,2] ) ) ]
  }

  ################################### Output "Y_hat" ##################################################

  aux <- inn - l1 * rho
  aux[1,1] <- sqrt( 1 - z * rho^2 )
  omg1 <- t( aux %*% d ) %*% ( aux %*% d )
  detr <- ( C %*% solve( omg1, tol=1e-30 ) %*% t(C) )
  bc <- solve( t( xc ) %*% solve( detr, tol=1e-30 ) %*% xc, tol=1e-30 ) %*% ( t( xc ) %*% solve( detr, tol=1e-30 ) %*% y )
  ea <- y - xc %*% bc
  yc <- x %*% bc + solve( omg1, tol=1e-30 ) %*% t(C) %*% solve( detr, tol=1e-30 ) %*% ea

  ssr_ols <- t(ea) %*% (ea)
  k <- nrow(bc)
  sigmaols <- sqrt(ssr_ols/(a-k))
  ssr_g <- ( t(ea) %*% solve(detr, tol=1e-30) %*% ea )
  sigma2 <- as.numeric( ssr_g / ( a - k ) )
  sigmagls <- sqrt( sigma2 )
  cov_bc <- sigma2 * solve( t(xc) %*% solve(detr, tol=1e-30) %*% xc, tol=1e-30 )
  se_bc <- as.matrix( sqrt( diag(cov_bc) ), ncol=1, nrow=k )
  smo <- t(yc) %*% t(ss) %*% t(ss) %*% ss %*% ss %*% yc / n
  z_statistic <- bc / se_bc
  ll <- -0.5 * ( log( det(detr) ) + a * log(a) ) - (a/2) * (1 + log(2*pi) + log(ssr_g))
  p_value <- c()

  for ( i in 1:k ) {
    p_value[i] <- c( ifelse( z_statistic[i]<0, 2 * pnorm( z_statistic[i], lower.tail=T ), 2 * pnorm( z_statistic[i], lower.tail=F ) ) )
  }

  tt_table <- cbind( bc, se_bc, z_statistic, p_value )
  colnames(tt_table) <- c( "estimate", "standard error", "z-statistic", "p-value" )
  rownames(tt_table) <- if ( length(ct_p)==0 ) { x.names } else { if ( length(data.x)!=0 ) { c( colnames(ct)[ct_p], x.names ) } else { colnames(ct)[ct_p] }  }

  yc_df1 <- c( rep(da:du, each=s), rep(du+1, times=tri) )

  if (tri!=0) {
    yc_df2 <- c( rep(1:s, times=length(da:du)) , 1:tri )
  }

  else {
    yc_df2 <- rep( 1:s, times=length(da:du) )
  }


  yc_df <- data.frame( yc_df1, yc_df2, yc )
  colnames( yc_df ) <- c( "year", "period", "y-hat" )

  ################################### Full output ##################################################


  lst <- list( RHO=rho, LOGLIK=ll, SIGMA_GLS=sigmagls, SIGMA_OLS=sigmaols, SMOOTH=sqrt(smo), BETA_ESTIMATION=tt_table, OBJECTIVE_FUNCTION=if (search==1) {of} else {NA}, Y_HAT=yc_df )

  tit2 <- c( "Maximum Likelihood 'rho' estimation", "Generalised Least Squares 'rho' estimation", "Previously fixed 'rho'" )
  tit2 <- if ( search==1 ) { if (ML==1) { tit2[1] } else { tit2[2] } } else { tit2[3] }

  lst <- lst[ c( 1:6, if ( search==1 ) { 7 }, 8 ) ]


  if ( length(data.x)!=0 ) {
    lst$X_DATA <- data.x1
  }

  lst$Y_DATA <- y1
  lst$RESIDUALS <- ea

  print.lst <- function(obj) {
    cat("\n")
    cat("METHOD: ", tit1, "\n")
    cat(paste0(rep("_",90), collapse = ""), "\n\n\n\n")
    cat("PARAMETERS ESTIMATION: ", "\n")
    cat(paste0(rep("_",90), collapse = ""), "\n\n")
    cat(paste0(tit2,": "), obj$RHO, "\n")
    cat("(Loglik: ", obj$LOGLIK, ")\n\n")
    cat("Sigma GLS: ", obj$SIGMA_GLS, "\n")
    cat("Sigma OLS: ", obj$SIGMA_OLS, "\n\n")
    cat("(Smooth: ", obj$SMOOTH, ")\n\n\n")
    cat("Model coefficients: \n")
    print( obj$BETA_ESTIMATION )
    cat("\n\n")
    cat(paste0("Estimated values (", y.name, "): \n"))
    print( obj$Y_HAT, row.names = FALSE)
    cat("\n\n\n")

    if ( "OBJECTIVE_FUNCTION" %in% names(lst) ) {
      cat("OBJECTIVE FUNCTION: \n")
      cat(paste0(rep("_",90), collapse = ""), "\n\n")
      print( obj$OBJECTIVE_FUNCTION, row.names = FALSE )
      cat("\n")
    }
  }

  if ( plots==1 ) {
    par( mfrow=c( search+plots, 1 ) )

    if (search==1) {
      label=paste0("(Rho = ",rho,")")
      plot(of[,1], of[,2], xlab=" Rho ", ylab="Objective function", main=tit2, type="l", col="lightblue4", lwd=2)
      points( rho, max(of$value), pch=21, col="lightblue4" )
      grid( col="grey" )
      legend("bottom", label, pch=21, col="lightblue4", bty="n", pt.lwd=2, pt.cex=0.9, cex = 0.6 )
    }

    x.plot <- seq( da, du+2, by=(1/s) )[1:n]
    plot( x.plot, yc, xlab="Year", ylab=y.name, main="Series estimation", type="l", col="chocolate", lwd=2)
    points( x.plot, yc, col="chocolate", pch=21 )
    grid( col="grey" )
  }

  print.lst(lst)
  return( invisible(lst) )
}
