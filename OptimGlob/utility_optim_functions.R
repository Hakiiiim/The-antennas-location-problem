### utility functions for plotting and handling optim related data structures 
### R. Le Riche

plotfmat <- function(fmat,xlab="",ylab="", xlim=NULL, ylim=NULL,add=FALSE,icol=1){
  m <- apply(fmat, 2, mean)
  if (dim(fmat)[1]>1) {
    std <- apply(fmat, 2, sd)
  } else {std <- 0}
  low68 <- m-std
  up68 <- m+std
  x <- seq(1,length(m))
  if(is.null(xlim)) {
    xlim <- range(x)
  }
  if(is.null(ylim)) {
    ylim <- range(low68-0.5,up68+0.5)
  }

  if(!add){
    plot(x, m, type="n", xlab=xlab,ylab=ylab, xlim=xlim, ylim=ylim, cex.axis=1.5,cex.lab=2)
  }
  # polygon(c(x,rev(x)),c(up68,rev(low68)),border=NA,col=gray.colors(12)[icol])
  # polygon(c(x,rev(x)),c(up68,rev(low68)),border=NA,col="transparent")
  lines(x,m,col=icol,lwd=3)
  lines(x,low68,col=icol)
  lines(x,up68,col=icol)
}


hist2hmat <- function(h){
  nrows <- length(h)
  ncols <- 0
  for (i in 1:length(h)) {
    ncols <- max(length(h[[i]]$fhist),ncols)
  }
  hmat <- matrix(data=NA,nrow=nrows,ncol=ncols)
  for (i in 1:length(h)) {
    run_lgth <- length(h[[i]]$fhist)
    hmat[i,1:run_lgth] <- h[[i]]$fhist
    if (run_lgth<ncols) {
      hmat[i,(run_lgth+1):ncols] <- h[[i]]$fhist[run_lgth]
    }
  }
  return(hmat)
}
