#Derived from parcoord function

#Basic function declaration, color=black, line type=1, default is False for label
perf_analysis<-function (x, col = 1, lty = 1, xlab, ylab,main, ...){
  matplot(1L:ncol(x), t(x), type = "l", lwd=2, col = col, lty = lty, 
          xlab = xlab, ylab = ylab, axes = FALSE, ...)
  
  #X axis settings
  axis(1, at = 1L:ncol(x), labels = colnames(x),)
  axis(2)
  
  for (i in 1L:ncol(x)) {
    lines(c(i, i), c(0, max(x)+0.1*max(x)), col = "grey80")
    text(c(i,i),c(0,1), labels = format(x[1,i], digits = 3), 
           xpd = NA, offset = 3, pos = c(1, 1), cex = 0.8)}
  title(main)
  legend(locator(1),
           +        fill<-c("blue","pink"),        # the colors to fill legend boxes
           +        bty<-"n")
  invisible()
}