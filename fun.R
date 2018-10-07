


lags <- function(maxlag, x){
   xtem = list()
   nn = vector()
   if (is.null(dim(x))) {
    for (i in 1:maxlag) {
      x.lag = x[((maxlag + 1 - i) : (length(x) - i))]
      xtem[[i]] = x.lag
      nn[i] = paste0( deparse(substitute(x)) , ".lag", i)
    }
     names(xtem) = nn
    lagged.x = do.call(cbind, xtem)
   }
  else {
    for (i in 1:maxlag) {
      x.lag = x[((maxlag + 1 - i) : (dim(x)[1] - i)) , ]
      colnames(x.lag) = paste0( colnames(x.lag) , ".lag", i)
      xtem[[i]] = x.lag
    }
    lagged.x = do.call(cbind, xtem)
  }
  return(lagged.x)
}


lag0 <- function(maxlag, y){
  if (is.null(dim(y))) {
    yt = y[((maxlag + 1) : length(y))]
    names(yt) = paste0( deparse(substitute(y)) )
  }
  else {
    yt = y[((maxlag + 1) : (dim(y)[1])) , ]
    colnames(yt) = colnames(y)
  }
  return(yt)
}


lag1<- function(maxlag, y){
  if (is.null(dim(y))) {
    yt = y[( (maxlag) : (length(y)-1) )]
    names(yt) = paste0( deparse(substitute(y)) , ".lag1")
  }
  else {
    yt = y[( (maxlag) : (dim(y)[1]-1) ) , ]
    colnames(yt) = paste0( colnames(y) , ".lag1")
  }
  return(yt)
}












trace.plot.table <- function(model, lambda.seq, lambda.fix) {
  plot.lasso <- rbind(lambda.seq, coef(model)[-1,]) %>% t() %>%
    as.matrix() %>% as.data.frame(stringsAsFactors = F)
  colnames(plot.lasso)[1] <- "lambda"
  melt.lasso <- melt(plot.lasso, id = "lambda")
  coef.fix <- subset(melt.lasso, lambda == 
                       lambda[which.min(abs(lambda - lambda.fix))] &
                       value != 0)
  rownames(coef.fix) <- c()
  plot <- ggplot(melt.lasso, aes(x= lambda, y= value, color= variable))+
    geom_line() + theme(legend.position = "none") +
    labs(x='lambda', y='coefficients') + 
    geom_text_repel(data=subset(melt.lasso, lambda.seq == 
                                  lambda.seq[which.min(abs(lambda.seq-lambda.fix))] &
                                  value != 0), aes(label= variable, color=variable)) + 
    ggtitle("Coeffcients trace plot") 
  table = kable(coef.fix[,(2:3)],
                caption = "Non-zero Coefficients with fixed lambda")
  list = list("plot"=plot, "table"= table, "coef" = coef.fix)
  return(list)
}








