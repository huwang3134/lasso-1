


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

out.mse <- function(x, y, lambda.seq, p){
  x.train = x[1:ceiling(0.8*dim(x)[1]) , ]
  y.train = y[1:ceiling(0.8*dim(x)[1])]
  x.test = x[ (ceiling(0.8*dim(x)[1])+1):dim(x)[1] , ]
  y.test = y[(ceiling(0.8*dim(x)[1])+1):dim(x)[1]]
  p = p
  n = dim(x.train)[1]
  lambda = sqrt(log( p ) / n )
  lasso = glmnet(x.train, y.train, alpha=1, thresh=1E-5, lambda= lambda.seq, maxit = 10^9)
  y.fitted = predict(lasso , s=lambda.seq[which.min(abs(lambda.seq-lambda))] , newx = x.test)
  mse = mean((y.test-y.fitted)^2)
  return(mse)
}


adf.auto <- function(data, criteria) {
  tb = data.frame(Series = character(), Conclusion = character(), Type = character() , 
                  Lags = numeric() , 
                  stringsAsFactors = FALSE)
  if (is.null(dim(data))) { # if it is a single series
    nm = deparse(substitute(data))
    data = na.omit(data)
    t = length(data)
    x = 1:t
    pmax = ceiling(12 * (t/100)^0.25 )
    llmm = lm(data ~ x) %>% summary.lm()
    if (llmm$coefficients[8] < 0.05) {
      type = "trend"
      rr = 3
    } else {
      type = "drift"
      rr = 2
    }
    adf = ur.df(data, selectlags = criteria, lags=pmax, type=type)
    calc = adf@teststat[1,1]
    if (calc < adf@cval[1,2]) {
      con = "I(0)"
    } else {
      con = "I(1)"
    }
    ll = dim(adf@testreg$coefficients)[1]-rr
    tb[1,] = c(nm, con, type, ll)
  } # if it is a data frame
  else {
    for (i in 1:dim(data)[2]){
      s = data[,i]
      s = na.omit(s)
      t = length(s)
      x = 1:t
      pmax = ceiling(12 * (t/100)^0.25 )
      llmm = lm(s ~ x) %>% summary.lm()
      if (llmm$coefficients[8] < 0.05) {
        type = "trend"
        rr = 3
      } else {
        type = "drift"
        rr = 2
      }
      adf = ur.df(s, selectlags = criteria, lags=pmax, type=type)
      calc = adf@teststat[1,1]
      if (calc < adf@cval[1,2]) {
        con = "I(0)"
      } else {
        con = "I(1)"
      }
      ll = dim(adf@testreg$coefficients)[1]-rr
      tb[i,] = c(colnames(data)[i], con, type, ll)
    }
  }
 return(tb)
}


lasso111 <- function(data, tcode){
  data.sw = data
  t.sw = tcode
  gdp <- data[,1]
  D.gdp <- diff(gdp) %>% as.numeric() # (2-224=223)
  y <- lag0(maxlag , D.gdp) %>% scale()# dependent variable (3-224=222), done
  D.gdp.lags <- lags(maxlag, D.gdp)
  gdp.lag1 <- lag1(maxlag, gdp[-1])
  # create x 
  x.sw = data.sw[,-1]
  x.t.sw = t.sw[-1]
  xtem.sw = x.sw[-1,]
  xtem.sw[ , x.t.sw=="2"] <- apply( x.sw[ , x.t.sw=="2"], 2, diff )
  colnames(xtem.sw)[x.t.sw=="2"] = 
    paste0("D.", colnames(xtem.sw)[x.t.sw=="2"])
  xtem.sw = lags(maxlag , xtem.sw)  # explanatory variables (2-223=222), done
 
  x.sw = cbind(xtem.sw, D.gdp.lags, gdp.lag1) %>% scale()
  p.sw = dim(x.sw)[2]
  n.sw = dim(x.sw)[1]
  lambda.fix.sw = sqrt(log( p.sw ) / n.sw )
  # fit model, collect data
  lasso.sw = glmnet(x.sw, y, alpha=1, thresh=1E-5, lambda= lambda.seq, maxit = 10^9)
  fitted.sw = predict(lasso.sw , s=lambda.seq[which.min(abs(lambda.seq-lambda.fix.sw))] , 
                      newx = x.sw)
  mse.sw = mean((y-fitted.sw)^2)
  mse.out.sw = out.mse(x.sw, y, lambda.seq, p.sw)
  list = list(lasso = lasso.sw, mse = mse.sw, mse.out = mse.out.sw, lambda = lambda.fix.sw)
  return(list)
}

lasso222 <- function(data, tcode){
  # create data
  i0 = data[-(1:2) , tcode=="0"]
  di1 = apply(data[,tcode=="1"], 2, diff)
  di1 = di1[-1,]
  colnames(di1) = paste0("D.", colnames(di1))
  di2 = apply(data[,tcode=="2"], 2, diff)
  d2i2 = apply(di2, 2, diff)
  colnames(d2i2) = paste0("D2.", colnames(d2i2))
  x = cbind(i0, di1, d2i2)
  x = lags(maxlag , x) # x lost 2 obs
  gdp = bpraw[, 1]
  D.gdp = diff(gdp)
  D.gdp = D.gdp[-1]
  gdp = gdp[-(1:2)]
  y.lag1 = lag1(maxlag , gdp)
  y.5 = lag0(maxlag, D.gdp) %>% scale()
 
  x.5 = cbind( y.lag1 , x) %>% scale()
  p.5 = dim(x.5)[2]
  n.5 = dim(x.5)[1]
  lambda.fix = sqrt(log( p.5 ) / n.5 )
  # fit model, collect data
  lasso.5 = glmnet(x.5, y.5, alpha=1, thresh=1E-5, lambda= lambda.seq, maxit = 10^9)
  fitted.5 = predict(lasso.5 , s=lambda.seq[which.min(abs(lambda.seq-lambda.fix))] , newx = x.5)
  mse.5 = mean((y.5-fitted.5)^2)
  mse.out5 = out.mse(x.5, y.5, lambda.seq, p.5)
  list = list(lasso = lasso.5, mse = mse.5, mse.out = mse.out5, lambda = lambda.fix)
  return(list)
}

lasso333 <- function(data, tcode){
  # create data
  gdp = bpraw[, 1]
  D.gdp = diff(gdp)[-1]
  y = lag0(maxlag, D.gdp) %>% scale()
  i0 = data[ , tcode=="0"] # *
  i1 = data[,tcode=="1"] # *
  i2 = data[,tcode=="2"]
  di1 = apply(i1, 2, diff) # *
  di2 = apply(i2, 2, diff) # *
  d2i2 = apply(di2, 2, diff) # *
  colnames(di1) = paste0("D.", colnames(di1))
  colnames(di2) = paste0("D.", colnames(di2))
  colnames(d2i2) = paste0("D2.", colnames(d2i2))
  xtem = cbind(i0[-(1:2),] , di1[-1,] , d2i2)
  xtem = lags(maxlag, xtem)
  xtem2 = cbind( i1[-(1:2),] , di2[-1,])
  xtem2 = lag1(maxlag, xtem2)
 
  x = cbind( xtem, xtem2) %>% scale()
  p = dim(x)[2]
  n = dim(x)[1]
  lambda.fix = sqrt(log( p ) / n )
  # fit model, collect data
  lasso.5 = glmnet(x, y, alpha=1, thresh=1E-5, lambda= lambda.seq, maxit = 10^9)
  fitted.5 = predict(lasso.5 , s=lambda.seq[which.min(abs(lambda.seq-lambda.fix))] , newx = x)
  mse.5 = mean((y-fitted.5)^2)
  mse.out5 = out.mse(x, y, lambda.seq, p)
  list = list(lasso = lasso.5, mse = mse.5, mse.out = mse.out5, lambda = lambda.fix)
  return(list)
}
 
sumtable <- function(t1, t2, t3){
  variable = c(t1$coef[,2] %>% as.character(), 
               t2$coef[,2] %>% as.character(), 
               t3$coef[,2] %>% as.character()) %>% unique() %>% 
    as.data.frame(stringsAsFactors = FALSE)
  colnames(variable) = "variable"
  sum.table = left_join(variable, t1$coef[,2:3]) %>% 
    left_join(., t2$coef[,2:3], by = "variable") %>%
    left_join(., t3$coef[,2:3], by = "variable") %>%
    as.data.frame(stringsAsFactors = FALSE)
  sum.table[ , -1] = round(sum.table[ , -1] , digits = 6)
  return(sum.table)
}


lasso111ng <- function(data,  yr, tcode){
  x.sw = data
  x.t.sw = tcode
  ytem = yr[-1]
  y <- lag0(maxlag , ytem) %>% scale() %>% as.numeric()# dependent variable (3-224=222), done
  #y.lags <- lags(maxlag, ytem)
  xtem.sw = x.sw[-1,]
  xtem.sw[ , x.t.sw=="2"] <- apply( x.sw[ , x.t.sw=="2"], 2, diff )
  colnames(xtem.sw)[x.t.sw=="2"] = 
    paste0("D.", colnames(xtem.sw)[x.t.sw=="2"])
  xtem.sw = lags(maxlag , xtem.sw)  # explanatory variables (2-223=222), done
  #tt = 1:dim(xtem.sw)[1] # remove time trend in lasso model 13 Nov
  #x.sw = cbind(tt, xtem.sw) %>% scale()
  x.sw = xtem.sw %>% scale()
  p.sw = dim(x.sw)[2]
  n.sw = dim(x.sw)[1]
  lambda.fix.sw = sqrt(log( p.sw ) / n.sw )
  # fit model, collect data
  lasso.sw = glmnet(x.sw, y, alpha=1, thresh=1E-5, lambda= lambda.seq, maxit = 10^9)
  fitted.sw = predict(lasso.sw , s=lambda.seq[which.min(abs(lambda.seq-lambda.fix.sw))] , 
                      newx = x.sw)
  mse.sw = mean((y-fitted.sw)^2)
  mse.out.sw = out.mse(x.sw, y, lambda.seq, p.sw)
  list = list(lasso = lasso.sw, mse = mse.sw, mse.out = mse.out.sw, lambda = lambda.fix.sw)
  return(list)
}

lasso222ng <- function(data, yr, tcode){
  # create data
  i0 = data[-(1:2) , tcode=="0"]
  di1 = apply(data[,tcode=="1"], 2, diff)
  di1 = di1[-1,]
  colnames(di1) = paste0("D.", colnames(di1))
  di2 = apply(data[,tcode=="2"], 2, diff)
  d2i2 = apply(di2, 2, diff)
  colnames(d2i2) = paste0("D2.", colnames(d2i2))
  x = cbind(i0, di1, d2i2)
  x = lags(maxlag , x) # x lost 2 obs
  ytem = yr[-(1:2)]
  y.5 = lag0(maxlag, ytem) %>% scale()
  #y.lags = lags(maxlag, ytem)
  #tt = 1:dim(x)[1] # remove time trend
  #x.5 = cbind( tt , x) %>% scale()
  x.5 = scale(x)
  p.5 = dim(x.5)[2]
  n.5 = dim(x.5)[1]
  lambda.fix = sqrt(log( p.5 ) / n.5 )
  # fit model, collect data
  lasso.5 = glmnet(x.5, y.5, alpha=1, thresh=1E-5, lambda= lambda.seq, maxit = 10^9)
  fitted.5 = predict(lasso.5 , s=lambda.seq[which.min(abs(lambda.seq-lambda.fix))] , newx = x.5)
  mse.5 = mean((y.5-fitted.5)^2)
  mse.out5 = out.mse(x.5, y.5, lambda.seq, p.5)
  list = list(lasso = lasso.5, mse = mse.5, mse.out = mse.out5, lambda = lambda.fix)
  return(list)
}


lasso333ng <- function(data, tcode, y){
  # create data
  ytem = y[-(1:2)]
  y = lag0(maxlag, ytem) %>% scale()
  i0 = data[ , tcode=="0"] # *
  i1 = data[,tcode=="1"] # *
  i2 = data[,tcode=="2"]
  di1 = apply(i1, 2, diff) # *
  di2 = apply(i2, 2, diff) # *
  d2i2 = apply(di2, 2, diff) # *
  colnames(di1) = paste0("D.", colnames(di1))
  colnames(di2) = paste0("D.", colnames(di2))
  colnames(d2i2) = paste0("D2.", colnames(d2i2))
  xtem = cbind(i0[-(1:2),] , di1[-1,] , d2i2)
  xtem = lags(maxlag, xtem)
  xtem2 = cbind( i1[-(1:2),] , di2[-1,])
  xtem2 = lag1(maxlag, xtem2)
  #tt = 1:dim(xtem)[1] # remove time trend
  #x = cbind(tt, xtem, xtem2) %>% scale()
  x = cbind(xtem, xtem2) %>% scale()
  p = dim(x)[2]
  n = dim(x)[1]
  lambda.fix = sqrt(log( p ) / n )
  # fit model, collect data
  lasso.5 = glmnet(x, y, alpha=1, thresh=1E-5, lambda= lambda.seq, maxit = 10^9)
  fitted.5 = predict(lasso.5 , s=lambda.seq[which.min(abs(lambda.seq-lambda.fix))] , newx = x)
  mse.5 = mean((y-fitted.5)^2)
  mse.out5 = out.mse(x, y, lambda.seq, p)
  list = list(lasso = lasso.5, mse = mse.5, mse.out = mse.out5, lambda = lambda.fix)
  return(list)
}





