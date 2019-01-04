


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
  else if(dim(x)[2]==1) {
    for (i in 1:maxlag) {
      x.lag = x[((maxlag + 1 - i) : (length(x) - i))]
      xtem[[i]] = x.lag
      nn[i] = paste0( deparse(substitute(x)) , ".lag", i)
    }
    names(xtem) = nn
    lagged.x = do.call(cbind, xtem)
  } else {
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
  else if (dim(y)[2]==1) {
    yt = y[((maxlag + 1) : length(y))]
    names(yt) = paste0( deparse(substitute(y)) )
  } else {
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

adf.drift <- function(data, criteria) {
  tb = data.frame(Series = character(), Conclusion = character(), Type = character() , 
                  Lags = numeric() , 
                  stringsAsFactors = FALSE)
  if (is.null(dim(data))) { # if it is a single series
    nm = deparse(substitute(data))
    data = na.omit(data)
    t = length(data)
    #x = 1:t
    pmax = ceiling(12 * (t/100)^0.25 )
    type = "drift"
    adf = ur.df(data, selectlags = criteria, lags=pmax, type=type)
    calc = adf@teststat[1,1]
    if (calc < adf@cval[1,2]) {
      con = "I(0)"
    } else {
      con = "I(1)"
    }
    ll = dim(adf@testreg$coefficients)[1]-2
    tb[1,] = c(nm, con, type, ll)
  } # if it is a data frame
  else {
    for (i in 1:dim(data)[2]){
      s = data[,i]
      s = na.omit(s)
      t = length(s)
      #x = 1:t
      pmax = ceiling(12 * (t/100)^0.25 )
      type = "drift"
      adf = ur.df(s, selectlags = criteria, lags=pmax, type=type)
      calc = adf@teststat[1,1]
      if (calc < adf@cval[1,2]) {
        con = "I(0)"
      } else {
        con = "I(1)"
      }
      ll = dim(adf@testreg$coefficients)[1]-2
      tb[i,] = c(colnames(data)[i], con, type, ll)
    }
  }
  return(tb)
}
 
sumtable <- function(t1, t2, t3){# t1, t2, t3 are the outputs from trace function
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



ng1train <- function(data,  yr, tcode){
  yr = yr[-1] %>% scale()
  y <- lag0(maxlag , yr) %>% as.numeric()# dependent variable (3-224=222), done
  y.lags <- lags(maxlag, yr)
  xtem = data[-1,]
  xtem[ , tcode=="2"] <- apply( data[ , tcode=="2"], 2, diff )
  colnames(xtem)[tcode=="2"] = 
    paste0("D.", colnames(xtem)[tcode=="2"])
  xtem = lags(maxlag , xtem) %>% scale() 
  x = cbind(xtem, y.lags)  # explanatory variables (2-223=222), done
  # train
  l = ceiling(dim(x)[1]*0.8)
  x.train = x[1:l,] 
  x.test = x[(l+1):dim(x)[1] , ]
  y.train = y[1:l]
  y.test = y[(l+1):length(y)]
  p.train = dim(x.train)[2] 
  n.train = dim(x.train)[1] 
  lambda.train = sqrt(log( p.train ) / n.train ) 
  # fit model, collect data
  lasso = glmnet(x.train, y.train, alpha = 1, thresh = 1E-7, lambda = lambda.seq, maxit = 10^9)
  fitted = predict(lasso , s=lambda.seq[which.min(abs(lambda.seq-lambda.train))] , newx = x.test) #10
  out.mse = mean((y.test-fitted)^2)
  list = list(lasso=lasso, mse=out.mse, lambda = lambda.train)
  return(list)
}

ng2train <- function(data, yr, tcode){
  # create data
  i0 = data[-(1:2) , tcode=="0"]
  di1 = apply(data[,tcode=="1"], 2, diff)
  di1 = di1[-1,]
  colnames(di1) = paste0("D.", colnames(di1))
  di2 = apply(data[,tcode=="2"], 2, diff)
  d2i2 = apply(di2, 2, diff)
  colnames(d2i2) = paste0("D2.", colnames(d2i2))
  xtem = cbind(i0, di1, d2i2) %>% scale()
  yr = yr[-(1:2)] %>% scale()
  y.lags <- lags(maxlag, yr)
  y = lag0(maxlag, yr) %>% as.numeric()
  xtem = lags(maxlag , xtem)  # x lost 2 obs
  x = cbind(y.lags, xtem)
  # train
  l = ceiling(dim(x)[1]*0.8)
  x.train = x[1:l,] 
  x.test = x[(l+1):dim(x)[1] , ]
  y.train = y[1:l]
  y.test = y[(l+1):length(y)]
  p.train = dim(x.train)[2]
  n.train = dim(x.train)[1]
  lambda.train = sqrt(log( p.train ) / n.train )
  # fit model, collect data
  lasso = glmnet(x.train, y.train, alpha=1, thresh=1E-7, lambda= lambda.seq, maxit = 10^9)
  fitted = predict(lasso , s=lambda.seq[which.min(abs(lambda.seq-lambda.train))] , newx = x.test)
  out.mse = mean((y.test-fitted)^2)
  list = list(lasso=lasso, mse=out.mse, lambda = lambda.train)
  return(list)
}


ng3train <- function(data, tcode, yr){
  # create data
  yr = yr[-(1:2)] %>% scale() %>% as.matrix()
  y = lag0(maxlag, yr) 
  y.lags = lags(maxlag, yr)
  i0 = data[ , tcode=="0"] # *
  i1 = data[,tcode=="1"] # *
  i2 = data[,tcode=="2"]
  di1 = apply(i1, 2, diff) # *
  di2 = apply(i2, 2, diff) # *
  d2i2 = apply(di2, 2, diff) # *
  colnames(di1) = paste0("D.", colnames(di1))
  colnames(di2) = paste0("D.", colnames(di2))
  colnames(d2i2) = paste0("D2.", colnames(d2i2))
  xtem = cbind(i0[-(1:2),] , di1[-1,] , d2i2) %>% scale()
  xtem = lags(maxlag, xtem)
  xtem2 = cbind( i1[-(1:2),] , di2[-1,]) %>% scale()
  xtem2 = lag1(maxlag, xtem2)
  x = cbind(xtem, xtem2, y.lags)
  # train
  l = ceiling(dim(x)[1]*0.8)
  x.train = x[1:l,] 
  x.test = x[(l+1):dim(x)[1] , ]
  y.train = y[1:l]
  y.test = y[(l+1):length(y)]
  p.train = dim(x.train)[2]
  n.train = dim(x.train)[1]
  lambda.train = sqrt(log( p.train ) / n.train )
  # fit model, collect data
  lasso = glmnet(x.train, y.train, alpha=1, thresh=1E-7, lambda= lambda.seq, maxit = 10^9)
  fitted = predict(lasso , s=lambda.seq[which.min(abs(lambda.seq-lambda.train))] , newx = x.test)
  out.mse = mean((y.test-fitted)^2)
  list = list(lasso=lasso, mse=out.mse, lambda = lambda.train)
  return(list)
}

############### only sig lags #############################################################
#x = x[, which( substr(colnames(x), nchar(colnames(x))-3, nchar(colnames(x))) == "lag1" |
#                 substr(colnames(x), nchar(colnames(x))-3, nchar(colnames(x))) =="lag2" |
#                 substr(colnames(x), nchar(colnames(x))-3, nchar(colnames(x))) =="lag3" |
#                 substr(colnames(x), nchar(colnames(x))-4, nchar(colnames(x))) =="lag14" |
#                 substr(colnames(x), nchar(colnames(x))-4, nchar(colnames(x))) =="lag12" |
#                 substr(colnames(x), nchar(colnames(x))-4, nchar(colnames(x))) =="lag13" |
#                 substr(colnames(x), nchar(colnames(x))-4, nchar(colnames(x))) =="lag24" |
#                 substr(colnames(x), nchar(colnames(x))-4, nchar(colnames(x))) =="lag25")]
############### only sig lags #############################################################

#### GDP #### out-of-sample MSE ######## rolling window estimation #######
GDP1train <- function(data, tcode){
  GDP = data[ , colnames(data)=="GDP"]
  D.GDP = diff(GDP) %>% as.numeric()
  GDP.lag1 = lag1(maxlag, GDP[-1])
  y = lag0(maxlag, D.GDP) %>% scale() %>% as.matrix()
  temp = data[-1, ]
  temp[, tcode=="2"] <- apply( data[ , tcode=="2"], 2, diff )
  colnames(temp)[tcode=="2"] = 
    paste0("D.", colnames(temp)[tcode=="2"])
  temp = lags(maxlag, temp)
  x = cbind(temp, GDP.lag1) %>% scale() %>% as.matrix()
  # train
  l = ceiling(dim(x)[1]*0.8)
  x.train = x[1:l,] 
  x.test = x[(l+1):dim(x)[1] , ]
  y.train = y[1:l]
  y.test = y[(l+1):length(y)]
  p.train = dim(x.train)[2]
  n.train = dim(x.train)[1]
  lambda.train = sqrt(log( p.train ) / n.train )
  # fit model, collect data
  lasso.first = glmnet(x.train, y.train, alpha=1, thresh=1E-5, lambda= lambda.seq, maxit = 10^9)
  # exp 23rd Nov
  tau=1
  first.step.coef=coef(lasso.first)[-1]
  penalty.factor=abs(first.step.coef+1/sqrt(nrow(x.train)))^(-tau)
  lasso = glmnet(x.train, y.train, alpha=1, thresh=1E-5, lambda= lambda.seq, 
                       maxit = 10^9, penalty.factor=penalty.factor)
  ## exp 23rd Nov
  fitted = predict(lasso , s=lambda.seq[which.min(abs(lambda.seq-lambda.train))] , newx = x.test)
  out.mse = mean((y.test-fitted)^2)
  list = list(lasso=lasso, mse=out.mse, lambda = lambda.train)
}

GDP2train <- function(data, tcode){
  i0 = data[-(1:2) , tcode=="0"]
  di1 = apply(data[,tcode=="1"], 2, diff)
  di1 = di1[-1,]
  colnames(di1) = paste0("D.", colnames(di1))
  di2 = apply(data[,tcode=="2"], 2, diff)
  d2i2 = apply(di2, 2, diff)
  colnames(d2i2) = paste0("D2.", colnames(d2i2))
  x = cbind(i0, di1, d2i2)
  x = lags(maxlag , x) %>% scale() %>% as.matrix()# x lost 2 obs
  GDP = data[, colnames(data)=="GDP"]
  D.GDP = diff(GDP)
  D.GDP = D.GDP[-1]
  #GDP = GDP[-(1:2)]
  #GDP.lag1 = lag1(maxlag , GDP)
  y = lag0(maxlag, D.GDP) %>% scale() %>% as.matrix()
  # train
  l = ceiling(dim(x)[1]*0.8)
  x.train = x[1:l,] 
  x.test = x[(l+1):dim(x)[1] , ]
  y.train = y[1:l]
  y.test = y[(l+1):length(y)]
  p.train = dim(x.train)[2]
  n.train = dim(x.train)[1]
  lambda.train = sqrt(log( p.train ) / n.train )
  # fit
  lasso = glmnet(x.train, y.train, alpha=1, thresh=1E-5, lambda= lambda.seq, maxit = 10^9)
  fitted = predict(lasso , s=lambda.seq[which.min(abs(lambda.seq-lambda.train))] , newx = x.test)
  out.mse = mean((y.test-fitted)^2)
  list = list(lasso=lasso, mse=out.mse, lambda = lambda.train)
}

GDP3train <- function(data, tcode){
  # create data
  GDP = data[, colnames(data)=="GDP"]
  D.GDP = diff(GDP)[-1]
  y = lag0(maxlag, D.GDP) %>% scale()
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
  # train
  l = ceiling(dim(x)[1]*0.8)
  x.train = x[1:l,] 
  x.test = x[(l+1):dim(x)[1] , ]
  y.train = y[1:l]
  y.test = y[(l+1):length(y)]
  p.train = dim(x.train)[2]
  n.train = dim(x.train)[1]
  lambda.train = sqrt(log( p.train ) / n.train )
  # fit
  lasso = glmnet(x.train, y.train, alpha=1, thresh=1E-5, lambda= lambda.seq, maxit = 10^9)
  fitted = predict(lasso , s=lambda.seq[which.min(abs(lambda.seq-lambda.train))] , newx = x.test)
  out.mse = mean((y.test-fitted)^2)
  list = list(lasso=lasso, mse=out.mse, lambda = lambda.train)
}

#### Inflation #### out-of-sample MSE ######## rolling window estimation #######
inf1train <- function(data, tcode){
  CPI = data[ , colnames(data)=="CPI"]
  Inflation = diff(CPI) %>% as.numeric()
  CPI.lag1 = lag1(maxlag, CPI[-1])
  y = lag0(maxlag, Inflation) %>% scale() %>% as.matrix()
  temp = data[-1, ]
  temp[, tcode=="2"] <- apply( data[ , tcode=="2"], 2, diff )
  colnames(temp)[tcode=="2"] = 
    paste0("D.", colnames(temp)[tcode=="2"])
  temp = lags(maxlag, temp)
  x = cbind(temp, CPI.lag1) %>% scale() %>% as.matrix()
  # train
  l = ceiling(dim(x)[1]*0.8)
  x.train = x[1:l,] 
  x.test = x[(l+1):dim(x)[1] , ]
  y.train = y[1:l]
  y.test = y[(l+1):length(y)]
  p.train = dim(x.train)[2]
  n.train = dim(x.train)[1]
  lambda.train = sqrt(log( p.train ) / n.train )
  # fit model, collect data
  lasso = glmnet(x.train, y.train, alpha=1, thresh=1E-5, lambda= lambda.seq, maxit = 10^9)
  fitted = predict(lasso , s=lambda.seq[which.min(abs(lambda.seq-lambda.train))] , newx = x.test)
  out.mse = mean((y.test-fitted)^2)
  list = list(lasso=lasso, mse=out.mse, lambda = lambda.train)
}

inf2train <- function(data, tcode){
  i0 = data[-(1:2) , tcode=="0"]
  di1 = apply(data[,tcode=="1"], 2, diff)
  di1 = di1[-1,]
  colnames(di1) = paste0("D.", colnames(di1))
  colnames(di1)[colnames(di1)=="D.CPI"] = "Inflation"
  di2 = apply(data[,tcode=="2"], 2, diff)
  d2i2 = apply(di2, 2, diff)
  colnames(d2i2) = paste0("D2.", colnames(d2i2))
  x = cbind(i0, di1, d2i2)
  x = lags(maxlag , x) %>% scale() %>% as.matrix()# x lost 2 obs
  CPI = data[, colnames(data)=="CPI"]
  Inflation = diff(CPI)
  Inflation = Inflation[-1]
  #CPI = CPI[-(1:2)]
  #CPI.lag1 = lag1(maxlag , CPI)
  y = lag0(maxlag, Inflation) %>% scale() %>% as.matrix()
  # train
  l = ceiling(dim(x)[1]*0.8)
  x.train = x[1:l,] 
  x.test = x[(l+1):dim(x)[1] , ]
  y.train = y[1:l]
  y.test = y[(l+1):length(y)]
  p.train = dim(x.train)[2]
  n.train = dim(x.train)[1]
  lambda.train = sqrt(log( p.train ) / n.train )
  # fit
  lasso = glmnet(x.train, y.train, alpha=1, thresh=1E-5, lambda= lambda.seq, maxit = 10^9)
  fitted = predict(lasso , s=lambda.seq[which.min(abs(lambda.seq-lambda.train))] , newx = x.test)
  out.mse = mean((y.test-fitted)^2)
  list = list(lasso=lasso, mse=out.mse, lambda = lambda.train)
}

inf3train <- function(data, tcode){
  # create data
  CPI = data[, colnames(data)=="CPI"]
  Inflation = diff(CPI)[-1]
  y = lag0(maxlag, Inflation) %>% scale()
  i0 = data[ , tcode=="0"] # *
  i1 = data[,tcode=="1"] # *
  i2 = data[,tcode=="2"]
  di1 = apply(i1, 2, diff) # *
  di2 = apply(i2, 2, diff) # *
  d2i2 = apply(di2, 2, diff) # *
  colnames(di1) = paste0("D.", colnames(di1))
  colnames(di1)[colnames(di1)=="D.CPI"] = "Inflation"
  colnames(di2) = paste0("D.", colnames(di2))
  colnames(d2i2) = paste0("D2.", colnames(d2i2))
  xtem = cbind(i0[-(1:2),] , di1[-1,] , d2i2)
  xtem = lags(maxlag, xtem)
  xtem2 = cbind( i1[-(1:2),] , di2[-1,])
  xtem2 = lag1(maxlag, xtem2)
  x = cbind( xtem, xtem2) %>% scale()
  # train
  l = ceiling(dim(x)[1]*0.8)
  x.train = x[1:l,] 
  x.test = x[(l+1):dim(x)[1] , ]
  y.train = y[1:l]
  y.test = y[(l+1):length(y)]
  p.train = dim(x.train)[2]
  n.train = dim(x.train)[1]
  lambda.train = sqrt(log( p.train ) / n.train )
  # fit
  lasso = glmnet(x.train, y.train, alpha=1, thresh=1E-5, lambda= lambda.seq, maxit = 10^9)
  fitted = predict(lasso , s=lambda.seq[which.min(abs(lambda.seq-lambda.train))] , newx = x.test)
  out.mse = mean((y.test-fitted)^2)
  list = list(lasso=lasso, mse=out.mse, lambda = lambda.train)
}



inf1.2train <- function(data, tcode){
  CPI <- data[ , colnames(data)=="CPI"]
  Inflation <- diff(CPI) %>% as.numeric() # (2-224=223)
  D.Inflation <- diff(Inflation) %>% as.numeric() #(3-224=222)
  y <- lag0(maxlag , D.Inflation) %>% scale()
  D.Inflation.lags <- lags(maxlag, D.Inflation)
  Inflation.lag1 <- lag1(maxlag, Inflation[-1])
  # create x 
  x = data[ , colnames(data)!="CPI"]
  x.tcode = tcode[colnames(data)!="CPI"]
  xtem = x[-1,]
  xtem[ , x.tcode=="2"] <- apply( x[ , x.tcode=="2"], 2, diff )
  colnames(xtem)[x.tcode=="2"] = 
    paste0("D.", colnames(xtem)[x.tcode=="2"])
  xtem = lags(maxlag , xtem[-1,])  # explanatory variables (2-223=222), done
  x = cbind(xtem, D.Inflation.lags, Inflation.lag1) %>% scale()
  # train
  l = ceiling(dim(x)[1]*0.8)
  x.train = x[1:l,] 
  x.test = x[(l+1):dim(x)[1] , ]
  y.train = y[1:l]
  y.test = y[(l+1):length(y)]
  p.train = dim(x.train)[2]
  n.train = dim(x.train)[1]
  lambda.train = sqrt(log( p.train ) / n.train )
  # fit model, collect data
  lasso = glmnet(x.train, y.train, alpha=1, thresh=1E-5, lambda= lambda.seq, maxit = 10^9)
  fitted = predict(lasso , s=lambda.seq[which.min(abs(lambda.seq-lambda.train))] , newx = x.test)
  out.mse = mean((y.test-fitted)^2)
  list = list(lasso=lasso, mse=out.mse, lambda = lambda.train)
}

inf2.2train <- function(data, tcode){
  CPI = data[, colnames(data)=="CPI"]
  Inflation = diff(CPI)
  D.Inflation = diff(Inflation) #(3-224=222)
  y = lag0(maxlag, D.Inflation) %>% scale()
  i0 = data[-(1:2) , tcode=="0"]
  di1 = apply(data[,tcode=="1"], 2, diff)
  di1 = di1[-1,]
  colnames(di1) = paste0("D.", colnames(di1))
  di2 = apply(data[,tcode=="2"], 2, diff)
  d2i2 = apply(di2, 2, diff)
  colnames(d2i2) = paste0("D2.", colnames(d2i2))
  colnames(d2i2)[colnames(d2i2)=="D2.CPI"] = "D.Inflation"
  x = cbind(i0, di1, d2i2)
  x = lags(maxlag , x) # x lost 2 obs
  Inflation.lag1 = lag1(maxlag , Inflation[-1])
  x = cbind( Inflation.lag1 , x) %>% scale()
  # train
  l = ceiling(dim(x)[1]*0.8)
  x.train = x[1:l,] 
  x.test = x[(l+1):dim(x)[1] , ]
  y.train = y[1:l]
  y.test = y[(l+1):length(y)]
  p.train = dim(x.train)[2]
  n.train = dim(x.train)[1]
  lambda.train = sqrt(log( p.train ) / n.train )
  # fit model, collect data
  lasso = glmnet(x.train, y.train, alpha=1, thresh=1E-5, lambda= lambda.seq, maxit = 10^9)
  fitted = predict(lasso , s=lambda.seq[which.min(abs(lambda.seq-lambda.train))] , newx = x.test)
  out.mse = mean((y.test-fitted)^2)
  list = list(lasso=lasso, mse=out.mse, lambda = lambda.train)
}

inf3.2train <- function(data, tcode){
  CPI = data[, colnames(data)=="CPI"]
  Inflation = diff(CPI)
  D.Inflation = diff(Inflation) # (3-224=222)
  y = lag0(maxlag, D.Inflation) %>% scale()
  i0 = data[ , tcode=="0"] # *
  i1 = data[,tcode=="1"] # *
  i2 = data[,tcode=="2"]
  di1 = apply(i1, 2, diff) # *
  di2 = apply(i2, 2, diff) # *
  d2i2 = apply(di2, 2, diff) # *
  colnames(di1) = paste0("D.", colnames(di1))
  colnames(di2) = paste0("D.", colnames(di2))
  colnames(di2)[colnames(di2)=="D.CPI"] = "Inflation"
  colnames(d2i2) = paste0("D2.", colnames(d2i2))
  colnames(d2i2)[colnames(d2i2)=="D2.CPI"] = "D.Inflation"
  xtem = cbind(i0[-(1:2),] , di1[-1,] , d2i2)
  xtem = lags(maxlag, xtem)
  xtem2 = cbind( i1[-(1:2),] , di2[-1,])
  xtem2 = lag1(maxlag, xtem2)
  x = cbind( xtem, xtem2) %>% scale()
  # train
  l = ceiling(dim(x)[1]*0.8)
  x.train = x[1:l,] 
  x.test = x[(l+1):dim(x)[1] , ]
  y.train = y[1:l]
  y.test = y[(l+1):length(y)]
  p.train = dim(x.train)[2]
  n.train = dim(x.train)[1]
  lambda.train = sqrt(log( p.train ) / n.train )
  # fit model, collect data
  lasso = glmnet(x.train, y.train, alpha=1, thresh=1E-5, lambda= lambda.seq, maxit = 10^9)
  fitted = predict(lasso , s=lambda.seq[which.min(abs(lambda.seq-lambda.train))] , newx = x.test)
  out.mse = mean((y.test-fitted)^2)
  list = list(lasso=lasso, mse=out.mse, lambda = lambda.train)
}



trace.interest <- function(model, lambda.seq) {
  plot.lasso <- rbind(lambda.seq, coef(model)[-1,]) %>% t() %>%
    as.matrix() %>% as.data.frame(stringsAsFactors = F)
  colnames(plot.lasso)[1] <- "lambda"
  plot.lasso <- plot.lasso[, substr(colnames(plot.lasso), 1,
                                    nchar(colnames(plot.lasso))-5) == "TB-3Mth" |
                             substr(colnames(plot.lasso), 1,
                                    nchar(colnames(plot.lasso))-5) == "Com Paper" | 
                             colnames(plot.lasso)=="lambda"]
  melt.lasso <- melt(plot.lasso, id = "lambda")
  plot <- ggplot(melt.lasso, aes(x= lambda, y= value, color= variable))+
    geom_line() + theme(legend.position = "none") +
    labs(x='lambda', y='coefficients') + 
    ggtitle("Trace plot of two interest rates") + 
    geom_text(data=subset(melt.lasso, lambda==0),
                    aes(label= variable, color=variable)) +
    xlim(0, 0.02)
  return(plot)
}

ddriven.l = function(x, y){
  l = ceiling(nrow(x)*0.8)
  x.train = x[1:l,] 
  x.test = x[(l+1):nrow(x) , ]
  y.train = y[1:l]
  y.test = y[(l+1):length(y)]
  lasso = glmnet(x.train, y.train, alpha = 1, thresh = 1E-7, lambda = lambda.seq, maxit = 10^9)
  fitted.in = list()
  fitted.out = list()
  mse.in = vector(length = length(lambda.seq))
  mse.out = vector(length = length(lambda.seq))
  for (i in 1:length(lambda.seq)) {
    fitted.in[[i]] = predict(lasso, s=lambda.seq[i] , newx = x.train)
    fitted.out[[i]] = predict(lasso, s=lambda.seq[i] , newx = x.test)
    mse.in[i] = mean((y.train-fitted.in[[i]])^2)
    mse.out[i] = mean((y.test-fitted.out[[i]])^2)
  }
  idx.in = which.min(mse.in)
  idx.out = which.min(mse.out)
  lambda.in = lambda.seq[idx.in]
  lambda.out = lambda.seq[idx.out]
  mse.plot = cbind(lambda.seq[-201], mse.in[-201], mse.out[-201]) %>% as.data.frame()
  colnames(mse.plot) = c("lambda.seq", "mse.in", "mse.out")
  mse.plot = melt(mse.plot, id = "lambda.seq")
  plot = ggplot(mse.plot, aes(x=lambda.seq, y=value)) + 
    geom_point(aes(color=variable), size = 1) +
    geom_vline(xintercept = lambda.out)
  list = list(plot = plot, lambda.in = lambda.in, lambda.out = lambda.out, 
              mse.in=mse.in[idx.in], mse.out=mse.out[idx.out], lasso = lasso)
}

ddriven.l1l2 = function(x, y){
  l1 = ceiling(nrow(x)*0.6)
  l2 = ceiling(nrow(x)*0.8)
  x.train = x[1:l1,] 
  y.train = y[1:l1]
  x.valid = x[(l1+1):l2 , ]
  y.valid = y[(l1+1):l2]
  x.test = x[(l2+1):nrow(x) , ]
  y.test = y[(l2+1):length(y)]
  lasso = glmnet(x.train, y.train, alpha = 1, thresh = 1E-7, lambda = lambda.seq, maxit = 10^9)
  fitted.train = list()
  fitted.valid = list()
  fitted.test = list()
  mse.train = vector(length = length(lambda.seq))
  mse.valid = vector(length = length(lambda.seq))
  mse.test = vector(length = length(lambda.seq))
  for (i in 1:length(lambda.seq)) {
    fitted.train[[i]] = predict(lasso, s=lambda.seq[i] , newx = x.train)
    fitted.valid[[i]] = predict(lasso, s=lambda.seq[i] , newx = x.valid)
    fitted.test[[i]] = predict(lasso, s=lambda.seq[i] , newx = x.test)
    mse.train[i] = mean((y.train-fitted.train[[i]])^2)
    mse.valid[i] = mean((y.valid-fitted.valid[[i]])^2)
    mse.test[i] = mean((y.test-fitted.test[[i]])^2)
  }
  idx.valid = which.min(mse.valid)
  lambda.valid = lambda.seq[idx.valid]
  mse.plot = cbind(lambda.seq[-201], mse.train[-201], mse.valid[-201], mse.test[-201]) %>% as.data.frame()
  colnames(mse.plot) = c("lambda.seq", "mse.train", "mse.valid", "mse.test")
  mse.plot = melt(mse.plot, id = "lambda.seq")
  plot = ggplot(mse.plot, aes(x=lambda.seq, y=value)) + 
    geom_point(aes(color=variable), size = 1) +
    geom_vline(xintercept = lambda.valid)
  list = list(plot = plot, lambda.valid = lambda.valid, 
              mse.out=mse.test[idx.valid], lasso = lasso)
}




