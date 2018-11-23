library("glmnet")
library("reshape2")
library("cellranger")
library("tidyverse")
library("knitr")
library("ggrepel")
library("pander")
library("kableExtra")
library("tseries")
library("urca")
source("/Users/stanza/documents/github/lasso/fun.R")

maxlag = 4
ngraw <- read.csv("/Users/stanza/documents/github/lasso/ng/ngraw.csv")  
tcode <- read.csv("/Users/stanza/documents/github/lasso/ng/tcode.csv") 
short <- read.delim2("/Users/stanza/documents/github/lasso/ng/series.txt", header = FALSE,
                     stringsAsFactors=FALSE) %>% as.data.frame()
colnames(ngraw) = str_trim(short$V1)
ydata <- readxl::read_xls("/Users/stanza/documents/github/lasso/ng/yvar.xls")
ydata <- select(ydata, yr2, yr3, yr4, yr5)
ydata = ydata %>% mutate(time=(1:480))
dt.yr2 = lm(data = ydata, yr2 ~ time)$residuals
dt.yr3 = lm(data = ydata, yr3 ~ time)$residuals
dt.yr4 = lm(data = ydata, yr4 ~ time)$residuals
dt.yr5 = lm(data = ydata, yr5 ~ time)$residuals
fd.yr2 = diff(ydata$yr2)
fd.yr3 = diff(ydata$yr3)
fd.yr4 = diff(ydata$yr4)
fd.yr5 = diff(ydata$yr5)
yr2 = ydata$yr2
yr3 = ydata$yr3
yr4 = ydata$yr4
yr5 = ydata$yr5
ngraw = ngraw[1:dim(ydata)[1] , ]
ngraw[, tcode==4 | tcode==5 | tcode==6] = log(ngraw[, tcode==4 | tcode==5 | tcode==6])
# keep full set of data before remove "spread"
full.raw <- ngraw
full.tcode <- tcode
full.short <- short
# remove spread
idx2 = which(colnames(ngraw) == "scp90" | 
               colnames(ngraw) == "sfygm3" | 
               colnames(ngraw) == "sfygm6" | 
               colnames(ngraw) == "sfygt1" | 
               colnames(ngraw) == "sfygt5" |
               colnames(ngraw) == "sfygt10" |
               colnames(ngraw) == "sfyaaac" |
               colnames(ngraw) == "sfybaac")
tcode <- tcode[-idx2, ]
short <- short[-idx2, ]
ngraw <- ngraw[ , -idx2]

t.aic = data.frame(integration = character(length = 131), stringsAsFactors = FALSE)
t.bic = data.frame(integration = character(length = 131), stringsAsFactors = FALSE)
t.ng = data.frame(integration = character(length = 131), stringsAsFactors = FALSE)
adf.aic = adf.auto(full.raw, criteria = "AIC")
adf.bic = adf.auto(full.raw, criteria = "BIC")
t.aic[adf.aic$Conclusion=="I(0)",] = "0"
t.bic[adf.bic$Conclusion=="I(0)",] = "0"
diff.full.raw = apply(full.raw , 2, diff)
adf.aic.diff = adf.auto(diff.full.raw, criteria = "AIC")
adf.aic = cbind(adf.aic, adf.aic.diff)
adf.bic.diff = adf.auto(diff.full.raw, criteria = "BIC")
adf.bic = cbind(adf.bic, adf.bic.diff)
write.csv(adf.aic, "/Users/stanza/documents/github/lasso/ng/adf_aic_ng.csv")
write.csv(adf.bic, "/Users/stanza/documents/github/lasso/ng/adf_bic_ng.csv")
t.aic[adf.aic[,6]=="I(1)",] = "2"
t.bic[adf.bic[,6]=="I(1)",] = "2"
t.aic[t.aic=="",] = "1"
t.bic[t.bic=="",] = "1"
t.ng[full.tcode==1,] = "0"
t.ng[full.tcode==2,] = "1"
t.ng[full.tcode==4,] = "0"
t.ng[full.tcode==5,] = "1"
t.ng[full.tcode==6,] = "2"
full.t.aic = t.aic
full.t.bic = t.bic
full.t.ng = t.ng
t.aic = t.aic[-idx2,]
t.bic = t.bic[-idx2,]
t.ng = t.ng[-idx2,]

# Investigation of the all variables
inves.raw2 = full.raw
inves.tcode2 = full.tcode
ng2 = inves.raw2[,inves.tcode2==6]
#adf.auto(ng2, criteria = "AIC")
temp = apply(ng2, 2, diff)
temp = rbind(NA, temp)
inves.raw2[,inves.tcode2==6] = temp
adf.aic = adf.auto(inves.raw2, criteria = "AIC")
adf.bic = adf.auto(inves.raw2, criteria = "BIC")
ng = data.frame(Series = colnames(inves.raw2), ng = character(length = 131), 
                stringsAsFactors = FALSE)
ng$ng[inves.tcode2==1] = "I(0)"
ng$ng[inves.tcode2==2] = "I(1)"
ng$ng[inves.tcode2==4] = "I(0)"
ng$ng[inves.tcode2==5] = "I(1)"
ng$ng[inves.tcode2==6] = "I(2)"
aabbss <- left_join(ng, adf.aic, by = "Series") %>% left_join(., adf.bic, by = "Series")
aabbss <- aabbss[,-4]
aabbss <- aabbss[,c(1,2,3,5,6,4,7)]
colnames(aabbss) = c("Series", "ng", "AIC", "BIC", "Type", "AIC lags", "BIC lags")
tt = aabbss[inves.tcode2==6,]
tt$AIC[tt$AIC=="I(1)"] <- "I(2)"
tt$BIC[tt$BIC=="I(1)"] <- "I(2)"
tt$AIC[tt$AIC=="I(0)"] <- "I(1)"
tt$BIC[tt$BIC=="I(0)"] <- "I(1)"
aabbss[inves.tcode2==6,] = tt
#saveRDS(aabbss, file = "/Users/stanza/documents/github/lasso/ng/compare.rds")
maxlag = 4
lambda.seq <- seq(0.5, 0, -0.0025)
lambda.seq.fd <- seq(0.2, 0, -0.001)
save.image(file = "/Users/stanza/documents/github/lasso/ng/ngdata.RData")







