####### GDP data
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
lambda.seq <- seq(0.5, 0, -0.0025)
bpraw <- read.csv("/Users/stanza/documents/github/lasso/new/bpraw.csv") 
bptrans <- read.csv("/Users/stanza/documents/github/lasso/new/bptrans.csv", na.strings = NaN)
tcode <- read.csv("/Users/stanza/documents/github/lasso/new/tcode.csv", stringsAsFactors = F) 
short <- read.csv("/Users/stanza/documents/github/lasso/new/short.csv", stringsAsFactors = F) %>% unname()
long <- read.csv("/Users/stanza/documents/github/lasso/new/long.csv", stringsAsFactors = F) %>% unname()
# remove NA
idx <- apply(bpraw, 2, anyNA) %>% which() %>% unname()
colnames(bpraw) = short
colnames(bptrans) = short
bpraw <- bpraw[ , -idx]
bptrans <- bptrans[, -idx] %>% na.omit()
tcode <- tcode[-idx,] 
short <- short[, -idx] %>% t() %>% as.vector()
long <- long[, -idx] %>% t() %>% as.vector()
# log()
bpraw[, tcode==5 | tcode==6] = log(bpraw[, tcode==5 | tcode==6])
# keep full set of data before remove "spread"
full.raw <- bpraw
full.tran <- bptrans
full.tcode <- tcode
full.short <- short
full.long <- long
# remove "BAA_GS10", "tb6m_tb3m, GS1_tb3m, GS10_tb3m, CP_Tbill Spread"
idx2 = which(colnames(bpraw) == "BAA_GS10" | 
               colnames(bpraw) == "tb6m_tb3m" | 
               colnames(bpraw) == "GS1_tb3m" | 
               colnames(bpraw) == "GS10_tb3m" | 
               colnames(bpraw) == "CP_Tbill Spread")
tcode <- tcode[-idx2]
short <- short[-idx2]
long <- long[-idx2]
bpraw <- bpraw[ , -idx2]
bptrans <- bptrans[ , -idx2]
t.aic = data.frame(integration = character(length = 146), stringsAsFactors = FALSE)
t.bic = data.frame(integration = character(length = 146), stringsAsFactors = FALSE)
t.sw = data.frame(integration = character(length = 146), stringsAsFactors = FALSE)
adf.aic = adf.auto(full.raw, criteria = "AIC")
adf.bic = adf.auto(full.raw, criteria = "BIC")
t.aic[adf.aic$Conclusion=="I(0)",] = "0"
t.bic[adf.bic$Conclusion=="I(0)",] = "0"
diff.full.raw = apply(full.raw , 2, diff)
adf.aic.diff = adf.auto(diff.full.raw, criteria = "AIC")
adf.aic = cbind(adf.aic, adf.aic.diff)
adf.bic.diff = adf.auto(diff.full.raw, criteria = "BIC")
adf.bic = cbind(adf.bic, adf.bic.diff)
#write.csv(adf.aic, "/Users/stanza/documents/github/lasso/adf_aic.csv")
#write.csv(adf.bic, "/Users/stanza/documents/github/lasso/adf_bic.csv")
t.aic[adf.aic[,6]=="I(1)",] = "2"
t.bic[adf.bic[,6]=="I(1)",] = "2"
t.aic[t.aic=="",] = "1"
t.bic[t.bic=="",] = "1"
t.sw[full.tcode==1,] = "0"
t.sw[full.tcode==2 | full.tcode==5,] = "1"
t.sw[full.tcode==6,] = "2"
full.t.aic = t.aic
full.t.bic = t.bic
full.t.sw = t.sw
t.aic = t.aic[-idx2,]
t.bic = t.bic[-idx2,]
t.sw = t.sw[-idx2,]
compare = cbind(full.short, full.t.sw, full.t.aic, full.t.bic)
save.image(file = "/Users/stanza/documents/github/lasso/sab_gdp.RData")