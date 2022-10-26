#####################
# Download daily settlement price of futures from Quandl
#####################

rm(list = ls())
# ------------------------------------------------------------------------------------------

# install.packages("Quandl")
library(Quandl)
# BZ – Brent Crude Oil Futures (Oil Futures used in Pricing International Oil) 
BZ = Quandl("CHRIS/CME_BZ1")

# CL – West Texas Intermediate Crude Oil Futures (Oil Futures used in Pricing North American Oil)
CL = Quandl("CHRIS/CME_CL1")

# ES - E-Mini S&P 500 Futures
ES = Quandl("CHRIS/CME_ES1")

# GC – Gold Futures
GC = Quandl("CHRIS/CME_GC1")

# GE – EuroDollar Futures (providing information about market belief of future interest rates)
GE = Quandl("CHRIS/CME_ED1")

# GF – Feeder Cattle Futures
GF = Quandl("CHRIS/CME_FC1")

# HE – Lean Hog Futures
HE = Quandl("CHRIS/CME_LN1")

# HG – Copper Futures
HG = Quandl("CHRIS/CME_HG1")

# HO – Heating Oil Futures
HO = Quandl("CHRIS/ICE_O1")

# LE – Live Cattle Futures
LE = Quandl("CHRIS/CME_LC1")

# NG – Natural Gas Futures
NG = Quandl("CHRIS/CME_QG1")

# NQ – E-Mini Nasdaq-100 Futures
NQ = Quandl("CHRIS/CME_NQ1")

# PA – Palladium Futures
PA = Quandl("CHRIS/CME_PA1")

# PL – Platinum Futures
PL = Quandl("CHRIS/CME_PL1")

# RB – RBOB Gasoline Futures
RB = Quandl("CHRIS/CME_RB1")

# RTY – E-Mini Russell 2000 Futures
RTY = Quandl("CHRIS/ICE_TF1")

# SI – Silver Futures
SI = Quandl("CHRIS/CME_SI1")

# YM – E-Mini Dow Jones Industrial Average Futures
YM = Quandl("CHRIS/OSE_DJIA1")

# ZC – Corn Futures
ZC = Quandl("CHRIS/CME_C1")

# ZF – Treasury Futures on 5-Year Maturity T-Note
ZF = Quandl("CHRIS/CME_FV1")

# ZL – Soybean Oil Futures
ZL = Quandl("CHRIS/CME_BO1")

# ZM – Soybean Meal Futures
ZM = Quandl("CHRIS/CME_SM1")

# ZN – Treasury Futures on 10-Year Maturity Bond
ZN = Quandl("CHRIS/CME_TY1")

# ZS – Soybean Futures
ZS = Quandl("CHRIS/CME_S1")

# ZT – Treasury Futures on 2-Year Maturity T-Note
ZT = Quandl("CHRIS/CME_TU1")

cocoa = Quandl("CHRIS/LIFFE_C1")

wheat = Quandl("CHRIS/CME_W1")

oats = Quandl("CHRIS/CME_O1")

ethanol = Quandl("CHRIS/CME_EH1")


#######################
var_names=c('BZ','CL','ES','GC','GE','GF','HE','HG','HO','LE','NG','NQ',
            'PA','PL','RB','RTY','SI','YM','ZC','ZF','ZL','ZM','ZN',
            'ZS','ZT','cocoa','wheat','oats','ethanol')
p=length(var_names);

## First look at the colomn names
col_ind = c();
for (i in 1:p){
  temp=get(var_names[i]);
  col_ind[i] = which(colnames(temp) == "Settle" |
                       colnames(temp) == "Settlement Price")
}
for (i in 1:p){
  temp=get(var_names[i]);
  print(i)
  print(colnames(temp)[col_ind[i]])
  print(min(temp[,1]))
  print(dim(temp))
}

## Outer join futures price by date
df=get(var_names[1])[,c(1, col_ind[1])];
colnames(df)[2]=var_names[1];
for (i in 2:p){
  temp=get(var_names[i])[,c(1, col_ind[i])];
  colnames(temp)[2]=var_names[i];
  df=merge(df,temp,by='Date',all=T);
}

write.csv(df,'futures_price_raw.csv',row.names = F)





