#calculating IACT using CODA package
library(coda)
theta<-read_excel("PMMH_SmallScale_v4_nonfixed_USED.xlsx", col_names = FALSE)
ESS <-effectiveSize(theta)
length.theta=dim(theta)[1]
IACT <- length.theta/ESS
IACT
mean(IACT)
max(IACT)
