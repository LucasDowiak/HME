setwd("~/Dropbox/HME/")
library(data.table)

# XUMAGBD - monthly avg; XUMLGBD - end of month?
dtfFX <- fread("data/monthly_fx.csv", skip=4)
dtfFX[, DATE := as.Date(DATE, format="%d %b %Y")]
lubridate::day(dtfFX$DATE) <- 1
dtfFX[, DATE := as.IDate(DATE)]


dtfIndProd <- fread("oecd_ind_prod.csv")
dtfIndProd[, DATE := as.IDate(paste0(TIME, "-01"))]
dtfIndProd <- dcast(dtfIndProd[LOCATION %in% c("USA", "GBR")], "DATE ~ LOCATION",
                    value.var="Value")
setnames(dtfIndProd, c("GBR", "USA"), c("IndGBR", "IndUSA"))

dtfCoreCPI <- fread("data/core_cpi.csv")
dtfCoreCPI[, DATE := as.IDate(paste0(TIME, "-01"))]
dtfCoreCPI <- dcast(dtfCoreCPI[LOCATION %in% c("USA", "GBR")], "DATE ~ LOCATION",
                    value.var="Value")
setnames(dtfCoreCPI, c("GBR", "USA"), c("CpiGBR", "CpiUSA"))


dtfLTBond <- fread("data/tenyrbond.csv")
dtfLTBond[, DATE := as.IDate(paste0(TIME, "-01"))]
dtfLTBond <- dcast(dtfLTBond[LOCATION %in% c("USA", "GBR")], "DATE ~ LOCATION",
                   value.var="Value")
setnames(dtfLTBond, c("GBR", "USA"), c("IntGBR", "IntUSA"))

dtf <- merge(dtfFX, dtfIndProd, all.x=TRUE, by="DATE")
dtf <- merge(dtf, dtfCoreCPI, all.x=TRUE, by="DATE")
dtf <- merge(dtf, dtfLTBond, all.x=TRUE, by="DATE")
dtf <- na.omit(dtf)



detrend <- function(x)
{
  N <- length(x)
  x - c(x[1], x[1] + (x[N] - x[1]) / (N - 1) * 1:(N-1))
}
par(mfrow=c(4,1))
plot(dtf[,.(DATE, 1/XUMAGBD)], type="l", main="Avg. FX ($/L)")
grid()
plot(dtf[,.(DATE, detrend(IndUSA))],
     type="l", col="orange", main="Orange=IndUSA, Sky Blue=IndGBR")
lines(dtf[,.(DATE, detrend(IndGBR))], col="sky blue")
grid()
plot(dtf[,.(DATE, CpiGBR)], col="orange", type="l", main="Orange=CpiUSA, Sky Blue=CpiGBR")
lines(dtf[,.(DATE, CpiUSA)], col="sky blue")
grid()
plot(dtf[,.(DATE, IntGBR)], col="orange", type="l", main="Orange=IntUSA, Sky Blue=IntGBR")
lines(dtf[,.(DATE, IntUSA)], col="sky blue")
grid()


dtf[, `:=`(dXUMAGBD=c(NA_real_, diff(1/XUMAGBD)), dXUMLGBD=c(NA_real_, diff(1/XUMLGBD)))]
dtf[, `:=`(realUSA=IntUSA-CpiUSA, realGBR=IntGBR-CpiGBR)]
dtf[, realDiff := realUSA - realGBR]
dtf[, DDInd := c(NA_real_, diff(log(IndUSA)) - diff(log(IndGBR)))]
tmp <- tsHME("dXUMAGBD ~ .", DATA=dtf, DEPTH=2, NEXPERTS=2, TOL=1e-3,
             AR=list(list(1,FALSE), list(1,FALSE)), MAXIT=10)
respecify(tmp, criterion="BIC")
tmp2 <- tsHME("dXUMAGBD ~ .", DATA=dtf, DEPTH=2, NEXPERTS=2, TOL=1e-3,
              AR=list(list(1, FALSE), list(1:2, FALSE)), MAXIT=10)


debugonce(tsHME)
tmp3 <- tsHME("XUMAGBD ~ (IntUSA - IntGBR)", DATA=dtf, DEPTH=2, NEXPERTS=2, TOL=1e-3,
              AR=list(list(0,TRUE, FALSE, TRUE),
                      list(0,TRUE, FALSE,TRUE)),
              MAXIT=10, MAXEPOCHS=250)


respecify(tmp3, criterion = "BIC")
tmp4 <- tsHME("dXUMAGBD ~ realDiff + DDInd", DATA=dtf, DEPTH=2, NEXPERTS=3, TOL=1e-3,
              AR=list(list(1,FALSE), list(1:2,FALSE), list(1,FALSE)), MAXIT=10, MAXEPOCHS=700,
              PRIOR=tmp3$weights$prior)





