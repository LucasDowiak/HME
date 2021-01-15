library(data.table)

dtf <- fread("data/wage_equation.csv")
setnames(dtf, "occup", "occ80")

# Cross-Walk from Sarah Porter SAS code
occ_files <- dir("data", pattern="occ", full.names = T)[-3]
occ_csvs <- napply(occ_files, fread)

dtf <- merge(dtf, occ_csvs[["data/occ80_occ90.csv"]], by="occ80", all.x=T)
dtf[, occ90 := ifelse(is.na(occ90), occ80, occ90)]
dtf <- merge(dtf, occ_csvs[["data/occ90_occ00.csv"]], by="occ90", all.x=T)
dtf <- merge(dtf, occ_csvs[["data/occ_soccoda.csv"]], by="occ00", all.x=T)
dtf <- merge(dtf, occ_csvs[["data/occ_soccodew.csv"]], by="occ00", all.x=T)


onet_files <- c("data/Analytical_Thinking.csv", "data/Programming.csv", 
                "data/Analyzing_Data_or_Information.csv", "data/Design.csv",
                "data/Thinking_Creatively.csv", "data/Social_Perceptiveness.csv",
                "data/Social_Orientation.csv")

g_ <- function(x)
{
  x[x == "Not relevant"] <- "0"
  x[x == "Not available"] <- NA_character_
  return(as.integer(x))
}

process_onet_files <- function(x)
{
  tit <- gsub(".*\\/ *(.*?) *\\.csv*", "\\1", x)
  csv <- fread(x)
  
  for (v in intersect(names(csv), c("Importance", "Level"))) {
    newname <- paste(tit, v, sep="_")
    setnames(csv, v, newname)
    
    if (v == "Level") {
      xx <- csv[, eval(parse(text=newname))]
      csv[, (newname) := eval(parse(text=paste("g_(", newname, ")")))]
    }
  }
  if (tit != "Social_Orientation")
    csv[, Occupation := NULL]
  return(csv)
}

cobb_douglass <- function(K, L, a, b)
{
  stopifnot(all(c(a, b) > 0),
            all(c(a, b) < 1),
            a + b == 1)
  
  return(K**(a) * L**(b))
}

onet_csvs <- napply(onet_files, process_onet_files)
dtfONet <- Reduce(function(x, y) merge(x, y, by="Code"), onet_csvs)
dtfONet[, code := as.integer(substring(gsub("-", "", Code), 1, 6))]

a <- 2/3; b <- 1/3
dtfONet[, Analytics := cobb_douglass(Analyzing_Data_or_Information_Importance,
                                     Analyzing_Data_or_Information_Level, a, b)]
dtfONet[, Creativity := cobb_douglass(Thinking_Creatively_Importance,
                                      Thinking_Creatively_Level, a, b)]
dtfONet[, Design := cobb_douglass(Design_Importance, Design_Level, a, b)]
dtfONet[, Perceptive := cobb_douglass(Social_Perceptiveness_Importance,
                                      Social_Perceptiveness_Level, a, b)]

dtfONet[, Design := scale(Design)]
dtfONet[, Perceptive := scale(Perceptive)]
dtfONet[, Creativity := scale(Creativity)]
dtfONet[, Analytics := scale(Analytics)]

dtfONetMeans <- dtfONet[, .(code, Creativity, Design, Analytics, Perceptive)]
dtfONetMeans <- dtfONetMeans[, lapply(.SD, mean, na.rm=T), by=code]

dtf <- merge(dtf, dtfONetMeans, by.x="soccoda", by.y="code", all.x=T)

rm(g_, onet_files, onet_csvs, occ_files, occ_csvs, dtfONet, dtfONetMeans, 
   cobb_douglass, a, b, process_onet_files)



par(mfrow=c(2,2))


dtf[sex==1 & !is.na(Design), plot(density(Design),
                                  main="Design", xlab="", ylab="",
                                  col="goldenrod2", lty=4)]
dtf[sex==0 & !is.na(Design), lines(density(Design), col="skyblue", lty=1)]
grid()

dtf[sex==0 & !is.na(Analytics), plot(density(Analytics),
                                     main="Analytics", xlab="", ylab="",
                                     col="skyblue", lty=1)]
dtf[sex==1 & !is.na(Analytics), lines(density(Analytics), col="goldenrod2", lty=4)]
grid()

dtf[sex==0 & !is.na(Creativity), plot(density(Creativity),
                                      main="Creativity", xlab="", ylab="",
                                      col="skyblue", lty=1)]
dtf[sex==1 & !is.na(Creativity), lines(density(Creativity), col="goldenrod2", lty=4)]
grid()

dtf[sex==0 & !is.na(Perceptive), plot(density(Perceptive),
                                      main = "Perceptive", xlab="", ylab="",
                                      col="skyblue", lty=1)]
dtf[sex==1 & !is.na(Perceptive), lines(density(Perceptive), col="goldenrod2", lty=4)]
grid()

legend(-6, 1.35, legend = c("Male", "Female"), col=c("skyblue", "goldenrod2"),
       lwd=2, lty=c(1,4), xpd="NA")

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "center",inset = 0,
       legend = c("Male", "Female"), 
       col=c("skyblue", "goldenrod2"), lwd=5, cex=.5, horiz = TRUE)

summary(dtf[sex==0, .(Creativity, Design, Analytics, Perceptive)])
sum(complete.cases(dtf[sex==0, .(Design, Creativity, Analytics, Perceptive)]))
sum(complete.cases(dtf[sex==0, .(Creativity, Analytics, Perceptive)]))


