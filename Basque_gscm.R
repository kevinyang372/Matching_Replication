install.packages("gsynth")
library(gsynth)
?gsynth

basquegscm <- basque
basquegscm$treat <- rep(0, 774)
basquegscm[704:731,"treat"] <- 1
#test
data(gsynth)
out <- gsynth(gdpcap ~ treat 
              # + sec.agriculture + sec.energy + sec.industry + sec.construction + sec.services.venta +sec.services.nonventa + school.illit + school.prim
              # + school.med + school.high + school.post.high + popdens + invest
              ,data = basquegscm, Y = gdpcap, D = treat, 
              # X = c(sec.agriculture, sec.energy, sec.industry , sec.construction , sec.services.venta , sec.services.nonventa , school.illit , school.prim
                                                               #, school.med ,school.high, school.post.high, popdens,invest)
              , parallel = FALSE, na.rm = FALSE,
              index = c("regionno","year"), force = "two-way",
              CV = TRUE, r = c(0, 5), se = FALSE) 
print(out) 
plot(out, type = "counterfactual")

?gsynth

# dataprep: prepare data for synth
dataprep.out <-
  dataprep(
    foo = basque
    ,predictors= c("school.illit",
                   "school.prim",
                   "school.med",
                   "school.high",
                   "school.post.high"
                   ,"invest"
    )
    ,predictors.op = c("mean")
    ,dependent     = c("gdpcap")
    ,unit.variable = c("regionno")
    ,time.variable = c("year")
    ,special.predictors = list(
      list("gdpcap",1960:1969,c("mean")),                            
      list("sec.agriculture",seq(1961,1969,2),c("mean")),
      list("sec.energy",seq(1961,1969,2),c("mean")),
      list("sec.industry",seq(1961,1969,2),c("mean")),
      list("sec.construction",seq(1961,1969,2),c("mean")),
      list("sec.services.venta",seq(1961,1969,2),c("mean")),
      list("sec.services.nonventa",seq(1961,1969,2),c("mean")),
      list("popdens",1969,c("mean")))
    ,treatment.identifier  = 17
    ,controls.identifier   = c(2:16,18)
    ,time.predictors.prior = c(1964:1969)
    ,time.optimize.ssr     = c(1960:1969)
    ,unit.names.variable   = c("regionname")
    ,time.plot            = c(1955:1997) 
  )

# 1. combine highest and second highest 
# schooling category and eliminate highest category
dataprep.out$X1["school.high",] <- 
  dataprep.out$X1["school.high",] + 
  dataprep.out$X1["school.post.high",]
dataprep.out$X1                 <- 
  as.matrix(dataprep.out$X1[
    -which(rownames(dataprep.out$X1)=="school.post.high"),])
dataprep.out$X0["school.high",] <- 
  dataprep.out$X0["school.high",] + 
  dataprep.out$X0["school.post.high",]
dataprep.out$X0                 <- 
  dataprep.out$X0[
    -which(rownames(dataprep.out$X0)=="school.post.high"),]

# 2. make total and compute shares for the schooling catgeories
lowest  <- which(rownames(dataprep.out$X0)=="school.illit")
highest <- which(rownames(dataprep.out$X0)=="school.high")

dataprep.out$X1[lowest:highest,] <- 
  (100 * dataprep.out$X1[lowest:highest,]) /
  sum(dataprep.out$X1[lowest:highest,])
dataprep.out$X0[lowest:highest,] <-  
  100 * scale(dataprep.out$X0[lowest:highest,],
              center=FALSE,
              scale=colSums(dataprep.out$X0[lowest:highest,])
  )

# run synth
synth.out <- synth(data.prep.obj = dataprep.out)

# Get result tables
synth.tables <- synth.tab(
  dataprep.res = dataprep.out,
  synth.res = synth.out
) 

path.plot(synth.res = synth.out,
          dataprep.res = dataprep.out,
          Ylab = c("real per-capita GDP (1986 USD, thousand)"),
          Xlab = c("year"), 
          Ylim = c(0,13), 
          Legend = c("Basque country","synthetic Basque country"),
) 
gaps.plot(synth.res = synth.out,
          dataprep.res = dataprep.out, 
          Ylab = c("gap in real per-capita GDP (1986 USD, thousand)"),
          Xlab = c("year"), 
          Ylim = c(-1.5,1.5), 
)