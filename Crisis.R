#====================================================
# Crisis Data -- Replication of Main Results
# Date: 7/8/04
# Data file: "crisis4d.sav"
# Functions: "fns.R"
#
# Tested in R 1.9.0 (http://www.r-project.org/) with
# MatchIt Version 0.1-8 (http://gking.harvard.edu/matchit/)
# 
# Codebook:
#
# (1) The following variables are the same as in the Spaeth dataset:
#     dir, analu, value, term, lctdir
#
# (2) The following variables are coded as described in Lee Epstein,
#     Daniel E. Ho, Gary King & Jeff Segal, The Supreme Court During
#     Crisis: How War Affects Only Nonwar Cases, N.Y.U. L. Rev (2005)
#
#      war2 = 1 if during the day of the oral argument, the US was
#             involved in an interstate militarized dispute; 0 otherwise
#
#      warcri = 1 if during the day of the oral argument, the US was
#               involved in an interstate militarized dispute or international
#               crisis; 0 otherwise
#
#      bef75 = indicator variable for whether the case was heard in a
#              term before 1975, during which the Nixon appointees
#              substantially changed the ideological composition of
#              the court
#
#      warcase = 1 if a case is a result of the war; 0 otherwise
#
#      scm = a rescaled measure of the median Segal-Cover score (+0.5)
#
#      scm2 = scm^2
#
#      nyt = 1 if the case was reported on the front page of the New
#            York Times; 0 otherwise
#
#      rally2 = 1 if there was a rally effect; 0 otherwise
#
#      crisis = 1 if there was an international crisis; 0 otherwise
#
#====================================================

#clearing memory
remove(list=ls())

#loading libraries
library(Matching)
library(lattice)
library(mvtnorm)
library(rbounds)

court <- read.csv("/Users/kevin/desktop/Github/Matching_Replication/Court.csv",header=TRUE)
#reading in data
dta0 <- read.table("/Users/johannes/Desktop/Minerva/CS112/Matching_Replication/Crisis.txt",header=TRUE)
dta <- subset.data.frame(dta0, select= c(dir, war2, warcri, bef75, warcase, scm, scm2, lctdir, analu, value, term, nyt, rally2, crisis,scmedian))
dta <- na.omit(dta)

#reading in functions
source("/Users/kevin/desktop/Github/Matching_Replication/fns.R")

#=============================================================================
#  Main Table, Model 1 (Logit)

# initializing summary table
# columns: Non-warcases (ATE, SE, N), Warcases (ATE, SE, N)
st <- matrix(0,5,6)

# Non-warcases
logitm0 <- glm(dir~scm+lctdir+war2+bef75+scm2,family=binomial(link=logit),
             data=dta[dta$warcase==0,])
st[1,1:3] <- as.numeric(sim.logit(logitm0))

#Warcases
logitm1 <- glm(dir~scm+lctdir+war2+bef75+scm2,family=binomial(link=logit),
             data=dta[dta$warcase!=0,])
st[1,4:6] <- as.numeric(sim.logit(logitm1))

#=============================================================================
# Main Table, Model 2 (Exact matching except for Term)

#Non-warcases

psm.nonwar <- glm(war2 ~ lctdir + scm + nyt + bef75, data = dta[dta$warcase==0,], family = "binomial" )
X.nonwar <- psm.nonwar$fitted
Tr <- dta[dta$warcase==0,]$war2
Y <- dta[dta$warcase==0,]$dir

nonwar.match <- Match(Y = Y, Tr = Tr, X = X.nonwar, exact = FALSE)
mb.nonwar  <- MatchBalance(war2 ~ lctdir + scm + nyt + bef75, data=dta[dta$warcase==0,], match.out=nonwar.match, nboots=100)

nonwar.match <- Match(Y=Y, Tr = Tr, X = X.nonwar, exact = TRUE)
summary(nonwar.match)

st[2,1:3] <- c(nonwar.match$est[1], nonwar.match$se, nonwar.match$orig.treated.nobs - length(nonwar.match$index.dropped))

#Warcases
psm.war <- glm(war2 ~ lctdir + scm + nyt + bef75, data = dta[dta$warcase!=0,], family = "binomial" )
X.war <- psm.war$fitted
Tr <- dta[dta$warcase!=0,]$war2
Y <- dta[dta$warcase!=0,]$dir

war.match <- Match(Tr = Tr, X = X.war, exact = TRUE)
mb.war  <- MatchBalance(war2 ~ lctdir + scm + nyt + bef75, data=dta[dta$warcase!=0,], match.out=war.match, nboots=100)

war.match <- Match(Y=Y, Tr = Tr, X = X.war, exact = TRUE)
summary(war.match)

st[2,4:6] <- c(war.match$est[1], war.match$se, war.match$orig.treated.nobs - length(war.match$index.dropped))

#=============================================================================
# Main Table, Model 3 (Exact matching on everything)

#Non-warcases
psm.nonwar <- glm(war2 ~ lctdir + scm + nyt + bef75 + term, data = dta[dta$warcase==0,], family = "binomial" )
X.nonwar <- psm.nonwar$fitted
Tr <- dta[dta$warcase==0,]$war2
Y <- dta[dta$warcase==0,]$dir

nonwar.match <- Match(Tr = Tr, X = X.nonwar, exact = TRUE)
mb.nonwar  <- MatchBalance(war2 ~ lctdir + scm + nyt + bef75 + term, data=dta[dta$warcase==0,], match.out=nonwar.match, nboots=100)

nonwar.match <- Match(Y=Y, Tr = Tr, X = X.nonwar, exact = TRUE)
summary(nonwar.match)

st[3,1:3] <- c(nonwar.match$est[1], nonwar.match$se, nonwar.match$orig.treated.nobs - length(nonwar.match$index.dropped))

#Warcases
psm.war <- glm(war2 ~ lctdir + scm + nyt + bef75 + term, data = dta[dta$warcase!=0,], family = "binomial" )
X.war <- psm.war$fitted
Tr <- dta[dta$warcase!=0,]$war2
Y <- dta[dta$warcase!=0,]$dir

war.match <- Match(Tr = Tr, X = X.war, exact = TRUE)
mb.war  <- MatchBalance(war2 ~ lctdir + scm + nyt + bef75 + term, data=dta[dta$warcase!=0,], match.out=war.match, nboots=100)

war.match <- Match(Y=Y, Tr = Tr, X = X.war, exact = TRUE)
summary(war.match)

st[3,4:6] <- c(NA, NA, NA)

#=============================================================================
# Main Table, Models 4 & 5 (Pscore matching (without and with logistic
# adjustment))

# Calculating variance and point estimate by averaging
# over variability of matching estimates (due to ties in propensity
# score, we use multiple imputation principles here)

# Failed to match a good result
psm.nonwar <- glm(war2 ~ lctdir + scm + term + nyt + bef75 + 
                  I(term*term) + I(lctdir*lctdir) + I(scm*scm) + 
                  I(nyt*nyt) + I(bef75*bef75) + I(lctdir*scm) +
                  I(lctdir*term) + I(lctdir*nyt) + I(lctdir*bef75)
                  + I(scm*term) + I(scm*nyt) + I(scm*bef75)
                  + I(term*nyt) + I(term*bef75) + I(nyt*bef75), data = dta[dta$warcase==0,], family = "binomial")

X.nonwar <- psm.nonwar$fitted
Tr <- dta[dta$warcase==0,]$war2
Y <- dta[dta$warcase==0,]$dir

nonwar.match <- Match(Y= Y,Tr = Tr, X = X.nonwar, replace = FALSE)
summary(nonwar.match)
mb.nonwar  <- MatchBalance(war2 ~ lctdir + scm + nyt + bef75 + term, data=dta[dta$warcase==0,], match.out=nonwar.match, nboots=100)

nonwar.match <- Match(Y=Y, Tr = Tr, X = X.nonwar, exact = TRUE)
summary(nonwar.match)
mb.nonwar  <- MatchBalance(war2 ~ lctdir + scm + nyt + bef75 + term, data=dta[dta$warcase==0,], match.out=nonwar.match, nboots=100)

# Using GenMatch Instead
dta.nonwar <- dta[dta$warcase==0,]
X.nonwar.gen <- cbind(dta.nonwar$lctdir, dta.nonwar$scm, dta.nonwar$term, dta.nonwar$nyt, dta.nonwar$bef75)
Tr <- dta.nonwar$war2
Y <- dta.nonwar$dir



genout <- GenMatch(Tr=Tr, X=X.nonwar.gen, pop.size=20, max.generations=100, wait.generations=25, exact=c(FALSE, TRUE, TRUE, FALSE, FALSE))

# genout <- GenMatch(Tr=Tr, X=X.nonwar.gen, pop.size=200, max.generations=100, wait.generations=10, caliper = c(FALSE, .09, .09, FALSE, FALSE))

mout  <- Match(Tr=Tr, X=X.nonwar.gen, Weight.matrix = genout$Weight.matrix, exact=c(FALSE, TRUE, TRUE, FALSE, FALSE))
summary(mout)

mb  <- MatchBalance(war2 ~ lctdir + scm + term + nyt + bef75, data=dta.nonwar, match.out=mout, nboots=100)

mout.nonwar.nb  <- Match(Y=Y, Tr=Tr, X=X.nonwar.gen, Weight.matrix = genout$Weight.matrix)
summary(mout.nonwar.nb)

st[4,1:3] <- c(mout.nonwar.nb$est[1], mout.nonwar.nb$se, mout.nonwar.nb$orig.treated.nobs - length(mout.nonwar.nb$index.dropped))

sensout <- psens(mout.nonwar.nb, Gamma=1.1, GammaInc=.01)
plot(sensout$bounds$Gamma, sensout$bounds$`Upper bound`, xlab = "Gamma",  ylab = "Upper bound p-value")
grid(nx = NULL ,ny = F)
abline(h = 0.05, col = "red")

hlsens(mout.nonwar.nb, Gamma=1.5, GammaInc=.1, .1)

mout.nonwar.b  <- Match(Y=Y, Tr=Tr, X=X.nonwar.gen, BiasAdjust = TRUE, Weight.matrix = genout$Weight.matrix)
summary(mout.nonwar.b)

st[5,1:3] <- c(mout.nonwar.b$est[1], mout.nonwar.b$se, mout.nonwar.b$orig.treated.nobs - length(mout.nonwar.b$index.dropped))

psens(mout.nonwar.b, Gamma=1.5, GammaInc=.1)

# War Case

dta.war <- dta[dta$warcase!=0,]

dta.war$term[which(dta.war$term >= 1980)] <- 5
dta.war$term[which(dta.war$term >= 1970)] <- 4
dta.war$term[which(dta.war$term >= 1960)] <- 3
dta.war$term[which(dta.war$term >= 1950)] <- 2
dta.war$term[which(dta.war$term >= 1940)] <- 1

X.war.gen <- cbind(dta.war$lctdir, dta.war$scm, dta.war$term, dta.war$nyt, dta.war$bef75)
Tr <- dta.war$war2
Y <- dta.war$dir

genout.war <- GenMatch(Tr=Tr, X=X.war.gen, pop.size=200, max.generations=100, wait.generations=50, caliper=c(1, 0.5, 1, 1, 1))

mout  <- Match(Tr=Tr, X=X.war.gen, Weight.matrix = genout.war$Weight.matrix, caliper=c(1, 0.5, 1, 1, 1))
summary(mout)

mb  <- MatchBalance(war2 ~ lctdir + scm + term + nyt + bef75, data=dta.war, match.out=mout.war, nboots=100)

mout.war.nb  <- Match(Y=Y, Tr=Tr, X=X.war.gen, Weight.matrix = genout.war$Weight.matrix)
summary(mout.war.nb)

st[4,4:6] <- c(mout.war.nb$est[1], mout.war.nb$se, mout.war.nb$orig.treated.nobs - length(mout.war.nb$index.dropped))

psens(mout.war.nb, Gamma=1.5, GammaInc=.1)

# Bias Adjust

mout.war.b  <- Match(Y=Y, Tr=Tr, X=X.war.gen, BiasAdjust = TRUE, Weight.matrix = genout.war$Weight.matrix)
summary(mout.war.b)

st[5,4:6] <- c(mout.war.b$est[1], mout.war.b$se, mout.war.b$orig.treated.nobs - length(mout.war.b$index.dropped))

psens(mout.war.b, Gamma=1.5, GammaInc=.1)

#formatting table
st <- as.data.frame(st)
names(st) <- rep(c("ATE","SE","N"),2)
row.names(st) <- c("Logistic Model", "Exact Matching Except for Term",
                   "Exact Matching", "Genetic Matching",
                   "Genetic Matching (Bias Adjusted)")

#=============================================================================
# Main graphs
#

#------- SC Densities conditional on war --------------
trellis.device(device="pdf",file="denssc.pdf",color=FALSE,width=6.5,height=4)
plot.ideol(dta)
dev.off()

#------- empirical distn of cases --------------
trellis.device(device="pdf",file="fullv3.pdf",color=FALSE,width=6.5,height=4)
plot.dist(dta)
dev.off()

#------- Liberalness across terms --------------
trellis.device(device="pdf",file="termlib.pdf",color=FALSE,width=6.5,height=4)
plot.lib(dta)
dev.off()

#------- SC over terms --------------
trellis.device(device="pdf",file="termsc.pdf",color=FALSE,width=6.5,height=4)
plot.sc(dta)
dev.off()

#Comparison of propensity scores
#Density plots
foo <- matchit(war2 ~ lctdir + scm + term + nyt + bef75,
               data = dta, replace=FALSE, discard=1)
trellis.device(device="pdf",file="pscores.pdf",color=FALSE,width=6.5,height=4)
doverlay(foo)
dev.off()
#Jitter plot
trellis.device(device="pdf",file="pjit.pdf",color=FALSE,width=6.5,height=4)
jit.plot(foo)
dev.off()

#------- Histogram of ATE --------------
g0 <- matchit(war2 ~ lctdir + scm + term + nyt + bef75,
              data = dta[dta$warcase==0,], replace=FALSE, discard=1,
              reestimate=TRUE)
mg0 <- mlogit(g0)
g1 <- matchit(war2 ~ lctdir + scm + term + nyt + bef75,
              data = dta[dta$warcase!=0,], replace=FALSE, discard=1,
              reestimate=TRUE)
mg1 <- mlogit(g1)
trellis.device(device="pdf",file="atehist.pdf",color=FALSE,width=6.5,height=4)
plot.h(mg0,mg1)
dev.off()

#------- ATE by Term Pairs  --------------
g0 <- matchit(war2 ~ lctdir + scm + term + scm2 + nyt + bef75 + I(bef75*nyt)
        + I(term*nyt) + I(term*scm2) + I(term^2) + I(lctdir*nyt) +
        I(scm*nyt) + I(scm*lctdir) + I(lctdir*term^2) + I(nyt*term^2),
        data = dta[dta$warcase==0,], replace=TRUE, discard=0, reestimate=TRUE)
trellis.device(device="pdf",file="matchedATEall-tmp.pdf",color=FALSE,width=6.5,height=4)
dterm <- plot.term(g0,sims=20)
dev.off()

#------- By Issue area  --------------
foo <- match.by.issue(dta[dta$warcase==0,])
trellis.device(device="pdf",file="issues.pdf",color=FALSE,width=6.5,height=4)
xyplot(foo$issue ~ foo$ate.data,cex=1.5,ylab="",xlab="Treatment Effect",
         aspect=0.5, main = "War Effect by Issue Area",
         panel = function(x,y,...){ 
           panel.abline(v = 0, col = "grey", lwd = 1,...)
           panel.bwplot(x,y,...)})
dev.off()

#=============================================================================
# Post-treatment models

#Excluding nyt
dan <- matchit(war2 ~ lctdir + scm + term + bef75, exact=T,
               data=dta[dta$warcase==0,])
calc.ate(dan)

#Excluding lctdir
dan <- matchit(war2 ~ scm + term + nyt + bef75, exact=T,
               data=dta[dta$warcase==0,])
calc.ate(dan)

#Excluding both lctdir and nyt
dan <- matchit(war2 ~ scm + term + bef75, exact=T,
               data=dta[dta$warcase==0,])
calc.ate(dan)

#Sensitivity to SC score
dan <- matchit(war2 ~ lctdir + term + nyt + bef75, exact=T,
               data=dta[dta$warcase==0,])
calc.ate(dan)

#--------------------------------------------
# Now examining effects on government claims
#--------------------------------------------

# A. Effect on US Govt Position 
# Re-Reading in data
dta <- subset.data.frame(dta0, select = c(dir, war2, warcri, bef75,
                                 warcase, scm, scm2, lctdir, analu,
                                 value, term, nyt, rally2, crisis,
                                 warcase, uswin))
dta <- na.omit(dta)
g0 <- matchit(war2 ~ lctdir + scm + term + nyt + bef75,
              data = dta[dta$warcase==0,], replace=FALSE, discard=1,
              reestimate=TRUE)
neyman(uswin,g0)
g1 <- matchit(war2 ~ lctdir + scm + term + nyt + bef75,
              data = dta[dta$warcase!=0,], replace=FALSE, discard=1,
              reestimate=TRUE)
neyman(uswin,g0)

# B. Effect on Solicitor General's Position
dta <- subset.data.frame(dta0, select = c(dir, war2, warcri, bef75,
                                 warcase, scm, scm2, lctdir, analu,
                                 value, term, nyt, rally2, crisis,
                                 warcase, sgwin))
dta <- na.omit(dta)
g0 <- matchit(war2 ~ lctdir + scm + term + nyt + bef75,
              data = dta[dta$warcase==0,], replace=FALSE, discard=1,
              reestimate=TRUE)
neyman(sgwin,g0)
#not enough warcases
g1 <- matchit(war2 ~ lctdir + scm + term + nyt + bef75,
              data = dta[dta$warcase!=0,], replace=FALSE, discard=1,
              reestimate=TRUE)

