# Functions for replication of main crisis results
# Date: 7/8/04

#Function to calculate exact match ATE from matchit output
calc.ate <- function(x){
  dta <- eval(x$data, sys.parent())
  weighted.var <- function(x, w) sum(w * (x - weighted.mean(x,w))^2)/(sum(w) - 1)
  n1 <- sum(x$psweights!=0 & dta$war2==1)
  n0 <- sum(x$psweights!=0 & dta$war2==0)
  wy1 <- weighted.mean(dta$dir[dta$war2==1],x$psweights[dta$war2==1])
  wy0 <- weighted.mean(dta$dir[dta$war2==0],x$psweights[dta$war2==0])
  ate <- wy1-wy0
  vary1 <- weighted.var(dta$dir[dta$war2==1],x$psweights[dta$war2==1])/n1
  vary0 <- weighted.var(dta$dir[dta$war2==0],x$psweights[dta$war2==0])/n0
  sdate <- sqrt(vary1+vary0)
  cat("ATE is:",round(ate,2))
  cat("\nSD(ATE) is:",round(sdate,2))
  cat("\n")
  xx <- c(ate,sdate,n1+n0)
}

#function to simulate logit results
sim.logit <- function(logitm){
  dta <- model.frame(logitm)
  cov <- summary(logitm)$cov.scaled
  coef <- logitm$coef
  library(mvtnorm)
  sims <- 1000
  simco <- rmvnorm(sims,coef,cov)
  steps <- rep(0,nrow(dta))
  diff3 <- matrix(0,length(steps),sims)
  ate <- matrix(0,nrow(dta),sims)
  for (i in 1:length(steps)){
    X <- c(1,dta$scm[i],dta$lctdir[i],0,dta$bef75[i],dta$scm2[i])
    p0 <- 1/(1+exp(-(simco%*%X)))
    X[4] <- 1
    p1 <- 1/(1+exp(-(simco%*%X)))
    p <- p1-p0
    diff3[i,] <- p
  }
  diff <- apply(diff3,2,mean)
  ate <- mean(diff)
  se <- sd(diff)
  list(ate = ate, se = se, n = nrow(dta))
}

#function for logistic correction to matched sample
mlogit <- function(tt){
  require(mvtnorm)
  dta <- eval(tt$data, sys.parent())
  xx <- dta[tt$matched,]
  fml <- as.formula(paste("dir ~", paste(c(tt$treat,"scm","bef75","nyt","lctdir"),collapse=" + ")))
  ff <- glm(fml,data=xx,family=binomial(link=logit))
  bb <- summary(ff)$coefficients[,1]
  ss <- summary(ff)$cov.un
  simbs <- rmvnorm(1000,mean=bb,sigma=ss)
  XX1 <- XX0 <- XX <- model.matrix(ff)[,names(ff$coefficients)[!is.na(ff$coefficients)]]
  ttt <- XX[,2]==1
  XX1[,2] <- rep(1,nrow(XX1))
  XX0[,2] <- rep(0,nrow(XX0))
  y0 <- cbind(1/(1+exp(-simbs%*%t(XX0))))
  y1 <- cbind(1/(1+exp(-simbs%*%t(XX1))))
  te <- apply(y1-y0,1,mean)
  te
}

#drawing circles
circles <- function(x, y, radius, col=NA, border=par("fg")) 
    { 
        nmax <- max(length(x), length(y)) 
        if (length(x) < nmax) x <- rep(x, length=nmax) 
        if (length(y) < nmax) y <- rep(y, length=nmax) 
        if (length(col) < nmax) col <- rep(col, length=nmax) 
        if (length(border) < nmax) border <- rep(border, length=nmax) 
        if (length(radius) < nmax) radius <- rep(radius, length=nmax) 
        theta <- 2* pi * seq(0, 355, by=5) / 360 
        ct <- cos(theta) 
        st <- sin(theta) 
        for(i in 1:nmax) 
            polygon(x[i] + ct * xinch(radius[i]), y[i] + st * yinch(radius[i]), 
                col=col[i], border=border[i],lty=1,lwd=1) 
#                col=col[i], border=border[i],lty=1,lwd=3) 
      }

#function to plot ATE by term
plot.term <- function(g0,sims=100){
  tt <- sort(unique(g0$data$term[g0$matched & g0$data$war2==1]))
  nt <- length(tt)
  dir1 <- matrix(0,nt,sims)
  dir0 <- matrix(0,nt,sims)
  ww <- matrix(0,nt,sims)
  for(j in 1:sims){
    x <- eval(g0$call, envir = parent.frame())
    mm <- na.omit(x$match.matrix)
    pdta <- x$data[as.character(mm[,1]),]
    wdta <- x$data[row.names(mm),]
    for(i in 1:length(tt)){
      idx <- wdta$term==tt[i]
      dir1[i,j] <- mean(wdta$dir[idx])
      dir0[i,j] <- mean(pdta$dir[idx])
      ww[i,j] <- sqrt(sum(idx)/sum(x$psweights))/(2*pi)
    } # end terms
  } #end sims
  dd <- matrix(0,nt,2)
  dd[,1] <- apply(dir1,1,mean)
  dd[,2] <- apply(dir0,1,mean)
  dd <- as.data.frame(dd)
  row.names(dd) <- tt
  names(dd) <- c("war","peace")
  ww <- apply(ww,1,mean)
  names(ww) <- tt
  #plotting
  par(mar=c(2, 2, 2, 2) + 0.1,cex.lab=0.6,cex.axis=0.6,mgp=c(1,0.5,0),cex.main=0.8)
  mtitle <- "War Cases and Matched Peace Cases"
  ytitle <- "Proportion of Decisions Supporting the Rights Claim"
  plot(tt,tt,type="n",ylim=c(0,1),main="",ylab=ytitle,xlab="Term")
  a <- (dd[,1]<dd[,2])
  arrows(tt[a],dd[,2][a],tt[a],(dd[,1][a]+ww[a]/4),lty=1,length=0.08,lwd=2,col="black")
  arrows(tt[a==0],dd[,2][a==0],tt[a==0],(dd[,1][a==0]-ww[a==0]/4),lty=2,length=0.08,lwd=2,col="grey")
  circles(tt,dd[,2],ww,border=0,col="grey")
  circles(tt,dd[,1],ww,border="black")
  legend(1980,1,pch=c(21,16),col=c("black","grey"),legend=c("War","Peace"))
  return(list(dirs = dd, weights = ww))
}


# Graphing the empirical distn of the cases
plot.dist <- function(dta){
  wdata <- dta[dta$war2==1,]
  pdata <- dta[dta$war2==0,]
  termsc <- unique(data.frame(cbind(dta$term,dta$scmedian)))
  levels <- termsc[,1]
  courts <- length(levels)
  #dirwar,numberwar,dirpeace,numberpeace
  store <- matrix(0,courts,4)
  for (i in 1:courts){
    wa <- wdata$term==levels[i]
    pa <- pdata$term==levels[i]
    store[i,2] <- sum(wa)
    if(store[i,2]!=0){store[i,1] <- mean(wdata$dir[wa])}
    store[i,4] <- sum(pa)
    if(store[i,4]!=0){store[i,3] <- mean(pdata$dir[pa])}
  }
  storeall <- cbind(termsc,store)
  storeall <- storeall[sort(termsc[,1],index.return=T)$ix,]
  all <- (store[,2]*store[,1]+store[,4]*store[,3])/(store[,2]+store[,4])

  #Reading in the circle function (contributed by Ross Ihaka)
  circles <- function(x, y, radius, col=NA, border=par("fg")) 
    { 
      nmax <- max(length(x), length(y)) 
      if (length(x) < nmax) x <- rep(x, length=nmax) 
      if (length(y) < nmax) y <- rep(y, length=nmax) 
      if (length(col) < nmax) col <- rep(col, length=nmax) 
      if (length(border) < nmax) border <- rep(border, length=nmax) 
      if (length(radius) < nmax) radius <- rep(radius, length=nmax) 
      theta <- 2* pi * seq(0, 355, by=5) / 360 
      ct <- cos(theta) 
      st <- sin(theta) 
      for(i in 1:nmax) 
        polygon(x[i] + ct * xinch(radius[i]), y[i] + st * yinch(radius[i]), 
                col=col[i], border=border[i],lty=1,lwd=3) 
    }
  pweight <- store[,4]/sum(c(store[,4],store[,2]))
  pweight <- pweight[store[,3]!=0]
  wweight <- store[,2]/sum(c(store[,4],store[,2]))
  wweight <- wweight[store[,1]!=0]
  aweight <- (store[,2]+store[,4])/sum(c(store[,4],store[,2]))
  # Now converting from numbers to area, A=pi*r^2, r=sqrt(A/pi)
  prweight <- sqrt(pweight/pi)
  wrweight <- sqrt(wweight/pi)
  arweight <- sqrt(aweight/pi)
  par(mar=c(2, 2, 2, 2) + 0.1,cex.lab=0.6,cex.axis=0.6,mgp=c(1,0.5,0),cex.main=0.8)
  ytitle <- "Proportion of Decisions Supporting the Rights Claim"
  xtitle <- "Political Ideology of the Court"
  plot(levels,store[,3],ylim=c(0,1),xlim=rev(c(-0.5,0.75)),type="n",ylab=ytitle,xlab=xtitle,main="")
  circles(termsc[,2][store[,3]!=0],store[,3][store[,3]!=0],prweight,border=0,col="grey")
  circles(termsc[,2][store[,1]!=0],store[,1][store[,1]!=0],wrweight,border=1)
  legend(-0.2,1,legend=c("War","No War"),pch=c(1,19),col=c(1,"grey"),bty="n",cex=0.8)
  rect(-0.2,0.85,-0.42,0.99)
  text(termsc[,2][store[,1]!=0]-0.08,jitter(store[,1][store[,1]!=0],
                         factor=60),termsc[,1][store[,1]!=0],cex=0.6,col="black")
}

# plot liberalness of cases across time (w/SEs)
plot.lib <- function(dta){
  terms <- sort(unique(dta$term))
  tlib <- matrix(0,length(unique(term)),4)
  for (i in 1941:2001){
    a <- dta$dir[dta$term==i]
    tlib[i-1940,1] <- mean(a)
    tlib[i-1940,2] <- length(a)
    ss <- sqrt(mean(a)*(1-mean(a))/length(a))
    tlib[i-1940,3:4] <- c(mean(a)-1.96*ss,mean(a)+1.96*ss)
  }
  par(mar=c(2, 2, 2, 2) + 0.1,cex.lab=0.6,cex.axis=0.6,mgp=c(1,0.5,0),cex.main=0.8)
  mtitle <- c("Civil Rights and Liberties Claims over Time")
  ytitle <- c("Proportion of Cases Supporting Rights Claim")
  xtitle <- c("Term")
  plot(terms,tlib[,1],type="n",xlab=xtitle,ylab=ytitle,col="black",lwd=3,main="",ylim=c(0,1))
  rect(1941.9,0,1945,1,col="lightgrey",lty=0)
  rect(1950.5,0,1953.6,1,col="lightgrey",lty=0)
  rect(1965,0,1973,1,col="lightgrey",lty=0)
  rect(1991,0,1991.3,1,col="lightgrey",lty=0)
  rect(2001.8,0,2001,1,col="lightgrey",lty=0)
  lines(terms,tlib[,1],col="black",lwd=3)
  lines(terms,tlib[,3],col="black",lty=2)
  lines(terms,tlib[,4],col="black",lty=2)
  box(col="black")
}

#plotting histogram of ATE
plot.h <- function(x,y,both=T){
  par(mar=c(2, 2, 2, 2) + 0.1,cex.lab=0.6,cex.axis=0.6,mgp=c(1,0.5,0),cex.main=0.8)
  xtit <- "Effect of War on Probability of Supporting the Rights Claim"
  if(both){
    hist(x,main="",freq=F,xlab=xtit,col="grey",
         xlim=range(c(x,y)), breaks=20)
    hh <- hist(y, plot=F, freq=F,add=T,breaks=30)
    lines(c(hh$breaks[1],hh$breaks), c(0,hh$density,0), type="s", lwd=3)
    legend(0.2,max(density(mg0)$y),lty=2:3,
           legend=c("90% CI","95% CI"))
    text(-0.017,8.43,"Nonwarcases",pos=4)
    text(0.14,2.2,"Warcases", pos=4)    
  } else {
    hist(x,main="",freq=F,xlab=xtit,col="grey",
         xlim=c(-0.25,0.05), breaks=30)
    legend(-0.4,max(density(mg0)$y),lty=2:3,
           legend=c("90% CI","95% CI"))
  }
  d <- quantile(x,probs=c(0.05,0.95))
  d2 <- quantile(x,probs=c(0.025,0.975))
  box()
  abline(v=d[1],lty=2)
  abline(v=d[2],lty=2)
  abline(v=d2[1],lty=3)
  abline(v=d2[2],lty=3)
}

# Condition density estimates of ideology
plot.ideol <- function(dta){
  par(mar=c(2, 2, 2, 2) + 0.1,cex.lab=0.6,cex.axis=0.6,mgp=c(1,0.5,0),cex.main=0.8)
  x1 <- dta$scmedian[dta$war2==1]
  x0 <- dta$scmedian[dta$war2==0]
  minobs <- min(c(x1,x0))
  maxobs <- max(c(x1,x0))
  dx <- density(x1,from=minobs,to=maxobs)
  dy <- density(x0,from=minobs,to=maxobs)
  xtitle <- "Political Ideology of the Court"
  plot(dx$x,dx$y,type="l",xlab=xtitle,ylab="Density",lty=1,col=1,lwd=3,xlim=rev(range(dx$x)))
  lines(dx$x,dy$y,lty=2,lwd=2)
  text(-0.34,2.4,"No War")
  text(0.35,3,"War")
}

#distribution of propensity scores
doverlay <- function(foo){
  par(mar=c(2, 2, 2, 2) + 0.1,cex.lab=0.6,cex.axis=0.6,mgp=c(1,0.5,0),cex.main=0.8)
  weights <- foo$data$psweights
  x <- foo$data$pscore
  treat <- foo$data$war2
  xmiss <- !is.na(x)
  x <- x[xmiss]
  treat <- treat[xmiss]
  weights <- weights[xmiss]
  minobs <- min(x)
  maxobs <- max(x)
  dx1 <- density(x[treat==1],from=minobs,to=maxobs)
  dx0 <- density(x[treat==0],from=minobs,to=maxobs)
  x1 <- x[treat==1&weights!=0]
                                        # sample from C distn b/c of weights
  x0 <- sample(x[treat==0], size=min(10000,(100*length(x[treat==0]))),
               replace=TRUE,prob=weights[treat==0]/sum(weights[treat==0]))
  d1 <- density(x1,from=minobs,to=maxobs)
  bw <- d1$bw #need to store bandwidth here so that sampled
                                        # density plot bw doesn't automatically become too fine
  d0 <- density(x0,from=minobs,to=maxobs,bw=bw)
  par(mfrow=c(1,2))
  
  plot(dx0$x,dx1$y,type="n",ylim=range(c(dx0$y,dx1$y,d1$y,d0$y)),
       ylab="Density",xlab="Propensity Score",main="All Cases")
  polygon(c(min(dx1$x), dx1$x, max(dx1$x)), c(0, dx1$y,0), col = "grey",border=0)
  lines(dx0$x,dx0$y,type="l",lwd=2)  
  legend(minobs,max(c(d1$y,d0$y,dx1$y,dx0$y)), col=c("black","grey"), lwd=2,
         legend=c("Peace","War"))
  plot(dx0$x,d0$y,type="n",ylim=range(c(dx0$y,dx1$y,d1$y,d0$y)),
          ylab="Density",xlab="Propensity Score",main="Matched Cases")
  polygon(c(min(d1$x), d1$x, max(d1$x)), c(0, d1$y,0), col = "grey",border=0)
  lines(dx0$x,d0$y,type="l",lwd=2)
  legend(minobs,max(c(d1$y,d0$y,dx1$y,dx0$y)), col=c("black","grey"), lwd=2,
         legend=c("Peace","War"))
  par(mfrow=c(1,1))
}


#Jitter plot
jit.plot <- function(foo){
  weights <- foo$data$psweights
  pscore <- foo$data$pscore
  treat <- foo$data$war2
  par(mar=c(2, 2, 2, 2) + 0.1,cex.lab=0.6,cex.axis=0.6,mgp=c(1,0.5,0),cex.main=0.8)
  jitp <- jitter(rep(1,length(treat)),factor=10)-(treat==0)
  cwt <- sqrt(weights)
  plot(pscore,xlim=range(na.omit(pscore)),ylim=c(-1,2),type="n",ylab="",xlab="Propensity Score",axes=F,main="Distribution of Propensity Scores")
  points(pscore[treat==1&weights!=0],jitp[treat==1&weights!=0],pch=18,cex=cwt[treat==1&weights!=0])
  points(pscore[treat==1&weights==0],jitp[treat==1&weights==0],pch=5,col="grey",cex=0.5)
  points(pscore[treat==0&weights==0],jitp[treat==0&weights==0],pch=5,col="grey",cex=0.5)
  points(pscore[treat==0&weights!=0],jitp[treat==0&weights!=0],pch=18,cex=cwt[treat==0&weights!=0])
  axis(1)
  text(sum(range(na.omit(pscore)))/2,1.5,"War Cases")
  text(sum(range(na.omit(pscore)))/2,-0.5,"Peace Cases")
  box()
}

# Plotting ideology over time
plot.sc <- function(dta){
  par(mar=c(2, 2, 2, 2) + 0.1,cex.lab=0.6,cex.axis=0.6,mgp=c(1,0.5,0),cex.main=0.8)
  terms <- sort(unique(dta$term))
  tscmed <- matrix(0,length(unique(dta$term)),2)
  for (i in 1941:2001){
  a <- dta$scmedian[dta$term==i]
  tscmed[i-1940,1] <- mean(a)}
  ytitle <- c("Political Ideology of the Court")
  xtitle <- c("Term")
  plot(terms,tscmed[,1],type="n",xlab=xtitle,ylab=ytitle,main="",col="black",lwd=3,ylim=c(-1,1))
  rect(1941.9,-1,1945,1,col="lightgrey",lty=0)
  rect(1950.5,-1,1953.6,1,col="lightgrey",lty=0)
  rect(1965,-1,1973,1,col="lightgrey",lty=0)
  rect(1991,-1,1991.3,1,col="lightgrey",lty=0)
  rect(2001.8,-1,2001,1,col="lightgrey",lty=0)
  lines(terms,tscmed[,1],col="black",lwd=3)
  box(col="black")
}

#by issue area 
match.by.issue <- function(dta){
  require(mvtnorm)
  require(lattice)
  sims <- 1000
  values <- 4
  ate <- matrix(0,sims,values)
  ss <- matrix(0,values,4)
  for(i in 1:values){
    dan <- matchit(war2 ~ lctdir + scm + term + nyt + bef75,
                   data = dta[dta$value==i,], replace=FALSE,
                   discard=1, reestimate=TRUE, maxit=500) 
    foo <- glm(dir~war2+scm+scm2+bef75+nyt+lctdir,
               data=dta[dta$value==i,],family=binomial(link=logit),
               subset=(dan$matched),maxit=500)
    ss[i,4] <- sum(dan$matched)
    cov <- summary(foo)$cov.scaled
    coef <- foo$coef
    simco <- t(rmvnorm(sims,coef,cov))
    X <- model.matrix(foo)
    X1 <- X
    X1[,2] <- rep(1,nrow(X))
    X0 <- X
    X0[,2] <- rep(0,nrow(X))
    p1 <- 1/(1+exp(-(X1%*%simco)))
    p0 <- 1/(1+exp(-(X0%*%simco)))
    w <- model.frame(foo)$war2
    d <- model.frame(foo)$dir
    for (j in 1:sims){
      ysim1 <- c(d[w==1],p1[w!=1,j])
      ysim0 <- c(p0[w==1,j],d[w!=1])
      ate[j,i] <- mean(ysim1-ysim0)}}
  ate.data <- as.vector(ate)
  issue2<- rep(c("Criminal Procedure","Civil Rights","First Amendment","Due Process"),each = sims)
  issue2 <- as.factor(issue2)
  ss[,1] <- apply(ate,2,mean)
  ss[,2] <- apply(ate,2,sd)
  ss[,3] <- apply(ate>0,2,mean)
  ss <- as.data.frame(ss)
  row.names(ss) <- c("Criminal Procedure","Civil Rights","First Amendment","Due Process")
  names(ss) <- c("ATE","SE","p-value","N")
  list(ate.data = ate.data, issue = issue2, out=ss)
}
