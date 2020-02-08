# Stat 624 Final Project Code
# Cason Wight

# This code uses the EM algorithm to estimate the Markov-chain parameters,
# p_11 and p_22. After showing that this model does not fit the data, an 
# initial approach to a better model is proposed.



set.seed(6)

# True infected data and number of students
influenza <- c(3,8,28,75,221,291,255,235,190,125,70,28,12,5)
n <- 743

# The other two time series in the SIR data set are unknown
susceptible <- c(n,rep(NA,length(influenza)))
infected <- c(0,influenza)
recovered <- c(0,rep(NA,length(influenza)))
cbind(susceptible,infected,recovered)

# Plotting the infected students
pdf("NumInfected.pdf")
par(mar = c(5,4,4,4) + 0.1)
plot(1:14,influenza, pch = 19, xlab = "Day", ylab = "Students with Influenza", main = "Number of Students with Influenza over Time")
axis(side=4,at=seq(0,300,by=74.3), labels = sapply(seq(0,40,by=10), function(x) paste0(x,"%")))
mtext("As a Percent of Population", side = 4, line = 2)
dev.off()

# EM Algorithm
get.ests <- function(influenza, new.n = n) {
  infected <- c(0,influenza)
  
  # Starting parameters
  ps <- matrix(c(.5,.5),nrow = 1)
  convergence <- FALSE
  j <- 1
  while(convergence == FALSE) {
    
    #EXPECTATION STEP
    for(i in 2:length(infected)){
      # E(St) = S{t-1} - It - p22 * I{t-1}
      susceptible[i] <- max(susceptible[i-1] - ceiling(max(infected[i]-ps[j,2]*infected[i-1],0)),0)
      # E(Rt)= n - St - It
      recovered[i] <- new.n - susceptible[i] - infected[i]
    }
    
    cbind(susceptible,infected,recovered)
    
    # MAXIMIZATION STEP
    likelihood <- function(curr.ps){
      log.like <- 0
      # Likelihood for St
      log.like <- log.like + sum(dbinom(susceptible[-1],susceptible[-length(susceptible)],curr.ps[1], log = TRUE))
      # Likelihood for It
      log.like <- log.like + sum(dbinom(infected[-1]+diff(susceptible), infected[-length(infected)],curr.ps[2], log = TRUE))
      -log.like
    }
    
    curr.ps <- ps[j,]
    
    # MLE for the probabilities
    ps <- rbind(ps,optim(curr.ps,likelihood, method = "L-BFGS-B" , lower = c(0.00001,0.00001), upper = c(.99999,.99999))$par)
    
    # If the combined movement of both parameters exceeds .000001
    if(sum(abs(ps[j+1,]-ps[j,]))<.000001) convergence <- TRUE
    j <- j+1
    
  }
  
  est.ps <- ps[nrow(ps),]
  # Converged Parameters
  est.ps
}


### Simulation study

simulate <- function(transProbs, new.n = n, newmodel = FALSE){
  # Probability matrix determined by given probs p11 and p22
  pMatrix <- matrix(c(transProbs[1],1-transProbs[1],0,
                      0,transProbs[2],1-transProbs[2],
                      0,0,1), nrow = 3, ncol = 3, byrow = TRUE)
  # starting susceptible = n, starting infected and recovered = 0
  susceptible <- new.n
  infected <- 0
  recovered <- 0
  for(i in 2:(length(influenza)+1)){
    if(newmodel) { # If using the new model, p11 changes based on the number infected
      pMatrix[1,1] <- min(max(transProbs[1] - (infected[i-1])^transProbs[4]/(743*transProbs[3]),0),1)
      pMatrix[1,2] <- 1-pMatrix[1,1]
    }
    # A binomial draw to get how many stay susceptible
    susceptible[i] <- rbinom(1,susceptible[i-1],pMatrix[1,1])
    # A binomial draw to get how many stay infected
    infected[i] <- susceptible[i-1]-susceptible[i]+rbinom(1,infected[i-1],pMatrix[2,2])
    # Rt = n - St - It 
    recovered[i] <- new.n-sum(infected[i],susceptible[i])
  }
  # Exclude starting values
  cbind(susceptible[-1],infected[-1],recovered[-1])
}


set.seed(6)
# True parameters that we will try to estimate
true.ps <- c(.7,.4)
ests <- c()

# Simulate 100 data sets off of these true pars
for(i in 1:100){
  samp.influenza <- simulate(true.ps)[,2]
  # Use the EM algorithm on each data set to get 100 parameter estimate pairs
  ests <- rbind(ests,tryCatch({get.ests(samp.influenza)}, error = function(cond) c(NA,NA)))
}
# Get rid of those with NAs (due to impossible likelihood in some cases)
ests <- if(length(which(is.na(ests[,1])))==0) {ests} else ests[-which(is.na(ests[,1])),]

# Bias in the estimates
apply(ests,2,mean)-true.ps


# Graph the samples of the par estimates for both parameters
pdf("ParEstimates.pdf")
layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4, 1))
par(mai=rep(0.5, 4))


plot(density(ests[,1]), xlim = c(.6,.9), main = expression(paste("Density of ",p[11]," estimates")), xlab = "Estimate")
abline(v=mean(ests[,1]), lty = "dashed")
abline(v=.7, col = 2)


plot(density(ests[,2]), xlim = c(.3,.55), main = expression(paste("Density of ",p[22]," estimates")), xlab = "Estimate")
abline(v=mean(ests[,2]), lty = "dashed")
abline(v=.4, col = 2)

par(mai=c(0,0,0,0))
plot.new()
legend(x="center",nco=3,legend=c("Density of Samples", "Mean of Samples", "True Probability"), col = c(1,1,2), lty = c("solid","dashed","solid"),cex=.9)
dev.off()



### Results 

## EM for this data set

# Use EM algorithm to get pars
est.ps <- get.ests(influenza)
est.ps

# Get simulations of number infected from EM pars
samples <- sapply(1:100000, function(x) simulate(est.ps)[,2])
est.inf <- apply(samples,1,mean)
low <- apply(samples,1,quantile,.005)
upp <- apply(samples,1,quantile,.995)

# SSE
sum((est.inf-influenza)^2)

# Plot the estimates to the truth
pdf("EMtoActual.pdf")
par(mar = c(5,4,4,4) + 0.1)
plot(1:length(influenza),influenza, pch = 19, xlab = "Day", ylab = "Students with Influenza", main = "Expected Infected from EM to Actual")
points(1:length(influenza), est.inf, pch = 19, col = 2)
arrows(x0=1:length(influenza),y0=low,y1=upp,col=2, angle = 90, length = .03)
arrows(x0=1:length(influenza),y0=upp,y1=low,col=2, angle = 90, length = .03)
axis(side=4,at=seq(0,300,by=74.3), labels = sapply(seq(0,40,by=10), function(x) paste0(x,"%")))
mtext("As a Percent of Population", side = 4, line = 2)
legend("topright", legend = c("Observed","Expected", "99% PI"), col=c(1,2,2), pch=c(19,19,NA), lty = c("blank","blank","solid"), cex = .8)
dev.off()

## Alternative model

# Get lots of combinations of potential pars
xs <- t(sapply(1:100000, function(x) c(runif(3),rlnorm(1,0,.3))))

result <- c()
# Get the SSE for one sample of each par combo
for(i in 1:nrow(xs)) result[i] <- sum((simulate(xs[i,], newmodel = TRUE)[,2]-influenza)^2)

# Use the best fitting combination to simulate with those pars
samps <- sapply(1:10000, function(x) simulate(xs[which.min(result),], newmodel = TRUE)[,2])
means <- apply(samps, 1, mean)
low <- apply(samps, 1, quantile, .005)
upp <- apply(samps, 1, quantile, .995)

# SSE
sum((means-influenza)^2)

# Plot the estimates using those pars to the truth
pdf("AltModeltoActual.pdf")
par(mar = c(5,4,4,4) + 0.1)
plot(1:14,means, pch = 19, col = 2, ylim = c(0,max(c(upp,influenza))), 
     xlab = "Day", ylab = "Students with Influenza", main = "Expected Infected from Alternative Model to Actual")
arrows(x0=1:14, y0=low, y1 = upp, angle = 90, length = .03, col = 2)
arrows(x0=1:14, y1=low, y0 = upp, angle = 90, length = .03, col = 2)
points(1:14, influenza, col = 1, pch = 19)
axis(side=4,at=seq(0,300,by=74.3), labels = sapply(seq(0,40,by=10), function(x) paste0(x,"%")))
mtext("As a Percent of Population", side = 4, line = 2)
legend("topright", legend = c("Observed","Expected", "99% PI"), col=c(1,2,2), pch=c(19,19,NA), lty = c("blank","blank","solid"), cex = .8)
dev.off()
