###############################################
# Simple Panel Data Example                   #
###############################################

library(tictoc)
library(WeibullR)
library(optimr)

tic()

nobs <- 100
delta1 <- .1
delta2 <- .1


# Model visualization
plot(c(1,5,9),c(5,5,5), pch = 19, cex = 15, col = 2:4, ylim = c(3,7), xlim = c(-1,11), xlab = "", ylab = "")
arrows(x0 = c(2,4,6), x1 = c(4,2,8), y0 = c(4.5,5.5,5), y1 = c(4.5,5.5,5), length = .1, lwd = 2)
text(c(5,5,5,0) ~c(1,5,9,8), labels=c("Healthy", "Sick", "Dead"), cex=0.9, font=2, col = c(1,1,"white"))

# True Parameters for the model
#All transitional distributions ~gamma(a,b)
p21 <- .8
a12 <- 10
b12 <- 1.5
a21 <- 3.5
b21 <- 1
a23 <- 1.5
b23 <- .8

ogPars <- c(p21, a12, b12, a21, b21, a23, b23)

#####################################
# Generating Data from this process #
#####################################

# Assume all start in state 1 at time 0, 
# All ended in state 3 or 4, or were censored at time 300 in state 1 or 2

data <- list()
panel_data <- list()

# Function that generates next transition and cumulative time, given current state and time
transition <- function(curr_state, curr_time, pars){
  if(curr_state==1){
    # State 1 can only transition into state 2
    next_state <- 2
    # Simulate the time until the transition, given the current and randomly selected next state
    cum_time <- curr_time + rgamma(1,a12,b12)
  }else if(curr_state==2){
    # Randomly decide next transition state, given transitional probabilities
    next_state <- ifelse(rbinom(1,1,p21)==1,1,3)
    # Simulate the time until the transition, given the current and randomly selected next state
    cum_time <- curr_time + rgamma(1,ifelse(next_state==1,a21,a23),ifelse(next_state==1,b21,b23))
  } else{
    # No more transitions can occur if the current state is 3 or 4 (both are absorbing states)
    print("error: trying to transition from 3")
  }
  rbind(next_state,cum_time)
}

###############################
# Generate actual transitions #
###############################
for (i in 1:nobs) {
  # Each obs starts in state 1 at time 0
  state <- 1
  time <- 0
  # Transitions occur until state 3 or 4 is reached. Until then, the "transistion" function is called
  transitioning <- "Not done"
  out <- rbind(state, time)
  while(transitioning!="Done"){
    out <- cbind(out,transition(state,time,ogPars))
    state <- out[1,ncol(out)]
    time <- out[2,ncol(out)]
    # When the state is 3 or 4, transitioning stops
    if(state == 3) {
      transitioning <- "Done"
    }
  }
  out <- as.matrix(out, nrow = 2)
  rownames(out) <- c()
  colnames(out) <- c()
  data[[i]] <- out
}


###########################################
# Collecting Panel Data from this process #
###########################################

# Assume that we observe at 0, and at exponentially distributed intervals (mean = 50) until t=100
cutoff.time <- 150
avg.obs.time <- 10

for (i in 1:nobs) {
  actual.data <- as.matrix(data[[i]],nrow = 2)
  last.time <- 0
  ob_times <- numeric()
  # Observation times are determined at intervals, randomly picked from the exponential (1/50) until the total time is greater than 300
  while(last.time < cutoff.time){
    ob_times <- c(ob_times, sum(c(rexp(1,1/avg.obs.time),ob_times)))
    last.time <- ob_times[length(ob_times)]
  }
  # The last observation is always after 300, so we rid of that observation
  ob_times <- ob_times[-length(ob_times)]
  # "num" is the total number of observations on the process
  num <- length(ob_times)
  # "out" is a matrix that will have the panel data. For each observation, there are "num" observations and each have a time and 
  # an observed state at that time.
  out <- matrix(NA,nrow=2,ncol=num)
  
  for (j in 1:num) {
    # This line takes the actual state at the time of a given observation 
    out[1,j] <- actual.data[1,max(which(actual.data[2,]<ob_times[j]))]
  }
  # "out" now has all of the states at observation, and this line adds the times of observation
  out[2,] <- ob_times
  # If the state is 3 or 4 (greater than 2), then additional observations tell us nothing, so we delete these observations
  index <- which(out[1,]>2)
  if (length(index)>1) {
    out <- out[,-index[-1],drop=F]
  }
  # Reformats as matrix
  panel_data[[i]] <- as.matrix(out, nrow = 2)  
}


##############################################
# Produce pseudo-data based on panel data   #
#############################################

stochastic_data <- function(data_num, pars, num_samps){
  
  # Get the original panel data for the obs
  panel_dats <- as.matrix(panel_data[[data_num]], nrow = 2)
  final.out <- numeric()
  
  # Each sample should have its own set of transition paths to be used in the likelihood function
  for(i in 1:num_samps){
    matches <- FALSE
    # This loop generates paths for each sample until it makes paths that would have yielded 
    # the panel data that we have for that obs.
    while(!matches){
      # Starts at state 1 and time 0
      state <- 1
      time <- 0
      # Transitions will occur until reaching absorbing states or time > 300
      transitioning <- "Not done"
      out <- numeric()
      while(transitioning!="Done"){
        # Appends new transition to previous transitions (including first)
        out <- cbind(out,transition(state,time,pars))
        # Updates current state and time
        state <- out[1,ncol(out)]
        time <- out[2,ncol(out)]
        # If time is ever over 300, it is beyond the panel data window, so data past this point are worthless
        if(time > cutoff.time){
          transitioning <- "Done"
          # Transition after 300 must be ignored
          out <- out[,-ncol(out)]
          # Transitions will occur until reaching state 3  (absorbing state)
        }else if(state == 3) {
          transitioning <- "Done"
        }
      }
      # Appends the starting state and time to the transitions
      out <- cbind(c(1,0),out)
      rownames(out) <- c("State","Time")
      colnames(out) <- c(1:ncol(out))
      
      matches <- TRUE
      # This loop goes through each panel data point for the obs, and checks to see if the 
      # transition state immediately before the panel data time matches the panel data state
      # at that time. If there are any conflicts, "matches" will return FALSE.
      for(j in 1:ncol(panel_dats)){
        # Gets a panel data timestamp and state
        panel_time <- panel_dats[2,j]
        panel_state <- panel_dats[1,j]
        # Sees which of the pseudo transition times happen before the panel data point
        times.line.up <- out[2,] < panel_time
        # Picks out the last of the transitions that happen before the panel data point
        most.recent.time.col <- which.max(out[2,][times.line.up])
        # Gets the state at that transition before the panel data point
        state.at.recent.time <- out[1,most.recent.time.col]
        # Verifies that the two states (panel data and generated transition data) match
        col.state.matches <- state.at.recent.time == panel_state
        # If this or any other does not match, "matches" will return FALSE
        matches <- matches * col.state.matches
      }
    }
    
    #########################################################
    # This line throws out any simulated transitions that occur after the last observed panel data point. Is this valid?
    out.new <- as.matrix(out[,out[2,] < panel_dats[2,ncol(panel_dats)]],nrow = 2)
    #########################################################
    
    
    out.new <- as.matrix(rbind(out.new,NA),nrow=3)
    # Time of censoring is time of last simulated transition
    last.time <- max(out.new[2,])
    # Time of last panel data obs
    last.panel.time <- panel_dats[2,ncol(panel_dats)]
    
    # If there is at least one transition (from 1 at time 0)
    if(ncol(out.new)!=1){
      # Marks which transition is the last transition "TRUE" non-censored transitions are marked "FALSE"
      out.new[3,] <- c(rep(FALSE,ncol(out.new)-1),TRUE) #ifelse(out.new[1,ncol(out.new)] %in% c(1,2),TRUE,FALSE))
      # Takes out first transition, which is always 1 at time 0 and gets the times between transitions instead of cumulative times
      out.new <- rbind(out.new[1,-1],diff(out.new[2,]),out.new[3,-1])
      
      out.new[2,] <- ifelse(out.new[1,] == 3, out.new[2,], # If the state is 3, no need for censoring time
                            ifelse(out.new[3,]==TRUE, last.panel.time-last.time,# Else, if it is the last transition,
                                   # get how much time is between the last 
                                   # transition and the panel data time
                                   out.new[2,]))# otherwise, put the regular time
      # If there is only the one point (state 1 at time 0 with no other transitions)
    } else{
      out.new[2,] <- last.panel.time # Get censoring time
      out.new[3,] <- TRUE 
      
    }
    # Appends each sample of each observation into one wide matrix that can be run through the likelihood function
    final.out <- cbind(final.out,out.new) 
  }
  final.out
}


###########################################
#    Likelihood contribution function     #
#     for each "observed" transition      #
###########################################

# Still need to comment better
likelihood <- function(pars){
  my.data <- dat_gen.lik
  state <- my.data[,1]
  time <- my.data[,2]
  last <- my.data[,3]
  
  P21 <- pars[1]
  A12 <- pars[2]
  B12 <- pars[3]
  A21 <- pars[4]
  B21 <- pars[5]
  A23 <- pars[6]  
  B23 <- pars[7]  
  
  likelihood.contribution <- numeric()
  
  for(i in 1:nrow(my.data)) {
    if(last[i]==FALSE){
      
      likelihood.contribution[i] <- ifelse(state[i]==1, 
                                           log(P21) + dgamma(time[i],A21,B21,log=TRUE),
                                           dgamma(time[i],A12,B12,log=TRUE))
      
    } else{
      likelihood.contribution[i] <- ifelse(state[i]==1,  
                                           pgamma(time[i],A12,B12,lower.tail=FALSE,log=TRUE),
                                           ifelse(state[i]==2, 
                                                  log(P21 * pgamma(time[i],A21,B21,lower.tail=FALSE) +
                                                        (1-P21) * pgamma(time[i],A23,B23,lower.tail=FALSE)),
                                                  log(1-P21) + dgamma(time[i],A23,B23,log=TRUE)))
    }
  }
  -sum(likelihood.contribution)
}




##################################################
# SEM Algorithm: Iterate until convergence holds #
##################################################

# Assign arbitrary prior parameters
newPars <- rep(.5,7)
#newPars <- ogPars
# COnvergence must hold for three iterations, so must track 3 iterations of convergence
convergence    <- FALSE
curr.converge  <- FALSE
prev1.converge <- FALSE
prev2.converge <- FALSE
iteration <- 1

allPars <- numeric()

# Convergence is required to end, which is described at the top of this script
while(!convergence){
  # Updates parameter estimates with each iteration
  oldPars <- newPars
  
  # m_1=1, m_2=2, ..., m_n=n according to Aralis and Brookmeyer, which is why num_samps=iteration
  # For each obs, samples must be simulated that match the panel data for the obs
  generated_data <- lapply(1:nobs, function(x) stochastic_data(x, newPars, iteration))
  
  # Formats the data for the likelihhod function
  dat_gen.lik <- c(numeric(),numeric(),numeric())
  
  # Combines all of the generated data into a wide matrix
  for(i in 1:nobs){
    dat_gen.lik <- cbind(dat_gen.lik,generated_data[[i]])
  }
  
  # Transposes for convenience
  dat_gen.lik <- t(dat_gen.lik)
  
  # Optimizes the likelihood function to get the most likely parameters, given the generated 
  # data that matched the panel data
  newPars <- optimr(par = oldPars, fn = likelihood, upper = c(.99999,rep(50,6)), lower = rep(0.0001,7), method = "L-BFGS-B")$par
  
  # Reassigns convergence based on equation described at beginning
  prev2.converge <- prev1.converge
  prev1.converge <- curr.converge
  curr.converge <- max(abs(newPars-oldPars)/(abs(oldPars)+delta1)) < delta2
  
  # Checks to see that this iteration and the two prior meet the convergence threshold
  converge <- prod(prev2.converge,prev1.converge,curr.converge)
  
  # Tracks all parameters throughout process
  allPars <- rbind(allPars,newPars)
  
  # For convenience, prints parameters and iteration
  print(list(iteration,ogPars,newPars))
  
  iteration <- iteration + 1
}



allPars

list(ogPars,newPars)

toc()

