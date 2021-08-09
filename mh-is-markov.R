#R code for Assignment EBS2043
#Authors:
#Rebecca Walter, Student Number: i6186657
#Meet Shah, Student Number: i6196781
#This code needs the addition of truncated normal from the R library.
library(truncnorm)

#My seed is set to 2097
set.seed(2097)

#Function to calculate steady state probabilities in a first order, 2 state Markov Chain
#Input:
#p11: Scalar, Probability of transitioning from State 1 to State 1
#P21: Scalar, Probability of transitioning from State 2 to State 1
#Output:
#Steady state probability of State 1 and 2.
steady_state_calculator <- function(p11,p21) {
  temp = 1 + (1-p11)/p21
  p1 = 1/temp
  p2 = 1- p1
  return(list(p1,p2))
}

#Function to Simulate a First order, 2 state Markov Chain
#Input: 
#p11: Scalar, Probability of transitioning from State 1 to State 1
#P21: Scalar, Probability of transitioning from State 2 to State 1
#Observations: Number of Observed States
#Output:
#Markov Chain: Array containing the sequence of States
#Scalar with an estimator for steady state probability of state 1
#Scalar with an estimator for steady state probability of state 2
#Assumption that State 1 is started at
sim_markov_chain <- function(p11, p21, observations=1000) {
  p1 = steady_state_calculator(p11,p21)[[1]]
  markov_chain = rep(0, observations)
  starter <- runif(1)
  markov_chain[1][starter <= p1] = 1
  markov_chain[1][starter >= p1] = 2
  count_1 = 0
  count_2 = 0
  draws = rep(0,observations)
  for (i in 2:observations) {
    draw <- runif(1)
    draws[i] = draw
    if (markov_chain[i-1] == 1) {
      markov_chain[i][draw <= p11] = 1
      markov_chain[i][draw > p11] = 2
      count_1 = count_1 + 1
    }
    else {
      markov_chain[i][draw <= p21] = 1
      markov_chain[i][draw > p21] = 2
      count_2 = count_2 + 1
    }
  }
  if (markov_chain[observations]== 1) {
    count_1 = count_1 + 1
  }
  else {
    count_2 =  count_2 + 1
  }
  return(list(markov_chain,count_1/observations,count_2/observations))
}

#Function which counts the number of times a particular transition has occurred.
#Input:
#Markov Chain: Vector containing sequence of states
#Output:
#List containing the number of times state transitions have occurred (Eg 1 to 2)
#And first state in chain
transition_counter <- function(markov_chain){
  count_11 = 0
  count_12 = 0
  count_21 = 0
  count_22 = 0
  for (k in 2:length(markov_chain)){
    count_11[markov_chain[k-1] == 1 & markov_chain[k] == 1] = count_11 + 1
    count_12[markov_chain[k-1] == 1 & markov_chain[k] == 2] = count_12 + 1
    count_21[markov_chain[k-1] == 2 & markov_chain[k] == 1] = count_21 + 1
    count_22[markov_chain[k-1] == 2 & markov_chain[k] == 2] = count_22 + 1
  }
  return(list(count_11,count_12,count_21,count_22,markov_chain[1]))
}

#Function which gives the log-likelihood of a particular markov chain.
#Input:
#p11: Scalar, Probability of transitioning from State 1 to State 1
#P21: Scalar, Probability of transitioning from State 2 to State 1
#Counts: List of the Count of the Number of State Transitions in the chain, and
#first state in the chain
#Output:
#Log-likelihood: Scalar, log likelihood of probabilities given sequence.
log_likelihood <- function(p11,p21, counts) {
  transition_likelihood = counts[[1]] * log(p11) + counts[[2]] * log(1-p11) + counts[[3]] * log(p21) + counts[[4]] * log(1-p21)
  first_state_likelihood = 0
  first_state_likelihood[counts[[5]]==1] = steady_state_calculator(p11,p21)[[1]]
  first_state_likelihood[counts[[5]]==2] = steady_state_calculator(p11,p21)[[2]]
  log_likelihood = log(first_state_likelihood) + transition_likelihood
  return (log_likelihood)
}

#Function which gives the target density (posterior density) of a particular markov chain.
#Input:
#p11: Scalar, Probability of transitioning from State 1 to State 1
#P21: Scalar, Probability of transitioning from State 2 to State 1
#Counts: List of the Count of the Number of State Transitions in the chain, and
#first state in the chain
#Output:
#target density: which is proportional to the likelihood times the prior 
#(one in our case because we use the uniform dist.)
target_density <- function(p11,p21,counts) {
  first_state_likelihood = 0
  first_state_likelihood[counts[[5]]==1] = steady_state_calculator(p11,p21)[[1]]
  first_state_likelihood[counts[[5]]==2] = steady_state_calculator(p11,p21)[[2]]
  transition_likelihood = p11 ^ counts[[1]] * (1-p11) ^ counts[[2]] * p21 ^counts[[3]] * (1-p21) ^ counts[[4]]
  prior_likelihood = dunif(p11) * dunif(p21)
  return(first_state_likelihood*transition_likelihood*prior_likelihood)
}

#Function which draws a tuple of candidates for the Metropolis-Hastings algorithm.
#Uses a truncated normal function, which is dependent on the previously accepted draw
#Input:
#p11_old: Scalar, Previously drawn p11
#p21_old: Scalar, Previously drawn p21
#Sigma: Scalar, Standard deviation of truncated normal distribution
#Output:
#A draw containing two new probabilities
mh_candidate_draw <- function(p11_old,p21_old,sigma) {
  p11_new <- rtruncnorm(1,0,1,p11_old,sigma)
  p21_new <- rtruncnorm(1,0,1,p21_old,sigma)
  return(c(p11_new,p21_new))
}

#Function which gives the candidate density for the Metropolis-Hastings algorithm
#Uses a truncated normal density, dependent on previous draw
#Input:
#p11_new, p21_new: scalars, newly drawn probabilities
#p11_old, p21_old: scalars, previously accepted probabilities
#sigma: Scalar, Standard Deviation of truncated normal distribution
#Output:
#Density of the new draw (Assuming independence)
mh_candidate_density<-function(p11_new,p21_new, p11_old, p21_old,sigma) {
  return(dtruncnorm(p11_new, 0,1, p11_old,sigma)*dtruncnorm(p21_new,0,1,p21_old,sigma))
}

#Function calculating acceptance probability for Metropolis Hastings algorithm
#Input:
#Current Draw of Probabilities, Previous draw of Probabilities, (-2) draw of Probabilities.
#Standard deviation of candidate distribution
#Number of Transition counts in Markov chain.
#Output:
#Acceptance Probability
#Logs used to provide numerical stability
#Ultimately could use log_likelihood function
#But that prevents generalization to other priors
mh_acc_prob <- function(p11_new,p21_new, p11_old,p21_old,p11_old2, p21_old2, counts,sigma) {
  new_prob_nom <- log(target_density(p11_new,p21_new,counts))
  new_prob_denom <- log(mh_candidate_density(p11_new,p21_new,p11_old,p21_old,sigma))
  old_prob_nom <- log(target_density(p11_old,p21_old,counts))
  old_prob_denom <- log(mh_candidate_density(p11_old,p21_old,p11_old2,p21_old2,sigma))
  a =exp(new_prob_nom - old_prob_nom + old_prob_denom -new_prob_denom)
  return(min(1,a))
}

#Function to run the Metropolis Hastings Algorithm using a Truncated Normal Candidate Distribution
#Input:
#Markov Chain sample,
#Number of draws/iterations to perform
#Standard deviation of Candidate Truncated Normal
#Output:
#List Containing Accepted/Repeated Draws in Metropolis Hastings and
#Number of Acceptances
mh_algorithm <- function(markovchain, n_sim,sigma){
  draws <- matrix(ncol = 2, nrow= n_sim)
  counts <- transition_counter(markovchain)
  counter = 0
  draws[1,] <- c(0.5,0.5)
  draws[2,] <- mh_candidate_draw(draws[1,1], draws[1,2], sigma) #First Draw
  for (i in 3:n_sim) {
    current_draw = mh_candidate_draw(draws[i-1,1], draws[i-1,2],sigma) #Draw a new candidate
    alpha = mh_acc_prob(current_draw[1], current_draw[2], draws[i-1,1], draws[i-1,2],draws[i-2,1], draws[i-2,2],counts,sigma)
    #Get acceptance probability
    if (runif(1) <= alpha) { 
      draws[i,] = current_draw
      counter = counter +1
    }
    else {
      draws[i,] = draws[i-1,]
    }
  }
return(list(draws,counter))
}

#Function to perform Importance Sampling Algorithm
#Input:
#Number of Draws to Take
#Markov Chain
#Standard Deviation of Truncated Normal
#Output:
#List containing estimate for Expected Value of P11, P21
#And Weights used to calculate them.
importance_sampling <- function (n, markov_chain,sigma) {
  counts <- transition_counter(markov_chain)
  
  p11_estimate <- counts[[1]]/(counts[[1]]+counts[[2]])
  p21_estimate <- counts[[3]]/(counts[[3]]+counts[[4]])
  draws_p11 <- rtruncnorm(n,0,1,p11_estimate,sigma)
  draws_p21 <- rtruncnorm(n,0,1,p21_estimate,sigma)
  
  #Uncomment the following two lines and comment the former two when trying to use beta distribution as sampler
  #draws_p11 <- rbeta(n, counts[[1]]/10, counts[[2]]/10)
  #draws_p21 <- rbeta(n, counts[[3]]/10, counts[[4]]/10)
  
  weights <- rep(0,n)
  for (i in 1:n) {
  weights[i] = target_density(draws_p11[i], draws_p21[i],counts)/ dtruncnorm(draws_p11[i],0,1,p11_estimate,sigma)/ dtruncnorm(draws_p21[i],0,1,p21_estimate,sigma)
  
  #Un-comment following line to see results using log_likelihood as target density and truncated normal as sampling density).
  #weights[i] = log_likelihood(draws_p11[i], draws_p21[i],counts) - log(dtruncnorm(draws_p11[i],0,1,p11_estimate,sigma)) -log(dtruncnorm(draws_p21[i],0,1,p21_estimate,sigma))
  
  #Un-comment followling line to see results using beta density as sampling density.
  #weights[i] = target_density(draws_p11[i], draws_p21[i],counts)/ dbeta(draws_p11[i],counts[[1]]/10, counts[[2]]/10)/ dbeta(draws_p21[i],counts[[3]]/10, counts[[4]]/10)
  }
  
  #Scaling to increase variance and range, comment next line when using log-likelihood
  weights = weights* 10^274
  
  p11 <- sum(weights*draws_p11)/sum(weights)
  p21 <- sum(weights*draws_p21)/sum(weights)
  return(list(p11,p21,weights, draws_p11, draws_p21))
}

#Simulate Markov Chain
p11 <- 0.2
p21 <- 0.45
output_sim_MC <- sim_markov_chain(p11,p21,1000)
markovchain <- output_sim_MC[[1]]
counts <- transition_counter(markovchain)
steady_states <- steady_state_calculator(p11, p21)


#Metropolis-Hastings 
n<-10000
metropolis_hastings <- mh_algorithm(markovchain, n, 0.05)
burn_in <- 1000
draws_p11 <- metropolis_hastings[[1]][burn_in:n,1] #Burnt-in
draws_p21 <- metropolis_hastings[[1]][burn_in:n,2] #Burnt-in
acceptance_rate <- metropolis_hastings[[2]]/n
rmse_mh <- (((mean(draws_p11) - p11)^2 + (mean(draws_p21) - p21)^2)/2)^0.5

#MH - plotting the draws of p11 and p21 from MH algorithm. 
#the red line giving the mean of the draws; the green line giving the initial p11 and p21 from the Markov Chain
plot(1:n,metropolis_hastings[[1]][,1], type= 'l', main="Accepted draws of p11 in MH algorithm", ylab="p11",xlab="n")
abline(h=p11,col=2,lty=2)
abline(h=mean(metropolis_hastings[[1]][,1]),col=3,lty=2)
legend("bottomright",c("p11_MC","estimated_p11"),fill=c(2 ,3))

plot(1:n,metropolis_hastings[[1]][,2], type= 'l', main="Accepted draws of p21 in MH algorithm", ylab="p21",xlab="n")
abline(h=p21,col=2,lty=2)
abline(h=mean(metropolis_hastings[[1]][,2]),col=3,lty=2)
legend("topright",c("p21_MC","estimated_p21"),fill=c(2 ,3))

#Plot of burned-in draws from p11 and p21 resulted from MH algorithm 
plot(burn_in:n,draws_p11, type= 'l', main="Burned-in draws of p11 in MH algorithm", ylab="p11",xlab="n")
abline(h=p11,col=2,lty=2)
abline(h=mean(draws_p11),col=3,lty=2)
legend("bottomright",c("p11_MC","estimated_p11"),fill=c(2 ,3))

plot(burn_in:n,draws_p21, type= 'l', main="Burned-in draws of p21 in MH algorithm", ylab="p21",xlab="n")
abline(h=p21,col=2,lty=2)
abline(h=mean(draws_p21),col=3,lty=2)
legend("topright",c("p21_MC","estimated_p21"),fill=c(2 ,3))

#Plot of the Autocorrelation between lags of draws from p11 and lags of draws from p21
acf(metropolis_hastings[[1]][,1],10, main =  "ACF of draws for p11")
acf(metropolis_hastings[[1]][,2],10, main =  "ACF of draws for p21")



#Importance Sampling
n_is <- 10000
sigma_IS = 0.05
output_IS <- importance_sampling(n_is,markovchain, sigma_IS)
plot(1:n_is, output_IS[[3]], main="Importance Sampling calculated Weights ", xlab='Draw', ylab = 'Weight')
plot(1:n_is, output_IS[[4]], main="Importance Sampling draws for p11", xlab='Draw', ylab = 'Draw P11')
plot(1:n_is, output_IS[[5]], main="Importance Sampling draws for p21", xlab='Draw', ylab = 'Draw P21')
rmse_is <- (((output_IS[[1]]-p11) ^2 + (output_IS[[2]]-p21) ^2)/2)^0.5

# running importance sampling for a different number of n.
# then we calculate the different p11 and p21 for that specific n.
n_draws <- seq(10,10000,length=50)
e_p11  <- rep(0,50)
e_p21  <- rep(0,50)
for(i in 1:50){
  etheta = importance_sampling(n_draws[i],markovchain,sigma_IS)
  e_p11[i] = etheta[[1]]
  e_p21[i] = etheta[[2]]
}

# plot of p11 for different number of obersavtions n.
plot(n_draws,e_p11,xlab='n',ylab=expression(E(P11)),type='l', ylim=c(p11-0.1, p11+0.1), main="Importance Sampling p11 for different n")
abline(h=p11,col=2,lty=2)

# plot of p21 for different number of observations n.
plot(n_draws,e_p21,xlab='n',ylab=expression(E(P21)),type='l', ylim=c(p21-0.1, p21+0.1), main="Importance Sampling p21 for different n")
abline(h=p21,col=2,lty=2)



{
cat('Markov Chain Simulated with Probabilities: P11=' , p11, 'and P21=', p21, fill=TRUE)
cat('Probability of state 1 ', output_sim_MC[[2]],'Probability of state 2' , output_sim_MC[[3]], fill = TRUE)  
cat('Transitions from 1-1:', counts[[1]],'Transitions from 1-2:' , counts[[2]], fill = TRUE)
cat('Transitions from 2-1:' , counts[[3]],  'Transitions from 2-2:' , counts[[4]], fill= TRUE)
cat('steady state probabilitly state 1:' , steady_states[[1]] , fill= TRUE) 
cat('steady state probabilitly state 2:' , steady_states[[2]] , fill= TRUE) 
cat('\nFOR METROPOLIS HASTINGS ALGORITHM:', fill = TRUE)
cat('Number of Simulations:', n, fill = TRUE)
cat('Acceptance Rate:', acceptance_rate, fill= TRUE)
cat("Number of Burned-in draws:", 1000, fill=TRUE)
cat('Estimated P11 Mean:', mean(draws_p11), fill = TRUE)
cat('Estimated P11 Standard Deviation:', sd(draws_p11), fill = TRUE)
cat('Estimated P21 Mean:', mean(draws_p21), fill = TRUE)
cat('Estimated P21 Standard Deviation:', sd(draws_p21), fill = TRUE)
cat('Root Mean squared error:', rmse_mh, fill = TRUE)
cat('View plots for pre and post-burn in draws and autocorrelation function', fill = TRUE)
cat('\nFOR IMPORTANCE SAMPLING ALGORITHM with', n_is, 'draws', fill = TRUE)
cat('Estimated P11 Mean:', output_IS[[1]], fill = TRUE)
cat('Estimated P21 Mean:', output_IS[[2]], fill = TRUE)
cat('Root Mean Squared Error:', rmse_is, fill=TRUE)
cat('Check Plots for weights and draws and different mean estimates with different number of draws', fill = TRUE)
}
