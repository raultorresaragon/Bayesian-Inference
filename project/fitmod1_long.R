library(tidyverse)
library(rethinking)
library(rjags)
rm(list = ls())
fd <- read_csv("project/gun_deaths full.csv")
ad <- read_csv("project/an_data.csv")


#~~~~~~~~~#
# MODEL 1 #
#~~~~~~~~~#

m1s1_code <- "
model{
    #---------#
    # STAGE 1 #
    #---------#
    for(i in 1:n) {
                    
       # Likelihood
         X[i] ~ dnorm(mu_x[i], tau) #tau[i])
         
       # Posterior predictive
         Xnew[i] ~ dnorm(mu_x[i], tau) #tau[i])  
         
       # Sigma
       # sigma[i]) ~ dexp()
    }
         
    # Linear model
    mu_x = b_xz0 + b_xz1*Z 
    # log(sigma[i]) =  I %*% gamma
    
    # Priors
    b_xz0 ~ dnorm(0, pow(5000, -2))
    b_xz1 ~ dnorm(1, pow(1, -2))
    sigma ~ dexp(lambda)
    tau <- pow(sigma, -2)
}"

m1s2_code <- "
model{
    #---------#
    # STAGE 2 #
    #---------#
    for(i in 1:n) {
    
       # Likelihood stage 2
         Y[i] ~ dnorm(mu_y[i], tau) #tau[i])
    }
    
    # linear model
    mu_y = b_yxhat0 + b_yxhat1*logXnew
    
    # Priors
    b_yxhat0 ~ dnorm(2, pow(1, -2))
    b_yxhat1 ~ dnorm(0.2, pow(0.1, -2))   
}"

time_start <- Sys.time()
end_time <- Sys.time()
time._taken <- end_time - start_time
print(time.taken)

# fit model
m1s1 <- jags.model(file = textConnection(m1s1_code), 
                   data = list(#Y = ad$log_gun_deaths - mean(ad$log_gun_deaths),
                               X = ad$gun_permits - mean(ad$gun_permits),
                               Z = ad$hunt_lic - mean(ad$hunt_lic),
                               #years = ydummies,
                               #I = sdummies,
                               lambda = 1,
                               n = nrow(ad)))

nsim <- 1e3

m1s1_samples <- coda.samples(model = m1s1, 
                             variable.names = c("Xnew"), 
                             n.iter = nsim)


Xhat <- as.data.frame(m1s1_samples[[1]]) |> colMeans() |> unlist() |> `names<-`(NULL)

m1s2_code <- "
model{
    #---------#
    # STAGE 2 #
    #---------#
    for(i in 1:n) {
    
       # Likelihood stage 2
         Y[i] ~ dnorm(mu_y[i], tau) #tau[i])
    }
    
    # linear model
    mu_y = b_yxhat0 + b_yxhat1*logXhat
    
    # Priors
    b_yxhat0 ~ dnorm(2, pow(1, -2))
    b_yxhat1 ~ dnorm(0.2, pow(0.1, -2))   
    sigma ~ dexp(lambda)
    tau <- pow(sigma, -2)
}"

# fit model
m1s2 <- jags.model(file = textConnection(m1s2_code), 
                   data = list(Y = ad$log_gun_deaths,
                               logXhat = log(Xhat),
                              #years = ydummies,
                              #I = sdummies,
                               lambda = 1,
                               n = nrow(ad)),
                   n.chains = 3, 
                   n.adapt  = 1e3
)

nsim <- 1e3

m1s2_samples <- coda.samples(model = m1s2, 
                             variable.names = c("b_yxhat0","b_yxhat1","sigma"), 
                             n.iter = nsim,
                             n.chains = 3,
                             n.adapt = 1e3)


# MCMC DIAGNOSTICS
# ----------------
par(mar=c(1,1,1,1))
m1_mcmc_checks <- plot(m1_samples)
m1_gelman <- gelman.diag(m1_samples)
m1_ess <- effectiveSize(m1_samples)

# ESTIMATES
#----------
m1_precis <- precis(m1_df, precis = 0.95, digits = 3)

save.image(file = "m1_workspace.RData")


# POSTERIOR CHECKS
#-----------------

# obtain expected value of y given x and credible interval
post <- m1_df
mu_link <- function(x) post$b_yxhat0 + post$b_yxhat*x
gunp_seq <- seq(0, 2500000, length.out = 1e4)
mu <- sapply(gunp_seq, mu_link)
mu_median <- apply(mu, 2, median)
mu_ci <- apply(mu, 2, PI, prob = 0.95)

# obtain predicted interval
ysim_fun <- function(x) {
  rnorm(n = nrow(post),
        mean = post$b_yxhat0 + post$b_yxhat*x,
        sd = post$sigma)
}
ysim <- sapply(gunp_seq, ysim_fun)
ysim_pi <- apply(ysim, 2, PI, prob = 0.95)

# plot raw data
plot(gun_deaths ~ gun_permits, data=ad, 
     ylab = "gun deaths per 10K", 
     xlab = "gun purchases per 10K",
     main = "gun deaths and gun purchases \n per 10,000")

# plot MAP line
lines(gunp_seq, mu_median, lwd=2)

# shaded region for credible interval
shade(mu_ci, gunp_seq, col = col.alpha("lightsalmon", 0.55))

# shaded region for predictive interval
shade(ysim_pi, gunp_seq, col = col.alpha("lightskyblue1", 0.35))

# overlaying the observed data density 
# against the density of 200 simulated datasets
plot(density(ad$gun_deaths), lwd = 3, ylim = c(0, 0.03))
for(i in 1:200) {
  y_new <- sapply(ad$gun_permits, 
                  function(x) rnorm(1, 
                                    post$byxhat0[1:200] + post$byxhat[1:200]*x, 
                                    post$sigma[1:200]))
  lines(density(y_new), col = col.alpha("dodgerblue1", 0.15))
}







