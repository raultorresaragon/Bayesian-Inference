library(tidyverse)
library(rethinking)
library(rjags)
rm(list = ls())
fd <- read_csv("project/gun_deaths full.csv")
ad <- read_csv("project/an_data.csv")



model1_code <- "
model{
    #---------#
    # STAGE 1 #
    #---------#
    for(i in 1:n) {
                    
       # Likelihood
         X[i] ~ dnorm(mu_x[i], pow(sigma[i], -2)) #tau)
         
       # Posterior predictive
         Xnew[i] ~ dnorm(mu_x[i], pow(sigma[i], -2)) #tau)
    }
         
    # Linear model
    mu_x = b_xz0 + b_xz1*Z 
    #log(sigma) =  I %*% gamma
    
    # Because log is not vectorized we put this here
    for(i in 1:n) {
       log(sigma[i]) =  I %*% gamma
    }
    
    # Priors
    b_xz0 ~ dnorm(0, pow(5000, -2))
    b_xz1 ~ dnorm(1, pow(1, -2))
    for(s in 1:ns) {
      gamma[s] ~ dnorm(0, pow(1, -2))
    }
    #sigma ~ dexp(lambda)
    #tau <- pow(sigma, -2)


    #---------#
    # STAGE 2 #
    #---------#
    for(i in 1:n) {
    
       # Likelihood stage 2
         Y[i] ~ dnorm(mu_y[i], pow(sigma[i], -2)) #tau)
    }
    
    # linear model
    mu_y = b_yxhat0 + b_yxhat1*Xnew
    
    # Priors
    b_yxhat0 ~ dnorm(2, pow(1, -2))
    b_yxhat1 ~ dnorm(0.2, pow(0.1, -2))   
}"

ad$state <- factor(ad$state)
sdummies <- model.matrix(~ . -1, 
                         data=ad[, "state"],
                         contrasts.arg = lapply(ad[, "state"],
                                                contrasts, contrasts = FALSE))


time_start <- Sys.time()

# fit model
m1 <- jags.model(file = textConnection(model1_code), 
                 data = list(Y = ad$log_gun_deaths,
                             X = 1e4*(ad$gun_permits/ad$population),
                             Z = 1e4*(ad$hunt_lic/ad$population) - 
                               mean(1e4*(ad$hunt_lic/ad$population)),
                             #years = ydummies,
                             I = sdummies,
                             #lambda = 1,
                             ns = 50,
                             n = nrow(ad))#,
                 #n.chains = 3, 
                 #n.adapt  = 1e3
)
nsim <- 1e3

#update(m1, nsim)

m1_samples <- coda.samples(model = m1, 
                           variable.names = c("b_yxhat1","b_yxhat0","sigma"), 
                           n.iter = nsim)
time_end <- Sys.time()
time_taken <- time_end - time_start
print(time_taken)

m1_df <- as.data.frame(m1_samples[[1]])


# MCMC DIAGNOSTICS
# ----------------
par(mar=c(1,1,1,1))
m1_mcmc_checks <- plot(m1_samples)
m1_gelman <- gelman.diag(m1_samples)
m1_ess <- effectiveSize(m1_samples)

# ESTIMATES
#----------
precis(m1_df, precis = 0.95, digits = 3)

save.image(file = "m1_workspace.RData")


# POSTERIOR CHECKS
#-----------------

# obtain expected value of y given x and credible interval
post <- m1_df
mu_link <- function(x) {
  exp(post$b_yxhat0) + exp(post$b_yxhat1)*x
}
gunp_seq <- seq(min(1e4*(ad$gun_permits/ad$population)), 
                max(1e4*(ad$gun_permits/ad$population)), length.out = 1e4)

mu <- sapply(gunp_seq, mu_link)
mu_median <- apply(mu, 2, median)
mu_mean <- apply(mu, 2, mean)
mu_ci <- apply(mu, 2, PI, prob = 0.95)

summary(mu_median)
summary(mu_mean)
mu_ci[1,] |> summary()
mu_ci[2,] |> summary()


# obtain predicted interval
ysim_fun <- function(x) {
  rnorm(n = nrow(post),
        mean = exp(post$b_yxhat0) + exp(post$b_yxhat1)*x,
        sd = log(post$sigma))
}
ysim <- sapply(gunp_seq, ysim_fun)
ysim_pi <- apply(ysim, 2, PI, prob = 0.95)

ysim_pi[1, ] |> summary()
ysim_pi[2, ] |> summary()
ysim_pi[ ,sample(1:800, 10)]
par(mfrow=c(1,2))

# plot raw data
ad$gun_permits_p10 <- 1e4*(ad$gun_permits/ad$population)
plot(gun_deaths ~ gun_permits_p10, data=ad[sample(1:800, 300),], 
     ylab = "gun deaths per 10,000 people", 
     xlab = "gun purchases per 10,000 people",
     main = "Median line and\ncredible intervals from M1",
     xaxt = "n",
     xlim = c(0, 1500),
     ylim = c(0, 4000))
axis(side = 1,
     at   =   seq(0, 3000, by=500), 
     labels = seq(0, 3000, by=500))

# plot MAP line
lines(gunp_seq, mu_median, lwd=2, col="darkgreen")

# shaded region for credible interval
shade(mu_ci, gunp_seq, col = col.alpha("lightsalmon", 0.45))

# shaded region for predictive interval
shade(ysim_pi, gunp_seq, col = col.alpha("lightskyblue1", 0.35))

# overlaying the observed data density 
# against the density of 200 simulated datasets
plot(density(ad$gun_deaths), lwd = 3, ylim = c(0, 0.0025), 
     main = "Gun deaths observed\nand in 200 simulated datasets\n from M1")
for(i in 1:200) {
  y_new <- sapply(ad$gun_permits_p10, 
                  function(x) rnorm(1, 
                                    exp(post$b_yxhat0) + exp(post$b_yxhat1)*x, 
                                    log(post$sigma)))
  lines(density(y_new), col = col.alpha("dodgerblue1", 0.15))
}







