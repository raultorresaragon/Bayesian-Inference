library(tidyverse)
library(rethinking)
library(rjags)
rm(list = ls())
fd <- read_csv("project/gun_deaths full.csv")
ad <- read_csv("project/an_data.csv")


#~~~~~~~~~#
# MODEL 2 #
#~~~~~~~~~#

# create within state means for Y, X, and Z
ad <- ad |>
      group_by(state) |>
      mutate(mean_logY = mean(log_gun_deaths),
             mean_X = mean(gun_permits),
             mean_Z = mean(hunt_lic)) |> 
      mutate(D_Y = log_gun_deaths - mean_logY,
             D_X = gun_permits - mean_X,
             D_Z = hunt_lic - mean_Z)


model2_code <- "
model{
    #---------#
    # STAGE 1 #
    #---------#
    for(i in 1:n) {
                    
       # Likelihood
         DX[i] ~ dnorm(mu_x[i], pow(sigma[i], -2)) #tau)
         
       # Posterior predictive
         DXnew[i] ~ dnorm(mu_x[i], pow(sigma[i], -2))  #tau)

    }
         
    # Linear model
    mu_x = b_xz0 + b_xz1*DZ + years %*% pi
    #log(sigma) =  I %*% gamma
    
    # Because log is not vectorized we put this here
    for(i in 1:n) {
       log(sigma[i]) =  I %*% gamma
    }
    
    # Priors
    for(j in 1:ns){
        gamma[j] ~ dnorm(0, pow(1, -2))
    }    
    for(j in 1:nt){
        pi[j] ~ dnorm(0, pow(1, -2))
    }
    b_xz0 ~ dnorm(0, pow(1000, -2))
    b_xz1 ~ dnorm(1, pow(1, -2))
    
    # HOMOSKEDASTIC B/C CAN'T GET HETEROSKEDASTIC TO WORK
    # ---------------------------------------------------
    #sigma ~ dexp(lambda)
    #tau <- pow(sigma, -2)
    
    #---------#
    # STAGE 2 #
    #---------#
    for(i in 1:n) {
    
       # Likelihood stage 2
         DY[i] ~ dnorm(mu_y[i], pow(sigma[i], -2)) #tau)
    }
    
    # linear model
    mu_y = b_yxhat0 + b_yxhat1*DXnew + years %*% pis
    
    # Priors
    for(j in 1:nt){
        pis[j] ~ dnorm(0, pow(1000, -2))
    }
    b_yxhat0 ~ dnorm(0, pow(5, -2))
    b_yxhat1 ~ dnorm(0.2, pow(0.1, -2))   
}"

# construct fixed effects matrices
ad$year <- factor(ad$year)
ydummies <- model.matrix(~ . -1, 
                  data=ad[, "year"],  
                  contrasts.arg = lapply(ad[, "year"], 
                                         contrasts, contrasts = FALSE))
ad$state <- factor(ad$state)
sdummies <- model.matrix(~ . -1, 
                  data=ad[, "state"],
                  contrasts.arg = lapply(ad[, "state"],
                                         contrasts, contrasts = FALSE))

time_start <- Sys.time()

# fit model
m2 <- jags.model(file = textConnection(model2_code), 
                 data = list(DY = ad$D_Y,
                             DX = ad$D_X,
                             DZ = ad$D_Z,
                             years = ydummies,
                             I = sdummies,
                             #lambda = 500,
                             n = nrow(ad),
                             nt = length(unique(ad$year)),
                             ns = length(unique(ad$state)))#,
                 #n.chains = 3, 
                 #n.adapt  = 1e3
                 )
nsim <- 1e3

#update(m2, nsim)

m2_samples <- coda.samples(model = m2, 
                           variable.names = c("b_yxhat1","b_yxhat0","sigma"), 
                           n.iter = nsim)

m2_samples_full <- coda.samples(model = m2, 
                                variable.names = c("b_yxhat1","b_yxhat0","sigma","pis"), 
                                n.iter = nsim)

time_end <- Sys.time()
time_taken <- time_end - time_start
print(time_taken)

m2_df <- as.data.frame(m2_samples[[1]])

m2_df_full <- as.data.frame(m2_samples_full[[1]])


# MCMC diagnostics
par(mar=c(1,1,1,1))
m2_mcmc_checks <- plot(m2_samples)
m2_gelman <- gelman.diag(m2_samples)
m2_ess <- effectiveSize(m2_samples)

# Estimates
m2_precis <- precis(exp(m2_df[,1:2]), prob = 0.95, digits = 3)

save.image(file = "m2_workspace.RData")


# POSTERIOR CHECKS
#-----------------

# obtain expected value of y given x and credible interval for the AVERAGE year
post <- m2_df_full
mu_link <- function(x) {
    exp(post$b_yxhat0)  + 
    exp(post$b_yxhat1) * x + 
    post$`pis[1]` *(1/16) + post$`pis[2]` *(1/16) + 
    post$`pis[3]` *(1/16) + post$`pis[4]` *(1/16) +
    post$`pis[5]` *(1/16) + post$`pis[6]` *(1/16) + 
    post$`pis[7]` *(1/16) + post$`pis[8]` *(1/16) +
    post$`pis[9]` *(1/16) + post$`pis[10]`*(1/16) + 
    post$`pis[11]`*(1/16) + post$`pis[12]`*(1/16) +
    post$`pis[12]`*(1/16) + post$`pis[13]`*(1/16) + 
    post$`pis[14]`*(1/16) + post$`pis[15]`*(1/16) +
    post$`pis[16]`*(1/16)
}

gunp_seq <- seq(min(1e4*(ad$gun_permits/ad$population)), 
                max(1e4*(ad$gun_permits/ad$population)), length.out = 1e4)
mu <- sapply(gunp_seq, mu_link)
mu_median <- apply(mu, 2, median)
mu_ci <- apply(mu, 2, PI, prob = 0.95)

mu[1:5] 
mu_median[1:5] 
mu_ci[1:2,1:5] 

# obtain predicted interval
#ysim_fun <- function(x) {
#  rnorm(n = nrow(post),
#        mean = exp(post$b_yxhat0) + post$b_yxhat1*x,
#        sd = exp(post$sigma))
#}
#ysim <- sapply(gunp_seq, ysim_fun)
#ysim_pi <- apply(ysim, 2, PI, prob = 0.95)

# plot raw data
ad$gun_permits_p10 <- 1e4*(ad$gun_permits/ad$population)

par(mfrow = c(1,2))
plot(gun_deaths ~ gun_permits_p10, data=ad[sample(1:800, 500), ], 
     ylab = "gun deaths per 10,000 people", 
     xlab = "gun purchases per 10,000 people",
     main = "Median line and\ncredible intervals from M2",
     xaxt = "n",
     xlim = c(0, 1500),
     ylim = c(0, 4000))
axis(side = 1,
     at = seq(0, 3000, by=500), 
     labels = seq(0, 3000, by=500))

# plot MAP line
lines(gunp_seq, mu_median, lwd=2, col = "darkgreen")

# shaded region for credible interval
shade(mu_ci, gunp_seq, col = col.alpha("lightsalmon", 0.15))

# shaded region for predictive interval
#shade(ysim_pi, gunp_seq, col = col.alpha("lightskyblue1", 0.35))

# overlaying the observed data density 
# against the density of 200 simulated datasets
plot(density(ad$gun_deaths), lwd = 3, ylim = c(0, 0.0021),
     main = "Gun deaths observed\nand in 200 simulated datasets\n from M2")
for(i in 1:100) {
  y_new <- sapply(ad$gun_permits_p10, 
                  function(x) rnorm(1, 
                    exp(post$b_yxhat0 )+ exp(post$b_yxhat1)*x + 
                    post$`pis[1]` *(1/16) + post$`pis[2]` *(1/16) + 
                    post$`pis[3]` *(1/16) + post$`pis[4]` *(1/16) +
                    post$`pis[5]` *(1/16) + post$`pis[6]` *(1/16) + 
                    post$`pis[7]` *(1/16) + post$`pis[8]` *(1/16) +
                    post$`pis[9]` *(1/16) + post$`pis[10]`*(1/16) + 
                    post$`pis[11]`*(1/16) + post$`pis[12]`*(1/16) +
                    post$`pis[12]`*(1/16) + post$`pis[13]`*(1/16) + 
                    post$`pis[14]`*(1/16) + post$`pis[15]`*(1/16) +
                    post$`pis[16]`*(1/16), 
                    log(post$sigma))
                  )
  lines(density(y_new), col = col.alpha("dodgerblue1", 0.15))
}

