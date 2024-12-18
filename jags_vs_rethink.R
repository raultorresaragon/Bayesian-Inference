library(rjags)

# ~~~~ #
# JAGS #
# ~~~~ #
linear_model_code <- "
  data{
    D <- dim(x)
    n <- D[1]
    p <- D[2]
  }
  
  model{
    
   for(i in 1:n){
      # likelihood
      y[i] ~ dnorm(mu[i], tau)
      
      # posterior predictive
      ynew[i] ~ dnorm(mu[i], tau)
   }
    
    # conditional mean using matrix algebra
    mu <- x %*% beta
    
    # priors
    for(j in 1:p){
      beta[j] ~ dnorm(bm[j], pow(bs[j], -2))
    }
    sigma ~ dexp(lambda)
    tau <- pow(sigma, -2)
    
  }
"  
X <- model.matrix(~ . -1, 
                  data=d[, c("wine","judge")], 
                  contrasts.arg = lapply(d[,c("wine","judge")], contrasts, contrasts = FALSE))
y <- d$score_s

lm_jags <- jags.model(file = textConnection(linear_model_code), 
                      data = list(x = X,
                                  y = y,
                                  bm = c(rep(0, 20), rep(0, 9)),
                                  bs = c(rep(.5, 20), rep(1,9)),
                                  lambda = 1))

nsim <- 1e4
lmjags_samples <- coda.samples(model = lm_jags, variable.names = c("beta", "sigma"), n.iter = nsim)
lmjags_df <- as.data.frame(lmjags_samples[[1]])
names(lmjags_df) <- c(paste0("alpha", 1:20), paste0("beta", 1:9), "sigma")
precis(lmjags_df, prob = 0.95)



# ~~~~~~~~~~ #
# RETHINKING #
# ~~~~~~~~~~ #
d$judge_id <- coerce_index(d$judge)
d$wine_id  <- coerce_index(d$wine)
lm <- quap(
          alist(
             score_s ~ dnorm(mu, sigma),
             mu <- a[wine_id] + b[judge_id],
             a[wine_id] ~ dnorm(0,1),
             b[judge_id] ~ dnorm(0,1),
             sigma ~ dexp(1)
          ), data = d
      )
precis(lm, depth = 2)

d$flight_id <- coerc_index(d$flight)
d$wine.amer_id <- d$wine.amer + 1
d$judge.amer_id <- d$wine.amer + 1
lm2 <- quap(
   alist(
      score_s ~ dnorm(mu, sigma),
      mu <- a[flight_id] + b[wine.amer_id] + c[judge.amer_id],
      a[flight_id] ~ dnorm(0,.5),
      b[wine.amer_id] ~ dnorm(0,.5),
      c[judge.amer_id] ~ dnorm(0,.5),
      sigma ~ exp(1)
   ), data = d
)
precis(lm2, depth = 2)


x <- model.matrix(~ -1 + judge + wine,
                  data = d,
                  contrasts.arg = list(wine = contrasts(d$wine, contrasts = FALSE)))










