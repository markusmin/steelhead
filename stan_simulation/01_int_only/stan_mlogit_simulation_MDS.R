# multinomial/categorical logit simulation

# In this script, we will simulate data from a categorical logit similar in structure
# to our detection history data and fit a simple version of the model.

##### Simulate data

# here, let's simulate data with three possible outcomes (1, 3, and 5),
# and two outcomes that are not possible (2 and 4).
# State 5 will function as the loss state
# This mimics the actual transitions in that most states are not reachable and therefore
# there will always be some zeros.

# P_vec using mlogit
# beta1 = beta3 = 1
p_vec <- c(exp(1)/(1 + exp(1) + exp(1)), 0, exp(1)/(1 + exp(1) + exp(1)), 0, 1 - (exp(1)/(1 + exp(1) + exp(1)) + exp(1)/(1 + exp(1) + exp(1))))

# 100 trials

set.seed(123)
dat <- rmultinom(1000, 1, prob = p_vec)

dat_2 <- vector(length = 1000)

# Now fill it in
for (i in 1:1000){
  det_hist <- dat[,i]
  dat_2[i] <- which(det_hist == 1, arr.ind = TRUE)
}

# A new dat2 with high probability of moving to 3 or loss, but low to 1
dat2 <- c(rep(1, 450), rep(3, 100), rep(5, 450))


library(cmdstanr)

# Run stan model
data <- list(y = dat2, N = 1000)

mod <- cmdstan_model("stan_mlogit_simulation_MDS.stan", compile = FALSE)

# Step 2: Compile the model, set up to run in parallel
mod$compile(cpp_options = list(stan_threads = TRUE))

# Step 3: Run MCMC (HMC)
fit <- mod$sample(
  data = data, 
  seed = 123, 
  chains = 1, 
  parallel_chains = 1,
  refresh = 100, # print update every 100 iters
  iter_sampling = 1000,
  iter_warmup = 1000,
  threads_per_chain = 7,
  init = 1,
)


fit$summary(variables = c("p_vec"))
fit$summary()
