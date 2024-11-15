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

# p_vec <- c(exp(2.7)/(1 + exp(2.7) + exp(2.7)), 0, exp(2.7)/(1 + exp(2.7) + exp(2.7)), 0, 1 - (exp(2.7)/(1 + exp(2.7) + exp(2.7)) + exp(2.7)/(1 + exp(2.7) + exp(2.7))))

# 100 trials

set.seed(123)
dat <- rmultinom(1000, 1, prob = p_vec)

dat_2 <- vector(length = 1000)

# Now fill it in
for (i in 1:1000){
  det_hist <- dat[,i]
  dat_2[i] <- which(det_hist == 1, arr.ind = TRUE)
}


# Run stan model
data <- list(y = dat_2, N = 1000)

mod <- cmdstan_model("stan_mlogit_simulation.stan", compile = FALSE)

# Step 2: Compile the model, set up to run in parallel
mod$compile(cpp_options = list(stan_threads = TRUE))

# Step 3: Run MCMC (HMC)
fit <- mod$sample(
  data = data, 
  seed = 123, 
  chains = 1, 
  parallel_chains = 1,
  refresh = 10, # print update every 10 iters
  iter_sampling = 500,
  iter_warmup = 500,
  threads_per_chain = 7,
  init = 1,
)

# Run stan model again, this time with different struture
data2 <- list(y = dat_2, 
             N = 1000, # number of trials
             K = 5,
             D = 1,
             x = matrix(1, nrow = 100, ncol = 1))

mod <- cmdstan_model("stan_mlogit_simulation2.stan", compile = FALSE)

mod$compile(cpp_options = list(stan_threads = TRUE))

fit <- mod$sample(
  data = data2, 
  seed = 123, 
  chains = 1, 
  parallel_chains = 1,
  refresh = 10, # print update every 10 iters
  iter_sampling = 500,
  iter_warmup = 500,
  threads_per_chain = 7,
  init = 1,
)

fit

