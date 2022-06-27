# multinomial/categorical logit simulation

# In this script, we will simulate data from a categorical logit similar in structure
# to our detection history data and fit a simple version of the model.
# 
# K <- 10
# N <- 11
# D <- 1
# y <- c(2,7,2,1,2,7,2,7,2,6,10)
# x <- matrix(1, nrow = 11, ncol = 1)
# beta <- matrix(1, nrow = D, ncol = K)
# 
# x_beta <- x %*% beta
# 
# 
# for (n in 1:N){
#   print(y[n])
#   t(x_beta[n])
# }


##### Simulate data

# here, let's simulate data with three possible outcomes (1, 3, and 5),
# and two outcomes that are not possible (2 and 4)

# P_vec using mlogit
p_vec <- c(exp(1)/(1 + exp(1) + exp(1)), 0, exp(1)/(1 + exp(1) + exp(1)), 0, 1 - (exp(1)/(1 + exp(1) + exp(1)) + exp(1)/(1 + exp(1) + exp(1))))

# 100 trials

set.seed(123)
dat <- rmultinom(100, 1, prob = p_vec)

dat_2 <- vector(length = 100)

# Now fill it in
for (i in 1:100){
  det_hist <- dat[,i]
  dat_2[i] <- which(det_hist == 1, arr.ind = TRUE)
}


# Run stan model
data <- list(y = dat_2, N = 100)

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

  