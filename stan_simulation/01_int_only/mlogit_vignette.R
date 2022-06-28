# source: https://vasishth.github.io/bayescogsci/book/multinomial-processing-tree-mpt-models.html

# Probabilities of different answers
Pr_NR <- function(a, t, f, c)
  1 - a
Pr_Neologism <- function(a, t, f, c)
  a * (1 - t) * (1 - f) * (1 - c) + a * t * (1 - f) * (1 - c)
Pr_Formal <- function(a, t, f, c)
  a * (1 - t) * (1 - f) * c +  a * t * (1 - f) * c
Pr_Mixed <- function(a, t, f, c)
  a * (1 - t) * f
Pr_Correct <- function(a, t, f, c)
  a * t * f
# true underlying values for simulated data
a_true <- .75
t_true <- .9
f_true <- .8
c_true <- .1
# Probability of the different answers:
Theta <- tibble(NR = Pr_NR(a_true, t_true, f_true, c_true),
                Neologism = Pr_Neologism(a_true, t_true, f_true, c_true),
                Formal = Pr_Formal(a_true, t_true, f_true, c_true),
                Mixed = Pr_Mixed(a_true, t_true, f_true, c_true),
                Correct = Pr_Correct(a_true, t_true, f_true, c_true))
N_trials <- 200
(ans <- rmultinom(1, N_trials, c(Theta)))

N_obs <- 50
complexity <- rnorm(N_obs, mean = 0, sd = 2)
## choose some hypothetical values:
alpha_f <- .3
# the negative sign indicates that
# increased complexity will lead to a reduced value of f
beta_f <- -.3
# f' as a linear function of complexity
f_prime <- alpha_f + complexity * beta_f
head(f_prime)

## probabilities f for each item
f_true <- plogis(f_prime)
head(f_true)

theta_NR_v <- rep(Pr_NR(a_true, t_true, f_true, c_true), N_obs)
theta_Neologism_v <- Pr_Neologism(a_true, t_true, f_true, c_true)
theta_Formal_v <- Pr_Formal(a_true, t_true, f_true, c_true)
theta_Mixed_v <- Pr_Mixed(a_true, t_true, f_true, c_true)
theta_Correct_v <- Pr_Correct(a_true, t_true, f_true, c_true)
theta_item <- matrix(
  c(theta_NR_v,
    theta_Neologism_v,
    theta_Formal_v,
    theta_Mixed_v,
    theta_Correct_v),
  ncol = 5)
dim(theta_item)

library(hesim)
sim_data_cx <- tibble(item = 1:N_obs,
                      complexity = complexity,
                      w_ans = c(rcat(N_obs,theta_item)))
sim_data_cx



##### 19.2.1.5 A hierarchical MPT in Stan #####

# Data:
N_item <- 20
N_subj <- 30
N_obs <- N_item * N_subj 

subj <- rep(1:N_subj, each = N_item)
item <- rep(1:N_item, time = N_subj)

complexity <- rep(rnorm(N_item, 0, 2), times = N_subj)

(exp_sim <- tibble(subj = subj,
                   item = item,
                   complexity = complexity))

# New parameters, in log-odds space:
tau_u_a <- 1.1
## generate subject adjustments in log-odds space:
u_a <- rnorm(N_subj, 0, tau_u_a)
str(u_a)

## convert the intercept to log-odds space:
alpha_a <- qlogis(a_true)
## a_h' in log-odds space:
a_h_prime <-  alpha_a + u_a[subj]
## convert back to probability space
a_true_h <- plogis(a_h_prime)
str(a_true_h)

f_true <- plogis(alpha_f + complexity * beta_f)

# Aux. parameters that define the probabilities:
theta_NR_v_h <- Pr_NR(a_true_h, t_true, f_true, c_true) 
theta_Neologism_v_h <- Pr_Neologism(a_true_h, t_true, f_true, c_true)
theta_Formal_v_h <- Pr_Formal(a_true_h, t_true, f_true, c_true)
theta_Mixed_v_h <- Pr_Mixed(a_true_h, t_true, f_true, c_true)
theta_Correct_v_h <- Pr_Correct(a_true_h, t_true, f_true, c_true)
theta_h <- matrix(
  c(theta_NR_v_h,
    theta_Neologism_v_h,
    theta_Formal_v_h,
    theta_Mixed_v_h,
    theta_Correct_v_h),
  ncol = 5)
dim(theta_h)

# simulated data:
(sim_data_h <- mutate(exp_sim,
                      w_ans = rcat(N_obs,theta_h)))

# Run the stan model

setwd("/Users/markusmin/Documents/CBR/steelhead/stan_simulation/01_int_only")

sim_list_h <-  list(N_obs = nrow(sim_data_h),
                    w_ans = sim_data_h$w_ans,
                    N_subj = max(sim_data_h$subj),
                    subj = sim_data_h$subj,
                    complexity = sim_data_h$complexity)

fit_sh <- stan(file = "mlogit_vignette.stan", data = sim_list_h, chains = 1, iter = 100)  









