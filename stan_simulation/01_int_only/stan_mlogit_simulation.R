# multinomial/categorical logit simulation

# In this script, we will simulate data from a categorical logit similar in structure
# to our detection history data and fit a simple version of the model.

K <- 10
N <- 11
D <- 1
y <- c(2,7,2,1,2,7,2,7,2,6,10)
x <- matrix(1, nrow = 11, ncol = 1)
beta <- matrix(1, nrow = D, ncol = K)

x_beta <- x %*% beta


for (n in 1:N){
  print(y[n])
  t(x_beta[n])
}


