// stan_mlogit_simulation

data {
  int N;
  array[N] int y;
}


parameters {
  real beta1;
  real beta3;
}


transformed parameters {
  // derived proportions
  simplex[5] p_vec;
  // vector for logits
  vector[5] logits;
  // first prob
  logits[1] = beta1;
  // effectively a 0 in logit space
  logits[2] = -100000;
  // second prob
  logits[3] = beta3;
  // effectively a 0 in logit space
  logits[4] = -100000;
  // loss param
  logits[5] = 0;
  // proper proportion vector
  p_vec = softmax(logits);
}


model {

  for (n in 1:N){
    // print("y[n]: ", y[n]);
    // print("p_vec, final: ", p_vec);
    target += categorical_lpmf(y[n] | p_vec);
  }
  
// target += categorical_lpmf(y | p_vec);  
}
