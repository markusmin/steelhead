// detection_efficiency_script
// this script applies to all tributaries, such that only one of these models must be run 
// and the outputs can be used as inputs for each of the ESU-level primary models.
// This stan model must be run prior to the primary stan model, in order to get posteriors for 
// detection efficiencies in tributaries that are then passed as data to the primary stan model.

data {
  // Fish detections at upstream and river mouth arrays
  // number of observations
  int N;
  // response (detected or not detected)
  int y[N];
  
  // total number of eras
  int J;
  
  // total number of tributaries
  int K;
  
  // The design matrix, X; populated with discharge and intercept terms for eras
  matrix[N, J+K] X;

}


parameters {
  // Initialize era (categorical site configuration) terms - alphas
  vector[J] alpha;
  // real tuc_era1;
  // real tuc_era2;
  // real aso_era1;
  // real aso_era2;
  // real imn_era1;
  
  // Initialize discharge slope terms (betas)
  vector[K] beta;
  // real beta_tuc;
  // real beta_aso;
  // real beta_imn;
  

}

// this is where we estimate eta, our linear predictor
transformed parameters {
  // Make one vector that contains all parameters
  vector[J + K] params;
  
  // concatenate alphas and betas into params vector
  params[1:J] = alpha;
  params[J+1:J+K] = beta;
  
  // calculate linear predictor
  vector[N] eta;
  eta = X * params;
  
}


model {
  // priors
  alpha ~ normal(0, 5); // this is the prior on the intercepts, aka the inverse logit of the detection probability at zero discharge (a biologically nonsensical parameter)
  // note that a value of 5 gives you 99%; a value of 0 gives you 0.5, and a value of -5 gives you 0.6%. So I think normal(0,5) is reasonable.
  beta ~ normal(0, 1); // this is the prior on the slopes; we would expect it to be around 0. All of our slopes in the Rmd are really close to 1, like under 0.1. So let's go normal (0,1),
  // which is honestly still too much. Hopefully it runs
  
  
  
  y ~ bernoulli_logit(eta);
}
