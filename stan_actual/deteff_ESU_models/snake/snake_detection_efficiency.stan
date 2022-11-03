// 03_snake_detection_efficiency
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
  
  // Make one vector that contains all parameters
  vector[N + K] params;
}

// this is where we estimate eta, our linear predictor
transformed parameters {
  
  // concatenate alphas and betas into params vector
  params[1:J] = alpha;
  params[J+1:J+K] = beta;
  
  // calculate linear predictor
  eta = params * X;
  
}


model {
  // priors
  
  
  
  y ~ bernoulli_logit(eta);
}
