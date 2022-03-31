# 05_multistate_covariates

# we would need to get a list of site + time
# state matrix, that has ascend (a), fallback (f), stray (s), home (h), loss (l) options
# each covariate is a matrix that has rows = time, columns = state
# But then how do you account for not all options being available at each time step/each state?


# one model structure - looping by individual
cat("
model {

#state-space likelihood 
for(i in 1:n.ind){ # We have covariates that are individual-based: natal origin, run year
# we have covariates that are based on the site + time
  for(j in 1:n.locations){
    for(k in 1:n.times){
    
  a2[i,j,k] <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[j, k] + # temperature (continuous)
    bFlow * flow[j, k] + # flow (continuous)
    bSpill * spill[j, k] # spill (continuous)
    
  logit(f[i,j]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[j, k] + # temperature (continuous)
    bFlow * flow[j, k] + # flow (continuous)
    bSpill * spill[j, k] # spill (continuous)
    
  logit(s[i,j]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[j, k] + # temperature (continuous)
    bFlow * flow[j, k] + # flow (continuous)
    bSpill * spill[j, k] # spill (continuous)
    
  logit(h[i,j]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[j, k] + # temperature (continuous)
    bFlow * flow[j, k] + # flow (continuous)
    bSpill * spill[j, k] # spill (continuous)
    
  logit(l[i,j]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[j, k] + # temperature (continuous)
    bFlow * flow[j, k] + # flow (continuous)
    bSpill * spill[j, k] # spill (continuous)
    }
    
    # Put the movements into a vector
    p <- c(a, f, s, h, l)
    # Evaluate the multinomial likelihood for the counts of detection probabilities
    y[1:N] ~ dmulti(p[1:N], 2121)
    z[i,j] ~ dbern(z[i,j-1]*phi[i,j-1]) #state process 
    y[i,j] ~ dbern(z[i,j]*p) #observation process 
  }

}

}
", fill=TRUE, file=here::here("model_files", "full_model.txt"))


# A second model structure - looping through movements/state transitions

cat("
model {

#state-space likelihood 
for(i in 1:n.ind){ # We have covariates that are individual-based: natal origin, run year
# we have covariates that are based on the site + time
  for(j in 1:n.locations){
    for(k in 1:n.times){
    
  a2[i,j,k] <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[j, k] + # temperature (continuous)
    bFlow * flow[j, k] + # flow (continuous)
    bSpill * spill[j, k] # spill (continuous)
    
  logit(f[i,j]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[j, k] + # temperature (continuous)
    bFlow * flow[j, k] + # flow (continuous)
    bSpill * spill[j, k] # spill (continuous)
    
  logit(s[i,j]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[j, k] + # temperature (continuous)
    bFlow * flow[j, k] + # flow (continuous)
    bSpill * spill[j, k] # spill (continuous)
    
  logit(h[i,j]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[j, k] + # temperature (continuous)
    bFlow * flow[j, k] + # flow (continuous)
    bSpill * spill[j, k] # spill (continuous)
    
  logit(l[i,j]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[j, k] + # temperature (continuous)
    bFlow * flow[j, k] + # flow (continuous)
    bSpill * spill[j, k] # spill (continuous)
    }
    
    # Put the movements into a vector
    p <- c(a, f, s, h, l)
    # Evaluate the multinomial likelihood for the counts of detection probabilities
    y[1:N] ~ dmulti(p[1:N], 2121)
    z[i,j] ~ dbern(z[i,j-1]*phi[i,j-1]) #state process 
    y[i,j] ~ dbern(z[i,j]*p) #observation process 
  }

}

}
", fill=TRUE, file=here::here("model_files", "full_model.txt"))