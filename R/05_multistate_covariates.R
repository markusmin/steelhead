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

# Use as an example area: BON to MCN
# possible_next_states: mouth to BON, MCN to ICH or PRA, Hood River, Fifteenmile Creek, Deschutes River, John Day River, Umatilla River

# Loop through individuals - each individual is a matrix, with rows as possible sites, and columns as site visits


all_states = c("mainstem, mouth to BON", "mainstem, BON to MCN", 
                                     "mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS",
                                     "mainstem, RIS to RRE", "mainstem, RRE to WEL", 
                                     "mainstem, upstream of WEL", "mainstem, ICH to LGR",
                                     "mainstem, upstream of LGR", "Asotin Creek", 
                                     "Clearwater River",
                                     "Deschutes River", 
                                     "Entiat River", 
                                     "Fifteenmile Creek", 
                                     "Grande Ronde River", 
                                     "Imnaha River",
                                     "John Day River", 
                                     "Methow River", 
                                     "Okanogan River", 
                                     "Salmon River", 
                                     "Tucannon River", 
                                     "Umatilla River",
                                     "Walla Walla River",
                                     "Wenatchee River", 
                                     "Yakima River",
                                     "Upstream LGR other tributaries",
                                     "ICH to LGR other tributaries",
                                     "BON to MCN other tributaries",
                                     "WEL_upstream_other_trib_sites")


example_matrix <- matrix(rep(0, 29*4), nrow = 29, ncol = 4)
rownames(example_matrix) <- all_states

# populate with some example movements
example_matrix["mainstem, BON to MCN", 1] <- 1
example_matrix["mainstem, MCN to ICH or PRA", 2] <- 1
example_matrix["mainstem, BON to MCN", 3] <- 1
example_matrix["John Day River", 4] <- 1

cat("
model {

#state-space likelihood 
for(i in 1:n.ind){ # We have covariates that are individual-based: natal origin, run year
# we have covariates that are based on the site + time
  for(l in 1:n.locations){
    for(t in 1:n.times){
    
  a2[i,l,t] <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[l, t] + # temperature (continuous)
    bFlow * flow[l, t] + # flow (continuous)
    bSpill * spill[l, t] # spill (continuous)
    
  logit(f[i,l]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[l, t] + # temperature (continuous)
    bFlow * flow[l, t] + # flow (continuous)
    bSpill * spill[l, t] # spill (continuous)
    
  logit(s[i,l]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[l, t] + # temperature (continuous)
    bFlow * flow[l, t] + # flow (continuous)
    bSpill * spill[l, t] # spill (continuous)
    
  logit(h[i,l]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[l, t] + # temperature (continuous)
    bFlow * flow[l, t] + # flow (continuous)
    bSpill * spill[l, t] # spill (continuous)
    
  logit(l[i,l]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[l, t] + # temperature (continuous)
    bFlow * flow[l, t] + # flow (continuous)
    bSpill * spill[l, t] # spill (continuous)
    }
    
    # Put the movements into a vector
    p <- c(a, f, s, h, l)
    # Evaluate the multinomial litelihood for the counts of detection probabilities
    y[1:N] ~ dmulti(p[1:N], 2121)
    z[i,l] ~ dbern(z[i,l-1]*phi[i,l-1]) #state process 
    y[i,l] ~ dbern(z[i,l]*p) #observation process 
  }

}

}
", fill=TRUE, file=here::here("model_files", "full_model.txt"))