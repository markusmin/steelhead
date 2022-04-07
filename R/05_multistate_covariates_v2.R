# 05_multistate_covariates

# Structure: Each individual detection history is a matrix, with rows as all possible states,
# and columns as site occupancy

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
                                     "WEL_upstream_other_trib_sites",
               "loss")

nstates <- length(all_states)

# Put together a matrix that shows the relationships between various states
state_relationships <- data.frame(state = all_states, ascend = NA, descend = NA,
                                  loss = "loss")


# Create an example detection matrix - an individual that ascended McNary, fell back,
# and then went to the John Day River
example_matrix <- matrix(rep(0, nstates*4), nrow = nstates, ncol = 4)
rownames(example_matrix) <- all_states

# populate with some example movements
example_matrix["mainstem, BON to MCN", 1] <- 1
example_matrix["mainstem, MCN to ICH or PRA", 2] <- 1
example_matrix["mainstem, BON to MCN", 3] <- 1
example_matrix["John Day River", 4] <- 1


# Use dmulti, but have one p vector for each of the l locations
# in effect, we would then have a p matrix, basically with from and to columns
# same as the transition matrix we previously constructed
# p will be mostly zeros




cat("
model {

#state-space likelihood 
for(i in 1:n.ind){ # Loop through the detection matrices for each individual
  for(l in 1:n.locations){
    for(t in 1:n.times){
    
    
    # Probability of ascending a dam, for an individual i at location l at time t
    # This way, each ascension probability can be generalized to a[l], rather than having
    # to write out separate ones
    
    # Ascend dam
    logit(a[i,l,t]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) 
    bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[l, t] + # temperature (continuous)
    bFlow * flow[l, t] + # flow (continuous)
    bSpill * spill[l, t] # spill (continuous)
    
    # Descend dam
    logit(d[i,l,t]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) 
    bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[l, t] + # temperature (continuous)
    bFlow * flow[l, t] + # flow (continuous)
    bSpill * spill[l, t] # spill (continuous)
    
    # Enter tributary
    # May need one line here for each tributary
    logit([i,l,t]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) 
    bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[l, t] + # temperature (continuous)
    bFlow * flow[l, t] + # flow (continuous)
    bSpill * spill[l, t] # spill (continuous)
    
    # Loss (end - use e to not confuse with location l index)
    logit(e[i,l,t]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) 
    bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[l, t] + # temperature (continuous)
    bFlow * flow[l, t] + # flow (continuous)
    bSpill * spill[l, t] # spill (continuous)
  
    
    # Put the movement probabilities into a vector
    # All other movements probabilities are zero
    p <- rep(0, nstates)
    
    # Get the index of the ascend, descend, etc. states
    p_ascend_index <- state_relationships[l,ascend]
    p_descend_index <- state_relationships[l,descend]
    
    p[p_ascend_index] <- a
    p[p_descend_index] <- d
    # put in the loss values
    p[nstates] <- e
    # Evaluate the multinomial likelihood for the counts of detection probabilities
    y[1:nstates] ~ dmulti(p[1:nstates], 1)
  }

}

}
", fill=TRUE, file=here::here("model_files", "full_model.txt"))