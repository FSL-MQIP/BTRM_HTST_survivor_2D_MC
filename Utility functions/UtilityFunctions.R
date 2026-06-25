## Utility Functions
# Last edited: 120820

## (1) buchanan_log10N, (2) gompertz_log10N, (3) baranyi_log10N
# Source: All equations for the 3 models below are copied from the Nlms package in R (https://rdrr.io/cran/nlsMicrobio/src/R/growthmodels.R)

# Purpose: Calculate log10N using respective growth model (either buchanan, gompertz, or barayani)

# Parameters: Same for (1) buchanan, (2) gompertz, and (3) barayani
# (i) t: time in hours
# (ii) lag: length of lag phase
# (iii) mumax: growth rate (ln/hour)
# (iv) LOG10N0: initial microbial concentration 
# (v) LOG10Nmax: carrying capacity

# Functions:
buchanan_nl_log10 = function(t,mumax,LOG10N0,LOG10Nmax){
  LOG10N <- LOG10N0 + (t <= ((LOG10Nmax - LOG10N0) /mumax)) * mumax * t + (t > ((LOG10Nmax - LOG10N0) /mumax)) * (LOG10Nmax - LOG10N0)
  return(LOG10N)
}


buchanan_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax){
  ans <- LOG10N0 + (t >= lag) * (t <= (lag + (LOG10Nmax - LOG10N0) * log(10)/mumax)) * mumax * (t - lag)/log(10) + (t >= lag) * (t > (lag + (LOG10Nmax - LOG10N0) * log(10)/mumax)) * (LOG10Nmax -     LOG10N0)
  return(ans)
}

gompertz_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10N0 + (LOG10Nmax - LOG10N0) * exp(-exp(mumax * exp(1) * (lag - t)/((LOG10Nmax - LOG10N0) * log(10)) + 1))
  return(ans)
}

baranyi_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10Nmax + log10((-1 + exp(mumax * lag) + exp(mumax * t))/(exp(mumax * t) - 1 + exp(mumax * lag) * 10^(LOG10Nmax - LOG10N0)))
  return(ans)
}


## (4) log10N_func

#Purpose: Implement the appropriate growth model based on provided model_name, in order to calculate log10N 

#Parameters:
# (i) t, (ii) lag, (iii) mumax, (iv) LOG10N0, & (v) LOG10Nmax: Same as for above functions for the 3 growth models
# (vi) model_name: Model to use for calcuating log10N

# Function:
log10N_func <- function(t, lag, mumax, LOG10N0, LOG10Nmax, model_name="buchanan") {
  if (model_name == "buchanan") {
    return(buchanan_log10N(t, lag, mumax, LOG10N0, LOG10Nmax) )
  }
  else if(model_name == 'baranyi') {
    return(baranyi_log10N(t, lag, mumax, LOG10N0, LOG10Nmax) )
  }
  else if(model_name == 'gompertz') {
    return(gompertz_log10N(t, lag, mumax, LOG10N0, LOG10Nmax) )
  }
  else {
    stop(paste0(model_name, " is not a valid model name. Must be one of buchanan, baranyi, gompertz"))
  }
}


## (5) muAtNewTemp

# Source: this function uses Ratkowsky's square root model which describes the effect of temperature on the growth of microorganisms (https://www.ncbi.nlm.nih.gov/pubmed/22417595)

# Purpose: Calculate the new mu parameter at new temperature.

# Parameters: 
# (i) newTemp: New temperature for which we calculate mu
# (ii) oldMu: Previous mu value to adjust
# (iii) oldTemp: Temperature corresponding to previous mu; NOTE: If you don't specify "oldTemp", 
#then automatically oldTemp = 6C; this value (6C) was the temp at which growth curve experiments were originally performed at)
# (iv) T0: Parameter used to calculate new mu; NOTE: If you don't specify "T0", then automatically T0 = -3.62C; this value (-3.62C) 
#was determined using Ratkowsky's square root model and Paenibacillus ordorifer growth curves obtained at 4, 7, and 32C in BHI broth (N.H. Martin unpublished data)

# Function:
muAtNewTemp <- function(newTemp, oldMu, oldTemp = 6, T0 = -3.62) { 
  numerator <- newTemp - T0
  denom <- oldTemp - T0
  newMu <- ((numerator / denom)^2) * oldMu
  
  return(newMu)
}


## (6) adjustLag

# Purpose: Adjust the lag phase based on the Zwietering 1994 paper.

# Parameters:  
# (i) t: Current timestep (in days)
# (ii) oldLag: Lag time (in days) at the previous temperature 
# (iii) newLag: Lag time (in days) at the current temparature 
# (iv) restartExp: If true, then lag phase restarts even if already in exponential growth phase; NOTE: If you don't 
#specify "restartExp", then automatically restartExp = T;
# (v) adjustmentConstant: Amount to adjust lag; NOTE: If you don't provide "adjustmentConstant", 
#then automatically adjustmentConstant = 0.25; Zwietering 1994 paper suggests using adjustmentConstant = 0.25

# Function:
adjustLag <- function (t, oldLag, newLag, restartExp = T, adjustmentConstant = 0.25) {
  #determine the amount of lag phase completed
  remainingLag <- 1 - (t / oldLag)
  if(restartExp) {
    remainingLag <- ifelse(remainingLag < 0, 0, remainingLag)
  }
  else {
    adjustedLag <- ifelse(remainingLag <=0, oldLag,
                          t + remainingLag * newLag + adjustmentConstant*newLag)
  }
  
  adjustedLag <- t + remainingLag*newLag + adjustmentConstant*newLag
  return(adjustedLag)
}


## (7) lagAtNewTemp

# Purpose: Calculate the new lag parameter at new temperature.

# Parameters:
# (i) t: Current timestep (in days)
# (ii) newTemp: New temperature for which we calculate lag
# (iii) oldLag: Previous lag value to adjust
# (iv) oldTemp: Temperature corresponding to previous lag; NOTE: If you don't specify "oldTemp", 
#then automatically oldTemp = 6C; this value (6C) was the temp at which growth curve experiments were originally performed at)
# (v) T0: Parameter used to calculate new mu; NOTE: If you don't specify "T0", then automatically T0 = -3.62C; this value (-3.62C) 
#was determined using Ratkowsky's square root model and Paenibacillus ordorifer growth curves obtained at 4, 7, and 32C in BHI broth (N.H. Martin unpublished data)

# Function:
lagAtNewTemp <- function (newTemp, oldLag, oldTemp = 6, T0 = 1.15) {
  numerator <- oldTemp -T0
  denom <- newTemp - T0
  newLag <- ( (numerator / denom)^2) * oldLag
  return(newLag)
}

## (8) getPrevRow

getPrevRow <- function(df, sim_run, milk_unit, day) {
  old_temp <- df[df$lot_id == sim_run & df$milk_unit == milk_unit & df$day==day-1,] 
}

# Purpose: 

# Parameters:
# (i) df: dataframe
# (ii) sim_run:
# (iii) milk_unit:
# (iv) day:


## (9) SampleAT
SampleAT = function(len) {
  AT <- vector()
  for (i in 1:(length(len))){
    AT_samp <- sample(AT_freq, 1,replace = T)
    #while(AT_samp == "AT_23" || AT_samp == "AT_159"){
      #AT_samp <- sample(AT_freq, 1,replace = T)}
    AT<- c(AT,AT_samp)
  }
  return(AT)
}




## (10) sample category

#Purpose: function to sample categories for genera in each lot iteration
#Input: 
  # Data frame that contains lot_ids, unit_ids, and taxa
  # Data frame that defines the category proportions of each genus

assign_category <- function(ModelData, genus_cat) {
  # Ensure the Genus column in both data frames is treated consistently
  ModelData <- ModelData %>%
    mutate(Genus = as.character(Genus))
  
  genus_cat <- genus_cat %>%
    mutate(Genus = as.character(Genus))
  
  # Extract unique Genus-Iteration combinations
  unique_genus_iteration <- ModelData %>%
    dplyr::select(lot_id, Genus) %>%
    distinct()
  
  # Assign a abundance category to each Genus-Iteration combination
  sampled_categories <- unique_genus_iteration %>%
    rowwise() %>%
    mutate(
      category = {
        # Extract the current Genus value
        current_genus <- Genus
        
        # Extract probabilities for the current genus
        genus_row <- genus_cat %>%
          filter(Genus == current_genus)
        
        if (nrow(genus_row) == 0) {
          # If no matching genus is found, return NA
          NA
        } else {
          # Select the A, B, C, D columns as probabilities
          probabilities <- unlist(genus_row[, c("A", "B", "C", "D")])
          
          # Check if probabilities are valid
          if (all(is.na(probabilities)) || sum(probabilities) == 0) {
            NA  # If invalid, assign NA
          } else {
            # Normalize probabilities to ensure they sum to 1
            probabilities <- probabilities / sum(probabilities)
            
            # Sample a category based on probabilities
            sample(c("A", "B", "C", "D"), 1, prob = probabilities)
          }
        }
      }
    ) %>%
    ungroup()
  
  # Join sampled categories back to the full ModelData frame
  ModelData <- ModelData %>%
    left_join(sampled_categories, by = c("lot_id", "Genus"))
  
  return(ModelData)
}





fix_simulation <- function(sim) {
  if (nrow(sim) < row_run) {
    # Remove last row (assumed to be decimal time)
    sim_fixed <- sim[-nrow(sim), ]
    
    # Get last two rows
    last_two <- sim_fixed %>% tail(2)
    
    # Check if last row has NaNs
    if (is.na(last_two$N[2]) || is.na(last_two$logN[2])) {
      # Replace last row's values with second-to-last row's values
      sim_fixed[nrow(sim_fixed), c("Q","N", "logN")] <- last_two[1, c("Q","N", "logN")]
      source_row <- last_two[1, ]
    } else {
      source_row <- last_two[2, ]
    }
    
    # Create new rows from next integer time to 40
    new_times <- seq(last_two[2,]$time + 1, row_run-1)
    new_rows <- tibble(
      time = new_times,
      Q = source_row$Q,
      N = source_row$N,
      logN = source_row$logN
    )
    
    # Combine and return
    bind_rows(sim_fixed, new_rows)
  } else {
    sim
  }
}


