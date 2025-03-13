# ebv functions -------------------------------------------------------------------------------

# Calculate local contribution to Beta diversity
calc_lcbd <- function(comm_M, div_index = "BJ" ,  component = "repl" ){
  
  betadiv   <- beta.div.comp(comm_M , coef  = div_index , quant = FALSE) 
  lcbd_list <- LCBD.comp(betadiv[[component]] , sqrt.D =TRUE)
  data.frame(lcbd = lcbd_list$LCBD,trapID= rownames(comm_M))
  
}

# Function to calibrate reads against biological spike ins
calibrate_reads <- function(asv_counts , spike_in_ids, sample_ids){
  
  setkey(asv_counts, ASV_ID)
  spike_counts  <-  asv_counts[.(spike_in_ids), ..sample_ids, nomatch = 0]
  spike_sums    <- colSums(spike_counts)
  correction    <- log10(spike_sums) - mean(log10(spike_sums[spike_sums!=0]))
  correction[which(is.infinite(correction))] <- 1
  
  # array too big to vectorise this operation...
  calib_counts  <- lapply(seq_along(asv_counts)[-1] , 
                          function(i) ceiling(asv_counts[,..i] / 10^(correction[i]))) |> 
    bind_cols() 
  
}

# Calculate Shannon Diversity Index 
calc_shannon_index <- function(otu_abnd) {
  prp <- otu_abnd/ sum(otu_abnd)
  prp <- prp[prp > 0]
  H <- -sum(prp * log(prp))
  
}

# function to quickly check pcoa ordination for functional trait ordinations
plot_pcoa<- function(pcoa_obj,names=FALSE){
  
  # Step 3: Extract variance explained and scores
  varExp <- pcoa_obj$values$Eigenvalues / sum(pcoa_obj$values$Eigenvalues) * 100
  scores <- pcoa_obj$vectors
  
  plot(
    scores[, 1], scores[, 2],
    xlab = paste0("PCoA Axis 1 (", round(varExp[1], 2), "%)"),
    ylab = paste0("PCoA Axis 2 (", round(varExp[2], 2), "%)"),
    main = "PCoA of Categorical Traits (Gower Dissimilarity)",
    pch = 19
  )
  if(names) text(scores[, 1], scores[, 2], labels = attr(pcoa_obj$vectors , "dimnames")[[1]], pos = 4)
}


# loss functions ------------------------------------------------------------------------------

calc_rmse     <- function(preds, obs){sqrt(mean((preds-obs)^2))}
calc_mae      <- function(preds,obs){(mean(abs(preds-obs)))}
get_fold_loss <- function(mod,y,fit_data=NULL, stratified=FALSE){
  
  
  if(stratified){
    
    preds <- mod$pred[mod$folds[[1]]]
    fit_data$folds <- NA
    fit_data$folds[mod$folds[[1]]] <- 1
    fdat <- fit_data |> mutate(trapFold = as.numeric(interaction(trap_ID  , folds)))
    fInd <- unique(fdat$trapFold) |> na.omit()
    
    loss_df  <- map_dfr(fInd, function(x) {
      oInd <- which(fdat$trapFold==x)
      preds <- mod$pred[oInd]
      obs   <- y[oInd]
      
      fname <- unique(fdat$trap_ID[oInd])
      # Calculate loss 
     data.frame(
        fold = fname,
        RMSE = calc_rmse(preds = preds, obs = obs),
        MAE = calc_mae(preds = preds, obs = obs)
      )
     
  })
  }
  else{
  # map across folds
  loss_df  <- map_dfr(1:length(mod$folds), function(x) {
    preds <- mod$pred[mod$folds[[x]]]
    obs   <- y[mod$folds[[x]]]
    
    fname <- ifelse(is.null(names(mod$folds[x])),x,names(mod$folds[x]))
    # Calculate loss 
   data.frame(
      fold = fname,
      RMSE = calc_rmse(preds = preds, obs = obs),
      MAE = calc_mae(preds = preds, obs = obs)
    )
  }
  )
  }
  
  return(loss_df)
} # Function to get predictions from fitted model

