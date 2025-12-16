# libraries ----------------------------------------------------------------
  
library(tidyverse)
library(xgboost)
set.seed(127)  
# load feature & target matrices --------------------------------------------------------------
ft_M <- readRDS('../../../data/tidydata/full_M.rds')$madagascar$LCBD |> drop_na() 

### Model fitting code ###  
 # data set up ---------------------------------------------------------------------------------

# Set up X matrix, response family, & offsets

drop_cols <- c(colnames(ft_M)[2] , 'trap_ID' , 'week_year')
X         <- ft_M |> dplyr::select(-all_of(drop_cols)) |> as.matrix()
y         <- ft_M[,2]
dtrain    <- xgb.DMatrix(data = X , label = y) 

obj_response <-ifelse(colnames(ft_M)[2] == 'nOTU' , "count:poisson", "reg:squarederror")
 


# hyper parameter tuning ----------------------------------------------------------------------
cat("Tuning learning rate\n")

# Grid search for hyper-parameter tuning for eta first
eta_grid <- expand.grid(eta     = seq(0.01,0.3,l=20) , nrounds = c(300,500,1000,1500,2000))

# ---- setup
cv_resList <- vector(mode = 'list' , length = nrow(eta_grid)) # List to store

# ---- Loop over options 
# ----
for (i in 1:nrow(eta_grid)) {{
  # Make a new list of parameters for each combination in our grid search

  #cat(paste0("Running iteration " , i , " of " , nrow(eta_grid), "\n"))

  params <- list(
   objective         = obj_response,
    eta               = eta_grid$eta[i])
  
  # CV across options  
  cv_resList[[i]] <- xgb.cv(
    params                = params,
    data                  = dtrain,
    nrounds               = eta_grid$nrounds[i],
    nfold                 = 5,
    early_stopping_rounds = 10,
    metrics               = list('rmse'),
    verbose               = FALSE,
    tree_method='hist')
  
}}

# get optimal eta
eta_fit       <- lapply(1:length(cv_resList) , function(x) cv_resList[[x]]$evaluation_log) |> 
  bind_rows(.id='hypar_comb')
opt_eta_ind   <- as.numeric(eta_fit[which.min(eta_fit$test_rmse_mean),]$hypar_comb)
opt_eta       <- eta_grid[opt_eta_ind,]


# grid search for the rest of the parameters --------------------------------------------------
cat("Tuning remaining hyperparameters\n")

hypar_grid <- expand.grid(
  max_depth        = c(3,5,7,9,11),
  min_child_weight = c(1,3,5),
  eta              = opt_eta$eta,
  nrounds          = opt_eta$nrounds,
  subsample        = c(0.8 , 1.0),
  colsample_bytree = c(0.35,0.5, 0.75),
  lambda           = c(0.5, 1.5, 3.5),         # L2 regularization (Ridge regression)
  alpha            = c(0.5, 1.5, 3.5)          # L1 regularization (Lasso regression)
)

# ---- Loop over options 

cv_resList <- vector(mode = 'list' , length = nrow(hypar_grid)) # List to store
for (i in 1:nrow(hypar_grid)) {{
  # Make a new list of parameters for each combination in our grid search
  #cat(paste0("Running iteration " , i ," of ", nrow(hypar_grid), "\n"))
  params <- list(
   objective         = obj_response,
    max_depth         = hypar_grid$max_depth[i],
    eta               = hypar_grid$eta[i],
    nrounds           = hypar_grid$nrounds[i],
    colsample_bytree  = hypar_grid$colsample_bytree[i],
    lambda            = hypar_grid$lambda[i],
    alpha             = hypar_grid$alpha[i]
  )
  
  # CV across options  
  cv_resList[[i]] <- xgb.cv(
    params                = params,
    data                  = dtrain,
    nrounds               = hypar_grid$nrounds[i],
    nfold                 = 5,
    early_stopping_rounds = 10,
    metrics               = list('rmse'),
    tree_method = 'hist',
    verbose               = FALSE)
  
}}

hypar_fit <- lapply(1:length(cv_resList) , function(x) cv_resList[[x]]$evaluation_log) |> bind_rows(.id='hypar_comb')
opt_h_ind <- as.numeric(hypar_fit[which.min(hypar_fit$test_rmse_mean),]$hypar_comb)

opt_hypars <- list(
 objective         = obj_response,
  max_depth         = hypar_grid$max_depth[opt_h_ind],
  eta               = hypar_grid$eta[opt_h_ind],
  nrounds           = hypar_grid$nrounds[opt_h_ind],
  subsample         = hypar_grid$subsample[opt_h_ind],
  colsample_bytree  = hypar_grid$colsample_bytree[opt_h_ind],
  lambda            = hypar_grid$lambda[opt_h_ind],
  alpha             = hypar_grid$alpha[opt_h_ind]
)

# K-fold ------------------------------------------------------------------------------------
cat("K-fold cross validation \n")

tmod_kfold <- xgb.cv(
  params                = opt_hypars,
  data                  = dtrain,
  nrounds               = 5e3,
  nfold                 = 5,
  early_stopping_rounds = 10,
  metrics               = list('rmse','mae'), 
  tree_method = 'hist',
  prediction=TRUE,
  verbose=FALSE)

# stratified  test-train split -----------------------------------------------------------------
cat("Sample stratified test-train cross validation \n")

# Setup leave-one-site out folds
loo_id      <- ft_M |> mutate(loo_id = 1:n()) |> pull(loo_id)
samp        <- ft_M |> mutate(test_train = ifelse(runif(n()) < 0.3, 1, 2),.by=trap_ID)
strat_folds <- split(loo_id , samp$test_train)              

tmod_strat <- xgb.cv(
  params                = opt_hypars,
  data                  = dtrain,
  nrounds               = 5e3,
  folds                 = strat_folds,
  early_stopping_rounds = 10,
  metrics               = list('rmse','mae'), 
  tree_method = 'hist',
  prediction=TRUE,
  verbose=FALSE)
  
  
# Leave one site out ---------------------------------------------------------------------------
cat("LOSO cross validation \n")

# Setup leave-one-site out folds
loo_id <- ft_M |> mutate(loo_id = 1:n()) |> pull(loo_id) 
loo_folds <- split(loo_id , ft_M$trap_ID)              

tmod_loso <- xgb.cv(
  params                = opt_hypars,
  data                  = dtrain,
  nrounds               = 5e3,
  folds                 = loo_folds,
  early_stopping_rounds = 10,
  metrics               = list('rmse','mae'), 
  tree_method = 'hist',
  prediction=TRUE,
  verbose=FALSE)

# Normal model fit ---------------------------------------------------------------------------
cat("Final model fit")
tmod_fit <- xgboost(
  params                = opt_hypars,
  data                  = dtrain,
  nrounds               = 5e3,
  early_stopping_rounds = 10,
  tree_method = 'hist',
  verbose=FALSE)
  
# Save fit data ---------------------------------------------------------------------------
      saveList <- list(data       = list(X=X,y=y,full=ft_M),                                # Data
                       model_fits = list(full  = tmod_fit , 
                                         kfold = tmod_kfold, 
                                         strat = tmod_strat,
                                         loso  = tmod_loso ), # CV Results
                       hypars     = list(opt_hypars))                                           # Hyperparameters
      
  
#########################

# save ----------------------------------------------------------------------------------------
  
 saveRDS(saveList , '../../../results/results_madagascar_LCBD.rds')
  
