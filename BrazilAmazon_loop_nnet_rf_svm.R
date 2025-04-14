# ==============================================================================
# SCRIPT: Geospatial Data Preprocessing for Predictive Modeling
# Theresa R P Barbosa
# ==============================================================================

# =============================================================================
# CLEAN WORKSPACE
# =============================================================================
# Removes all objects from environment including hidden ones
# Best practice: Start scripts with clean workspace to avoid conflicts
rm(list = ls(all.names = TRUE))

# =============================================================================
# PACKAGE LOADING
# =============================================================================
# List of required packages for analysis
pacotes <- c("readr","dplyr","doParallel","caret",
             "terra", "sf", "sp", "tibble", "ggplot2", 
             "forcats", "corrplot", "ggcorrplot", "RSNNS", 
             "stats", "utils")

# Check if packages need installation
if(sum(as.numeric(!pacotes %in% installed.packages())) != 0){
  instalador <- pacotes[!pacotes %in% installed.packages()]
  # Installs missing packages (only first one due to break())
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = T)
    break()}
  sapply(pacotes, require, character = T) 
} else {
  sapply(pacotes, require, character = T) 
}
# Cleanup temporary variables
remove(i, instalador, pacotes)

# =============================================================================
# DATA WRANGLING
# =============================================================================

# Load training and validation datasets from CSV files
training <- readr::read_csv("data/training.csv")
validation <- readr::read_csv("data/validation.csv")

# Add dataset identifier column
training$conj <- rep("train", nrow(training))  # Marks all training rows
validation$conj <- rep("test", nrow(validation))  # Marks all validation rows

# Combine datasets into single dataframe
pt <- rbind(training, validation)  
# Remove original separate datasets to save memory
remove(training, validation)

# Remove specific columns (from "aspect" to "zonasalter")
pt <- pt |>  
  dplyr::select(-one_of(names(pt)[seq(which(names(pt) == "aspect"), 
                                      which(names(pt) == "zonasalter"), 
                                      by = 1)]))

# Convert to spatial object (sf)
sf_data <- sf::st_as_sf(pt, 
                        coords = c("point_x", "point_y"),  # Coordinate columns
                        crs = 32719)  # UTM zone 19S CRS

# Convert to terra's vector format
vect_data <- terra::vect(sf_data)

# Load raster stack
l = list.files('raster', glob2rx('*.tif'), full.names = TRUE)  # Finds all .tif files
st = terra::rast(l)  # Creates raster stack

# Extract raster values at point locations
df_extract <- terra::extract(st, sf_data, bind = TRUE)

# Convert back to regular dataframe
dfinicial <- as.data.frame(df_extract)

# Reattach original coordinates
dfinicial <- cbind(dfinicial, pt$point_x)
dfinicial <- cbind(dfinicial, pt$point_y)

# Rename coordinate columns
dfinicial <- dfinicial |> rename(point_x = `pt$point_x`)
dfinicial <- dfinicial |> rename(point_y = `pt$point_y`)

# Cleanup temporary objects
remove(df_extract, pt, sf_data, st, vect_data, l)

# Split into target variables (y) and features (x)
dfy = dfinicial |> dplyr::select(c(id:conj, point_x, point_y))  # Target + metadata
dfx = dfinicial |> dplyr::select(id, aspect:vv)  # Predictor variables

# =============================================================================
# FEATURE SELECTION
# =============================================================================

# Remove near-zero variance features
nzv = dfx |> nearZeroVar(names = TRUE)  # Identifies low-variance predictors
dfnz = dfx
if (length(nzv) > 0) {
  dfnz = dfx |> dplyr::select(-one_of(nzv))
  print(paste('Variable removed:', nzv))
} else {
  print('No near-zero variance features found')
}

# Remove highly correlated features (Spearman correlation)
limiar_correl = 0.95  # Correlation threshold
mcor = dfnz |> 
  dplyr::select_if(is.numeric) |>  # Only numeric variables
  select(-id) |>  # Exclude ID column
  scale() |>  # Standardize data
  cor(method = "spearman")  # Calculate correlation matrix

# Find correlated variables above threshold
vc = caret::findCorrelation(x = mcor, 
                            cutoff = limiar_correl, 
                            names = TRUE) 

# Visualize correlation matrix
corrplot(mcor, method = "circle", tl.col = "black", mar = c(0,0,5,0))

# Remove correlated features
dfcor = dfnz
if (length(vc) > 0) {
  dfcor = dfx |> dplyr::select(-one_of(vc))
  print(paste('Variable removed:', vc))
} else {
  print('No highly correlated features removed')
}

# =============================================================================
# FINAL DATA PREPARATION
# =============================================================================

# Merge target and selected features
dffinal = merge(dfy, dfcor, by = "id", all = TRUE)

# Save cleaned dataset
save(dffinal, file = './data/dados_limpos.RData')

# Final workspace cleanup
rm(list = ls(all.names = TRUE))

# //////////////////////////////////////////////////////////////////////////
# MODELING SECTION START --------------------------------------------------
# //////////////////////////////////////////////////////////////////////////

# Required packages -------------------------------------------------------
pacotes <- c("readr","dplyr","doParallel","caret",
             "terra", "sf", "sp", "tibble", "ggplot2", 
             "forcats", "corrplot", "ggcorrplot", 
             "stats", "utils", "randomForest", 
             "kernlab", "nnet", "glmnet", "Matrix",
             "RSNNS")

# Install missing packages and load all
if(sum(as.numeric(!pacotes %in% installed.packages())) != 0){
  instalador <- pacotes[!pacotes %in% installed.packages()]
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = T)
    break()} # Note: Only installs first missing package due to break()
  sapply(pacotes, require, character = T) 
} else {
  sapply(pacotes, require, character = T) 
}
# Cleanup temporary variables
remove(i, instalador, pacotes)

# Helper function for factor name padding ----------------------------------
pad3 <- function(s) {
  s = stringr::str_pad(s, 3, side = 'left', pad = '0') # Pads numbers with leading zeros (e.g., 5 -> "005")
  return(s)
}

# Load preprocessed data --------------------------------------------------
load(file = 'data/dados_limpos.RData') # Contains cleaned dataframe 'dffinal'
str(dffinal) # Check structure of loaded data

# Split into target (y) and predictor (x) variables -----------------------
dfy = dffinal |> dplyr::select(c(fe2o3_pct, mn_o_pct, nb_pct, ti_o2_pct, al2o3_pct)) # Target variables (chemical concentrations)
dfx = dffinal |> dplyr::select(id, conj, aspect:length(dffinal)) # Predictors (ID, dataset label, and all features from 'aspect' onward)

names(dfx) # Verify predictor variables
remove(dffinal) # Remove full dataset to save memory

# Model configuration -----------------------------------------------------
# Available models (commented examples):
#modelnames <- paste(names(getModelInfo())) # Lists all available caret models
#modelnames

# Hyperparameter lookup (commented examples):
#modelLookup("rf")       # Random Forest
#modelLookup("nnet")     # Neural Network
#modelLookup("svmRadial")# SVM with Radial Kernel
#modelLookup("glmnet")   # Elastic Net
#modelLookup("knn")      # K-Nearest Neighbors

# Model selection for RFE and training
modelos_rfe <- c('rf', 'nnet', 'svmRadial','glmnet', 'knn', 'mlp') # Models for Recursive Feature Elimination
modelos_train <- c('rf', 'nnet', 'svmRadial', 'glmnet', 'knn', 'mlp') # Models for final training

# Corresponding feature selection functions
funcs <- c('rfFuncs', 'caretFuncs', 'caretFuncs', 'caretFuncs', 'caretFuncs', 'caretFuncs') # Selection methods for each model

# Hyperparameter grids for tuning -----------------------------------------
rfGrid <- expand.grid(mtry = c(2, 4, 6)) # Random Forest: number of variables sampled at each split
nnetGrid <- expand.grid(size = c(5, 10, 15), decay = c(0.01, 0.001, 0.0001)) # Neural Net: hidden units and weight decay
svmRadialGrid <- expand.grid(sigma = c(0.1, 1, 10), C = c(1, 10, 100)) # SVM: kernel width and cost
glmnetGrid <- expand.grid(alpha = c(0, 0.5, 1), lambda = c(0.1, 1, 10)) # Elastic Net: mixing and regularization
knnGrid <- expand.grid(k = c(3, 5, 7)) # KNN: number of neighbors
mlpGrid <- expand.grid(size = c(5, 10, 15, 20)) # MLP: hidden layer neurons

# Combined grid list
grids <- list(
  rf = rfGrid,
  nnet = nnetGrid,
  svmRadial = svmRadialGrid,
  glmnet = glmnetGrid,
  knn = knnGrid,
  mlp = mlpGrid
)

# Experimental setup ------------------------------------------------------
nrep = 100 # Number of repetitions
set.seed(123) # Seed for reproducibility
vseed = sample(1:20000, nrep) # Vector of random seeds for each repetition

# Model and target dimensions
nmod = length(modelos_train) # Number of models (6)
ny = ncol(dfy) # Number of target variables (5)
nl = ny * nmod * nrep # Total model runs (5 targets × 6 models × 100 reps = 3000)

print(paste('Total model fits to perform =',nl))

# Initialize results dataframe --------------------------------------------
dfresult <- tibble(target = character(nl),         # Target variable name
                   repeticao = integer(nl),        # Repetition number
                   modelo = character(nl),         # Model name
                   rfe_numvar = integer(nl),       # Number of selected features
                   train_rmse = numeric(nl),       # Training RMSE
                   train_mae = numeric(nl),        # Training MAE
                   train_r2 = numeric(nl),         # Training R-squared
                   model_rmse = numeric(nl),       # Validation RMSE
                   model_mae = numeric(nl),        # Validation MAE
                   model_r2 = numeric(nl),         # Validation R-squared
                   model_bias = numeric(nl),       # Model bias
                   model_ccc = numeric(nl),        # Concordance correlation
                   null_model_rmse = numeric(nl),  # Null model RMSE
                   null_model_r2 = numeric(nl),    # Null model R-squared
                   null_model_mae = numeric(nl))   # Null model MAE

# Parallel processing setup -----------------------------------------------
nc = detectCores() # Detect available CPU cores
print(paste(nc, 'CPU cores available'))
cl <- makePSOCKcluster(12) # Create 12-worker cluster
doParallel::registerDoParallel(cl) # Register parallel backend

# //////////////////////////////////////////////////////////////////////////
# MAIN MODELING LOOP ------------------------------------------------------
# //////////////////////////////////////////////////////////////////////////

# Initialize counter for results storage
cont = 1  # Tracks position in results dataframe

# Triple nested loop structure:
# 1. Outer loop: Models (j)
# 2. Middle loop: Target variables (i)
# 3. Inner loop: Repetitions (k)
for (j in 1:length(modelos_train)) {
  for (i in 1:ny) {
    for (k in 1:nrep) {
      
      # Get current target variable name
      target = names(dfy)[i] 
      
      # Create directories for model outputs and prediction maps
      fp = file.path(getwd(), 'modelo', target)  # Model storage path
      if (dir.exists(fp) == FALSE) {
        dir.create(fp, recursive = TRUE)  # Creates nested directories if needed
      }
      fm = file.path(getwd(), 'mapa', target)  # Map storage path
      if (dir.exists(fm) == FALSE) {
        dir.create(fm, recursive = TRUE)
      }
      
      # Combine current target with predictors
      dfxy = data.frame(dfy[,i], dfx)
      names(dfxy)[1] = target  # Rename target column
      
      # Split into training/test sets
      treino = dfxy |> 
        dplyr::filter(!is.na(dfxy[1])) |>  # Remove NA targets
        dplyr::filter(conj == "train") |>  # Training set only
        dplyr::select(-c(id, conj))  # Remove metadata columns
      
      teste  = dfxy |> 
        dplyr::filter(!is.na(dfxy[1])) |> 
        dplyr::filter(conj == "test") |>  # Test set only
        dplyr::select(-c(id, conj))
      
      # Check for NA values
      any(is.na(treino))
      any(is.na(teste))
      
      # //////////////////////////////////////////////////////////////////////
      # RECURSIVE FEATURE ELIMINATION (RFE) --------------------------------
      # //////////////////////////////////////////////////////////////////////
      
      # Define feature subset sizes to test
      subsets <- c(2:25,30,35,40,46)  # Sequence of feature counts to evaluate
      form = as.formula(paste(target, '~ .'))  # Formula for modeling
      
      # Configure RFE control parameters
      ctrl_rfe <- rfeControl(functions = get(funcs[j]),  # Selection function for current model
                             method = "repeatedcv",      # Repeated cross-validation
                             repeats = 5,                # 5 repeats
                             number = 10,                # 10 folds
                             verbose = FALSE)            # No progress printing

      # Execute RFE
      rfe_fit <- rfe(form = form,
                     data = treino,
                     metric = 'MAE',         # Optimization metric
                     maximize = FALSE,       # Minimize MAE
                     method = modelos_rfe[j],# Current model type
                     rfeControl = ctrl_rfe,
                     size = subsets)        # Feature subset sizes to test
      
      # Display RFE results
      rfe_fit
      plot(rfe_fit)  # Visualize performance vs. feature count
      
      # Extract optimal feature set
      num_var = rfe_fit$bestSubset    # Best number of features
      var_sel = rfe_fit$optVariables  # Selected feature names
      
      # Filter training data to selected features
      vs = var_sel
      treino_sel = treino |> dplyr::select(all_of(target), one_of(vs))
      
      # Save RFE results
      fn = paste0(fp,'/','rfe_',modelos_rfe[j],'_',pad3(k),'.RData')
      save(rfe_fit, file = fn)
      
      # //////////////////////////////////////////////////////////////////////
      # MODEL TRAINING -----------------------------------------------------
      # //////////////////////////////////////////////////////////////////////
      
      # Configure training control
      ctrl <- trainControl(method = "repeatedcv", 
                           number = 10,  # 10-fold CV
                           repeats = 3)   # 3 repetitions
      
      # Train model with optimal features
      set.seed(123)  # Ensure reproducibility
      model_fit <- caret::train(form = form,
                                data = treino_sel,
                                method = modelos_train[j],  # Current algorithm
                                metric = 'MAE',            # Optimization metric
                                trControl = ctrl,
                                tuneGrid = grids[[modelos_train[j]]])  # Pre-defined parameter grid
      
      # Display model summary
      print(model_fit)
      
      # Generate test set predictions
      v = predict(model_fit, teste) 
      
      # Calculate test metrics
      vm = caret::postResample(v, teste[,1])  # Regression metrics (RMSE, R2, MAE)
      
      # Calculate additional metrics
      model_bias = mean(v - teste[,1])        # Prediction bias
      model_ccc = cor(v, teste[,1])          # Concordance correlation
      
      # Save trained model
      fn = paste0(fp,'/','model_',modelos_train[j],'_', pad3(k),'.RData')
      save(model_fit, file = fn)
      
      # //////////////////////////////////////////////////////////////////////
      # SPATIAL PREDICTION MAPPING ----------------------------------------
      # //////////////////////////////////////////////////////////////////////
      
      # Load raster data for selected features
      l = paste0('raster/', vs, '.tif')
      r = terra::rast(l)  # Create raster stack
      
      # Calculate null model metrics (mean prediction)
      vmn = rep(mean(teste[,1]), length(v))  # Vector of mean values
      vm_null = caret::postResample(vmn, teste[,1])  # Null model metrics
      
      # Prepare full raster data for prediction
      l = paste0('./raster/', var_sel, '.tif')
      r = terra::rast(l)
      
      # Convert raster to dataframe (excluding NA)
      dft = terra::as.data.frame(r, xy = TRUE) |> na.omit() 
      
      # Generate spatial predictions
      v = predict(model_fit, dft)  # Predict across entire raster
      
      # Create prediction dataframe with coordinates
      dat1 = data.frame(x = dft$x, y = dft$y, z = v)
      
      # Convert predictions to raster
      mapa = rast(dat1, type = "xyz", crs = crs(r))  # Maintain original CRS
      str_main = paste(modelos_train[j], '-', target, '_',  pad3(k))  # Plot title
      plot(mapa, main = str_main)  # Display map
      
      # Save prediction raster
      fn = paste0(fm,'/',target,'_', modelos_train[j],'_',pad3(k),'.tif')
      writeRaster(x = mapa, filename = fn, overwrite = TRUE )
      
      # //////////////////////////////////////////////////////////////////////
      # STORE RESULTS ------------------------------------------------------
      # //////////////////////////////////////////////////////////////////////
      
      # Populate results dataframe
      dfresult$target[cont] = target
      dfresult$repeticao[cont] = k
      dfresult$modelo[cont] = modelos_train[j]
      
      # RFE results
      dfresult$rfe_numvar[cont] = num_var
      
      # Training metrics
      dfresult$train_rmse[cont] = min(model_fit$results$RMSE)
      dfresult$train_mae[cont] = min(model_fit$results$MAE)
      dfresult$train_r2[cont] = max(model_fit$results$Rsquared)
      
      # Test metrics
      dfresult$model_rmse[cont] = vm[1]  # RMSE is first element
      dfresult$model_mae[cont] = vm[3]   # MAE is third element 
      dfresult$model_r2[cont] = vm[2]    # R2 is second element
      dfresult$model_bias[cont] = model_bias
      dfresult$model_ccc[cont] = model_ccc
      
      # Null model metrics
      dfresult$null_model_rmse[cont] = vm_null[1]
      dfresult$null_model_mae[cont] = vm_null[3]
      dfresult$null_model_r2[cont] = vm_null[2]
      
      # Increment counter and save progress
      cont = cont + 1
      readr::write_csv(dfresult, './modelo/resultados_samplying1.csv')  # Continuous saving
    } # End repetition loop (k)
  } # End target variable loop (i)
} # End model loop (j)

# Final print statements (debugging)
print(i)
print(j)
print(k)
