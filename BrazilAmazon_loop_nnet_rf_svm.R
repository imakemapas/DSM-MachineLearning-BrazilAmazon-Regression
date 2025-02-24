# Limpa environment -------------------------------------------------------
rm(list = ls(all.names = TRUE))


# Pacotes utilizados ------------------------------------------------------
pacotes <- c("readr","dplyr","doParallel","caret",
             "terra", "sf", "sp", "tibble", "ggplot2", 
             "forcats", "corrplot", "ggcorrplot", "RSNNS", 
             "stats", "utils")

if(sum(as.numeric(!pacotes %in% installed.packages())) != 0){
  instalador <- pacotes[!pacotes %in% installed.packages()]
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = T)
    break()}
  sapply(pacotes, require, character = T) 
} else {
  sapply(pacotes, require, character = T) 
}
remove(i, instalador, pacotes)


# Leitura dados -----------------------------------------------------------
training <- readr::read_csv("data/training.csv")
validation <- readr::read_csv("data/validation.csv")

training$conj <- rep("train", nrow(training))
validation$conj <- rep("test", nrow(validation))

pt <- rbind(training, validation) # Juntei os dois csv recebidos.
remove(training, validation)
names(pt)


# Data wrangling ----------------------------------------------------------

# Removi todas as colunas referentes a covariaveis para depois chamar elas de novo
pt <- pt  |>  
  dplyr::select(-one_of(names(pt)[seq(which(names(pt) == "aspect"), 
                                      which(names(pt) == "zonasalter"), 
                                      by = 1)]))

#names(pt)

# Convertendo o DataFrame em um objeto sf
sf_data <- sf::st_as_sf(pt, 
                        coords = c("point_x", "point_y"), 
                        crs = 32719)

# Convertendo sf para SpatVect (terra object)
vect_data <- terra::vect(sf_data)

# Stack raster
l = list.files('raster', glob2rx('*.tif'), full.names = TRUE)
st = terra::rast(l)

# Extract
df_extract <- terra::extract(st, sf_data, bind = TRUE)

# Retorna para data.frama
dfinicial <- as.data.frame(df_extract)

# Pega os valores de x e y de novo
dfinicial <- cbind(dfinicial, pt$point_x)
dfinicial <- cbind(dfinicial, pt$point_y)

dfinicial <- dfinicial |> rename(point_x = `pt$point_x`)
dfinicial <- dfinicial |> rename(point_y = `pt$point_y`)

names(dfinicial)
remove(df_extract, pt, sf_data, st, vect_data, l)

# Separa dataframe com variaveis alvo (y) e co-variaveis (x) --------------
dfy = dfinicial |> dplyr::select(c(id:conj, point_x, point_y))
dfx = dfinicial |> dplyr::select(id, aspect:vv) # peguei tudo do aspect adiante

# Detecta variáveis com variância zero ou quase zero ----------------------
nzv = dfx |> nearZeroVar(names = TRUE)
dfnz = dfx
if (length(nzv) > 0) {
  dfnz = dfx |> dplyr::select(-one_of(nzv))
  print(paste('varivavel', nzv, 'eliminada'))
} else (paste('Não há variáveis com variância zero ou quase zero'))

# "varivavel hill eliminada"
# "varivavel slope_idx eliminada"   
# "varivavel valley_idx eliminada"
# "varivavel hill_idx eliminada"
# "varivavel valley eliminada" 

# Detecta variaveis altamente correlacionadas -----------------------------

## Nao considera variaveis do tipo fator antes de calcular a matriz de correlacao
## Nao considera lat e long para calcular a matriz de correlacao
## Dados sao escalados para calcular a matriz de correlacao

limiar_correl = 0.95 
mcor = dfnz |> 
  dplyr::select_if(is.numeric) |>  # seleciona variaveis numericas
  select(-id) |> 
  scale() |>
  cor(method = "spearman")

vc = caret::findCorrelation(x = mcor, 
                            cutoff = limiar_correl, 
                            names = TRUE, ) 

# Cria a figura da matriz de correlação
corrplot(mcor, method = "circle", tl.col = "black", mar = c(0,0,5,0))

# Mostra as variáveis altamente correlacionadas
print(vc) 


dfcor = dfnz
if (length(vc) > 0) {
  dfcor = dfx |> dplyr::select(-one_of(vc))
  print(paste('varivavel', vc, 'eliminada'))
} else(paste('Nenhuma variável foi removida por alta correlação'))

# "varivavel ndvi eliminada"
# "varivavel savi eliminada"
# "varivavel B8 eliminada"
# "varivavel B7 eliminada"
# "varivavel tpi eliminada"
# "varivavel mrvbf eliminada"
# "varivavel mass_balan eliminada"
# "varivavel B5 eliminada"
# "varivavel terrain_ru eliminada"
# "varivavel real_surfa eliminada"
# "varivavel slope_degr eliminada"
# "varivavel clayminera eliminada"
# "varivavel FerrousSilic eliminada"
# "varivavel ferroussil eliminada"
# "varivavel ferricIron eliminada"
# "varivavel curv_longi eliminada"
# "varivavel gossan_ast eliminada"


dffinal = merge(dfy, dfcor, by = "id", all = TRUE)


# Salva dados limpos no diretorio para modelagem
save(dffinal, file = './data/dados_limpos.RData')
rm(list = ls(all.names = TRUE))

# //////////////////////////////////////////////////////////////////////////
# INCIO MODELAGEM ---------------------------------------------------------
# //////////////////////////////////////////////////////////////////////////

# Pacotes utilizados ------------------------------------------------------
pacotes <- c("readr","dplyr","doParallel","caret",
             "terra", "sf", "sp", "tibble", "ggplot2", 
             "forcats", "corrplot", "ggcorrplot", 
             "stats", "utils", "randomForest", 
             "kernlab", "nnet", "glmnet", "Matrix",
             "RSNNS")

if(sum(as.numeric(!pacotes %in% installed.packages())) != 0){
  instalador <- pacotes[!pacotes %in% installed.packages()]
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = T)
    break()}
  sapply(pacotes, require, character = T) 
} else {
  sapply(pacotes, require, character = T) 
}
remove(i, instalador, pacotes)

# Funcao para recuperar nomes de fatores --------------------------------------
pad3 <- function(s) {
  s = stringr::str_pad(s, 3, side = 'left', pad = '0')
  return(s)
}

# Leitura dados limpos  ---------------------------------------------------
load(file = 'data/dados_limpos.RData')
str(dffinal)

# separa dataframe com variaveis alvo (y) e co-variaveis (x) --------------
dfy = dffinal |> dplyr::select(c(fe2o3_pct, mn_o_pct, nb_pct, ti_o2_pct, al2o3_pct)) #, si_o2_pct, id, conj
dfx = dffinal |> dplyr::select(id, conj, aspect:length(dffinal)) # peguei tudo do aspect adiante

names(dfx)

remove(dffinal)

# Algoritmos - ajustes ----------------------------------------------------
#modelnames <- paste(names(getModelInfo()))   # consulta todos os modelos
#modelnames

#consulte os hiperparâmetros
#modelLookup("rf")
#modelLookup("nnet")
#modelLookup("svmRadial")
#modelLookup("glmnet")
#modelLookup("knn")

# modelos
modelos_rfe <- c('rf', 'nnet', 'svmRadial','glmnet', 'knn', 'mlp')
modelos_train <- c('rf', 'nnet', 'svmRadial', 'glmnet', 'knn', 'mlp')

# pacotes acessorios
funcs <- c('rfFuncs', 'caretFuncs', 'caretFuncs', 'caretFuncs', 'caretFuncs', 'caretFuncs')

#?rfFuncs
#?caretFuncs

# hiperparametros
# GridSearch
rfGrid <- expand.grid(mtry = c(2, 4, 6))
nnetGrid <- expand.grid(size = c(5, 10, 15), decay = c(0.01, 0.001, 0.0001))
svmRadialGrid <- expand.grid(sigma = c(0.1, 1, 10), C = c(1, 10, 100))
glmnetGrid <- expand.grid(alpha = c(0, 0.5, 1), lambda = c(0.1, 1, 10))
knnGrid <- expand.grid(k = c(3, 5, 7))
mlpGrid <- expand.grid(size = c(5, 10, 15, 20))

grids <- list(
  rf = rfGrid,
  nnet = nnetGrid,
  svmRadial = svmRadialGrid,
  glmnet = glmnetGrid,
  knn = knnGrid,
  mlp = mlpGrid
)


# # Obtendo os melhores hiperparâmetros para cada modelo
# bestHyperparameters <- list(
#   rf = rfTuned$bestTune,
#   nnet = nnetTuned$bestTune,
#   svmRadial = svmRadialTuned$bestTune,
#   glmnet = glmnetTuned$bestTune,
#   knn = knnTuned$bestTune
# )

# número de repetições
nrep = 2

# aleatoriedade 
set.seed(123)
vseed = sample(1:20000, nrep)

# número de modelos
nmod = length(modelos_train)

# número de variaveis alvo
ny = ncol(dfy) #- 2

# numero total execucoes
nl = ny * nmod * nrep 

print(paste('numero de ajustes de modelo =',nl))

dfresult <- tibble(target = character(nl),         #1
                   repeticao = integer(nl),        #2
                   modelo = character(nl),         #3
                   rfe_numvar = integer(nl),       #4
                   train_rmse = numeric(nl),       #5
                   train_mae = numeric(nl),        #6
                   train_r2 = numeric(nl),         #7
                   model_rmse = numeric(nl),       #8
                   model_mae = numeric(nl),        #9
                   model_r2 = numeric(nl),         #10
                   model_bias = numeric(nl),       #11
                   model_ccc = numeric(nl),        #12
                   null_model_rmse = numeric(nl),  #13
                   null_model_r2 = numeric(nl),    #14
                   null_model_mae = numeric(nl))   #17

# Paralelismo
nc = detectCores()
print(paste(nc, 'núcleos de cpu disponiveis'))
cl <- makePSOCKcluster(12)
doParallel::registerDoParallel(cl)

# loop for ----------------------------------------------------------------
j = 1; i = 1; k = 1
cont = 1

for (j in 1:length(modelos_train)) {
  for (i in 1:ny) {
    for (k in 1:nrep) {
      
      # variável a ser modelada 
      target = names(dfy)[i] 
      
      # cria uma pasta para cada variável y a ser modelada e para armazenar os modelos e estatísticas criadas
      fp = file.path(getwd(), 'modelo', target)
      if (dir.exists(fp) == FALSE) {
        dir.create(fp, recursive = TRUE)
      }
      fm = file.path(getwd(), 'mapa', target)
      if (dir.exists(fm) == FALSE) {
        dir.create(fm, recursive = TRUE)
      }
      
      # cria dataframe com a variável y a ser analisada e as co-variáveis
      dfxy = data.frame(dfy[,i], dfx)
      names(dfxy)[1] = target
      
      # Separa treino e teste para as etapas da modelagem
      treino = dfxy |> dplyr::filter(!is.na(dfxy[1])) |> dplyr::filter(conj == "train") |> dplyr::select(-c(id, conj)) 
      teste  = dfxy |> dplyr::filter(!is.na(dfxy[1])) |> dplyr::filter(conj == "test") |> dplyr::select(-c(id, conj))
      
      any(is.na(treino))
      any(is.na(teste))
      
      # RFE - Seleção de variáveis ----------------------------------------------
      
      # número de melhores variáveis que serao testadas
      #subsets <- c(2:25,30,35,40,46)  # número de melhores variáveis que serão testadas
      subsets <- c(10)  
      form = as.formula(paste(target, '~ .'))
      
      ctrl_rfe <- rfeControl(functions = get(funcs[j]),
                             method = "repeatedcv",
                             repeats = 5,
                             number = 10,
                             verbose = FALSE)

      rfe_fit <- rfe(form = form,
                     data = treino,
                     metric = 'MAE', 
                     maximize = FALSE,
                     method = modelos_rfe[j],
                     rfeControl = ctrl_rfe,
                     size = subsets)
      
      rfe_fit
      plot(rfe_fit)  
      num_var = rfe_fit$bestSubset
      var_sel = rfe_fit$optVariables
      
      vs = var_sel
      treino_sel = treino |> dplyr::select(all_of(target), one_of(vs))
      
      fn = paste0(fp,'/','rfe_',modelos_rfe[j],'_',pad3(k),'.RData')
      save(rfe_fit, file = fn)
      
      
      # Ajuste do Modelo --------------------------------------------------------
      ctrl <- trainControl(method = "repeatedcv", 
                           number = 10, # folds
                           repeats = 3) # repetições
      
      
      #grid_atual <- grids[[modelos_train[j]]]
      
      set.seed(123)
      model_fit <- caret::train(form = form,
                                data = treino_sel,
                                method = modelos_train[j],
                                metric = 'MAE',
                                trControl = ctrl,
                                #tuneLength = tuneLength[i])
                                tuneGrid = grids[[modelos_train[j]]]) # Seleciona o grid de ajuste com base no modelo atual
      
      
      print(model_fit)
      #table(model_fit$trainingData$.outcome)
      
      # prediz valores da variável target
      v = predict(model_fit, teste) 
      
      # métricas de ajuste - TESTE
      vm = caret::postResample(v, teste[,1])   # regress or two class
      #vm = caret::confusionMatrix(v, teste[,1]) # multi class
      
      #calculo do bias
      model_bias = mean(v - teste[,1])
      
      #calculo do coefficient correlation concordance
      model_ccc = cor(v, teste[,1])
      
      ### salva modelo ajustado
      fn = paste0(fp,'/','model_',modelos_train[j],'_', pad3(k),'.RData')
      save(model_fit, file = fn)
      
      l = paste0('raster/', vs, '.tif')
      #r = raster::stack(l) #
      r = terra::rast(l)
      #plot(r)  
      names(r)
      #raster::crs(r) = 32719
   
      
      ## calcula modelo nulo
      vmn = rep(mean(teste[,1]), length(v))
      vm_null = caret::postResample(vmn, teste[,1])
      model_bias_null = mean(vmn - teste[,1])
      model_ccc_null = cor(vmn, teste[,1])
      
      
      l = paste0('./raster/', var_sel, '.tif')
      r = terra::rast(l)
      names(r)
      
      dft = terra::as.data.frame(r, xy = TRUE) |> na.omit() # transforma pra dataframe. Tem as 5 variaveis que nos selecionamos , #dft é o raster inteiro
      #dft = dft |> janitor::clean_names() #limpar os nomes das colunas de um conjunto de dados, tornando-os mais amigáveis e padronizados.
      
      # USARQUANDO TIVER VARIAVEIS FACTOR
      # if (length(vf > 0)) {
      #   dft = dft |> mutate_at(vars(all_of(vf)), factor)
      # }
      
      # Previsoes
      v = predict(model_fit, dft)
      
      # Cria um data frame com as previsões (v) e as coordenadas x e y
      dat1 = data.frame(x = dft$x, y = dft$y, z = v)
      
      ## converte dataframe para raster
      mapa = rast(dat1, type = "xyz", crs = crs(r))
      str_main = paste(modelos_train[j], '-', target, '_',  pad3(k))
      plot(mapa, main = str_main)
      
      
      #ext(mapa) = ext(r)
      #mapa = terra::resample(x = mapa , y = terra::rast(r))
      #plot(mapa, main = str_main)
      
      
      fn = paste0(fm,'/',target,'_', modelos_train[j],'_',pad3(k),'.tif')
      writeRaster(x = mapa, filename = fn, overwrite = TRUE )
      
      # RODADA
      dfresult$target[cont] = target
      dfresult$repeticao[cont] = k
      dfresult$modelo[cont] = modelos_train[j]
      #RFE
      dfresult$rfe_numvar[cont] = num_var
      #TREINAMENTO
      dfresult$train_rmse[cont] = min(model_fit$results$RMSE)
      dfresult$train_mae[cont] = min(model_fit$results$MAE)
      dfresult$train_r2[cont] = max(model_fit$results$Rsquared)
      #TESTE
      dfresult$model_rmse[cont] = vm[1] # quando estiver com duvida, faz um 'print(vm) e vc vai ver q o vm nesse caso tem tres argumentos e o numero é a ordem que eles aparecem
      dfresult$model_mae[cont] = vm[3] 
      dfresult$model_r2[cont] = vm[2]
      dfresult$model_bias[cont] = model_bias # Aqui nao precisa, vc calculou separado, fora do vm
      dfresult$model_ccc[cont] = model_ccc
      #NULO
      dfresult$null_model_rmse[cont] = vm_null[1]
      dfresult$null_model_mae[cont] = vm_null[3]
      dfresult$null_model_r2[cont] = vm_null[2]
      
      cont = cont + 1
      
      readr::write_csv(dfresult, './modelo/resultados_samplying1.csv')
    } # for i
  } # for j
} #for k

print(i)
print(j)
print(k)

#readr::write_csv(dfresult, './modelo/resultadossamplying1.csv')
stopCluster(cl) ## desliga paralelismo e libera memória
