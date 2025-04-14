# 1. Carregar bibliotecas e definir variáveis --------------------------------
library(dplyr)
library(ggplot2)
library(gridExtra)

prop <- c("mn_o_pct", "nb_pct", "al2o3_pct", "si_o2_pct", "ti_o2_pct", "fe2o3_pct")
models = c("glmnet", "knn", "nnet", "rf", "svm") #"mlp"
scenario_folder <- 'modelo_proc_i/'
scenario_number <- 'I'
save_folder <- 'image_scenario_i/'

# 2. Leitura dos dados -------------------------------------------------------
load(file = 'data/dados_limpos.RData')

# 3. Inicialização dos Conjuntos de Teste para cada Prop deste design --------
testes <- list() # lista vazia para armazenar os conjuntos de teste.

for (var in prop) {

  testes[[var]] <- dffinal |> 
    dplyr::select(all_of(var), id, conj, aspect:length(dffinal)) %>%
    dplyr::filter(!is.na(.data[[var]])) |> 
    dplyr::filter(conj == "test") |> 
    dplyr::select(-c(id, conj))
  
  if (any(is.na(testes[[var]]))) {
    warning(paste("NAs found in the test data for", var))
  }
  
}
remove(dffinal, var)

# 4. Função load_and_predict ----------------------------------------------

# A função carrega modelos salvos em diretórios específicos e 
# faz previsões sobre novos dados (newdata ---> testes).

load_and_predict <- function(directory, pattern, newdata) {
  files <- list.files(directory, pattern = glob2rx(pattern))

  model_list <- lapply(files, function(file) {
    env <- new.env()
    load(file.path(directory, file), envir = env)
    env$model_fit
  })

  pred_matrix <- sapply(model_list, function(model) predict(model, newdata = newdata))
  mean_predictions <- rowMeans(pred_matrix)             # Media
  std_dev_predictions <- apply(pred_matrix, 1, sd)      # Desvio padrão
  upper_limit <- mean_predictions + std_dev_predictions # Limite superior
  lower_limit <- mean_predictions - std_dev_predictions # Limite inferior
  return(list(mean = mean_predictions, upper = upper_limit, lower = lower_limit))
}

predictions <- list() # Inicializar uma lista para armazenar os predict com teste

# Loop duplo para fazer previsões para cada variável e modelo
for (var in prop) {
  for (model in models) {
    directory <- paste0(scenario_folder, var, '/')
    pattern <- paste0('model_', model, '*.RData')
    predictions[[paste(var, model, sep = "_")]] <- load_and_predict(directory, pattern, testes[[var]])
  }
}


# 5. Função design_plot ------------------------------------------------------

design_plot <- function(data_plot, var, title) {
  x_label <- switch(var,
                    "fe2o3_pct" = expression("Observed Fe"[2]*"O"[3]*" wt(%)"),
                    "mn_o_pct"  = expression("Observed MnO wt(%)"),
                    "nb_pct"    = expression("Observed Nb wt(%)"),
                    "ti_o2_pct" = expression("Observed TiO"[2]*" wt(%)"),
                    "al2o3_pct" = expression("Observed Al"[2]*"O"[3]*" wt(%)"),
                    "si_o2_pct" = expression("Observed SiO"[2]*" wt(%)"),
                    "Unknown Variable"
  )
  y_label <- switch(var,
                    "fe2o3_pct" = expression("Predicted Fe"[2]*"O"[3]*" wt(%)"),
                    "mn_o_pct"  = expression("Predicted MnO wt(%)"),
                    "nb_pct"    = expression("Predicted Nb wt(%)"),
                    "ti_o2_pct" = expression("Predicted TiO"[2]*" wt(%)"),
                    "al2o3_pct" = expression("Predicted Al"[2]*"O"[3]*" wt(%)"),
                    "si_o2_pct" = expression("Predicted SiO"[2]*" wt(%)"),
                    "Unknown Variable"
  )

  ggplot(data_plot, aes(x = Observed, y = Predicted)) +
    geom_errorbar(aes(ymin = lower, ymax = upper, color = Predicted), linewidth = 1) +
    geom_point(aes(color = Predicted), size = 1.5) +
    #geom_abline(intercept = 0, slope = 1, color = "black", linetype = "solid", linewidth = 0.5) +
    
    #geom_smooth(method = "lm", color = "black", linetype = "solid", linewidth = 0.5) +
    geom_smooth(method = "loess", color = "black", linetype = "solid", linewidth = 0.5) +
    #geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "black", linetype = "solid", linewidth = 0.5) +
    
    scale_color_gradient2(low = "#527E87FF", mid ="#B8B69EFF", high = "#A37903FF",
                          midpoint = mean(range(data_plot$Predicted))) +
    labs(title = title,
         x = x_label,
         y = y_label) +
    theme_bw() +
    theme(
      text = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 8)
    ) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.5))
}


# 6. Função create_plots --------------------------------------------------

create_plots <- function(predictions, testes) {
  plots <- list()

  model_names <- c("glmnet" = "Generalized Linear Model via Elastic Net (glmnet)",
                   "knn" = "k-Nearest Neighbors (knn)",
                   #"mlp" = "Multilayer Perceptron (mlp)",
                   "nnet" = "Neural Network (nnet)",
                   "rf" = "Random Forest (rf)",
                   "svm" = "Support Vector Machine with Radial Basis Function Kernel (svm)")

  for (var in prop) {
    for (model in models) {
      data_plot <- data.frame(Observed = testes[[var]][[1]],
                              Predicted = predictions[[paste(var, model, sep = "_")]]$mean,
                              lower = predictions[[paste(var, model, sep = "_")]]$lower,
                              upper = predictions[[paste(var, model, sep = "_")]]$upper)

      plot_title <- paste(model_names[model])
      plots[[paste(var, model, sep = "_")]] <- design_plot(data_plot, var, plot_title)
    }
  }

  return(plots)
}

# Criar os gráficos
plots <- create_plots(predictions, testes)

for (plot_name in names(plots)) {
  print(plots[[plot_name]])
}


# 7. Plot combinado ----------------------------------------------------------
grid_plots_list <- list() # Lista para armazenar os grid_plots separados por prop

# prop-based
for (var in prop) {
  prop_plots <- lapply(models, function(model) {
    plot_name <- paste(var, model, sep = "_")
    plots[[plot_name]]
  })
  grid_plot <- grid.arrange(grobs = prop_plots, ncol = 1, top = paste("Scenario", scenario_number))
  grid_plots_list[[var]] <- grid_plot
}

for (var in prop) {
  filename <- file.path(save_folder, paste("Prop_Comp_", var, "_scen_", scenario_number, ".png", sep = ""))
  ggsave(filename, grid_plots_list[[var]], width = 6, height = 8, units = "in", dpi = 300)
}

#model_based
for (model in models) {
  model_plots <- lapply(prop, function(var) {
    plot_name <- paste(var, model, sep = "_")
    plots[[plot_name]]
  })
  grid_plot <- grid.arrange(grobs = model_plots, ncol = 1, top = paste("Scenario", scenario_number))
  grid_plots_list[[model]] <- grid_plot
}

for (model in models) {
  filename <- file.path(save_folder, paste("Model_Comp_", model, "_scen_", scenario_number, ".png", sep = ""))
  ggsave(filename, grid_plots_list[[model]], width = 6, height = 10, units = "in", dpi = 300)
}
