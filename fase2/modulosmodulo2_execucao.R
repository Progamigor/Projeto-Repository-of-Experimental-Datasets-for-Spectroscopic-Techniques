extrator_dados_brutos <- function(file_path, lista_var_modulo1) {
  
  if (lista_var_modulo1$extensao_ficheiro == "spc") {
    objeto_hyper <- tryCatch({
      withCallingHandlers(
        hyperSpec::read.spc(file_path),
        warning = function(w) invokeRestart("muffleWarning")
      )
    }, error = function(e) {
      cat("Aviso: ficheiro SPC com encoding incompatível.\n")
      return(NULL)
    })
    
    if (is.null(objeto_hyper)) return(NULL)
    
    # ISTO ESTAVA A FALTAR
    axis_x <- hyperSpec::wl(objeto_hyper)
    matrix_signal <- t(objeto_hyper$spc)
    
    if (is.null(colnames(matrix_signal))) {
      colnames(matrix_signal) <- paste0("Amostra_", 1:ncol(matrix_signal))
    }
    
    return(list(tipo = "spc", matriz = matrix_signal, eixo_x = axis_x))
    
  } else {
    
    raw_data <- read.table(
      file = file_path,
      skip = lista_var_modulo1$linhas_metadados_topo,
      sep = lista_var_modulo1$separador_colunas,
      dec = lista_var_modulo1$separador_decimal,
      header = lista_var_modulo1$cabecalho_presente,
      stringsAsFactors = FALSE
    )
    
    return(list(tipo = "texto", dados = raw_data))
  }
}
#-----------------------------------------------------------------------------~
extrair_metadados_globais <- function(file_path, lista_var_modulo1) {
  objeto_hyper <- tryCatch({
    withCallingHandlers(
      hyperSpec::read.spc(file_path),
      warning = function(w) invokeRestart("muffleWarning")
    )
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(objeto_hyper)) return(NULL)
  linhas_topo <- lista_var_modulo1$linhas_metadados_topo
  
  if (linhas_topo == 0) {
    return(NULL)
  }
  
  linhas <- readLines(file_path, n = linhas_topo, warn = FALSE)
  metadados_df <- data.frame(row.names = 1)
  contador_sem_chave <- 1
  
  for (linha in linhas) {
    linha <- trimws(linha)
    if (linha == "") next
    
    match_pos <- regexpr("[:=;\\t]", linha)
    
    if (match_pos > 0) {
      chave <- substr(linha, 1, match_pos - 1)
      valor <- substr(linha, match_pos + 1, nchar(linha))
      
      chave <- make.names(trimws(chave)) 
      valor <- trimws(valor)
      
      metadados_df[1, chave] <- valor
    } else {
      chave <- paste0("InfoExtra_", contador_sem_chave)
      metadados_df[1, chave] <- linha
      contador_sem_chave <- contador_sem_chave + 1
    }
  }
  return(metadados_df)
}

#------------------------------------------------------------------------------

tranformar_transposta <- function(objeto_dados, lista_var_modulo1) {
  if (objeto_dados$tipo == "spc") return(objeto_dados)
  
  if (lista_var_modulo1$matriz_transposta == TRUE) {
    objeto_dados$dados <- as.data.frame(t(objeto_dados$dados))
  }
  
  return(objeto_dados)
}
#------------------------------------------------------------------------------


definir_eixos <- function(objeto_dados) {

  if (objeto_dados$tipo == "spc") return(objeto_dados)
  
  axis_x <- as.numeric(objeto_dados$dados[, 1])
  matrix_signal <- as.matrix(objeto_dados$dados[, -1])
  
  return(list(tipo = "processado", eixo_x = axis_x, matriz = matrix_signal))
}
#------------------------------------------------------------------------------
inverter_eixo_x <- function(objeto_dados, lista_var_modulo1) {
  if (lista_var_modulo1$eixo_x_invertido) {
    objeto_dados$eixo_x <- rev(objeto_dados$eixo_x)
    objeto_dados$matriz <- objeto_dados$matriz[nrow(objeto_dados$matriz):1, , drop = FALSE]
  }
  return(objeto_dados)
}
#------------------------------------------------------------------------------
sanetizar_nomes <- function(objeto_dados) {
  nomes_limpos <- make.names(colnames(objeto_dados$matriz), unique = TRUE)
  colnames(objeto_dados$matriz) <- nomes_limpos
  rownames(objeto_dados$matriz) <- as.character(as.numeric(unlist(objeto_dados$eixo_x)))
  
  
  return(objeto_dados)
}

#------------------------------------------------------------------------------
empacotar_dados <- function(objeto_dados, tipo, metadados_globais = NULL) {
  if (is.null(colnames(objeto_dados$matriz))) {
    colnames(objeto_dados$matriz) <- paste0("amostra_", seq_len(ncol(objeto_dados$matriz)))
  }
  
  nomes_amostras <- colnames(objeto_dados$matriz)
  
  tabela_metadados <- data.frame(
    Amostra = nomes_amostras,
    row.names = nomes_amostras
  )
  

  matriz_bruta <- as.matrix(objeto_dados$matriz)
  matriz_pura <- matrix(as.numeric(matriz_bruta), nrow = nrow(matriz_bruta), ncol = ncol(matriz_bruta))
  colnames(matriz_pura) <- colnames(objeto_dados$matriz)
  rownames(matriz_pura) <- rownames(objeto_dados$matriz)
  
  nomes_amostras <- colnames(matriz_pura)

  tabela_metadados <- data.frame(
    Amostra = as.factor(nomes_amostras),
    row.names = nomes_amostras,
    stringsAsFactors = TRUE
  )
  
  if (!is.null(metadados_globais)) {
    for (coluna in colnames(metadados_globais)) {
      valor_limpo <- as.character(metadados_globais[1, coluna])
      tabela_metadados[[coluna]] <- as.factor(rep(valor_limpo, length(nomes_amostras)))
    }
    texto_descricao <- "Dataset harmonizado com metadados extraídos."
  } else {
    tabela_metadados$Grupo <- as.factor(rep("Amostra", length(nomes_amostras)))
    texto_descricao <- "Dataset harmonizado sem metadados extraídos."
  }
  
  eixo_x_limpo <- unname(as.numeric(as.character(unlist(objeto_dados$eixo_x))))
  
  tipo_str <- as.character(unlist(tipo))[1]
  nome_y <- ifelse(grepl("uvv", tipo_str), "Absorvância", "Intensidade")
  
  nome_x <- switch(tipo_str,
                   "uvv-spectra"   = "Comprimento de Onda (nm)",
                   "raman-spectra" = "Desvio de Raman (cm⁻¹)",
                   "ir-spectra"    = "Número de Onda (cm⁻¹)",
                   "Comprimento de Onda / Número de Onda"
  )

  validos <- !is.na(eixo_x_limpo) & complete.cases(matriz_pura)
  eixo_x_limpo <- eixo_x_limpo[validos]
  matriz_pura <- matriz_pura[validos, , drop = FALSE]

  dataset_final <- list(
    data = matriz_pura,
    metadata = tabela_metadados,
    x.labels = eixo_x_limpo,
    type = tipo_str,
    description = texto_descricao,
    x.label.text = nome_x,        
    y.label.text = nome_y         
  )
  
  class(dataset_final) <- "dataset"
  return(dataset_final)
}
#---------------------------------------------------------------
  
exportar_dataset <- function(dataset_specmine, file_path) {
  
  if (!dir.exists("output")) {
    cat("-> A criar a pasta 'output'\n")
    dir.create("output", recursive = TRUE)
  }
  
  nome_base <- tools::file_path_sans_ext(basename(file_path))
  
  nome_saida_dados <- paste0(nome_base, "_dados_harmonizados.csv")
  caminho_dados <- file.path("output", nome_saida_dados)
  
  tabela_dados <- data.frame(
    Eixo_X = dataset_specmine$x.labels,
    dataset_specmine$data,
    check.names = FALSE
  )
  
  write.csv(tabela_dados, file = caminho_dados, row.names = FALSE)
  cat(sprintf("-> Matriz numérica (CSV) guardada em: '%s'\n", caminho_dados))
  

  if (ncol(dataset_specmine$metadata) > 1) {
    nome_saida_meta <- paste0(nome_base, "_metadados.csv")
    caminho_meta <- file.path("output", nome_saida_meta)
    
 
    write.csv(dataset_specmine$metadata, file = caminho_meta, row.names = FALSE)
    cat(sprintf("-> Tabela de Metadados (CSV) guardada em: '%s'\n", caminho_meta))
  }

  nome_saida_rds <- paste0(nome_base, "_harmonizado.rds")
  caminho_rds <- file.path("output", nome_saida_rds)
  
  saveRDS(dataset_specmine, file = caminho_rds)
  cat(sprintf("-> O ficheiro nativo (RDS) guardado em: '%s'\n", caminho_rds))
  
  return(list(csv = caminho_dados, rds = caminho_rds))
}
