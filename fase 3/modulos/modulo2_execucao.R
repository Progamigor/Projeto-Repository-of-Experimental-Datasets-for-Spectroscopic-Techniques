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
    
    axis_x <- hyperSpec::wl(objeto_hyper)
    matrix_signal <- t(objeto_hyper$spc)
    
    if (is.null(colnames(matrix_signal))) {
      colnames(matrix_signal) <- paste0("Amostra_", 1:ncol(matrix_signal))
    }
    
    return(list(tipo = "spc", matriz = matrix_signal, eixo_x = axis_x))
    
  } else {
    

    linhas <- readLines(file_path, warn = FALSE)
    

    linha_cabecalho <- NULL
    if (lista_var_modulo1$cabecalho_presente && lista_var_modulo1$linhas_metadados_topo > 0) {
      linha_cabecalho <- linhas[lista_var_modulo1$linhas_metadados_topo]
    }
    

    linhas_matriz <- sapply(linhas, function(linha) {
      linha_limpa <- trimws(linha)

      if (grepl("^[-+]?[0-9]", linha_limpa) && grepl(lista_var_modulo1$separador_colunas, linha_limpa, fixed = TRUE)) {
        return(linha_limpa)
      }
      return(NA)
    })
    

    linhas_matriz <- linhas_matriz[!is.na(linhas_matriz)]
    

    if (!is.null(linha_cabecalho)) {
      bloco_limpo <- c(linha_cabecalho, linhas_matriz)
    } else {
      bloco_limpo <- linhas_matriz
    }
    
    texto_limpo <- paste(bloco_limpo, collapse = "\n")
    

    raw_data <- read.table(
      text = texto_limpo,
      sep = lista_var_modulo1$separador_colunas,
      dec = lista_var_modulo1$separador_decimal,
      header = lista_var_modulo1$cabecalho_presente,
      stringsAsFactors = FALSE,
      fill = TRUE
    )
    
    return(list(tipo = "texto", dados = raw_data))
  }
}
#-----------------------------------------------------------------------------~
extrair_metadados_globais <- function(file_path, lista_var_modulo1) {
  
  if (lista_var_modulo1$extensao_ficheiro == "spc") {

    objeto_hyper <- tryCatch({
      withCallingHandlers(
        hyperSpec::read.spc(file_path),
        warning = function(w) invokeRestart("muffleWarning")
      )
    }, error = function(e) return(NULL))
    
    if (is.null(objeto_hyper)) return(NULL)
    
    metadados_df <- as.data.frame(objeto_hyper@data)
    metadados_df$spc <- NULL # Remove a matriz espectral dos metadados
    
    if (ncol(metadados_df) > 0) return(metadados_df[1, , drop = FALSE])
    return(NULL)
  }
  
  
  linhas <- readLines(file_path, warn = FALSE)
  metadados_df <- data.frame(row.names = 1)
  contador_sem_chave <- 1
  bloco_atual <- 1                        # <-- novo
  ultimo_foi_numerico <- FALSE            # <-- novo
  
  for (linha in linhas) {
    linha_limpa <- trimws(linha)
    if (linha_limpa == "") next
    
    is_numeric_row <- grepl("^[-+]?[0-9]", linha_limpa) && 
      !grepl("[a-zA-Z]{4,}", linha_limpa)
    
    if (is_numeric_row) {
      ultimo_foi_numerico <- TRUE         # <-- novo
      next
    }
    
    
    if (ultimo_foi_numerico) {            # <-- novo bloco
      bloco_atual <- bloco_atual + 1
      ultimo_foi_numerico <- FALSE
    }
    
    match_pos <- regexpr("[:=;]", linha_limpa)
    
    if (match_pos > 0) {
      chave <- paste0("Bloco", bloco_atual, "_",          # <-- prefixo
                      make.names(trimws(substr(linha_limpa, 1, match_pos - 1))))
      valor <- trimws(substr(linha_limpa, match_pos + 1, nchar(linha_limpa)))
      metadados_df[1, chave] <- valor
    } else {
      chave <- paste0("Bloco", bloco_atual, "InfoExtra", # <-- prefixo
                      contador_sem_chave)
      metadados_df[1, chave] <- linha_limpa
      contador_sem_chave <- contador_sem_chave + 1
    }
  }
  
  if (ncol(metadados_df) == 0) return(NULL)
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
  
  nome_saida_dados <- paste0(nome_base, "_harmonizado.csv")
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