#------------------------------
variaveis_sinalizadoras <- function() {
  lista_var_sinal <- list(
    formato_suportado = FALSE,
    extensao_ficheiro = "",
    separador_colunas = "",
    separador_decimal = ".",         
    matriz_transposta = FALSE,
    cabecalho_presente = FALSE,
    escala_transmitancia = FALSE,
    eixo_x_invertido = FALSE,
    linhas_metadados_topo = 0 
  )
  return(lista_var_sinal)
}



#--------------------------------------------------------------------------------


escolha_tipo_dados <- function() {
  cat("Qual o tipo de dados espectroscópicos?\n")
  cat("  1 - Raman\n")
  cat("  2 - Infrared\n")
  cat("  3 - UV-Vis\n")
  cat("Introduza o número: ")
  
  opcao <- as.integer(readLines(con = stdin(), n = 1))
  
  tipos <- list(
    "1" = "raman-spectra",
    "2" = "ir-spectra",
    "3" = "uvv-spectra"
  )
  
  if (!as.character(opcao) %in% names(tipos)) {
    stop("Opção inválida. Introduza 1, 2 ou 3.")
  }
  
  tipo_selecionado <- tipos[[as.character(opcao)]]
  cat(sprintf("Tipo selecionado: %s\n", tipo_selecionado))
  return(tipo_selecionado)
}



#-----------------------------------------------------------------------------
escolha_de_ficheiro <- function() {
  cat("Por favor, selecione o ficheiro.\n")
  fi<-file.choose()
  return (fi)
}


#------------------------------

detetar_formato <- function(file_path, lista_var_sinal) {
  cat("A verificar extensão do ficheiro...\n")
  
  extensao_formato <- tolower(tools::file_ext(file_path))
  
  formatos_permitidos <- c("csv", "txt", "tsv","spc")
  
  if (!(extensao_formato %in% formatos_permitidos)) {
    cat(sprintf("Erro: O formato '.%s' não é suportado \n", extensao_formato))
    stop("Execução interrompida: formato inválido")
  }
  
  cat(sprintf("Formato '.%s' aprovado para leitura.\n", extensao_formato))
  lista_var_sinal$formato_suportado <- TRUE
  if (extensao_formato == "spc") {
    if (!requireNamespace("hyperSpec", quietly = TRUE)) {
      cat("-> A instalar o pacote essencial: 'hyperSpec'...\n")
      install.packages("hyperSpec", quiet = TRUE)
      cat("O pacote 'hyperSpec' foi instalado.\n")
    }
  }
  lista_var_sinal$formato_suportado <- TRUE
  lista_var_sinal$extensao_ficheiro <- extensao_formato
  
  return(lista_var_sinal)
}
#---------------------------------------------------------------------------------------
detetar_delimitador <- function(file_path, lista_var_sinal) {
  cat("A analisar delimitador e separador decimal...\n")
  
  if (lista_var_sinal$extensao_ficheiro == "spc") {
    #Ficheiro spc detetado. Passo ignorado
    return(lista_var_sinal)
  }
  
  amostra_linhas <- readLines(file_path, n = 5, warn = FALSE)
  texto_amostra <- paste(amostra_linhas, collapse = " ")
  
  num_pontos_virgula <- lengths(regmatches(texto_amostra, gregexpr(";", texto_amostra)))
  num_tabs <- lengths(regmatches(texto_amostra, gregexpr("\t", texto_amostra)))
  num_virgulas <- lengths(regmatches(texto_amostra, gregexpr(",", texto_amostra)))
  
  if (num_pontos_virgula > 0) {
    separador_colunas <- ";"
    separador_dec <- "," 
    
  } else if (num_tabs > 0) {
    separador_colunas <- "\t"
    # Procuramos o padrão de um número com vírgula
    if (grepl("[0-9],[0-9]", texto_amostra)) {
      separador_dec <- ","
    } else {
      separador_dec <- "."
    }
    
  } else if (num_virgulas > 0) {
    separador_colunas <- ","
    separador_dec <- "."
    
  } else {
    separador_colunas <- ","
    separador_dec <- "."
  }
  
  if (separador_colunas == "\t") {
    cat(sprintf("Separador_de_colunas: TAB | Separador_Decimal: '%s'\n", separador_dec))
  } else {
    cat(sprintf("separador_colunas: '%s' | Separador_Decimal: '%s'\n", separador_colunas, separador_dec))
  }
  
  lista_var_sinal$separador_colunas <- separador_colunas
  lista_var_sinal$separador_decimal <- separador_dec
  
  return(lista_var_sinal)
}




#-----------------------------------------------------------------------------------------
detetar_linhas_metadados <- function(file_path, lista_var_sinal) {
  linhas_raw <- readLines(file_path, n = 50)
  linhas_hash <- sum(grepl("^##", linhas_raw))
  
  if (linhas_hash > 0) {
    linha_seguinte <- linhas_raw[linhas_hash + 1]
    linha_vazia <- ifelse(grepl("^\\s*$", linha_seguinte), 1, 0)
    lista_var_sinal$linhas_metadados_topo <- linhas_hash + linha_vazia
    cat(sprintf("Padrão ## detetado: %d linhas de metadados.\n", lista_var_sinal$linhas_metadados_topo))
    return(lista_var_sinal)
  }
  
  delimitador_atual <- lista_var_sinal$separador_colunas
  
  amostra_linhas <- readLines(file_path, n = 150, warn = FALSE)
  

  contagem_colunas <- sapply(amostra_linhas, function(linha) {
    length(strsplit(linha, split = delimitador_atual, fixed = TRUE)[[1]])
  })
  

  tabela_frequencias <- table(contagem_colunas)
  colunas_reais <- as.numeric(names(tabela_frequencias)[which.max(tabela_frequencias)])
  

  primeira_linha_valida <- min(which(contagem_colunas == colunas_reais))
  

  linhas_de_metadados <- primeira_linha_valida - 1
  
  if (linhas_de_metadados > 0) {
    cat(sprintf("Foram detetadas %d linhas de metadados antes da tabela principal.\n", linhas_de_metadados))
  } else {
    cat("O ficheiro foi analisado quanto à presença de metadados.\n")
  }
  
  lista_var_sinal$linhas_metadados_topo <- linhas_de_metadados
  
  return(lista_var_sinal)
}

#------------------------------------------------------------------------------------

detetar_transposicao <- function(file_path, lista_var_sinal) {

  
  if (lista_var_sinal$extensao_ficheiro == "spc") {
    # Ficheiro spc detetado. Passo ignorado
    return(lista_var_sinal)
  }
  
  # Lemos uma amostra geométrica de 50 linhas
  tabela_amostra <- tryCatch({
    read.table(file_path, 
               sep = lista_var_sinal$separador_colunas, 
               skip = lista_var_sinal$linhas_metadados_topo, 
               nrows = 50, 
               header = FALSE, 
               stringsAsFactors = FALSE,
               fill = TRUE)
  }, error = function(e) return(NULL))
  
  if (is.null(tabela_amostra) || ncol(tabela_amostra) < 2) {
    cat("Não foi possível ler colunas suficientes para testar transposição.\n")
    return(lista_var_sinal)
  }
  

  # O TESTE DA PROGRESSÃO FÍSICA

  coluna_1 <- suppressWarnings(as.numeric(tabela_amostra[, 1]))
  coluna_1 <- coluna_1[!is.na(coluna_1)]
  
  linha_1 <- suppressWarnings(as.numeric(tabela_amostra[1, -1]))
  linha_1 <- linha_1[!is.na(linha_1)]
  

  is_col_ordered <- length(coluna_1) >= 3 && (all(diff(coluna_1) > 0) || all(diff(coluna_1) < 0))
  is_row_ordered <- length(linha_1) >= 3 && (all(diff(linha_1) > 0) || all(diff(linha_1) < 0))
  
  if (is_col_ordered && !is_row_ordered) {
    lista_var_sinal$matriz_transposta <- FALSE
    
  } else if (is_row_ordered && !is_col_ordered) {
    lista_var_sinal$matriz_transposta <- TRUE

    
  } else {
    if (ncol(tabela_amostra) > 100 && nrow(tabela_amostra) <= 50) {
      lista_var_sinal$matriz_transposta <- TRUE
    } else {
      lista_var_sinal$matriz_transposta <- FALSE
    }
  }
  
  return(lista_var_sinal)
}
#-------------------------------------------------------------------------
detetar_cabecalho <- function(file_path, lista_var_sinal) {
  
  if (lista_var_sinal$extensao_ficheiro == "spc") return(lista_var_sinal)
  
  tabela_amostra <- tryCatch({
    read.table(file_path, sep = lista_var_sinal$separador_colunas, dec = lista_var_sinal$separador_decimal, 
               skip = lista_var_sinal$linhas_metadados_topo, nrows = 50, header = FALSE, fill = TRUE)
  }, error = function(e) return(NULL))
  
  if (is.null(tabela_amostra) || ncol(tabela_amostra) < 2) return(lista_var_sinal)
  
  if (!lista_var_sinal$matriz_transposta) {
    linha_1_amostras <- suppressWarnings(as.numeric(tabela_amostra[1, -1]))
    # Se houver NAs (texto que não virou número), o cabeçalho ESTÁ presente
    lista_var_sinal$cabecalho_presente <- (sum(is.na(linha_1_amostras)) > 0)
  } else {
    coluna_1_amostras <- suppressWarnings(as.numeric(tabela_amostra[-1, 1]))
    lista_var_sinal$cabecalho_presente <- (sum(is.na(coluna_1_amostras)) > 0)
  }
  
  if (lista_var_sinal$cabecalho_presente) {
    cat("Cabeçalho Presente.\n")
  } else {
    cat("Cabeçalho Ausente.\n")
  }
  
  return(lista_var_sinal)
}

#---------------------------------------------------------------------

detetar_eixo_x <- function(file_path, lista_var_sinal) {
  
  if (lista_var_sinal$extensao_ficheiro == "spc") return(lista_var_sinal)
  
  tabela_amostra <- tryCatch({
    read.table(file_path, sep = lista_var_sinal$separador_colunas, dec = lista_var_sinal$separador_decimal, 
               skip = lista_var_sinal$linhas_metadados_topo, nrows = 50, header = FALSE, fill = TRUE)
  }, error = function(e) return(NULL))
  
  if (is.null(tabela_amostra) || ncol(tabela_amostra) < 2) return(lista_var_sinal)
  
  if (!lista_var_sinal$matriz_transposta) {
    eixo_x <- suppressWarnings(as.numeric(tabela_amostra[, 1]))
  } else {
    eixo_x <- suppressWarnings(as.numeric(tabela_amostra[1, ]))
  }
  
  eixo_x <- eixo_x[!is.na(eixo_x)] 
  
  if (length(eixo_x) >= 2) {
    if (eixo_x[1] > eixo_x[length(eixo_x)]) {
      lista_var_sinal$eixo_x_invertido <- TRUE
      cat("Os valores estão em ordem decrescente.\n")
    } else {
      lista_var_sinal$eixo_x_invertido <- FALSE
      cat("Os valores estão em ordem crescente.\n")
    }
  }
  return(lista_var_sinal)
}
