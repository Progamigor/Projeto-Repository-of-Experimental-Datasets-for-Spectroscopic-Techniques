source("modulos/modulo1_diagnostico.R")
source("modulos/modulo2_execucao.R")
cat("-------------------------------------------------------------------------\n")
cat("---------------------------Fase 1 - Diagnóstico--------------------------\n")
cat("-------------------------------------------------------------------------\n")

tipo_dados <- escolha_tipo_dados()

ficheiro_alvo <- escolha_de_ficheiro()

main_lista_var <- variaveis_sinalizadoras()

main_lista_var <- detetar_formato(file_path = ficheiro_alvo, lista_var_sinal = main_lista_var)

main_lista_var <- detetar_delimitador(file_path = ficheiro_alvo, lista_var_sinal = main_lista_var)

main_lista_var <- detetar_linhas_metadados(file_path = ficheiro_alvo, lista_var_sinal = main_lista_var)

main_lista_var <- detetar_transposicao(file_path = ficheiro_alvo, lista_var_sinal = main_lista_var)

main_lista_var <- detetar_cabecalho(file_path = ficheiro_alvo, lista_var_sinal = main_lista_var)

main_lista_var <- detetar_eixo_x(file_path = ficheiro_alvo, lista_var_sinal = main_lista_var)

print(main_lista_var)

cat("-------------------------------------------------------------------------\n")
cat("-----------------------------Fase 2 - Execução---------------------------\n")
cat("-------------------------------------------------------------------------\n")

metadados_extraidos <- extrair_metadados_globais(file_path = ficheiro_alvo, lista_var_modulo1 = main_lista_var)

objeto_execucao <- extrator_dados_brutos(file_path = ficheiro_alvo, lista_var_modulo1 = main_lista_var)

objeto_execucao <- tranformar_transposta(objeto_dados = objeto_execucao, lista_var_modulo1 = main_lista_var)

objeto_execucao <- definir_eixos(objeto_dados = objeto_execucao)

objeto_execucao <- inverter_eixo_x(objeto_dados = objeto_execucao, lista_var_modulo1 = main_lista_var)

objeto_execucao <- sanetizar_nomes(objeto_dados = objeto_execucao)

dataset_final <- empacotar_dados(objeto_dados = objeto_execucao, tipo = tipo_dados, metadados_globais = metadados_extraidos)

caminho_ficheiro_csv <- exportar_dataset(dataset_specmine = dataset_final, file_path = ficheiro_alvo)


