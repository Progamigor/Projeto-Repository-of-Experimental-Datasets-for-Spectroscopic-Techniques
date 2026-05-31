# pipeline-etl-espectroscopia
Pipeline de harmonização de datasets espectroscópicos experimentais para o ecossistema specmine / WebSpecmine
## O que este código faz:
* **Separação Inteligente:** Isola os metadados (texto) da matriz numérica (valores dos espetros).
* **Harmonização:** Estrutura os dados sem corromper a matriz original. **Não é feita qualquer interpolação ou extrapolação** nesta fase, garantindo a integridade dos dados brutos.
* **Exportação:** Gera automaticamente:
  * Ficheiros `.csv` (dados e metadados separados) para leitura no WebSpecmine.
  * Ficheiros `.rds` para leitura direta no pacote `specmine`.

## Estrutura dos Scripts:
* **Parsers de Harmonização:** Scripts responsáveis por diagnosticar e limpar os ficheiros originais.

## Autores
* **Pedro Vieira** - Universidade do Minho
