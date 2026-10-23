# ================================================================================
#             🧬 PIPELINE CHIP-SEQ: ANOTAÇÃO GENÔMICA COM CHIPSEEKER (SUMMIT V2)
# ================================================================================
cat("\n[BANCA] Inicializando o ecossistema do Bioconductor...\n")

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(clusterProfiler)
  library(ggplot2)
})

# 1. PARÂMETROS CRÍTICOS (O que a banca vai avaliar)
ARQUIVO_GFF <- "FungiDB-68_HcapsulatumG217B.gff"

# Janela Promotora Adaptada para Fungos (Compacta e Estrita)
PROMOTER_UPSTREAM   <- 1000
PROMOTER_DOWNSTREAM <- 50

# Lista das 6 condições biológicas geradas pelo Python
condicoes <- c("WT_Ryp1", "WT_Ryp2", "WT_Ryp3", "rpd3_Ryp1", "rpd3_Ryp2", "rpd3_Ryp3")

# ================================================================================
# 2. CONSTRUÇÃO DO BANCO DE DADOS DE ANOTAÇÃO (TxDb)
# ================================================================================
cat("[BANCA] Construindo base de dados genômica a partir do GFF...\n")
if (!file.exists(ARQUIVO_GFF)) {
  stop(paste("❌ Erro: Arquivo GFF não encontrado:", ARQUIVO_GFF))
}

txdb <- makeTxDbFromGFF(ARQUIVO_GFF, format = "gff")

# ================================================================================
# 3. CARREGAMENTO E CONVERSÃO DOS PICOS FILTRADOS (PRESERVAÇÃO DO SUMMIT)
# ================================================================================
cat("[BANCA] Carregando os picos filtrados com coordenadas de Summit...\n")

lista_granges <- list()

for (condicao in condicoes) {
  arquivo_csv <- paste0("Picos_Filtrados_", condicao, ".csv")
  
  if (!file.exists(arquivo_csv)) {
    cat(paste("⚠️ Aviso: Arquivo", arquivo_csv, "não encontrado. Pulando...\n"))
    next
  }
  
  # Ler tabela enriquecida do Python
  dados <- read.csv(arquivo_csv)
  
  if (nrow(dados) == 0) {
    cat(paste("⚠️ Aviso: Arquivo", arquivo_csv, "está vazio. Pulando...\n"))
    next
  }
  
  # MODELAGEM DE ALTA PRECISÃO: Centraliza o GRanges no abs_summit 
  # para que o ChIPseeker meça a distância TSS do ponto exato de ligação.
  # As fronteiras originais (start/end) e estatísticas são salvas como metadados.
  gr <- GRanges(
    seqnames = dados$chr,
    ranges = IRanges(start = dados$abs_summit, end = dados$abs_summit),
    strand = "*",
    pico_start = dados$start,
    pico_end = dados$end,
    fold_enrichment = dados$fold_enrichment,
    log10_qvalue = dados$log10_qvalue,
    peak_name = dados$name
  )
  
  lista_granges[[condicao]] <- gr
}

if (length(lista_granges) == 0) {
  stop("❌ Erro: Nenhum arquivo 'Picos_Filtrados_*.csv' foi carregado com sucesso.")
}

# ================================================================================
# 4. ANOTAÇÃO DOS PICOS VIA CHIPSEEKER
# ================================================================================
cat("[BANCA] Mapeando picos nos elementos genômicos via ChIPseeker...\n")

lista_anotada <- list()

for (nome in names(lista_granges)) {
  cat(paste(" -> Anotando de forma cirúrgica pelo Summit:", nome, "\n"))
  
  annot <- annotatePeak(
    lista_granges[[nome]],
    tssRegion = c(-PROMOTER_UPSTREAM, PROMOTER_DOWNSTREAM),
    TxDb = txdb,
    verbose = FALSE
  )
  
  lista_anotada[[nome]] <- annot
  
  # Salvar a tabela final detalhada em CSV (Pronta para o Excel e discussões)
  df_saida <- as.data.frame(annot)
  write.csv(df_saida, paste0("Anotacao_Final_Sestreito", nome, ".csv"), row.names = FALSE)
}

# Salva o ambiente completo para o script de gráficos ler depois
cat("[Banca] Salvando objetos processados em 'dados_anotados.RData'...\n")
save(lista_anotada, lista_granges, txdb, PROMOTER_UPSTREAM, PROMOTER_DOWNSTREAM, file = "dados_anotados.RData")
cat("🎉 Pronto! Os dados estão consolidados e prontos para a geração de gráficos.\n")

cat("\n================================================================================")
cat("\n 🚀 SCRIPT CONCLUÍDO: Dados processados com rigor científico máximo! \n")
cat("================================================================================\n")
