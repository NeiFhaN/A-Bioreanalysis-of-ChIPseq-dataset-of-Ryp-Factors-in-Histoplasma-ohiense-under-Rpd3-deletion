# ===============================================================================
# 📊 SCRIPT ATUALIZADO: VISUALIZAÇÃO GRÁFICA DIRETO DOS ARQUIVOS CSV
# ===============================================================================
cat("\n====================================================================\n")
cat("   🧬 INICIALIZANDO ECOSSISTEMA GRÁFICO (RECONSTRUÇÃO VIA CSV)\n")
cat("====================================================================\n\n")

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  library(ggplot2)
  library(GenomicRanges)
})

# ===============================================================================
# CONFIGURAÇÕES INICIAIS
# ===============================================================================
# Defina o tamanho da janela do promotor que você usou originalmente 
# (Ex: 3000bp pra cima e pra baixo do TSS). Ajuste se necessário.
PROMOTER_UPSTREAM <- 1000
PROMOTER_DOWNSTREAM <- 50

# Configuração do GFF local
ARQUIVO_GFF <- "FungiDB-68_HcapsulatumG217B.gff"
if (!file.exists(ARQUIVO_GFF)) {
  stop(paste("❌ Erro: O arquivo GFF '", ARQUIVO_GFF, "' não foi encontrado nesta pasta!"))
}

cat("🧬 Construindo ponteiros ativos do TxDb a partir do GFF...\n")
txdb <- makeTxDbFromGFF(ARQUIVO_GFF)
cat("✅ Conexão com o genoma estabelecida com sucesso!\n\n")

# ===============================================================================
# 1. RECONSTRUINDO OS DADOS A PARTIR DOS CSVS
# ===============================================================================
cat("⏳ Lendo arquivos CSV e recriando objetos genômicos (GRanges e csAnno)...\n")

condicoes <- c("WT_Ryp1", "WT_Ryp2", "WT_Ryp3", "rpd3_Ryp1", "rpd3_Ryp2", "rpd3_Ryp3")
lista_granges <- list()
lista_anotada <- list()

# Procura todos os arquivos .csv no diretório atual
arquivos_csv <- list.files(pattern = "\\.csv$")

for (cond in condicoes) {
  # Encontra qual csv pertence a essa condição (ignorando diferenças sutis de nome)
  arq_correspondente <- arquivos_csv[grepl(cond, arquivos_csv)]
  
  if (length(arq_correspondente) > 0) {
    arq <- arq_correspondente[1] # Pega o primeiro que der match
    cat("   -> Processando condição:", cond, "a partir de", arq, "\n")
    
    # Lê o CSV
    df <- read.csv(arq)
    
    # Reconstrói as coordenadas (usando pico_start e pico_end que possuem as larguras reais)
    if("pico_start" %in% colnames(df) && "pico_end" %in% colnames(df)) {
      gr <- GRanges(seqnames = as.character(df$seqnames),
                    ranges = IRanges(start = df$pico_start, end = df$pico_end))
    } else {
      # Fallback caso a coluna tenha nome diferente
      gr <- GRanges(seqnames = as.character(df$seqnames),
                    ranges = IRanges(start = df$start, end = df$end))
    }
    
    lista_granges[[cond]] <- gr
    
    # Recria o objeto csAnno anotando o pico novamente no genoma
    cat("      Anotando picos em relação ao genoma...\n")
    anno <- annotatePeak(gr, TxDb = txdb, 
                         tssRegion = c(-PROMOTER_UPSTREAM, PROMOTER_DOWNSTREAM), 
                         verbose = FALSE)
    lista_anotada[[cond]] <- anno
    
  } else {
    cat("   [AVISO] Arquivo CSV não encontrado para a condição:", cond, "\n")
  }
}

cat("✅ Dados processados e memórias reconstruídas com sucesso!\n\n")

# ===============================================================================
# CATEGORIA 1: VISUALIZAÇÃO DA COBERTURA GENÔMICA (WT PRETO / NOVOS EM VERMELHO)
# ===============================================================================
cat("[📊 CATEGORIA 1] Gerando Covplots com legenda em linhas simples...\n")

# Tema padronizado base
tema_covplot <- theme_bw() + theme(
  plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
  strip.text.y = element_text(angle = 0, size = 8, margin = margin(l = 8, r = 8)),
  axis.text.y = element_text(size = 8, margin = margin(r = 5)),
  panel.spacing = unit(1.2, "lines")
)

# Definimos as proteínas que vamos comparar
proteinas <- c("Ryp1", "Ryp2", "Ryp3")

for (prot in proteinas) {
  nome_wt <- paste0("WT_", prot)
  nome_mut <- paste0("rpd3_", prot)
  
  if (nome_wt %in% names(lista_granges) && nome_mut %in% names(lista_granges)) {
    cat("   -> Calculando matemática genômica (intersect) para:", prot, "\n")
    
    # 1. Isolar os dados genômicos
    gr_wt <- lista_granges[[nome_wt]]
    gr_mut <- lista_granges[[nome_mut]]
    
    # Simplificar os nomes dos cromossomos
    seqlevels(gr_wt) <- sub(".*0+([1-9][0-9]*)$", "Cromossomo \\1", seqlevels(gr_wt))
    seqlevels(gr_mut) <- sub(".*0+([1-9][0-9]*)$", "Cromossomo \\1", seqlevels(gr_mut))
    
    # 2. Matemática de Conjuntos
    picos_compartilhados <- intersect(gr_wt, gr_mut)
    exclusivos_wt        <- setdiff(gr_wt, gr_mut)
    exclusivos_mut       <- setdiff(gr_mut, gr_wt)
    
    # ========================================================================
    # 3. Gerar Covplot para WT (100% PRETO E ZERO LEGENDA)
    # ========================================================================
    pdf_wt <- paste0("Grafico_1_Covplot_Baseline_", nome_wt, ".pdf")
    cat("      Criando Baseline WT:", pdf_wt, "\n")
    pdf(pdf_wt, width = 12, height = 8)
    
    p_wt <- covplot(list("Picos WT" = exclusivos_wt,
                         "Picos WT " = picos_compartilhados),
                    title = paste("Cobertura Genômica -", nome_wt, "(Controle WT)")) +
            xlab("Tamanho do Cromossomo (pb)") + ylab("Cobertura") +
            scale_color_manual(values = c("Picos WT" = "black", "Picos WT " = "black")) +
            tema_covplot +
            theme(legend.position = "none") 
            
    print(p_wt)
    dev.off()
    
    # ========================================================================
    # 4. Gerar Covplot para o Mutante (CORES E BORDAS CORRIGIDAS)
    # ========================================================================
    pdf_mut <- paste0("Grafico_1_Covplot_Destacado_", nome_mut, ".pdf")
    cat("      Criando Mutante Destacado:", pdf_mut, "\n")
    pdf(pdf_mut, width = 12, height = 8)
    
    p_mut <- covplot(list("Picos Novos" = exclusivos_mut,
                          "Picos Iguais ao WT" = picos_compartilhados),
                     title = paste("Cobertura Genômica -", nome_mut, "(Destaque para Novos Picos)")) +
            xlab("Tamanho do Cromossomo (pb)") + ylab("Cobertura") +
            # TRAVAMOS TANTO A COR DA LINHA QUANTO O PREENCHIMENTO (FILL):
            scale_color_manual(values = c("Picos Novos" = "red", "Picos Iguais ao WT" = "black")) +
            scale_fill_manual(values = c("Picos Novos" = "red", "Picos Iguais ao WT" = "black")) +
            tema_covplot +
            theme(
              legend.position = "bottom",
              legend.title = element_blank(),
              legend.text = element_text(size = 10, face = "bold"),
              legend.key = element_blank(),
              legend.key.size = unit(0.5, "cm")
            ) +
            # A SUA SOLUÇÃO MANTIDA AQUI:
            guides(color = guide_legend(override.aes = list(linewidth = 0, linetype = 1)))
            
    print(p_mut)
    dev.off()
# ===============================================================================
# CATEGORIA 2: PERFIS DE LIGAÇÃO EM RELAÇÃO AO TSS (RESOLUÇÃO ESPACIAL)
# ===============================================================================
cat("[📊 CATEGORIA 2] Extraindo coordenadas reais dos promotores...\n")
promotores_granges <- getPromoters(TxDb = txdb, 
                                   upstream = PROMOTER_UPSTREAM, 
                                   downstream = PROMOTER_DOWNSTREAM)

cat("[📊 CATEGORIA 2] Calculando matriz de tags genômicas...\n")
tagMatrixList <- lapply(lista_granges, getTagMatrix, windows = promotores_granges)

# Limite máximo correto solicitado para o eixo Y
LIMITE_Y_MAX <- 0.0025

# Gráfico para WT
amostras_wt <- names(tagMatrixList)[names(tagMatrixList) %in% c("WT_Ryp1", "WT_Ryp2", "WT_Ryp3")]
if(length(amostras_wt) > 0) {
  cat("   -> Criando perfil individualizado para o grupo WT...\n")
  pdf("Grafico_2A_Perfil_Medio_TSS_Condicao_WT.pdf", width = 8, height = 5)
  print(
    plotAvgProf(tagMatrixList[amostras_wt], xlim = c(-PROMOTER_UPSTREAM, PROMOTER_DOWNSTREAM),
                xlab = "Distância em relação ao TSS (pb)", 
                ylab = "Frequência Média") + 
    theme_bw() +
    ggtitle("Assinatura de Ligação das Proteínas Ryp (Cenário WT)") +
    labs(color = "Amostra") + 
    # CORREÇÃO CRÍTICA: Limite fixado exatamente em 0.0025
    coord_cartesian(ylim = c(0, LIMITE_Y_MAX)) + 
    theme(
          plot.title = element_text(hjust = 0.5, face="bold", size=12),
          legend.position = "right",
          # AJUSTE FINO DO EIXO X: Inclinação de 90° e afastamento vertical (vjust) para não aglomerar
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 9),
          axis.title.x = element_text(margin = margin(t = 10)) # Dá um respiro para o título do eixo X
    )
  )
  dev.off()
}

# Gráfico para rpd3
amostras_rpd3 <- names(tagMatrixList)[names(tagMatrixList) %in% c("rpd3_Ryp1", "rpd3_Ryp2", "rpd3_Ryp3")]
if(length(amostras_rpd3) > 0) {
  cat("   -> Criando perfil individualizado para o grupo mutante rpd3...\n")
  pdf("Grafico_2B_Perfil_Medio_TSS_Condicao_rpd3.pdf", width = 8, height = 5)
  print(
    plotAvgProf(tagMatrixList[amostras_rpd3], xlim = c(-PROMOTER_UPSTREAM, PROMOTER_DOWNSTREAM),
                xlab = "Distância em relação ao TSS (pb)", 
                ylab = "Frequência Média") + 
    theme_bw() +
    ggtitle("Assinatura de Ligação das Proteínas Ryp (Cenário Mutante Δrpd3)") +
    labs(color = "Amostra") + 
    # CORREÇÃO CRÍTICA: Limite fixado exatamente em 0.0025
    coord_cartesian(ylim = c(0, LIMITE_Y_MAX)) + 
    theme(
          plot.title = element_text(hjust = 0.5, face="bold", size=12),
          legend.position = "right",
          # AJUSTE FINO DO EIXO X: Inclinação de 90° e afastamento vertical (vjust) para não aglomerar
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 9),
          axis.title.x = element_text(margin = margin(t = 10)) # Dá um respiro para o título do eixo X
    )
  )
  dev.off()
}

# ===============================================================================
# CATEGORIA 3: ANOTAÇÃO E CONTEXTO GENÔMICO (ESTRUTURA DA BARRA)
# ===============================================================================
cat("[📊 CATEGORIA 3] Gerando Distribuição Genômica (plotAnnoBar)...\n")
pdf("Grafico_3A_Distribuicao_Genomica_Barras.pdf", width = 10, height = 6)
print(
  plotAnnoBar(lista_anotada) + 
  ggtitle("Distribuição Genômica dos Picos (WT vs rpd3)") +
  ylab("Porcentagem (%)") +                    # Eixo Y original (que fica na horizontal)
  labs(fill = "Região") +                      
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title.y = element_text(margin = margin(r = 10))
  )
)
dev.off()

cat("[📊 CATEGORIA 3] Gerando Distribuição de Distância ao TSS (plotDistToTSS)...\n")
pdf("Grafico_3B_Distancia_Direta_TSS.pdf", width = 10, height = 6)
print(
  plotDistToTSS(lista_anotada) + 
  ggtitle("Proximidade dos Sítios de Ligação Ryp ao TSS") +
  ylab("Sítios de Ligação (%) (5' -> 3')") +   # Eixo Y original (que fica na horizontal)
  labs(fill = "Região") +                      
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title.y = element_text(margin = margin(r = 10))
  )
)
dev.off()

# ===============================================================================
# CATEGORIA 4: UPSET PLOT
# ===============================================================================
if("WT_Ryp1" %in% names(lista_anotada)) {
  cat("[📊 CATEGORIA 4] Gerando UpSet Plot de Anotação (Sobreposição de features)...\n")
  pdf("Grafico_4_UpSet_Intersecoes_WT_Ryp1.pdf", width = 10, height = 6)
  print(upsetplot(lista_anotada[["WT_Ryp1"]]))
  dev.off()
}

cat("\n====================================================================\n")
cat("🚀 SUCESSO: Todos os PDFs de nível de publicação foram restaurados!\n")
cat("====================================================================\n")
