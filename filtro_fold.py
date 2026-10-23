import pandas as pd
import glob
import os

print("================================================================================")
print("             🧬 PIPELINE CHIP-SEQ: FILTRAGEM DE ALTA ESTRITUDE (RIGOROSO)       ")
print("================================================================================\n")

# 1. Localizar os arquivos de consenso
arquivos_consenso = glob.glob("Picos_Consenso_*.csv")

if not arquivos_consenso:
    print("❌ Nenhum arquivo 'Picos_Consenso_*.csv' foi encontrado!")
    exit()

# ==============================================================================
# 🎛️ PAINEL DE CONTROLE DE RIGOR (ESCOLHA O SEU NÍVEL DE EXIGÊNCIA)
# ==============================================================================
# Opções disponíveis:
# 'q98'    -> Ultra estrito: o Top 2% dos picos mais fortes 
# 'q95'    -> Estrito: Mantém apenas a elite (Top 5% dos picos mais altos)
# 'q90' -> Top 10% mas obriga o pico a ter um Fold Mínimo Absoluto (ex: 2.0)
TIPO_FILTRO = 'q90' 
FOLD_MINIMO_ABSOLUTO = 2.0  # Usado apenas se escolher o modo 'hybrid'
# ==============================================================================

print(f"🧠 Modo do Filtro Ativado: [{TIPO_FILTRO.upper()}]")
print(f"{'Condição':<20} | {'Picos Originais':<15} | {'Corte Fold':<12} | {'Picos Filtrados':<15} | {'Status'}")
print("-" * 85)

for arquivo in sorted(arquivos_consenso):
    condicao = arquivo.replace("Picos_Consenso_", "").replace(".csv", "")
    
    df = pd.read_csv(arquivo)
    total_original = len(df)
    
    if total_original == 0:
        print(f"{condicao:<20} | {total_original:<15} | {'-':<12} | {'0':<15} | ⚠️ Vazio")
        continue
        
    # Aplicação da estatística rigorosa
    if TIPO_FILTRO == 'q98':
        corte_dinamico = max(df['fold_enrichment'].quantile(0.98), FOLD_MINIMO_ABSOLUTO)
        df_filtrado = df[df['fold_enrichment'] >= corte_dinamico].copy()
        
    elif TIPO_FILTRO == 'q95':
        corte_dinamico = max(df['fold_enrichment'].quantile(0.95), FOLD_MINIMO_ABSOLUTO)
        df_filtrado = df[df['fold_enrichment'] >= corte_dinamico].copy()
        
    elif TIPO_FILTRO == 'q90':
        corte_q90 = df['fold_enrichment'].quantile(0.90)
        corte_dinamico = max(corte_q90, FOLD_MINIMO_ABSOLUTO)
        df_filtrado = df[df['fold_enrichment'] >= corte_dinamico].copy()

    total_filtrado = len(df_filtrado)
    
    # Salvar os picos que passaram no teste de elite
    arquivo_saida = f"Picos_Filtrados_{condicao}.csv"
    df_filtrado.to_csv(arquivo_saida, index=False)
    
    print(f"{condicao:<20} | {total_original:<15} | {corte_dinamico:<12.2f} | {total_filtrado:<15} | ✅ Filtrado")

print("\n================================================================================")
print("  🚀 TRIAGEM DE ELITE CONCLUÍDA: Arquivos prontos para o ChIPseeker! ")
print("================================================================================")
