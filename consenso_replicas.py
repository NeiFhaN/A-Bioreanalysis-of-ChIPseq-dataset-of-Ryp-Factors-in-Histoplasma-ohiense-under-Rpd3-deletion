import pandas as pd
import glob
import os

print("================================================================================")
print("       🧬 PIPELINE CHIP-SEQ: CONSENSO RIGOROSO (REPLICAS A, B, C & SUMMIT)     ")
print("================================================================================\n")

# 1. MAPEAMENTO INTELIGENTE POR PADRÃO DE NOME (Utiliza as réplicas A, B, C)
experimentos = {
    "WT_Ryp1":   ["*WT_Ryp1_A*", "*WT_Ryp1_B*", "*WT_Ryp1_C*"],
    "WT_Ryp2":   ["*WT_Ryp2_A*", "*WT_Ryp2_B*", "*WT_Ryp2_C*"],
    "WT_Ryp3":   ["*WT_Ryp3_A*", "*WT_Ryp3_B*", "*WT_Ryp3_C*"],
    "rpd3_Ryp1": ["*rpd3_Ryp1_A*", "*rpd3_Ryp1_B*", "*rpd3_Ryp1_C*"],
    "rpd3_Ryp2": ["*rpd3_Ryp2_A*", "*rpd3_Ryp2_B*", "*rpd3_Ryp2_C*"],
    "rpd3_Ryp3": ["*rpd3_Ryp3_A*", "*rpd3_Ryp3_B*", "*rpd3_Ryp3_C*"]
}

def ler_macs2_tdt(padrao_busca, nome_replica):
    """Localiza o arquivo pelo padrão glob, detecta o cabeçalho e normaliza as colunas."""
    arquivos = glob.glob(padrao_busca)
    if not arquivos:
        return pd.DataFrame()
    
    caminho = arquivos[0] # Pega o arquivo correspondente encontrado
    
    with open(caminho, 'r') as f:
        linhas = f.readlines()
        
    pular = 0
    for i, linha in enumerate(linhas):
        if 'chr' in linha.lower() and 'start' in linha.lower():
            pular = i
            break
            
    df = pd.read_csv(caminho, sep='\t', skiprows=pular)
    df.columns = df.columns.str.strip().str.lower() # Normalização preventiva
    
    # Dicionário de tradução para garantir compatibilidade universal
    col_mapping = {}
    for col in df.columns:
        if col in ['chr', 'chromosome', 'chrom']: col_mapping[col] = 'chr'
        elif col == 'start': col_mapping[col] = 'start'
        elif col == 'end': col_mapping[col] = 'end'
        elif col in ['abs_summit', 'summit', 'abs_summit(bp)']: col_mapping[col] = 'abs_summit'
        elif col in ['fold_enrichment', 'foldenrichment', 'fc', 'signal value']: col_mapping[col] = 'fold_enrichment'
        elif col in ['-log10(qvalue)', 'qvalue', 'log10_qvalue', '-log10(q-value)']: col_mapping[col] = 'log10_qvalue'
        elif col in ['name', 'peak_name', 'id']: col_mapping[col] = 'name'
        
    df = df.rename(columns=col_mapping)
    df['replica'] = nome_replica
    
    # Manter apenas as colunas essenciais exigidas pela banca
    colunas_finais = ['chr', 'start', 'end', 'abs_summit', 'fold_enrichment', 'log10_qvalue', 'name', 'replica']
    df = df[[c for c in colunas_finais if c in df.columns]]
    return df

print(f"{'Condição':<15} | {'Réplicas Encontradas':<22} | {'Picos Consenso (>=2 Reps)':<25}")
print("-" * 75)

for condicao, padroes in experimentos.items():
    dfs = []
    reps_encontradas = []
    
    for i, padrao in enumerate(padroes):
        nome_rep = f"Rep_{chr(65+i)}" # Rep_A, Rep_B, Rep_C
        # Busca flexível aceitando variações de extensão (.tdt.txt, .txt)
        df = ler_macs2_tdt(padrao, nome_rep)
        if not df.empty:
            dfs.append(df)
            reps_encontradas.append(chr(65+i))
            
    if len(dfs) < 2:
        print(f"{condicao:<15} | {f'Apenas {len(dfs)} rep(s)':<22} | ❌ Pulado (Amostras insuficientes)")
        continue
        
    # Agrupar e ordenar geograficamente
    df_completo = pd.concat(dfs, ignore_index=True)
    df_completo['start'] = df_completo['start'].astype(int)
    df_completo['end'] = df_completo['end'].astype(int)
    df_completo = df_completo.sort_values(by=['chr', 'start'])
    
    # Algoritmo de Agrupamento por Sobreposição
    clusters = []
    cluster_atual = []
    fim_maximo = -1
    chr_atual = None
    
    for _, linha in df_completo.iterrows():
        if linha['chr'] != chr_atual or list(linha)[1] > fim_maximo: # list(linha)[1] é o start
            if cluster_atual:
                clusters.append(pd.DataFrame(cluster_atual))
            cluster_atual = [linha]
            fim_maximo = linha['end']
            chr_atual = linha['chr']
        else:
            cluster_atual.append(linha)
            fim_maximo = max(fim_maximo, linha['end'])
            
    if cluster_atual:
        clusters.append(pd.DataFrame(cluster_atual))
        
    # Filtragem por voto de maioria (Pelo menos 2 réplicas precisam conter o pico)
    picos_consenso = []
    for cl in clusters:
        if cl['replica'].nunique() >= 2:
            # Preservação física absoluta: elege a linha íntegra do pico mais forte
            pico_campeao = cl.loc[cl['fold_enrichment'].idxmax()]
            picos_consenso.append(pico_campeao)
            
    if picos_consenso:
        df_consenso_final = pd.DataFrame(picos_consenso).drop(columns=['replica'])
        df_consenso_final.to_csv(f"Picos_Consenso_{condicao}.csv", index=False)
        reps_str = ",".join(reps_encontradas)
        print(f"{condicao:<15} | {f'Reps ({reps_str})':<22} | {len(df_consenso_final):<25}")

print("\n================================================================================")
print(" 🚀 ETAPA 1 CONCLUÍDA: Arquivos de consenso gerados com todas as colunas MACS2! ")
print("================================================================================")
