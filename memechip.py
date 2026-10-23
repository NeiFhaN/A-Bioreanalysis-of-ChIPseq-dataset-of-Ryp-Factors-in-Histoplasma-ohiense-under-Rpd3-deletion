import os
import subprocess
from pathlib import Path

# =============================================================================
# CONFIGURAÇÕES DE DIRETÓRIOS E ARQUIVOS 
# =============================================================================
BASE_DIR = Path(".")  # Onde os seus arquivos estão salvos
OUTPUT_DIR = BASE_DIR / "meme_chip_results_500bp" # Novo diretório de output solicitado
THREADS = 4 

# Nome do genoma completo para criar o background
GENOME_FASTA = BASE_DIR / "FungiDB-64_HcapsulatumG217B_Genome.fasta"
BACKGROUND_FILE = BASE_DIR / "hcapsulatum_genome_order2.bg"

# Banco de dados de motivos
# (Mantenha o nome do arquivo .meme que você está usando na sua pasta)
MOTIF_DB = BASE_DIR / "A_B.meme"

# =============================================================================
# Fatores e condições (em minúsculo para bater com os arquivos gerados)
# =============================================================================
ryp_factors = ["ryp1", "ryp2", "ryp3"]
conditions = ["wt", "rpd3"]

# Criar pasta de resultados
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# TRAVA DE SEGURANÇA ANTIERRO HUMANO
# -----------------------------------------------------------------------------
if not MOTIF_DB.exists():
    print("\n" + "="*80)
    print(f"[ERRO CRÍTICO] O arquivo do banco de motivos não foi encontrado: {MOTIF_DB.name}")
    print(f"Por favor, certifique-se de que ele está salvo na mesma pasta que este script!")
    print("="*80 + "\n")
    raise FileNotFoundError(f"Arquivo ausente: {MOTIF_DB.name}")

# =============================================================================
# PASSO 0: CRIAÇÃO DO ARQUIVO DE BACKGROUND GENÔMICO
# =============================================================================
if not BACKGROUND_FILE.exists():
    print(f"-> Gerando modelo de background a partir do genoma completo...")
    if not GENOME_FASTA.exists():
        raise FileNotFoundError(f"Erro crítico: O arquivo do genoma {GENOME_FASTA} não foi encontrado.")
    
    cmd_bg = ["fasta-get-markov", "-m", "2", str(GENOME_FASTA), str(BACKGROUND_FILE)]
    
    try:
        subprocess.run(cmd_bg, check=True)
        print(f"[SUCESSO] Arquivo de background criado: {BACKGROUND_FILE}\n")
    except subprocess.CalledProcessError as e:
        print(f"[ERRO] Falha ao rodar fasta-get-markov.")
        exit(1)
else:
    print(f"-> Arquivo de background {BACKGROUND_FILE} já existe. Pulando etapa.\n")


# =============================================================================
# FASE 1: ANÁLISE GERAL (TODOS OS PICOS - 500pb)
# =============================================================================
print("=== FASE 1: Iniciando Análise Geral (Todos os Picos - 500pb) ===")
for ryp in ryp_factors:
    for cond in conditions:
        
        # ---------------------------------------------------------------------
        # NOVO PADRÃO DE NOME DE ARQUIVO:
        # Exatamente como: wt_ryp1_500pb.fasta, rpd3_ryp2_500pb.fasta
        # ---------------------------------------------------------------------
        fasta_name = f"{cond}_{ryp}_500pb.fasta"
        
        fasta_file = BASE_DIR / fasta_name
        out_dir = OUTPUT_DIR / f"{cond}_{ryp}_500bp"
        
        if not fasta_file.exists():
            print(f"[AVISO] Arquivo não encontrado: {fasta_name}. Pulando...")
            continue
        
        print(f"Processando análise global para: {cond}_{ryp}...")
        
        cmd = [
            "meme-chip",
            "-oc", str(out_dir),
            "-bfile", str(BACKGROUND_FILE),
            "-meme-p", str(THREADS),
            "-db", str(MOTIF_DB),
            str(fasta_file)
        ]
        
        try:
            subprocess.run(cmd, check=True)
            print(f"[OK] Concluído: {fasta_name}\n")
        except subprocess.CalledProcessError as e:
            print(f"[ERRO] Falha ao rodar o MEME-chip para a condição {cond}_{ryp}. Detalhes acima.")

print(f">>> Todas as análises foram concluídas! Verifique a pasta '{OUTPUT_DIR.name}'.")
