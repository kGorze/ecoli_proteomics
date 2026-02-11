# ==============================================================================
# 06_ai_analysis.py - AI-Driven Proteomic Noise Analysis (Interactive)
# ==============================================================================
# This script uses ESM-2 (Evolutionary Scale Modeling) embeddings to analyze
# the latent space of E. coli background proteins.
#
# Logic:
# 1. Load differential analysis results (CSV).
# 2. Identify "Noisy" E. coli proteins (False Positives) vs "Stable" ones.
# 3. Fetch sequences + METADATA (Gene, Description) from UniProt.
# 4. Generate embeddings using a lightweight ESM-2 model.
# 5. Visualize the latent space (UMAP) INTERACTIVELY (Plotly).
# ==============================================================================

import os
import pandas as pd
import numpy as np
import requests
import torch
from transformers import AutoTokenizer, AutoModel
import umap
import plotly.express as px
import plotly.io as pio
import time

# --- CONFIGURATION ---
RESULTS_DIR = "results"
TABLES_DIR = os.path.join(RESULTS_DIR, "tables")
FIGURES_DIR = os.path.join(RESULTS_DIR, "figures")
INPUT_CSV = os.path.join(TABLES_DIR, "differential_results_vsn.csv")
OUTPUT_HTML = os.path.join(FIGURES_DIR, "ai_latent_space_interactive.html")

# ESM-2 Model (Tiny: 6 layers, 8M params) - fast on CPU
MODEL_NAME = "facebook/esm2_t6_8M_UR50D"

def setup_directories():
    if not os.path.exists(FIGURES_DIR):
        os.makedirs(FIGURES_DIR)

def load_data():
    if not os.path.exists(INPUT_CSV):
        print(f"ERROR: {INPUT_CSV} not found. Run 03_differential_analysis.R first.")
        exit(1)
    
    df = pd.read_csv(INPUT_CSV)
    print(f"Loaded {len(df)} proteins from {INPUT_CSV}")
    return df

def fetch_uniprot_data(protein_ids):
    """Fetches sequences, gene names, and descriptions from UniProt."""
    print(f"Fetching metadata for {len(protein_ids)} proteins from UniProt...")
    
    base_url = "https://rest.uniprot.org/uniprotkb/accessions"
    
    # Cleaning IDs
    clean_ids_map = {pid: pid.split(';')[0].split('-')[0] for pid in protein_ids}
    unique_clean_ids = list(set(clean_ids_map.values()))
    
    chunk_size = 50
    results_map = {} # acc -> {seq, gene, name}
    
    for i in range(0, len(unique_clean_ids), chunk_size):
        chunk = unique_clean_ids[i:i+chunk_size]
        params = {"accessions": ",".join(chunk)}
        
        try:
            r = requests.get(base_url, params=params)
            r.raise_for_status()
            data = r.json()
            
            for res in data.get("results", []):
                acc = res["primaryAccession"]
                seq = res["sequence"]["value"]
                
                # Extract Gene Name
                gene = "Unknown"
                if "genes" in res and len(res["genes"]) > 0:
                    if "geneName" in res["genes"][0]:
                        gene = res["genes"][0]["geneName"]["value"]
                
                # Extract Protein Name (Description)
                name = "Unknown Protocol"
                if "proteinDescription" in res:
                    rec = res["proteinDescription"].get("recommendedName", {})
                    if "fullName" in rec:
                        name = rec["fullName"]["value"]
                    else:
                        sub = res["proteinDescription"].get("submissionNames", [{}])[0]
                        name = sub.get("fullName", {}).get("value", "Unknown")

                results_map[acc] = {"sequence": seq, "gene": gene, "name": name}
                
            print(f"  Fetched {len(results_map)}/{len(unique_clean_ids)}...")
            time.sleep(0.5)
            
        except Exception as e:
            print(f"  Error fetching chunk {i}: {e}")
            
    return results_map, clean_ids_map

def get_embeddings(sequences_dict, model_name=MODEL_NAME):
    """Generates mean embeddings."""
    print(f"Loading AI Model: {model_name}...")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name)
    model.eval()
    
    embeddings = []
    valid_accs = []
    
    print("Generating embeddings...")
    with torch.no_grad():
        for idx, (acc, seq) in enumerate(sequences_dict.items()):
            if idx % 50 == 0: print(f"  Processing AI embeddings {idx}/{len(sequences_dict)}...")
            
            # Truncate
            seq = seq[:1000] 
            
            inputs = tokenizer(seq, return_tensors="pt", padding=False, truncation=True)
            outputs = model(**inputs)
            
            # Mean pooling excluding special tokens
            emb = outputs.last_hidden_state[0, 1:-1, :].mean(dim=0).numpy()
            
            embeddings.append(emb)
            valid_accs.append(acc)
            
    return np.array(embeddings), valid_accs

def run_pipeline():
    setup_directories()
    df = load_data()
    
    # 1. infer Species
    if 'Species' not in df.columns:
        df['Species'] = df['Protein_ID'].apply(lambda x: 'ECOLI' if '_ECOLI' in str(x) else ('UPS' if 'UPS' in str(x) else 'OTHER'))
    
    df_ecoli = df[df['Species'] == 'ECOLI'].copy()

    # Define "Noisy" based on pure DEVIATION (logFC), not just statistical significance.
    # Since E. coli is constant, logFC should be 0. Any deviation is "Technical Noise".
    
    df_ecoli['Noise_Metric'] = df_ecoli['logFC'].abs()
    df_ecoli = df_ecoli.sort_values('Noise_Metric', ascending=False)
    
    # Top deviation = "Noisy" (Hard to measure consistently)
    noisy_subset = df_ecoli.head(200).copy()
    noisy_subset['Status'] = "High Noise (Unstable)"
    noisy_subset['Is_Noisy'] = True
    
    # Bottom deviation = "Stable" (Easy to measure)
    stable_subset = df_ecoli.tail(200).copy()
    stable_subset['Status'] = "Low Noise (Stable)" 
    stable_subset['Is_Noisy'] = False
    
    print(f"  Selected top {len(noisy_subset)} most deviating proteins as 'Noisy'.")
    print(f"  Selected bottom {len(stable_subset)} least deviating proteins as 'Stable'.")
    
    combined = pd.concat([noisy_subset, stable_subset]).copy()
    
    # 4. Fetch Data
    prot_ids = combined['Protein_ID'].tolist()
    uniprot_data, clean_map = fetch_uniprot_data(prot_ids)
    
    # 5. Map back to DataFrame
    combined['Clean_ID'] = combined['Protein_ID'].map(clean_map)
    combined = combined[combined['Clean_ID'].isin(uniprot_data.keys())]
    
    combined['Gene_Name'] = combined['Clean_ID'].apply(lambda x: uniprot_data[x]['gene'])
    combined['Protein_Name'] = combined['Clean_ID'].apply(lambda x: uniprot_data[x]['name'])
    combined['Sequence'] = combined['Clean_ID'].apply(lambda x: uniprot_data[x]['sequence'])
    # Status is already set in the selection step
    # combined['Status'] = combined['Is_False_Positive'].apply(lambda x: "Noisy (False Positive)" if x else "Stable (Background)")
    
    # 6. Embeddings
    seqs_to_embed = {row['Clean_ID']: row['Sequence'] for _, row in combined.iterrows()}
    embeds, valid_accs = get_embeddings(seqs_to_embed)
    
    # 7. UMAP
    print("Running UMAP dimensionality reduction...")
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42, n_jobs=1)
    embedding_2d = reducer.fit_transform(embeds)
    
    # Align DataFrame with embeddings (dict iteration order is preserved in Py3.7+, but safety first)
    # We iterated items() in get_embeddings, so valid_accs matches embeds order
    # Let's re-index DataFrame to match valid_accs
    combined = combined.set_index('Clean_ID').loc[valid_accs].reset_index()
    
    combined['UMAP1'] = embedding_2d[:, 0]
    combined['UMAP2'] = embedding_2d[:, 1]
    
    # 8. Interactive Plotting with Plotly
    print("Generating Interactive Plot...")
    
    # Highlight specific groups manually for the "Story"
    def highlight_biology(row):
        desc = str(row['Protein_Name']).lower() + " " + str(row['Gene_Name']).lower()
        if 'ribosom' in desc: return "Ribosome"
        if 'flagel' in desc: return "Flagella"
        if 'chaperon' in desc or 'groel' in desc or 'dnaj' in desc: return "Chaperone"
        if 'membran' in desc: return "Membrane"
        return "Other"

    # Ensure we don't have SettingWithCopyWarning
    combined = combined.copy()
    combined['Biology'] = combined.apply(highlight_biology, axis=1)
    
    fig = px.scatter(
        combined, 
        x='UMAP1', y='UMAP2',
        color='Status',
        symbol='Biology', # Different shapes for biological categories
        hover_data=['Protein_ID', 'Gene_Name', 'Protein_Name', 'logFC', 'Noise_Metric'],
        title='AI-Driven Proteomics: Latent Space of Technical Noise (High vs Low Deviation)',
        template='plotly_white',
        size_max=12
    )
    
    fig.update_traces(marker=dict(size=10, line=dict(width=1, color='DarkSlateGrey')))
    fig.write_html(OUTPUT_HTML)
    
    print(f"Interactive plot saved to: {OUTPUT_HTML}")

if __name__ == "__main__":
    run_pipeline()
