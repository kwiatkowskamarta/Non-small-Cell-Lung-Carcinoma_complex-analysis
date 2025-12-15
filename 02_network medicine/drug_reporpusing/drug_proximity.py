import pandas as pd
import networkx as nx
import separation  
import proximity

## EXCERSICE 3.2

# --- 1. Load Data ---
ppi = pd.read_csv('data/PPI.csv')
ppi = ppi[['Symbol_A', 'Symbol_B']].dropna()
G = nx.from_pandas_edgelist(ppi, 'Symbol_A', 'Symbol_B')
G.remove_edges_from(nx.selfloop_edges(G))
dga = pd.read_csv('data/disease_gene.tsv', sep='\t')

# --- 2. Define Helper to get genes ---
def get_disease_genes(df, search_term):
    # Search by name if IDs are not specific enough in the context
    subset = df[df['diseaseName'].str.contains(search_term, case=False, na=False)]
    genes = subset['geneSymbol'].unique()
    # Filter for genes only present in the graph G
    valid_genes = set(genes) & set(G.nodes())
    return valid_genes

# --- 3. Retrieve Gene Sets ---
genes_nsclc = get_disease_genes(dga, "Non-Small Cell Lung")

# Overlapping Candidate: Small Cell Lung Carcinoma
genes_sclc = get_disease_genes(dga, "Small Cell Lung")

# Distant Candidate: Schizophrenia
genes_schizo = get_disease_genes(dga, "Schizophrenia")

print(f"NSCLC Genes in Network: {len(genes_nsclc)}")
print(f"SCLC Genes in Network: {len(genes_sclc)}")
print(f"Schizophrenia Genes in Network: {len(genes_schizo)}")

# --- 4. Compute Separation ---

if len(genes_nsclc) > 0 and len(genes_sclc) > 0:
    # Compute separation for Overlapping candidate
    sep_overlap = separation.get_separation(G, genes_nsclc, genes_sclc)
    print(f"\nSeparation (NSCLC <-> Small Cell Lung Ca): {sep_overlap:.8f}")
    
    if sep_overlap < 0:
        print("-> Result: Negative separation. The modules overlap (common neighborhood).")
    else:
        print("-> Result: Positive separation. The modules are topologically distinct.")
else:
    print("\nInsufficient genes to calculate overlap separation.")

if len(genes_nsclc) > 0 and len(genes_schizo) > 0:
    # Compute separation for Distant candidate
    sep_distant = separation.get_separation(G, genes_nsclc, genes_schizo)
    print(f"\nSeparation (NSCLC <-> Schizophrenia): {sep_distant:.4f}")
    
    if sep_distant > 0:
        print("-> Result: Positive separation. The modules are topologically distinct.")
    else:
        print("-> Result: Negative separation. The modules overlap.")
else:
    print("\nInsufficient genes to calculate distant separation.")


## EXCERSICE 3.3

# --- 1. Load and Filter Drug Data ---
dt = pd.read_csv('data/drug_target.csv')
dt_human = dt[dt['organism'] == 'Humans']

# Create a dictionary mapping drugs to their target genes (list)
drug_to_targets = dt_human.groupby('Name')['Gene_Target'].apply(list).to_dict()

# --- 2. Identify Potential Candidates ---
candidates = []

for drug, targets in drug_to_targets.items():
    # Find overlap between drug targets and NSCLC disease genes
    overlap = set(targets) & set(genes_nsclc)
    
    # Filter: Keep drugs with at least 2 targets in the disease module
    if len(overlap) >= 2:
        candidates.append({
            'Drug': drug,
            'Targets': list(targets),
            'Overlap_Count': len(overlap),
            'Overlap_Genes': list(overlap)
        })

# Convert to DataFrame and sort by overlap count
df_candidates = pd.DataFrame(candidates).sort_values(by='Overlap_Count', ascending=False)

print(f"Found {len(df_candidates)} drugs with multiple targets in the NSCLC module.")
print("Top 5 candidates by target overlap:")
print(df_candidates.head(5))

# --- 3. Compute Proximity for Top Candidates ---
selected_drugs = df_candidates['Drug'].head(3).tolist() 

print(f"Calculating proximity for: {selected_drugs}...\n")

for drug_name in selected_drugs:
    targets = drug_to_targets[drug_name]
    valid_targets = set(targets) & set(G.nodes())
    
    if len(valid_targets) == 0:
        continue
        
    # Calculate proximity using the provided module 
    results = proximity.get_proximity(G, genes_nsclc, valid_targets, sims=1000)
    
    z = results['z_score']
    pval = results['p_value']
    d_obs = results['proximity']
    
    print(f"--- {drug_name} ---")
    print(f"Observed Distance: {d_obs:.4f}")
    print(f"Z-score: {z:.4f}")
    print(f"P-value: {pval:.4e}")
    
    if z < -1.5:
        print("RESULT: Significant Proximity (Potential Repurposing Candidate)\n")
    else:
        print("RESULT: Not Significant (Distant from Disease Module)\n")