import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import random
import numpy as np
from scipy.stats import norm

#data upload
ppi_data = pd.read_csv('C:/Users/kwiat/OneDrive/Pulpit/politechnika/magisterskie/sem2/complex data/network medicine/PPI.csv')
disease_data = pd.read_csv('/Users/kwiat/OneDrive/Pulpit/politechnika/magisterskie/sem2/complex data/network medicine/disease_gene.tsv', sep='\t')

# lowering disease names
disease_data['diseaseName'] = disease_data['diseaseName'].str.lower()

# filtering disease_gene associations
if 'DiseaseType' in disease_data.columns:
    disease_data = disease_data[~disease_data['DiseaseType'].isin(['group', 'phenotype'])]

# keeping only those that have at least 10 associated genes
disease_counts = disease_data['diseaseName'].value_counts()
valid_diseases = disease_counts[disease_counts >= 10].index
disease_data = disease_data[disease_data['diseaseName'].isin(valid_diseases)]

# our target disease
target_disease = 'non-small cell lung carcinoma'

#downloading list of genes for the NSCLC
disease_genes = disease_data[disease_data['diseaseName'] == target_disease]['geneSymbol'].unique()
print(f"Number of genes connected to the disease: {len(disease_genes)}")

# Interactome construction
G = nx.from_pandas_edgelist(ppi_data, 'Symbol_A', 'Symbol_B')
G.remove_edges_from(nx.selfloop_edges(G)) #removing self-loops

#filtering disease genes - keep only those present in our PPI network
valid_genes = [gene for gene in disease_genes if gene in G.nodes()]
print(f"Number of disease genes present in the PPI network: {len(valid_genes)}")

# helper function to calculate the size of the Largest Connected Component (LCC)
def calculate_lcc_size(graph, nodes):
    subgraph = graph.subgraph(nodes)
    if len(subgraph) > 0:
        # find the largest connected component
        largest_cc = max(nx.connected_components(subgraph), key=len)
        return len(largest_cc)
    return 0

#calculating LCC for our disease (observed value)
observed_lcc = calculate_lcc_size(G, valid_genes)
print(f"Size of the disease module (LCC): {observed_lcc}")

#random simulation
n_simulations = 1000
random_lcc_values = []
all_nodes = list(G.nodes())

print("Running random simulation...")
for i in range(n_simulations):
    # randomly select the same number of genes as the disease has
    random_genes = random.sample(all_nodes, len(valid_genes))
    rand_lcc = calculate_lcc_size(G, random_genes)
    random_lcc_values.append(rand_lcc)

#statistics: z-score and p-value
mean_random = np.mean(random_lcc_values)
std_random = np.std(random_lcc_values)
z_score = (observed_lcc - mean_random) / std_random

p_value = 1 - norm.cdf(z_score)
print(f"Random mean: {mean_random:.2f}")
print(f"Z-score: {z_score:.2f}")
print(f"P-value: {p_value:.10f}")

# histogram
plt.figure(figsize=(10, 6))
plt.hist(random_lcc_values, bins=20, color='skyblue', edgecolor='black', alpha=0.7, label='Random modules')
plt.axvline(observed_lcc, color='red', linestyle='dashed', linewidth=2, label=f'Our disease (LCC={observed_lcc})')
plt.title(f'Disease Module Significance (Z-score: {z_score:.2f})')
plt.xlabel('LCC Size')
plt.ylabel('Frequency')
plt.legend()
plt.grid(axis='y', alpha=0.5)
plt.savefig('disease_module_significance.png')
plt.show()

# saving file for Gephi
# creating a subgraph only from genes forming the disease module (LCC)
subgraph = G.subgraph(valid_genes)
if len(subgraph) > 0:
    lcc_nodes = max(nx.connected_components(subgraph), key=len)
    disease_module_graph = G.subgraph(lcc_nodes)
    nx.write_gexf(disease_module_graph, "nsclc_disease_module.gexf")
    print("File 'nsclc_disease_module.gexf' has been saved. Open it in Gephi.")
else:
    print("Disease module not found.")